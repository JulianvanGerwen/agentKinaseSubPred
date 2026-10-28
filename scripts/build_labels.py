#!/usr/bin/env python
"""Build kinase–substrate labels: cache/labels/ks_pairs.parquet (+ ks_unmapped.parquet) and
results/labels/{label_summary,label_funnel,kinase_label_counts}.tsv.

Requires results/labels/kinase_list.tsv + kinase_synonyms.tsv (scripts/build_kinase_list.py) and
cache/sites/universe.parquet (scripts/build_site_universe.py). Rules: src/ksp/labels.py, docs/labels_universe.md.
Usage: python -I scripts/build_labels.py [--min-identity 11] [--omnipath-include-aggregators]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import io_psp, kinases, labels, sites  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-identity", type=int, default=11, help="ortholog 15-mer fallback: min identical residues")
    ap.add_argument("--omnipath-include-aggregators", action="store_true",
                    help="also count dbPTM / Li2012 as curated OmniPath sources")
    a = ap.parse_args()
    t0 = time.time()
    out_c, out_r = ROOT / "cache" / "labels", ROOT / "results" / "labels"
    out_c.mkdir(parents=True, exist_ok=True)
    st = sites.default_store()
    res = kinases.load_resolver()
    kl = res.kl
    uni = pd.read_parquet(ROOT / "cache" / "sites" / "universe.parquet", columns=["site_id", "in_U_obs"])
    in_uobs = set(uni.site_id[uni.in_U_obs])
    funnel = []

    # ------------------------------------------------------------------ PSP
    raw = io_psp.load_psp()
    funnel.append(("PSP", "rows", len(raw)))
    funnel.append(("PSP", "rows_non_STY_dropped", int((~raw.is_sty).sum())))
    d = raw[raw.is_sty].copy()
    d = labels.resolve_psp_kinases(d, res)
    for h, n in d.kin_map_how.value_counts().items():
        funnel.append(("PSP", f"kinase_map:{h}", int(n)))
    ps = io_psp.load_psp_sites()
    ps = ps[ps.organism.eq("human") & ps.is_sty & ps.mod_type.eq("p")]
    sg = labels.sgid_to_human_sites(st, ps, raw)
    m = labels.map_psp_substrates(d, st, sg, a.min_identity)
    first = m.drop_duplicates("_row")
    for h, n in first.sub_map_method.value_counts().items():
        funnel.append(("PSP", f"substrate_map_rows:{h}", int(n)))
    for h, n in first.loc[first.sub_map_method.eq("unmapped"), "sub_map_reason"].value_counts().items():
        funnel.append(("PSP", f"substrate_unmapped_rows:{h}", int(n)))
    funnel.append(("PSP", "rows_exploding_to_>1_human_site", int((m.groupby("_row").site_id.nunique() > 1).sum())))
    bad = m.kinase.isna() | m.site_id.isna()
    b = m[bad].drop_duplicates("_row")
    unm = [pd.DataFrame({"provenance": "PSP",
                         "reason": np.where(b.kinase.isna(), "kinase_" + b.kin_map_how.astype(str),
                                            "site_" + b.sub_map_reason.astype(str)),
                         "kinase_name": b.kin_gene + "|" + b.kin_organism, "substrate": b.sub_acc_full + "|" + b.sub_organism,
                         "site": b.mod_rsd, "detail": b.window_raw})]
    pp = labels.psp_pairs(m[~bad], kl)

    # ------------------------------------------------------------------ external
    curated = labels.OMNIPATH_CURATED + (labels.OMNIPATH_OTHER_AGGREGATORS if a.omnipath_include_aggregators else ())
    op, olog = labels.omnipath_pairs(ROOT / "external/omnipath/omnipath_enzsub.tsv", st, res, curated)
    up, ulog = labels.uniprot_pairs(ROOT / "external/uniprot/ptm_tsv/UP000005640_reviewed_ptm.tsv.gz", st, res)
    unm += [olog, ulog]
    for lg in [olog, ulog]:
        for (pv, r), n in lg.groupby(["provenance", "reason"]).size().items():
            funnel.append((pv, f"dropped:{r}", int(n)))
    funnel.append(("OmniPath", "pairs_mapped", len(op)))
    funnel.append(("UniProt", "pairs_mapped", len(up)))

    # ------------------------------------------------------------------ mask-only rows (DESIGN §16 A7)
    tm = olog[olog.reason.eq("no_curated_non_psp_source")]
    mo = labels.mask_only_pairs(pd.concat([pp, op, up], ignore_index=True), olog, ulog, st, res, textmining=tm)
    for pv, n in mo.groupby(["provenance", "mask_sensitivity_only"]).size().items():
        funnel.append(("MASK", f"mask_only_pairs:{pv[0]}:sensitivity_only={pv[1]}", int(n)))

    # ------------------------------------------------------------------ combine + flags
    allp = pd.concat([pp, op, up, mo], ignore_index=True)
    allp["mask_only"] = allp.mask_only.astype("boolean").fillna(False).astype(bool)
    allp["mask_sensitivity_only"] = allp.mask_sensitivity_only.astype("boolean").fillna(False).astype(bool)
    cls = dict(zip(kl.gene, kl.residue_class))
    inc = dict(zip(kl.gene, kl.include))
    allp["kinase_acc"] = allp.kinase.map(dict(zip(kl.gene, kl.uniprot_acc)))
    allp["acc"] = allp.site_id.str.split("_").str[0]
    allp["res"] = allp.site_id.str.extract(r"_([STY])\d+$")[0]
    allp["pos"] = allp.site_id.str.extract(r"_[STY](\d+)$")[0].astype(int)
    allp["kinase_included"] = allp.kinase.map(inc).fillna(False).astype(bool)
    allp["kinase_residue_class"] = allp.kinase.map(cls)
    allp["residue_in_class"] = labels.residue_in_class(allp.res, allp.kinase_residue_class)
    allp["in_universe"] = allp.site_id.isin(in_uobs)
    for c in ["native_human", "hmr", "sub_native_human", "in_vivo", "in_vitro", "in_vivo_native", "in_vitro_native",
              "in_vivo_hmr", "in_vitro_hmr", "autophos", "res_changed", "sgid_multi", "omnipath_lists_psp"]:
        allp[c] = allp[c].astype("boolean").fillna(False).astype(bool)
    ext = allp.provenance.isin(["OmniPath", "UniProt"])
    nonpsp = allp.provenance.ne("PSP")
    allp.loc[nonpsp, "autophos"] = allp.loc[nonpsp, "kinase_acc"].eq(allp.loc[nonpsp, "acc"])
    allp["species_origin"] = allp.species_origin.fillna("")
    psp_raw_pairs = set(zip(pp.kinase, pp.site_id))
    allp["in_psp_ks"] = [(k, s) in psp_raw_pairs for k, s in zip(allp.kinase, allp.site_id)]
    allp["evidence_native"] = labels.evidence_class(allp.in_vivo_native, allp.in_vitro_native)
    allp["evidence_hmr"] = labels.evidence_class(allp.in_vivo_hmr, allp.in_vitro_hmr)
    ok = allp.kinase_included & allp.residue_in_class
    isp = allp.provenance.eq("PSP")
    allp["L1"] = isp & ok & allp.native_human
    allp["L2"] = isp & ok & allp.hmr
    allp["L3"] = isp & ok
    allp["L2_vitro"] = allp.L2 & allp.in_vitro_hmr
    allp["L2_noauto"] = allp.L2 & ~allp.autophos
    allp["EXT"] = ext & ok
    allp["EXT_not_in_psp"] = allp.EXT & ~allp.in_psp_ks
    cols = ["kinase", "kinase_acc", "site_id", "acc", "pos", "res", "provenance", "L1", "L2", "L3", "L2_vitro",
            "L2_noauto", "EXT", "EXT_not_in_psp", "native_human", "hmr", "sub_native_human", "species_origin",
            "in_vivo", "in_vitro", "in_vivo_native", "in_vitro_native", "in_vivo_hmr", "in_vitro_hmr",
            "evidence_native", "evidence_hmr", "autophos", "psp_autophos_flag", "kinase_included",
            "kinase_residue_class", "residue_in_class", "in_universe", "in_psp_ks", "sub_map_methods", "kin_map_hows",
            "res_changed", "sgid_multi", "n_rows", "psp_mod_rsd", "psp_sub_accs", "ext_sources", "ext_refs",
            "omnipath_lists_psp", "omnipath_all_sources", "uniprot_note", "psp_sgids", "mask_only",
            "mask_sensitivity_only", "mask_source_name"]
    allp = allp[cols].sort_values(["provenance", "kinase", "acc", "pos"]).reset_index(drop=True)
    allp["n_rows"] = allp.n_rows.astype("Int64")
    allp["psp_autophos_flag"] = allp.psp_autophos_flag.astype("boolean")
    allp.to_parquet(out_c / "ks_pairs.parquet", index=False)
    pd.concat(unm, ignore_index=True).astype(str).to_parquet(out_c / "ks_unmapped.parquet", index=False)

    # ------------------------------------------------------------------ summaries
    psp = allp[isp]
    funnel += [("PSP", "pairs_all_mapped", len(psp)),
               ("PSP", "pairs_kinase_excluded", int((~psp.kinase_included).sum())),
               ("PSP", "pairs_residue_out_of_class(kinase included)",
                int((psp.kinase_included & ~psp.residue_in_class).sum()))]
    oc = psp[psp.kinase_included & ~psp.residue_in_class]
    for c, n in oc.kinase_residue_class.value_counts().items():
        funnel.append(("PSP", f"pairs_residue_out_of_class:class={c}", int(n)))
    funnel.append(("PSP", "pairs_out_of_class_in_vivo_only_frac",
                   round(float((oc.in_vivo & ~oc.in_vitro).mean()), 4) if len(oc) else np.nan))
    for v in ["L1", "L2", "L3"]:
        funnel.append(("PSP", f"{v}_sites_not_in_U_obs", int(psp.loc[psp[v], "site_id"].pipe(
            lambda s: (~s.isin(in_uobs)).sum()))))
    for pv in ["OmniPath", "UniProt"]:
        x = allp[allp.provenance.eq(pv)]
        funnel += [(pv, "pairs_kinase_excluded", int((~x.kinase_included).sum())),
                   (pv, "pairs_residue_out_of_class", int((x.kinase_included & ~x.residue_in_class).sum())),
                   (pv, "EXT_pairs", int(x.EXT.sum())), (pv, "EXT_pairs_also_in_PSP", int((x.EXT & x.in_psp_ks).sum()))]
    pd.DataFrame(funnel, columns=["provenance", "metric", "value"]).to_csv(out_r / "label_funnel.tsv", sep="\t",
                                                                          index=False)
    V = []
    for v in ["L1", "L2", "L3", "L2_vitro", "L2_noauto"]:
        V.append(labels.variant_counts(allp, allp[v], v))
    for ev in ["in_vivo_only", "in_vitro_only", "both"]:
        V.append(labels.variant_counts(allp, allp.L1 & allp.evidence_native.eq(ev), f"L1[{ev}]"))
    V.append(labels.variant_counts(allp, allp.L2 & ~allp.native_human, "L2_minus_L1(rodent-only)"))
    for pv in ["OmniPath", "UniProt"]:
        V.append(labels.variant_counts(allp, allp.EXT & allp.provenance.eq(pv), f"EXT_{pv}"))
        V.append(labels.variant_counts(allp, allp.EXT_not_in_psp & allp.provenance.eq(pv), f"EXT_{pv}_not_in_PSP"))
    V.append(labels.variant_counts(allp, allp.EXT, "EXT_any"))
    V.append(labels.variant_counts(allp, allp.EXT_not_in_psp, "EXT_any_not_in_PSP"))
    V.append(labels.variant_counts(allp, allp.L2 | allp.EXT, "L2+EXT"))
    vs = pd.DataFrame(V)
    vs.to_csv(out_r / "label_summary.tsv", sep="\t", index=False, float_format="%.4f")
    pk = (allp.groupby("kinase").agg(**{v: (v, "sum") for v in ["L1", "L2", "L3", "L2_vitro", "L2_noauto", "EXT",
                                                                  "EXT_not_in_psp"]})
          .reset_index())
    pk = pk.merge(kl[["gene", "include", "residue_class", "manning_group", "manning_family", "is_pseudokinase"]],
                  left_on="kinase", right_on="gene", how="left").drop(columns="gene")
    pk.sort_values("L2", ascending=False).to_csv(out_r / "kinase_label_counts.tsv", sep="\t", index=False)
    print(pd.DataFrame(funnel, columns=["provenance", "metric", "value"]).to_string(index=False))
    print(vs.to_string(index=False))
    print(f"total {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
