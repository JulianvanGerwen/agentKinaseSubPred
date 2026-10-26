"""Kinase–substrate label construction (PSP Kinase_Substrate_Dataset + external curated pairs).

Output schema (cache/labels/ks_pairs.parquet): one row per (kinase, site_id, provenance), provenance in
{'PSP', 'OmniPath', 'UniProt'}. Rules are documented in docs/labels_universe.md; summary:

PSP rows (S/T/Y only)
  kinase   : human rows -> kinase list by UniProt accession, else symbol/synonym; non-human rows -> human kinase by
             (upper-case) symbol/synonym (ASSUMPTION: same symbol => ortholog); empty GENE -> PSP kinase name.
  substrate: human rows -> canonical site by 15-mer (ids.CanonicalStore.map_window: exact at stated position, else
             unique elsewhere in the same gene's canonical); non-human rows -> SITE_GRP_ID -> human sites of that
             group (PSP Phosphorylation_site_dataset human rows + KS human rows, each mapped to canonical by 15-mer);
             else ortholog sequence fallback (io_psp.map_sites_by_sequence: same symbol, ungapped 15-mer,
             >= min_identity/15 identical, same residue) on the Swiss-Prot canonical.
  flags    : native_human (>=1 row with kinase AND substrate organism human), sub_native_human, species_origin,
             in_vivo / in_vitro (any row) and their _native / _hmr versions, autophos (kinase accession == substrate
             accession), kinase_included, residue_in_class, in_universe.
  variants : L1 = native rows; L2 = kinase & substrate organism in {human, mouse, rat}; L3 = any organism;
             L2_vitro = L2 restricted to HMR rows with in vitro evidence; L2_noauto = L2 without autophos.
             Every variant additionally requires kinase_included & residue_in_class (& a mapped site).
External (never PSP variants): OmniPath enzsub rows with >= 1 curated non-PSP source; UniProt MOD_RES
  'Phospho...; by <KINASE>' notes with ECO:0000269 evidence. Variant column EXT = kinase_included & residue_in_class.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import io_psp
from .ids import CanonicalStore, map_table, site_ids
from .kinases import KinaseResolver

HMR = ("human", "mouse", "rat")
OMNIPATH_CURATED = ("SIGNOR", "HPRD", "phosphoELM", "NCI-PID_ProtMapper", "Reactome_ProtMapper", "SIGNOR_ProtMapper")
OMNIPATH_OTHER_AGGREGATORS = ("dbPTM", "Li2012")  # neither PSP/text-mining nor primary curation; off by default
PSP_TOKENS = ("PhosphoSite", "PhosphoSite_MIMP", "PhosphoSite_ProtMapper")
_MOD = re.compile(r'MOD_RES (\d+); /note="(Phospho(serine|threonine|tyrosine)[^"]*)"(?:; /evidence="([^"]*)")?')
_RES = {"serine": "S", "threonine": "T", "tyrosine": "Y"}


def residue_in_class(res: pd.Series, cls: pd.Series) -> pd.Series:
    return pd.Series([(isinstance(c, str) and r in c) if c != "ST" else r in ("S", "T")
                      for r, c in zip(res, cls)], index=res.index)


# ---------------------------------------------------------------------------------------------- PSP
def resolve_psp_kinases(d: pd.DataFrame, res: KinaseResolver) -> pd.DataFrame:
    """Add kinase (kinase-list gene or None) and kin_map_how."""
    hum_name2gene = (d[d.kin_organism.eq("human") & d.kin_gene_uc.ne("")]
                     .groupby(d.kin_name.str.upper()).kin_gene_uc.unique())
    name2gene = {k: v[0] for k, v in hum_name2gene.items() if len(v) == 1}
    out_g, out_h = [], []
    for org, g, acc, name in zip(d.kin_organism, d.kin_gene, d.kin_acc, d.kin_name):
        if org == "human":
            k, how = res.resolve(g, acc=acc)
        else:
            k, how = res.resolve(g)
            if k is None and not g:
                gg = name2gene.get(str(name).upper())
                if gg:
                    k, how = res.resolve(gg)
                    how = f"{how}(via_psp_name)"
            how = f"ortholog_{how}" if k is not None else how
        out_g.append(k)
        out_h.append(how)
    d = d.copy()
    d["kinase"], d["kin_map_how"] = out_g, out_h
    return d


def sgid_to_human_sites(store: CanonicalStore, psp_sites_human: pd.DataFrame, ks: pd.DataFrame) -> dict:
    """SITE_GRP_ID -> set of canonical human site_ids (PSP site dataset human rows + KS human rows, mapped by 15-mer)."""
    a = psp_sites_human[["acc_full", "position", "residue", "window", "gene", "site_grp_id"]].rename(
        columns={"acc_full": "acc", "position": "pos", "residue": "res"})
    h = ks[ks.sub_organism.eq("human") & ks.is_sty][["sub_acc_full", "position", "residue", "window", "sub_gene",
                                                      "site_grp_id"]].rename(
        columns={"sub_acc_full": "acc", "position": "pos", "residue": "res", "sub_gene": "gene"})
    x = pd.concat([a, h], ignore_index=True).dropna(subset=["site_grp_id"]).drop_duplicates()
    m = map_table(store, x, "acc", "pos", "window", "res", "gene")
    x = pd.concat([x.reset_index(drop=True), m.reset_index(drop=True)], axis=1).dropna(subset=["site_id"])
    return x.groupby("site_grp_id").site_id.agg(lambda s: sorted(set(s))).to_dict()


def map_psp_substrates(d: pd.DataFrame, store: CanonicalStore, sgid_map: dict, min_identity: int = 11) -> pd.DataFrame:
    """Return rows exploded over mapped canonical sites: original columns + site_id, sub_map_method, sub_map_reason,
    sgid_n_human_sites, seq_identity."""
    d = d.reset_index(drop=True).copy()
    d["_row"] = np.arange(len(d))
    is_h = d.sub_organism.eq("human")
    # human rows: window mapping
    hm = map_table(store, d[is_h].rename(columns={"sub_acc_full": "acc", "position": "pos", "window": "window"}),
                   "acc", "pos", "window", "residue", "sub_gene")
    hum = d[is_h].copy()
    hum["site_id"] = hm.site_id.values
    hum["sub_map_method"] = np.where(hm.site_id.notna(), "human_window_" + hm.map_method.values, "unmapped")
    hum["sub_map_reason"] = np.where(hm.site_id.notna(), "", "human_window_" + hm.map_reason.values)
    # non-human rows: SGID, then ortholog sequence fallback
    non = d[~is_h].copy()
    sg = non.site_grp_id.map(lambda g: sgid_map.get(g, []) if pd.notna(g) else [])
    non["_sg"] = sg.values
    has = non._sg.map(len) > 0
    a = non[has].explode("_sg").rename(columns={"_sg": "site_id"})
    a["sub_map_method"] = "sgid"
    a["sub_map_reason"] = ""
    a["sgid_n_human_sites"] = a.groupby("_row").site_id.transform("nunique")
    b = non[~has].drop(columns="_sg")
    if len(b):
        prot = [(g, acc, acc, store.seqs[acc]) for g, accs in store.gene_to_accs.items() for acc in accs]
        hp = pd.DataFrame(prot, columns=["gene_uc", "acc_full", "acc", "seq"])
        need = b[b.sub_gene_uc.isin(set(hp.gene_uc))]
        uniq = need.drop_duplicates(["sub_gene_uc", "window", "residue", "position"])[
            ["sub_gene_uc", "window", "residue", "position"]]
        hits = io_psp.map_sites_by_sequence(uniq, hp[hp.gene_uc.isin(set(uniq.sub_gene_uc))], min_identity)
        b = b.merge(hits, on=["sub_gene_uc", "window", "residue", "position"], how="left")
        ok = b.seq_h_sub_acc.notna()
        b["site_id"] = None
        b.loc[ok, "site_id"] = site_ids(b.loc[ok, "seq_h_sub_acc"], b.loc[ok, "residue"], b.loc[ok, "seq_h_position"])
        b["sub_map_method"] = np.where(ok, "ortholog_seq", "unmapped")
        b["sub_map_reason"] = np.where(ok, "", np.where(b.site_grp_id.notna(), "sgid_no_human_site;no_ortholog_hit",
                                                        "no_sgid;no_ortholog_hit"))
        b = b.drop(columns=["seq_h_sub_acc_full", "seq_h_sub_acc", "seq_h_position"])
    out = pd.concat([hum, a, b], ignore_index=True).sort_values("_row", kind="stable").reset_index(drop=True)
    out["sgid_n_human_sites"] = out.get("sgid_n_human_sites", pd.Series(np.nan, index=out.index))
    return out


def psp_pairs(d: pd.DataFrame, kl: pd.DataFrame, universe_ids: set | None = None) -> pd.DataFrame:
    """Aggregate mapped PSP rows (columns kinase, site_id, ...) to one row per (kinase, site_id)."""
    k = d.dropna(subset=["kinase", "site_id"]).copy()
    k["site_res"] = k.site_id.str.extract(r"_([STY])\d+$")[0]
    k["site_acc"] = k.site_id.str.split("_").str[0]
    acc = dict(zip(kl.gene, kl.uniprot_acc))
    k["kin_acc_h"] = k.kinase.map(acc)
    k["native"] = k.kin_organism.eq("human") & k.sub_organism.eq("human")
    k["hmr"] = k.kin_organism.isin(HMR) & k.sub_organism.isin(HMR)
    k["sub_native"] = k.sub_organism.eq("human")
    k["autophos_row"] = k.kin_acc_h.eq(k.site_acc)
    k["res_changed"] = k.residue.ne(k.site_res)
    k["origin"] = k.kin_organism.str.replace(" ", "_") + ">" + k.sub_organism.str.replace(" ", "_")
    for ev in ["in_vivo", "in_vitro"]:
        k[f"{ev}_native"] = k[ev] & k.native
        k[f"{ev}_hmr"] = k[ev] & k.hmr
    g = k.groupby(["kinase", "site_id"], sort=False)
    p = g.agg(native_human=("native", "any"), hmr=("hmr", "any"), sub_native_human=("sub_native", "any"),
              species_origin=("origin", lambda s: ";".join(sorted(set(s)))),
              in_vivo=("in_vivo", "any"), in_vitro=("in_vitro", "any"),
              in_vivo_native=("in_vivo_native", "any"), in_vitro_native=("in_vitro_native", "any"),
              in_vivo_hmr=("in_vivo_hmr", "any"), in_vitro_hmr=("in_vitro_hmr", "any"),
              autophos=("autophos_row", "any"), psp_autophos_flag=("autophos", "any"),
              sub_map_methods=("sub_map_method", lambda s: ";".join(sorted(set(s)))),
              kin_map_hows=("kin_map_how", lambda s: ";".join(sorted(set(s)))),
              res_changed=("res_changed", "any"), sgid_multi=("sgid_n_human_sites", lambda s: bool((s > 1).any())),
              n_rows=("_row", "nunique"), psp_mod_rsd=("mod_rsd", lambda s: ";".join(sorted(set(s)))),
              psp_sub_accs=("sub_acc_full", lambda s: ";".join(sorted(set(s))))).reset_index()
    p["provenance"] = "PSP"
    return p


def evidence_class(vivo: pd.Series, vitro: pd.Series) -> pd.Series:
    return pd.Series(np.select([vivo & vitro, vivo, vitro], ["both", "in_vivo_only", "in_vitro_only"], "none"),
                     index=vivo.index)


# ---------------------------------------------------------------------------------------------- OmniPath
def omnipath_pairs(path, store: CanonicalStore, res: KinaseResolver, curated=OMNIPATH_CURATED,
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Curated non-PSP OmniPath enzyme-substrate phosphorylation rows -> (pairs, unmapped log)."""
    o = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    o = o[o.modification.eq("phosphorylation")].copy()
    toks = o.sources.str.split(";")
    o["curated_sources"] = toks.map(lambda t: ";".join(x for x in t if x in curated))
    o["omnipath_lists_psp"] = toks.map(lambda t: any(x in PSP_TOKENS for x in t))
    log = []
    o_drop = o[o.curated_sources.eq("")]
    log.append(pd.DataFrame({"provenance": "OmniPath", "reason": "no_curated_non_psp_source",
                             "kinase_name": o_drop.enzyme_genesymbol, "substrate": o_drop.substrate,
                             "site": o_drop.residue_type + o_drop.residue_offset, "detail": o_drop.sources}))
    o = o[o.curated_sources.ne("")]
    cur_refs = []
    for refs, cs in zip(o.references, o.curated_sources):
        keep = set(c.replace("_ProtMapper", "") for c in cs.split(";")) | set(cs.split(";"))
        cur_refs.append(";".join(sorted({r.split(":", 1)[1] for r in refs.split(";") if ":" in r and
                                         r.split(":", 1)[0] in keep})))
    o["curated_refs"] = cur_refs
    kin, how, sid, reason = [], [], [], []
    for e, eg, s, rt, ro in zip(o.enzyme, o.enzyme_genesymbol, o.substrate, o.residue_type, o.residue_offset):
        k, h = res.resolve(eg, acc=e)
        m = store.map_position(s, ro, rt)
        kin.append(k), how.append(h), sid.append(m.site_id), reason.append(m.reason)
    o["kinase"], o["kin_map_how"], o["site_id"], o["sub_map_reason"] = kin, how, sid, reason
    bad = o.kinase.isna() | o.site_id.isna()
    b = o[bad]
    log.append(pd.DataFrame({"provenance": "OmniPath",
                             "reason": np.where(b.kinase.isna(), "kinase_" + b.kin_map_how,
                                                "site_" + b.sub_map_reason),
                             "kinase_name": b.enzyme_genesymbol, "substrate": b.substrate,
                             "site": b.residue_type + b.residue_offset, "detail": b.sources}))
    o = o[~bad]
    p = o.groupby(["kinase", "site_id"], as_index=False).agg(
        ext_sources=("curated_sources", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),
        ext_refs=("curated_refs", lambda s: ";".join(sorted(set(x for v in s for x in v.split(";") if x)))),
        omnipath_all_sources=("sources", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),
        omnipath_lists_psp=("omnipath_lists_psp", "any"), kin_map_hows=("kin_map_how", "first"))
    p["provenance"] = "OmniPath"
    return p, pd.concat(log, ignore_index=True)


# ---------------------------------------------------------------------------------------------- UniProt
def parse_uniprot_mod_res(path) -> pd.DataFrame:
    """Rows: acc, pos, res, note, evidence, kinase_name (one per 'by' name; 'autocatalysis' kept), eco269, pmids."""
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    recs = []
    for acc, gene, f in zip(t.Entry, t["Gene Names (primary)"], t["Modified residue"]):
        for pos, note, rname, ev in _MOD.findall(f):
            if "; by " not in note:
                continue
            by = note.split("; by ", 1)[1].split(";")[0].strip()
            pmids = ";".join(sorted(set(re.findall(r"ECO:0000269\|PubMed:(\d+)", ev or ""))))
            for name in [x.strip() for x in re.split(r",\s*|\s+and\s+|\s+or\s+", by) if x.strip()]:
                recs.append((acc, gene, int(pos), _RES[rname], note, ev, name, "ECO:0000269" in (ev or ""), pmids))
    return pd.DataFrame(recs, columns=["acc", "gene", "pos", "res", "note", "evidence", "kinase_name", "eco269",
                                       "pmids"])


def uniprot_pairs(path, store: CanonicalStore, res: KinaseResolver) -> tuple[pd.DataFrame, pd.DataFrame]:
    u = parse_uniprot_mod_res(path)
    log = [pd.DataFrame({"provenance": "UniProt", "reason": "no_ECO:0000269", "kinase_name": u.kinase_name[~u.eco269],
                         "substrate": u.acc[~u.eco269], "site": u.res[~u.eco269] + u.pos[~u.eco269].astype(str),
                         "detail": u.note[~u.eco269]})]
    u = u[u.eco269].copy()
    kin, how = [], []
    for name, acc in zip(u.kinase_name, u.acc):
        if name.lower() == "autocatalysis":
            k, h = res.resolve(acc=acc)
            h = "autocatalysis" if k is not None else "autocatalysis_not_kinase"
        else:
            k, h = res.resolve_note(name)
        kin.append(k), how.append(h)
    u["kinase"], u["kin_map_how"] = kin, how
    sid, reason = [], []
    for a, p, r in zip(u.acc, u.pos, u.res):
        m = store.map_position(a, p, r)
        sid.append(m.site_id), reason.append(m.reason)
    u["site_id"], u["sub_map_reason"] = sid, reason
    bad = u.kinase.isna() | u.site_id.isna()
    b = u[bad]
    log.append(pd.DataFrame({"provenance": "UniProt",
                             "reason": np.where(b.kinase.isna(), "kinase_" + b.kin_map_how, "site_" + b.sub_map_reason),
                             "kinase_name": b.kinase_name, "substrate": b.acc, "site": b.res + b.pos.astype(str),
                             "detail": b.note}))
    u = u[~bad]
    p = u.groupby(["kinase", "site_id"], as_index=False).agg(
        ext_sources=("kinase_name", lambda s: "UniProt:" + ";".join(sorted(set(s)))),
        ext_refs=("pmids", lambda s: ";".join(sorted(set(x for v in s for x in v.split(";") if x)))),
        uniprot_note=("note", "first"), kin_map_hows=("kin_map_how", lambda s: ";".join(sorted(set(s)))))
    p["provenance"] = "UniProt"
    return p, pd.concat(log, ignore_index=True)


# ---------------------------------------------------------------------------------------------- summaries
def variant_counts(df: pd.DataFrame, mask: pd.Series, name: str) -> dict:
    x = df[mask]
    per_k = x.groupby("kinase").site_id.nunique()
    return {"variant": name, "pairs": len(x.drop_duplicates(["kinase", "site_id"])), "sites": x.site_id.nunique(),
            "substrate_proteins": x.site_id.str.split("_").str[0].nunique(), "kinases": int((per_k > 0).sum()),
            "kinases_ge5": int((per_k >= 5).sum()), "kinases_ge20": int((per_k >= 20).sum()),
            "kinases_ge50": int((per_k >= 50).sum()),
            "top10_share": float(per_k.sort_values(ascending=False).head(10).sum() / max(per_k.sum(), 1))}
