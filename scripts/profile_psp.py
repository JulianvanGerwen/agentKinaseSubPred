#!/usr/bin/env python
"""Profile the PhosphoSitePlus Kinase_Substrate_Dataset -> results/data_profile/psp_*.tsv.

Only aggregate tables are written to results/ (PSP is CC BY-NC-SA; row-level human-mapped tables go
to cache/, which is gitignored).

If UniProt files from scripts/download_uniprot.py are present (external/uniprot), additionally:
  * non-human kinases mapped to human by symbol against UniProt human kinases (pkinfam / KW-0418),
  * non-human substrate sites without a human member of their SITE_GRP_ID in the KS dataset are mapped
    by ortholog sequence (same upper-case gene symbol, ungapped 15-mer, >= 11/15 identity, same residue),
  * kinase family (pkinfam group/family) and protein-kinase-ness flags.

Usage: PYTHONPATH=src python scripts/profile_psp.py [--out results/data_profile] [--uniprot external/uniprot]
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import io_psp  # noqa: E402

HMR = ["human", "mouse", "rat"]
THRESH = [1, 5, 10, 20, 50, 100]


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return np.nan
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def load_uniprot(udir):
    from ksp import uniprot as up
    udir = Path(udir)
    fa = up.read_fasta(udir / "proteome" / "UP000005640_all_isoforms.fasta.gz")
    ann = up.read_annotations(udir / "human_tsv" / "UP000005640_annotations.tsv.gz")
    pk = up.parse_pkinfam(udir / "pkinfam" / "pkinfam.txt")
    return fa, ann, pk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psp", default=str(io_psp.DEFAULT_PATH))
    ap.add_argument("--out", default="results/data_profile")
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--uniprot", default="external/uniprot")
    ap.add_argument("--min-identity", type=int, default=11)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    Path(a.cache).mkdir(exist_ok=True)
    S = {}  # key numbers

    def w(df, name, index=False):
        df.to_csv(out / f"psp_{name}.tsv", sep="\t", index=index)

    d = io_psp.load_psp(a.psp)
    S["release"] = io_psp.read_release(a.psp)
    S["n_rows"] = len(d)
    S["n_rows_window_ok"] = int(d.window_ok.sum())
    S["n_rows_centre_lowercase"] = int(d.center_lower.sum())

    # ---------------- 1. organisms ----------------
    ct = pd.crosstab(d.kin_organism, d.sub_organism, margins=True, margins_name="TOTAL")
    w(ct, "organism_crosstab", index=True)
    S["n_rows_human_human"] = int(((d.kin_organism == "human") & (d.sub_organism == "human")).sum())

    # ---------------- 2. residues ----------------
    res = d.groupby(["residue"]).size().rename("n_rows").reset_index()
    res["is_sty"] = res.residue.isin(list("STY"))
    w(res, "residue")
    nonsty = d[~d.is_sty][["kin_gene", "kin_organism", "sub_gene", "sub_organism", "mod_rsd", "window_raw",
                           "in_vivo", "in_vitro"]]
    w(nonsty, "nonsty_rows")  # small; documents the junk rows dropped
    S["n_rows_non_sty"] = int((~d.is_sty).sum())
    dd = d[d.is_sty].copy()

    # ---------------- 3. duplicates ----------------
    S["n_dup_rows_same_kinacc_subacc_rsd"] = int(dd.duplicated(["kin_acc_full", "sub_acc_full", "mod_rsd"]).sum())
    S["n_dup_rows_same_kingene_org_subacc_rsd"] = int(
        dd.duplicated(["kin_gene_uc", "kin_organism", "sub_acc_full", "mod_rsd"]).sum())

    # ---------------- 4. evidence ----------------
    rows = []
    for scope, sub in [("all_rows", dd), ("human_kin_human_sub", dd[(dd.kin_organism == "human") & (dd.sub_organism == "human")])]:
        for ev, n in sub.evidence.value_counts().items():
            rows.append((scope, "rows", ev, n))
        # unique (kinase gene, organisms, site) pairs: aggregate evidence over duplicate rows
        g = sub.groupby(["kin_gene_uc", "kin_organism", "sub_acc_full", "position"]).agg(v=("in_vivo", "any"), t=("in_vitro", "any"))
        ev = np.select([g.v & g.t, g.v, g.t], ["both", "in_vivo_only", "in_vitro_only"], "none")
        for e, n in pd.Series(ev).value_counts().items():
            rows.append((scope, "unique_pairs", e, n))
    w(pd.DataFrame(rows, columns=["scope", "unit", "evidence", "n"]), "evidence")

    # ---------------- UniProt ----------------
    have_up = (Path(a.uniprot) / "human_tsv" / "UP000005640_annotations.tsv.gz").exists()
    if have_up:
        fa, ann, pk = load_uniprot(a.uniprot)
        ann["gene_uc"] = ann["Gene Names (primary)"].str.upper()
        rev = ann[ann.Reviewed.eq("reviewed")]
        acc2gene = dict(zip(ann.Entry, ann.gene_uc))
        pk = pk.merge(ann[["Entry", "gene_uc"]], left_on="acc", right_on="Entry", how="left")
        pk_genes = set(pk.gene_uc.dropna()) - {""}
        kw = rev[rev["Keyword ID"].str.contains("KW-0418")]
        kw_pkdom = kw[kw["Domain [FT]"].str.contains("Protein kinase", case=False)]
        up_kinase_genes = pk_genes | set(kw_pkdom.gene_uc)
        S["uniprot_pkinfam_human_entries"] = int(len(pk))
        S["uniprot_kw0418_reviewed_human"] = int(len(kw))
        S["uniprot_kw0418_with_protein_kinase_domain"] = int(len(kw_pkdom))
        S["uniprot_human_protein_kinase_genes_union"] = len(up_kinase_genes)
    else:
        print("UniProt files not found; skipping UniProt-dependent parts", file=sys.stderr)
        acc2gene, pk, up_kinase_genes = {}, None, set()

    # ---------------- 5. kinase -> human ----------------
    psp_h_kin = io_psp.human_kinase_symbols(dd)
    S["n_human_kinase_symbols_in_psp"] = len(psp_h_kin)
    nh = dd[dd.kin_organism != "human"]
    kk = nh.groupby(["kin_organism", "kin_gene_uc"]).size().rename("n_rows").reset_index()
    kk["maps_psp_human_kinase"] = kk.kin_gene_uc.isin(psp_h_kin)
    if have_up:
        kk["maps_uniprot_human_kinase"] = kk.kin_gene_uc.isin(up_kinase_genes)
        kk["maps_uniprot_human_gene"] = kk.kin_gene_uc.isin(set(ann.gene_uc))
    w(kk.sort_values(["kin_organism", "n_rows"], ascending=[True, False]), "kinase_nonhuman_symbol_mapping")
    S["nonhuman_kinase_rows"] = int(len(nh))
    S["nonhuman_kinase_rows_mapped_psp_human_symbol"] = int(nh.kin_gene_uc.isin(psp_h_kin).sum())
    S["nonhuman_kinase_genes"] = int(len(kk))
    S["nonhuman_kinase_genes_mapped_psp_human_symbol"] = int(kk.maps_psp_human_kinase.sum())
    human_kin_set = psp_h_kin | up_kinase_genes
    if have_up:
        S["nonhuman_kinase_rows_mapped_psp_or_uniprot"] = int(nh.kin_gene_uc.isin(human_kin_set).sum())
        S["nonhuman_kinase_genes_mapped_psp_or_uniprot"] = int(kk.kin_gene_uc.isin(human_kin_set).sum())
        S["nonhuman_kinase_symbols_only_in_uniprot"] = sorted(set(kk.kin_gene_uc[kk.kin_gene_uc.isin(up_kinase_genes) & ~kk.maps_psp_human_kinase]))
    S["nonhuman_kinase_symbols_unmapped"] = sorted(set(kk.kin_gene_uc[~kk.kin_gene_uc.isin(human_kin_set)]))
    dk = io_psp.map_kinases_to_human(dd, human_symbols=human_kin_set)
    nhk = dk[dk.kin_organism != "human"]
    S["nonhuman_kinase_rows_mapped_final"] = int(nhk.kin_human_mapped.sum())
    S["nonhuman_kinase_rows_mapped_via_kinase_name"] = int((nhk.kin_human_mapped & nhk.kin_gene_from_name).sum())

    # ---------------- 6. substrate -> human ----------------
    m = io_psp.map_substrates_to_human(dk)
    m["map_method"] = np.where(m.sub_organism.eq("human"), "self", np.where(m.sub_human_mapped, "site_grp_id", "none"))
    if have_up:
        canon = fa[(fa.isoform == "") & fa.db.eq("sp")].copy()
        canon["gene_uc"] = canon.acc.map(acc2gene).fillna(canon.gene.str.upper())
        un = m[~m.sub_human_mapped][["sub_organism", "sub_acc_full", "position", "residue", "window", "sub_gene_uc"]].drop_duplicates(["sub_acc_full", "position"])
        mapped = io_psp.map_sites_by_sequence(un, canon[["gene_uc", "acc_full", "acc", "seq"]], a.min_identity)
        mapped = mapped.dropna(subset=["seq_h_sub_acc_full"])
        key = ["sub_acc_full", "position"]
        m = m.merge(mapped[key + ["seq_h_sub_acc_full", "seq_h_sub_acc", "seq_h_position", "seq_identity"]], on=key, how="left")
        fill = ~m.sub_human_mapped & m.seq_h_sub_acc_full.notna()
        m.loc[fill, "h_sub_acc_full"] = m.loc[fill, "seq_h_sub_acc_full"]
        m.loc[fill, "h_sub_acc"] = m.loc[fill, "seq_h_sub_acc"]
        m.loc[fill, "h_position"] = m.loc[fill, "seq_h_position"]
        m.loc[fill, "h_residue"] = m.loc[fill, "residue"]
        m.loc[fill, "map_method"] = "sequence"
        m["sub_human_mapped"] = m.h_sub_acc_full.notna()
        m = m.drop(columns=["seq_h_sub_acc_full", "seq_h_sub_acc", "seq_h_position"])
    m["h_position"] = m.h_position.astype("Int64")
    # residue consistency of SITE_GRP_ID mapping
    sg = m[m.map_method.eq("site_grp_id")]
    S["sgid_mapped_rows_residue_changed"] = int((sg.residue != sg.h_residue).sum())
    # site-level mapping table (unique non-human sites)
    sites = m[m.sub_organism != "human"].sort_values("map_method").groupby(["sub_organism", "sub_acc_full", "position"]).agg(
        method=("map_method", lambda x: "site_grp_id" if (x == "site_grp_id").any() else ("sequence" if (x == "sequence").any() else "none")),
        n_human_targets=("h_sub_acc_full", "nunique")).reset_index()
    tab = sites.groupby(["sub_organism", "method"]).size().unstack(fill_value=0)
    tab["total_sites"] = tab.sum(1)
    tab = tab.sort_values("total_sites", ascending=False)
    w(tab, "substrate_site_human_mapping", index=True)
    S["nonhuman_unique_sites"] = int(len(sites))
    S["nonhuman_unique_sites_mapped_sgid"] = int((sites.method == "site_grp_id").sum())
    S["nonhuman_unique_sites_mapped_sequence_fallback"] = int((sites.method == "sequence").sum())
    S["nonhuman_sites_sgid_mapping_to_gt1_human_site"] = int((sites.n_human_targets > 1).sum())
    # human gene of mapped substrate
    m["h_sub_gene"] = m.h_sub_acc.map(acc2gene) if have_up else np.nan
    m["h_sub_gene"] = m.h_sub_gene.where(m.h_sub_gene.notna() & (m.h_sub_gene != ""), np.where(m.sub_organism.eq("human"), m.sub_gene_uc, None))
    m.to_parquet(Path(a.cache) / "psp_human_mapped_rows.parquet", index=False)

    # ---------------- 7. pooling gain ----------------
    ok = m[m.kin_human_mapped & m.sub_human_mapped].copy()
    scen = {
        "H_kin__H_sub": (ok.kin_organism.eq("human")) & ok.sub_organism.eq("human"),
        "any_kin__H_sub": ok.sub_organism.eq("human"),
        "HMR_kin__HMR_sub": ok.kin_organism.isin(HMR) & ok.sub_organism.isin(HMR),
        "HMR_kin__HMR_sub_sgid_only": ok.kin_organism.isin(HMR) & ok.sub_organism.isin(HMR) & ok.map_method.isin(["self", "site_grp_id"]),
        "all_kin__all_sub": pd.Series(True, index=ok.index),
    }
    evs = {"any": None, "in_vivo": "in_vivo", "in_vitro": "in_vitro", "in_vivo_only_rows": "vivo_only", "in_vitro_only_rows": "vitro_only"}
    rows = []
    for sn, sm in scen.items():
        for en, ef in evs.items():
            x = ok[sm]
            if ef == "in_vivo":
                x = x[x.in_vivo]
            elif ef == "in_vitro":
                x = x[x.in_vitro]
            elif ef == "vivo_only":
                x = x[x.in_vivo & ~x.in_vitro]
            elif ef == "vitro_only":
                x = x[x.in_vitro & ~x.in_vivo]
            pairs = x.drop_duplicates(["kin_human_gene", "h_sub_acc_full", "h_position"])
            kc = pairs.kin_human_gene.value_counts()
            rows.append(dict(scenario=sn, evidence=en, n_rows=len(x),
                             n_sites=x[["h_sub_acc_full", "h_position"]].drop_duplicates().shape[0],
                             n_pairs=len(pairs), n_kinases=int(len(kc)), n_kinases_ge5=int((kc >= 5).sum()),
                             n_kinases_ge20=int((kc >= 20).sum()), n_substrate_proteins=x.h_sub_acc.nunique()))
    gain = pd.DataFrame(rows)
    w(gain, "pooling_gain")
    # evidence at pair level in pooled HMR set: is a pair supported in vivo in any species?
    hmr = ok[scen["HMR_kin__HMR_sub"]]
    pg = hmr.groupby(["kin_human_gene", "h_sub_acc_full", "h_position"]).agg(
        v=("in_vivo", "any"), t=("in_vitro", "any"), n_species=("sub_organism", "nunique"),
        human_row=("sub_organism", lambda s: (s == "human").any()))
    pg["evidence"] = np.select([pg.v & pg.t, pg.v, pg.t], ["both", "in_vivo_only", "in_vitro_only"], "none")
    w(pg.groupby(["evidence", "human_row"]).size().rename("n_pairs").reset_index(), "hmr_pair_evidence")
    S["hmr_pairs_supported_by_gt1_species"] = int((pg.n_species > 1).sum())

    # ---------------- 8. isoforms ----------------
    iso = dd.assign(is_iso=dd.sub_isoform.ne("")).groupby("sub_organism").agg(
        n_rows=("is_iso", "size"), n_rows_isoform_acc=("is_iso", "sum"),
        n_unique_isoform_accs=("sub_acc_full", lambda s: s[s.str.contains("-")].nunique()))
    iso = iso.sort_values("n_rows", ascending=False)
    w(iso, "isoforms_by_organism", index=True)
    w(dd[dd.sub_isoform.ne("")].sub_isoform.value_counts().rename_axis("isoform_suffix").rename("n_rows").reset_index(), "isoform_suffixes")
    S["sub_rows_isoform_acc"] = int(dd.sub_isoform.ne("").sum())
    S["sub_unique_isoform_accs"] = int(dd.sub_acc_full[dd.sub_isoform.ne("")].nunique())
    S["kin_rows_isoform_acc"] = int(dd.kin_isoform.ne("").sum())
    S["kin_isoform_accs"] = dd[dd.kin_isoform.ne("")].kin_acc_full.value_counts().to_dict()

    # ---------------- 9. kinase distribution ----------------
    def per_kin(x):
        p = x.drop_duplicates(["kin_human_gene", "h_sub_acc_full", "h_position"])
        return p.kin_human_gene.value_counts()

    hh = ok[scen["H_kin__H_sub"]]
    kin = pd.DataFrame({
        "n_pairs_human_only": per_kin(hh),
        "n_pairs_hmr": per_kin(hmr),
        "n_pairs_all_species": per_kin(ok),
        "n_pairs_hmr_in_vivo": per_kin(hmr[hmr.in_vivo]),
        "n_pairs_hmr_in_vitro": per_kin(hmr[hmr.in_vitro]),
        "n_substrate_proteins_hmr": hmr.groupby("kin_human_gene").h_sub_acc.nunique(),
    }).fillna(0).astype(int).sort_values("n_pairs_hmr", ascending=False)
    kin.index.name = "kin_human_gene"
    if have_up:
        g2pk = pk.dropna(subset=["gene_uc"]).drop_duplicates("gene_uc").set_index("gene_uc")
        kin["pkinfam_group"] = kin.index.map(g2pk["group"])
        kin["pkinfam_family"] = kin.index.map(g2pk["family"])
        kin["in_pkinfam"] = kin.index.isin(pk_genes)
        kin["uniprot_kw0418_pkdomain"] = kin.index.isin(set(kw_pkdom.gene_uc))
        kin["uniprot_kw0418"] = kin.index.isin(set(kw.gene_uc))
        kin["uniprot_human_gene"] = kin.index.isin(set(rev.gene_uc))
        kin["protein_kinase"] = kin.in_pkinfam | kin.uniprot_kw0418_pkdomain
    w(kin, "kinase_site_counts", index=True)
    thr = []
    for col in ["n_pairs_human_only", "n_pairs_hmr", "n_pairs_all_species", "n_pairs_hmr_in_vivo"]:
        v = kin[col]
        v = v[v > 0].sort_values(ascending=False)
        r = {"count_basis": col, "n_kinases": int(len(v))}
        for t in THRESH:
            r[f"ge{t}"] = int((v >= t).sum())
        r["total_pairs"] = int(v.sum())
        r["gini"] = round(gini(v.values), 3)
        r["frac_pairs_top10"] = round(v.iloc[:10].sum() / v.sum(), 3)
        r["frac_pairs_top20"] = round(v.iloc[:20].sum() / v.sum(), 3)
        r["frac_pairs_top50"] = round(v.iloc[:50].sum() / v.sum(), 3)
        r["top10"] = ",".join(v.index[:10])
        thr.append(r)
    w(pd.DataFrame(thr), "kinase_thresholds")

    # ---------------- 10. family distribution ----------------
    if have_up:
        fam = kin.assign(pkinfam_group=kin.pkinfam_group.fillna("(not in pkinfam)")).groupby("pkinfam_group").agg(
            n_kinases=("n_pairs_hmr", lambda s: int((s > 0).sum())), n_pairs_hmr=("n_pairs_hmr", "sum"))
        tot_pk = pk.drop_duplicates("acc").groupby("group").size().rename("n_human_kinases_in_pkinfam")
        fam = fam.join(tot_pk, how="outer").fillna(0).astype(int).sort_values("n_pairs_hmr", ascending=False)
        w(fam, "kinase_group_counts", index=True)
        S["uniprot_human_protein_kinases_with_ge1_hmr_pair"] = int(kin[kin.protein_kinase & (kin.n_pairs_hmr > 0)].shape[0])
        nonpk = kin[~kin.protein_kinase & (kin.n_pairs_all_species > 0)]
        w(nonpk.reset_index(), "non_protein_kinase_entries")
        S["n_kinase_genes_not_protein_kinase"] = int(len(nonpk))
        S["pairs_hmr_from_non_protein_kinases"] = int(nonpk.n_pairs_hmr.sum())

    # ---------------- 11. multi-site / multi-label ----------------
    pairs = hmr.drop_duplicates(["kin_human_gene", "h_sub_acc_full", "h_position"])
    site_k = pairs.groupby(["h_sub_acc_full", "h_position"]).kin_human_gene.nunique()
    w(site_k.value_counts().sort_index().rename_axis("n_kinases_per_site").rename("n_sites").reset_index(), "kinases_per_site")
    prot_s = pairs.drop_duplicates(["h_sub_acc_full", "h_position"]).groupby("h_sub_acc").size()
    w(prot_s.value_counts().sort_index().rename_axis("n_sites_per_protein").rename("n_proteins").reset_index(), "sites_per_protein")
    S["hmr_sites_with_gt1_kinase"] = int((site_k > 1).sum())
    S["hmr_sites_total"] = int(len(site_k))
    S["hmr_proteins_total"] = int(len(prot_s))
    S["hmr_proteins_with_gt1_site"] = int((prot_s > 1).sum())
    S["hmr_proteins_with_ge10_sites"] = int((prot_s >= 10).sum())
    S["hmr_sites_on_top_1pct_proteins_frac"] = round(float(prot_s.sort_values(ascending=False).iloc[:max(1, len(prot_s) // 100)].sum() / prot_s.sum()), 3)

    # ---------------- 12. autophosphorylation ----------------
    S["rows_autophos_same_org"] = int(dd.autophos.sum())
    auto_h = pairs.kin_human_gene == pairs.h_sub_gene
    S["hmr_pairs_autophos"] = int(auto_h.sum())
    S["hmr_pairs_autophos_frac"] = round(float(auto_h.mean()), 4)
    S["hmr_kinases_with_autophos_pair"] = int(pairs[auto_h].kin_human_gene.nunique())
    ap_tab = pairs[auto_h].kin_human_gene.value_counts().rename_axis("kin_human_gene").rename("n_autophos_pairs").reset_index()
    w(ap_tab, "autophos_by_kinase")

    # ---------------- 13. window: other phospho residues ----------------
    w(dd.n_other_phos_in_window.value_counts().sort_index().rename_axis("n_other_lowercase_in_window").rename("n_rows").reset_index(), "window_other_phos")

    json.dump(S, open(out / "psp_summary.json", "w"), indent=1, default=str)
    pd.Series({k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in S.items()}).rename_axis("metric").rename("value").reset_index().to_csv(out / "psp_summary.tsv", sep="\t", index=False)
    print(json.dumps(S, indent=1, default=str))
    print(gain.to_string())
    print(pd.DataFrame(thr).drop(columns="top10").to_string())


if __name__ == "__main__":
    main()
