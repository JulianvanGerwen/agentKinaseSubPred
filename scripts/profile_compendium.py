#!/usr/bin/env python
"""Profile the perturbation phosphoproteomics compendium (data/01-03), the Hijazi kinase-inhibitor matrix,
the overlap with PSP human sites, and the human proteome S/T/Y universe (UniProt).

Inputs : data/01_phosphoproteomics_conditions.csv, data/02_phosphosite_information.csv,
         cache/compendium_wide.parquet (from scripts/build_compendium_matrix.py),
         data/kinhib_inhib_comb__1.csv, cache/psp_human_mapped_rows.parquet (from scripts/profile_psp.py),
         external/uniprot/ (from scripts/download_uniprot.py; optional)
Outputs: results/data_profile/compendium_*.tsv, kinhib_*.tsv, proteome_*.tsv, compendium_summary.json

Perturbation types are assigned by a transparent keyword heuristic on the free-text Condition column
(see classify_condition); it is approximate and its 'other/unclear' share is reported.

Usage: PYTHONPATH=src python scripts/profile_compendium.py
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import io_psp  # noqa: E402

RULES = [  # (type, regex) - first match wins; case-insensitive
    ("genetic (KO/KD/si/sh/mutant/overexpression)",
     r"(^|[_\s\-(])(si|sh)[A-Z0-9\-]|\bKO\b|_KO|KO_|\bKD\b|knock|\bCRISPR|mutant|\(K\d+M\)|\bWT\b.*\b(vs|parkin)|"
     r"parkin|overexpress|\bvector\b|\bAID\b|-AID|degron|Doxycycline|\bS\d+A\b|\bC\d+A\b|HA-|FGFR4|Plk1\(|analog-sensitive|"
     r"3-MB-PP1|1NM-PP1|\bdel\b|deleted|siCtrl|siRNA|_kd\b"),
    ("infection", r"infect|virus|SARS|HIV|ZIKV|flexneri|aureus|trachomatis|salmonella|\bMOI\b|Vpr|influenza|bacteria"),
    ("physical/genotoxic stress",
     r"xray|x-ray|proton|\bIR\b|irradiat|\bUV|cold|heat|hypox|osmo|\bNCS\b|^HU_|HU_\d|hydroxyurea|bleomycin|"
     r"doxorubicin|carbo|cisplatin|etoposide|camptothecin|exercise|\bH2O2|oxidat|starv|SFM|glucose|3D\b|OIS|dispase|"
     r"temozolomide|TMZ|Cytarabine|Methotrexat|stretch|shear"),
    ("ligand/cytokine/growth-factor stimulation",
     r"EGF|FGF|TNF|\bIL-?\d|IFN|insulin|\bIGF|PDGF|HGF|NGF|BDNF|VEGF|TGF|PMA|TPA|serum|CXCL|CCL|PGE2|LPS|anti-CD3|"
     r"CD28|EDN|endothelin|thrombin|forskolin|glucagon|GLP|danuglipron|adrenal|isoproterenol|angiotensin|iloprost|"
     r"nicotine|caffeine|wnt|BMP|ionomycin|histamine|bradykinin|estradiol|\bE2\b|androgen|DHT|R1881|leptin|EPO|"
     r"\bAA\b|palmitic|oleic|propofol|24D|ATP\b|ADP\b|LPA|S1P|carbachol|insul"),
    ("drug/inhibitor (incl. kinase inhibitors, antibodies)",
     r"nib\b|nib_|nib\s|mab\b|mab_|mab\s|inhib|i_|i\b|i\s|\bi\(|[0-9.]+\s?[nuμµ]M|ng/mL|ug/ml|µg|AZD|GDC|PF\-?\d|"
     r"LY\d|BI\-?D?\d|MK\-?\d|SB\d|SP\d|U0126|PD\d|LRRK2-IN|IN-\d|JNK-IN|VX\-?\d|VE-\d|KU-?\d|Ku-794|CHIR|"
     r"rapamycin|torin|wortmannin|LY294002|staurosporine|H89|Go\d|Gö|GF109|MLN|ZM\d|AT\d|BX\d|AC220|PLX|"
     r"metformin|romidepsin|vorinostat|bortezomib|noscapine|flavopiridol|genistein|ginsenoside|SHP099|TPCA|RO3306|"
     r"A48[56]|CUDC|edelfosine|maraviroc|ibrutinib|TBCA|KN\d|NU\d|DNA_PK|RG\d|DG2|NMS|tofacitinib|BRAFi|AO\b|"
     r"DMSO"),
]
BASELINE = ("baseline / sample-type contrast (cell line, tissue, tumour)",
            r"cancer|tumou?r|normal|risk|^[A-Z0-9\-]{2,12}$|JHU|LNCaP|C4-2|SW\d|VMCUB|CAF|PazR|P31/Fuj|G\d$|"
            r"NCI-H|patient|grade|subtype|resistant|sensitive")


def classify_condition(txt: str) -> str:
    t = str(txt)
    for name, rx in RULES:
        if re.search(rx, t, flags=re.I):
            return name
    if re.search(BASELINE[1], t, flags=re.I if not BASELINE[1].startswith("^[A-Z") else 0):
        return BASELINE[0]
    return "other/unclear"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--uniprot", default="external/uniprot")
    ap.add_argument("--out", default="results/data_profile")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    S = {}

    def w(df, name, index=False):
        df.to_csv(out / f"{name}.tsv", sep="\t", index=index)

    # ---------------- 01 conditions ----------------
    c = pd.read_csv(Path(a.data) / "01_phosphoproteomics_conditions.csv")
    c = c.drop(columns=[x for x in c.columns if x.startswith("Unnamed")])
    S["n_conditions"] = len(c)
    S["n_unique_condition_text"] = int(c.Condition.nunique())
    S["n_pmids"] = int(c.PMID.nunique())
    S["n_samples_celltypes"] = int(c.Sample.nunique())
    S["source_counts"] = c.source.value_counts().to_dict()
    c["perturbation_type"] = c.Condition.map(classify_condition)
    for col in ["Labeling.method", "Enrichment.method", "Sample", "Timepoint", "source", "perturbation_type",
                "Identification.software"]:
        w(c[col].fillna("NA").value_counts().rename_axis(col).rename("n_conditions").reset_index(),
          "compendium_conditions_by_" + col.split(".")[0].lower())
    w(c.groupby("PMID").agg(n_conditions=("conditionID", "size"), source=("source", "first"),
                            n_samples=("Sample", "nunique"), labeling=("Labeling.method", "first"),
                            example=("Condition", "first")).sort_values("n_conditions", ascending=False).reset_index(),
      "compendium_conditions_by_pmid")
    S["pmid_36926954_share"] = round(float((c.PMID == 36926954).mean()), 3)

    # ---------------- 02 sites ----------------
    s = pd.read_csv(Path(a.data) / "02_phosphosite_information.csv")
    s = s.drop(columns=[x for x in s.columns if x.startswith("Unnamed")])
    x = s.peptide_UID.str.extract(r"^(?P<acc>[^_]+)_(?P<pos>\d+)_M(?P<mult>\d+)$")
    S["n_sites"] = len(s)
    S["n_site_ids_not_matching_ACC_POS_Mn"] = int(x.isna().any(axis=1).sum())
    S["multiplicity_suffix_counts"] = x["mult"].value_counts().to_dict()
    s["acc"], s["pos"] = x["acc"], x["pos"].astype(int)
    s["residue"] = s.sequence_window.str[7]
    S["n_proteins"] = int(s.acc.nunique())
    S["n_isoform_accs"] = int(s.acc.str.contains("-").sum())
    S["residue_counts"] = s.residue.value_counts().to_dict()
    S["windows_with_lowercase"] = int(s.sequence_window.str.contains("[a-z]").sum())
    S["windows_with_padding"] = int(s.sequence_window.str.contains("_").sum())
    S["windows_non_15"] = int((s.sequence_window.str.len() != 15).sum())
    S["n_acc_pos_duplicates"] = int(s.duplicated(["acc", "pos"]).sum())

    # ---------------- 03 quantification (wide) ----------------
    wide = pq.read_table(Path(a.cache) / "compendium_wide.parquet").to_pandas().set_index("peptide_UID")
    M = wide.to_numpy(np.float32)
    fin = np.isfinite(M)
    S["n_nonmissing"] = int(fin.sum())
    S["matrix_density"] = round(float(fin.mean()), 4)
    per_site = fin.sum(1)
    per_cond = fin.sum(0)
    qs = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]
    w(pd.DataFrame({"quantile": qs, "conditions_per_site": np.quantile(per_site, qs),
                    "sites_per_condition": np.quantile(per_cond, qs)}), "compendium_coverage_quantiles")
    bins = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
    w(pd.cut(pd.Series(per_site), bins=bins, right=False).value_counts().sort_index().rename_axis("conditions_per_site_bin")
      .rename("n_sites").reset_index().astype(str), "compendium_site_coverage_hist")
    pc = pd.DataFrame({"conditionID": wide.columns, "n_sites": per_cond}).merge(
        c[["conditionID", "PMID", "Condition", "Sample", "Labeling.method", "perturbation_type"]], on="conditionID")
    w(pc.sort_values("n_sites"), "compendium_per_condition_coverage")
    S["sites_with_0_values"] = int((per_site == 0).sum())
    for t in [10, 50, 100, 500]:
        S[f"sites_with_ge{t}_conditions"] = int((per_site >= t).sum())
    vals = M[fin]
    S["log2fc_quantiles"] = dict(zip(["q0.001", "q0.01", "q0.5", "q0.99", "q0.999"],
                                     [round(float(v), 3) for v in np.quantile(vals, [0.001, 0.01, 0.5, 0.99, 0.999])]))
    S["log2fc_exact_zero_frac"] = round(float((vals == 0).mean()), 4)
    S["log2fc_abs_gt1_frac"] = round(float((np.abs(vals) > 1).mean()), 4)
    # per-study median log2FC SD (scale heterogeneity)
    pmid_of = wide.columns.map(dict(zip(c.conditionID, c.PMID)))
    sd = pd.Series(np.nanstd(np.where(fin, M, np.nan), axis=0), index=wide.columns)
    w(pd.DataFrame({"PMID": pmid_of, "sd": sd.values}).groupby("PMID").sd.median().rename("median_condition_sd")
      .reset_index().sort_values("median_condition_sd"), "compendium_scale_by_pmid")
    S["condition_sd_quantiles"] = dict(zip(["q0.05", "q0.5", "q0.95"], [round(float(v), 3) for v in np.nanquantile(sd, [0.05, 0.5, 0.95])]))
    del M

    # ---------------- kinhib ----------------
    k = pd.read_csv(Path(a.data) / "kinhib_inhib_comb__1.csv", index_col=0)
    S["kinhib_shape_kinases_x_inhibitors"] = list(k.shape)
    S["kinhib_values"] = sorted(map(int, np.unique(k.values)))
    S["kinhib_kinases_with_ge1_inhibitor"] = int((k.sum(axis=1) > 0).sum())
    S["kinhib_inhibitors_with_0_targets"] = k.columns[k.sum(axis=0) == 0].tolist()
    w(k.sum(axis=0).sort_values(ascending=False).rename_axis("inhibitor").rename("n_target_kinases").reset_index(), "kinhib_targets_per_inhibitor")
    w(k.sum(axis=1).value_counts().sort_index().rename_axis("n_inhibitors").rename("n_kinases").reset_index(), "kinhib_inhibitors_per_kinase")
    # match inhibitor names in Condition (normalised alnum, substring); Hijazi = PMID 31959955
    rows = []
    for inh in k.columns:
        n = norm(inh)
        hit = c[c.Condition.map(norm).str.contains(n, regex=False)]
        for pm, g in hit.groupby("PMID"):
            rows.append((inh, pm, len(g), ";".join(g.Sample.astype(str).unique()), ";".join(g.Condition.head(3))))
    mt = pd.DataFrame(rows, columns=["inhibitor", "PMID", "n_conditions", "samples", "example_conditions"])
    w(mt, "kinhib_condition_matches")
    hij = c[c.PMID == 31959955].copy()
    hij["inhibitor"] = hij.Condition.str.split("_", n=1).str[1]
    hij["in_kinhib"] = hij.inhibitor.isin(k.columns)
    S["hijazi_conditions"] = int(len(hij))
    S["hijazi_cell_lines"] = hij.Sample.value_counts().to_dict()
    S["hijazi_unique_treatments"] = int(hij.inhibitor.nunique())
    S["hijazi_treatments_not_in_kinhib"] = sorted(set(hij.inhibitor[~hij.in_kinhib]))
    S["kinhib_inhibitors_not_in_hijazi_conditions"] = sorted(set(k.columns) - set(hij.inhibitor))
    S["kinhib_inhibitors_matched_in_other_studies"] = sorted(set(mt.inhibitor[mt.PMID != 31959955]))
    try:
        from ksp import uniprot as _up
        pkg = set(_up.parse_pkinfam(Path(a.uniprot) / "pkinfam" / "pkinfam.txt").gene_uc)
        kn = pd.Index(k.index.astype(str).str.upper())
        S["kinhib_row_names_not_in_pkinfam"] = sorted(set(kn) - pkg)
        S["kinhib_targeted_kinases_in_pkinfam"] = int(kn[(k.sum(axis=1) > 0).values].isin(pkg).sum())
        S["pkinfam_human_kinases_with_ge1_kinhib_inhibitor"] = int(len(set(kn[(k.sum(axis=1) > 0).values]) & pkg))
        S["pkinfam_human_kinases_total"] = len(pkg)
    except FileNotFoundError:
        pass
    w(hij[["conditionID", "Condition", "Sample", "inhibitor", "in_kinhib", "Control", "Timepoint"]], "kinhib_hijazi_conditions")

    # ---------------- PSP overlap ----------------
    psp = pd.read_parquet(Path(a.cache) / "psp_human_mapped_rows.parquet")
    psp = psp[psp.kin_human_mapped & psp.sub_human_mapped]
    site_cov = pd.Series(per_site, index=pd.MultiIndex.from_arrays([s.acc, s.pos]))
    win = pd.Series(s.sequence_window.values, index=site_cov.index)
    rows = []
    for name, sub in [("PSP human kinase x human substrate", psp[(psp.kin_organism == "human") & (psp.sub_organism == "human")]),
                      ("PSP HMR pooled (human-mapped)", psp[psp.kin_organism.isin(["human", "mouse", "rat"]) & psp.sub_organism.isin(["human", "mouse", "rat"])]),
                      ("PSP all species pooled (human-mapped)", psp)]:
        st = sub[["h_sub_acc_full", "h_position"]].drop_duplicates()
        idx = pd.MultiIndex.from_arrays([st.h_sub_acc_full, st.h_position.astype(int)])
        cov = site_cov.reindex(idx)
        n_in = int(cov.notna().sum())
        cv = cov.dropna()
        # pair-level coverage
        pr = sub.drop_duplicates(["kin_human_gene", "h_sub_acc_full", "h_position"])
        pidx = pd.MultiIndex.from_arrays([pr.h_sub_acc_full, pr.h_position.astype(int)])
        pin = site_cov.reindex(pidx).notna()
        # isoform-specific PSP sites never match (compendium has canonical accessions only)
        rows.append(dict(set=name, n_sites=len(st), n_sites_isoform_acc=int(st.h_sub_acc_full.str.contains("-").sum()),
                         n_in_compendium=n_in, frac=round(n_in / len(st), 3),
                         n_in_compendium_ge10_conditions=int((cv >= 10).sum()),
                         n_in_compendium_ge100_conditions=int((cv >= 100).sum()),
                         median_conditions_if_present=float(cv.median()) if len(cv) else np.nan,
                         n_pairs=len(pr), n_pairs_site_in_compendium=int(pin.sum()),
                         n_kinases_with_ge1_pair_in_compendium=int(pr[pin.values].kin_human_gene.nunique()),
                         n_kinases_with_ge5_pairs_in_compendium=int((pr[pin.values].kin_human_gene.value_counts() >= 5).sum())))
        if name.startswith("PSP HMR"):
            w(cv.astype(int).pipe(lambda v: pd.cut(v, bins=bins, right=False).value_counts().sort_index())
              .rename_axis("conditions_per_site_bin").rename("n_psp_sites").reset_index().astype(str), "compendium_psp_hmr_site_coverage_hist")
            # window agreement for matched sites (PSP 15-mer vs compendium 15-mer)
            hm = sub.drop_duplicates(["h_sub_acc_full", "h_position"])
            hm = hm[hm.sub_organism == "human"]
            ix = pd.MultiIndex.from_arrays([hm.h_sub_acc_full, hm.h_position.astype(int)])
            cw = win.reindex(ix)
            ok = cw.notna().values
            S["psp_human_sites_in_compendium_window_identical"] = int((cw.values[ok] == hm.window.values[ok]).sum())
            S["psp_human_sites_in_compendium_window_compared"] = int(ok.sum())
    ov = pd.DataFrame(rows)
    w(ov, "compendium_psp_overlap")

    # ---------------- PSP Phosphorylation_site_dataset overlap ----------------
    hs_path = Path(a.cache) / "psp_human_sites.parquet"
    if hs_path.exists():
        hs = pd.read_parquet(hs_path)
        hs["ev_default"] = io_psp.site_evidence_mask(hs)
        idx = pd.MultiIndex.from_arrays([hs.acc_full, hs.position.astype(int)])
        hs["comp_n"] = site_cov.reindex(idx).values
        hs["in_comp"] = hs.comp_n.notna()
        rows = []
        for name, msk in [("all human STY sites", hs.ev_default | True), ("LT>=1 or MS_total>=2", hs.ev_default),
                          ("LT>=1", hs.lt_lit >= 1), ("KS-annotated (HMR pooled)", hs.ks_hmr),
                          ("no KS annotation (dark)", ~hs.ks_hmr), ("dark & LT>=1 or MS_total>=2", ~hs.ks_hmr & hs.ev_default)]:
            x = hs[msk]
            rows.append(dict(psp_site_set=name, n_sites=len(x), n_in_compendium=int(x.in_comp.sum()),
                             frac_in_compendium=round(float(x.in_comp.mean()), 4),
                             n_in_compendium_ge10_conditions=int((x.comp_n >= 10).sum()),
                             n_in_compendium_ge100_conditions=int((x.comp_n >= 100).sum())))
        w(pd.DataFrame(rows), "compendium_psp_sitedb_overlap")
        comp_keys = pd.MultiIndex.from_arrays([s.acc, s.pos])
        known = pd.MultiIndex.from_arrays([hs.acc_full, hs.position.astype(int)])
        in_psp = comp_keys.isin(known)
        ks_keys = pd.MultiIndex.from_arrays([hs.acc_full[hs.ks_hmr], hs.position[hs.ks_hmr].astype(int)])
        S["compendium_sites_in_psp_sitedb"] = int(in_psp.sum())
        S["compendium_sites_in_psp_sitedb_frac"] = round(float(in_psp.mean()), 4)
        S["compendium_sites_with_ks_hmr"] = int(comp_keys.isin(ks_keys).sum())
        S["compendium_sites_not_in_psp"] = int((~in_psp).sum())
        win_known = s.sequence_window.isin(set(hs.window)).values
        S["compendium_sites_not_in_psp_but_window_in_psp"] = int((~in_psp & win_known).sum())
        S["compendium_sites_not_in_psp_median_conditions"] = float(np.median(per_site[~in_psp & ~win_known]))
        S["compendium_sites_in_psp_median_conditions"] = float(np.median(per_site[in_psp]))

    # ---------------- UniProt: compendium windows, proteome STY universe ----------------
    up_dir = Path(a.uniprot)
    if (up_dir / "proteome" / "UP000005640_all_isoforms.fasta.gz").exists():
        from ksp import uniprot as up
        fa = up.read_fasta(up_dir / "proteome" / "UP000005640_all_isoforms.fasta.gz")
        seqs = dict(zip(fa.acc_full, fa.seq))
        chk = io_psp.check_windows_against_sequences(
            s.rename(columns={"acc": "sub_acc_full", "pos": "position", "sequence_window": "window"}), seqs)
        w(chk.status.value_counts().rename_axis("status").rename("n_sites").reset_index(), "compendium_window_vs_uniprot")
        S["compendium_window_exact_frac"] = round(float((chk.status == "exact").mean()), 4)
        S["compendium_window_status"] = chk.status.value_counts().to_dict()
        S["compendium_accs_not_in_UP000005640"] = int(chk.loc[chk.status == "no_sequence", "sub_acc_full"].nunique())
        S["compendium_accs_swissprot"] = int(s.acc[s.acc.isin(set(fa.acc[fa.db == "sp"]))].nunique())
        S["compendium_accs_trembl"] = int(s.acc[s.acc.isin(set(fa.acc[fa.db == "tr"]))].nunique())
        # STY universe
        rows = []
        for name, sel in [("Swiss-Prot canonical", (fa.db == "sp") & (fa.isoform == "")),
                          ("Swiss-Prot canonical + isoforms", fa.db == "sp"),
                          ("UP000005640 all canonical (SP+TrEMBL)", fa.isoform == ""),
                          ("UP000005640 all incl. isoforms", pd.Series(True, index=fa.index))]:
            sq = fa.loc[sel, "seq"]
            cnt = {r: int(sq.str.count(r).sum()) for r in "STY"}
            rows.append(dict(set=name, n_sequences=int(sel.sum()), n_residues=int(sq.str.len().sum()), **cnt,
                             STY=sum(cnt.values())))
        uni = pd.DataFrame(rows)
        # distinct 15-mer-unique windows in SP canonical+isoforms (isoforms mostly duplicate canonical sites)
        sp = fa[fa.db == "sp"]
        wins = set()
        for sq in sp.seq:
            p = "_" * 7 + sq + "_" * 7
            for m_ in re.finditer("[STY]", sq):
                wins.add(p[m_.start(): m_.start() + 15])
        S["swissprot_distinct_STY_15mers_incl_isoforms"] = len(wins)
        sp_can = fa[(fa.db == "sp") & (fa.isoform == "")]
        can_acc = set(sp_can.acc)
        comp_sp = s[s.acc.isin(can_acc)]
        S["compendium_sites_on_swissprot_canonical"] = int(len(comp_sp))
        tot = int(uni.loc[uni.set == "Swiss-Prot canonical", "STY"].iloc[0])
        S["compendium_frac_of_swissprot_canonical_STY"] = round(len(comp_sp) / tot, 4)
        psp_h = psp[(psp.sub_organism.isin(["human", "mouse", "rat"])) & psp.kin_organism.isin(["human", "mouse", "rat"])][["h_sub_acc_full", "h_position"]].drop_duplicates()
        S["psp_hmr_sites_frac_of_swissprot_canonical_STY"] = round(len(psp_h) / tot, 5)
        w(uni, "proteome_sty_universe")
        S["proteome_sty_universe"] = uni.set_index("set")[["n_sequences", "STY"]].to_dict("index")
    json.dump(S, open(out / "compendium_summary.json", "w"), indent=1, default=str)
    print(json.dumps(S, indent=1, default=str))
    print(ov.to_string())


if __name__ == "__main__":
    main()
