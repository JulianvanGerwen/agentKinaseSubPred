#!/usr/bin/env python
"""Study lists for clean F6 variants.

(1) DESIGN §16 A12 — compendium study <-> KS-label-source list: which compendium PMIDs are cited as evidence for
    kinase-substrate records (OmniPath enzsub references incl. PhosphoSite-derived ones; iPTMnet phosphorylation rows
    with an enzyme) or for PSP Regulatory_sites rows (PSP KS rows themselves carry no PMIDs).
    -> results/features/phospho/label_source_studies.tsv
(2) DESIGN §6.7 — benchmarKIN perturbation studies (Hernandez-Armenta metadata PMIDs; Hijazi = PMID 31959955) and the
    CPTAC tumour benchmark cohorts, matched to compendium PMIDs (the benchmark metadata and the compendium have no
    PXD columns, so matching is by PMID only).
    -> results/features/phospho/benchmarkin_studies.tsv
(3) results/features/phospho/exclude_lists.json: {"clean_bench": {...}, "f6_clean": {...}} to pass as
    exclude_pmids / exclude_cptac to the builders.
Usage: PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_study_lists.py
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.features import phospho as ph  # noqa: E402

EXT = ROOT / "external"


def main():
    c = ph.load_conditions()
    comp = c.groupby("PMID").agg(n_conditions=("conditionID", "size"),
                                 blocks=("block_id", lambda x: ";".join(sorted(set(x)))),
                                 example=("Condition", "first")).reset_index()
    comp["PMID"] = comp.PMID.astype(int)
    pm = set(comp.PMID)
    # PSP regulatory sites
    reg = pd.read_csv(ROOT / "data" / "Regulatory_sites", sep="\t", skiprows=3, encoding="latin1", dtype=str)
    reg = reg[reg.MOD_RSD.str.endswith("-p", na=False)]
    reg_cnt = {}
    for s in reg.PMIDs.dropna():
        for p in re.findall(r"\d+", s):
            if int(p) in pm:
                reg_cnt[int(p)] = reg_cnt.get(int(p), 0) + 1
    # OmniPath enzsub (phosphorylation) references
    om = pd.read_csv(EXT / "omnipath" / "omnipath_enzsub.tsv", sep="\t", dtype=str)
    om = om[om.modification == "phosphorylation"]
    om_cnt, om_src = {}, {}
    for refs in om.references.dropna():
        for tok in refs.split(";"):
            if ":" not in tok:
                continue
            src, p = tok.rsplit(":", 1)
            if p.isdigit() and int(p) in pm:
                om_cnt[int(p)] = om_cnt.get(int(p), 0) + 1
                om_src.setdefault(int(p), set()).add(src)
    # iPTMnet phosphorylation rows with an enzyme
    ip = pd.read_csv(EXT / "iptmnet" / "ptm.txt", sep="\t", header=None, dtype=str,
                     names=["ptm_type", "source", "sub", "sub_gene", "org", "site", "enz", "enz_gene", "note", "pmid"])
    ip = ip[(ip.ptm_type == "PHOSPHORYLATION") & ip.enz.notna() & ip.org.str.contains("Human", na=False)]
    ip_cnt, ip_src = {}, {}
    for p, src in zip(ip.pmid.fillna(""), ip.source):
        for q in re.findall(r"\d+", p):
            if int(q) in pm:
                ip_cnt[int(q)] = ip_cnt.get(int(q), 0) + 1
                ip_src.setdefault(int(q), set()).add(src)
    comp["n_psp_regulatory_rows"] = comp.PMID.map(reg_cnt).fillna(0).astype(int)
    comp["n_omnipath_ks_refs"] = comp.PMID.map(om_cnt).fillna(0).astype(int)
    comp["omnipath_sources"] = comp.PMID.map(lambda p: ";".join(sorted(om_src.get(p, []))))
    comp["n_iptmnet_ks_rows"] = comp.PMID.map(ip_cnt).fillna(0).astype(int)
    comp["iptmnet_sources"] = comp.PMID.map(lambda p: ";".join(sorted(ip_src.get(p, []))))
    comp["ks_label_source"] = (comp.n_omnipath_ks_refs > 0) | (comp.n_iptmnet_ks_rows > 0)
    comp["regulatory_source"] = comp.n_psp_regulatory_rows > 0
    comp.sort_values(["ks_label_source", "n_conditions"], ascending=False).to_csv(
        ph.RESULTS / "label_source_studies.tsv", sep="\t", index=False)
    # benchmarKIN
    hm = pd.read_csv(EXT / "benchmarkin" / "benchmark_metadata_hernandez.csv")
    hj = pd.read_csv(EXT / "benchmarkin" / "benchmark_metadata_hijazi.csv")
    b = hm.groupby("PMID").agg(n_experiments=("id", "size"), targets=("target", lambda x: ";".join(sorted(set(
        map(str, x))))), example=("Description", "first")).reset_index()
    b["benchmark"] = "hernandez"
    b = pd.concat([b, pd.DataFrame([dict(PMID=31959955, n_experiments=len(hj), benchmark="hijazi",
                                         targets=";".join(sorted(set(hj.target.astype(str)))), example="")])])
    b["PMID"] = pd.to_numeric(b.PMID, errors="coerce").astype("Int64")
    b["in_compendium"] = b.PMID.isin(pm)
    b["compendium_n_conditions"] = b.PMID.map(dict(zip(comp.PMID, comp.n_conditions))).fillna(0).astype(int)
    cptac = [x.split("_")[0].upper() for x in
             ["brca", "ccrcc", "gbm", "hnscc", "lscc", "luad", "ucec"]]
    b = pd.concat([b, pd.DataFrame([dict(benchmark="cptac_tumour", example="cohort", targets=";".join(cptac))])])
    b.to_csv(ph.RESULTS / "benchmarkin_studies.tsv", sep="\t", index=False)
    bench_pm = sorted(int(p) for p in b.PMID.dropna() if int(p) in pm)
    ex = {"clean_bench": {"exclude_pmids": bench_pm, "exclude_cptac": True,
                          "note": "DESIGN 6.7: benchmarKIN perturbation PMIDs in the compendium + all CPTAC "
                                  "(the tumour benchmark uses 7 of our 10 cohorts; all are dropped)"},
          "clean_bench_perturbation": {"exclude_pmids": bench_pm, "exclude_cptac": False},
          "clean_bench_tumour": {"exclude_pmids": [], "exclude_cptac": True},
          "f6_clean": {"exclude_pmids": sorted(int(p) for p in comp.PMID[comp.ks_label_source]),
                       "exclude_cptac": False,
                       "note": "DESIGN 16 A12: compendium studies cited as evidence for KS records (OmniPath/iPTMnet)"}}
    (ph.RESULTS / "exclude_lists.json").write_text(json.dumps(ex, indent=1))
    print(f"compendium PMIDs {len(pm)}; KS-label-source {comp.ks_label_source.sum()} "
          f"({comp.n_conditions[comp.ks_label_source].sum()} conditions); regulatory-cited {comp.regulatory_source.sum()}")
    print(f"benchmarKIN Hernandez PMIDs {hm.PMID.nunique()}, in compendium {len(bench_pm)} "
          f"({b.compendium_n_conditions.sum()} conditions)")


if __name__ == "__main__":
    main()
