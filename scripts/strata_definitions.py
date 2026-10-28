#!/usr/bin/env python
"""Write the canonical stratum definitions (ksp.strata) to results/strata/ for transparency:
  stratum_definitions.tsv        rule + source per stratum
  kinase_strata.tsv              one row per index kinase with every kinase-level stratum column (+ K_cv/K_any flags)
  cantley_phospho_preference.tsv per S/T Cantley matrix: max pS/pT log2, qualifying offsets (priming_cantley rule)
  dev_hub_sites.tsv              DEV hub sites (A26): site, n non-paralogous kinase groups, kinases
  dev_site_strata_counts.tsv     DEV evaluable positives per site/pair stratum (CV-side counts only)
Usage: python -I scripts/strata_definitions.py
"""
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import index, labels, splits  # noqa: E402
from ksp import strata as ST  # noqa: E402

OUT = ROOT / "results" / "strata"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ST.definitions_table().to_csv(OUT / "stratum_definitions.tsv", sep="\t", index=False)
    T = ST.kinase_table()
    ks = splits.load_eval_kinase_sets().set_index("kinase_acc")
    T = T.join(ks[["in_K_cv", "in_K_eval", "in_K_any", "n_dev_pos"]], on="kinase_acc")
    T.to_csv(OUT / "kinase_strata.tsv", sep="\t", index=False)
    cp = ST.cantley_phospho_preference()
    K = index.load_kinases().set_index("uniprot_acc")
    cp["gene"] = cp.uniprot_acc.map(K.gene)
    cp["in_index"] = cp.uniprot_acc.isin(K.index)
    cp["offsets"] = cp.offsets.map(lambda o: ",".join(f"{x:+d}" for x in o))
    cp.sort_values("max_phospho_log2", ascending=False).to_csv(OUT / "cantley_phospho_preference.tsv", sep="\t",
                                                               index=False, float_format="%.4f")
    dev = labels.load_dev("all")
    ev = dev[labels.evaluable_positive(dev)].drop_duplicates(["kinase_acc", "site_id"])
    par = labels.paralog_sets()
    g = ST.site_kinase_groups(ev, par)
    hub = g[g >= ST.HUB_MIN_GROUPS]
    genes = ev.groupby("site_id").kinase.agg(lambda s: ",".join(sorted(set(s))))
    pd.DataFrame({"site_id": hub.index, "n_nonparalog_groups": hub.values,
                  "kinases": genes.reindex(hub.index).values}).sort_values(
        "n_nonparalog_groups", ascending=False).to_csv(OUT / "dev_hub_sites.tsv", sep="\t", index=False)
    # counts of DEV evaluable positives (U_ms) per site/pair stratum
    evu = ev[ev.in_U_ms.astype(bool)].copy()
    S = index.load_sites()
    evu["hub"] = evu.site_id.isin(hub.index)
    evu["evidence"] = ST.evidence_class(evu)
    evu["primed"] = ST.primed_mask(evu[["kinase_acc", "site_id"]], ST.priming_curated(), S)
    evu["primed_cantley"] = ST.primed_mask(evu[["kinase_acc", "site_id"]], ST.priming_cantley(), S)
    evu["is_priming"] = evu.kinase_acc.isin(ST.priming_curated())
    evu["is_priming_cantley"] = evu.kinase_acc.isin(ST.priming_cantley())
    rows = [{"stratum": "all", "value": "all", "n_pairs": len(evu), "n_sites": evu.site_id.nunique()}]
    for c in ["hub", "evidence"]:
        for v, x in evu.groupby(c):
            rows.append({"stratum": c, "value": str(v), "n_pairs": len(x), "n_sites": x.site_id.nunique()})
    for c, kc in [("primed", "is_priming"), ("primed_cantley", "is_priming_cantley")]:
        for v, x in evu[evu[kc]].groupby(c):
            rows.append({"stratum": c, "value": "primed" if v else "unprimed", "n_pairs": len(x),
                         "n_sites": x.site_id.nunique()})
    pd.DataFrame(rows).to_csv(OUT / "dev_site_strata_counts.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())
    print(T[T.priming | T.priming_cantley | T.pikk | T.secretory][
        ["gene", "priming", "priming_offsets", "priming_cantley", "priming_cantley_offsets", "pikk", "pikk_sq",
         "secretory", "in_K_cv", "n_dev_pos"]].to_string())


if __name__ == "__main__":
    main()
