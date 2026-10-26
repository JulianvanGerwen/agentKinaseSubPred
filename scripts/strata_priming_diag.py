#!/usr/bin/env python
"""Priming diagnostic (docs/strata_cv.md §4): for the curated priming kinases (ksp.strata.priming_curated, offsets
GSK3 +4, CK1 -3, CK2 +1/+3), split CV positives into
  primed              offset residue is S/T and an observed (U_ms) phosphosite      (= ksp.strata 'primed')
  offset_ST_unobs     offset residue is S/T but not in U_ms (potential priming site, never observed)
  offset_not_ST       offset residue is not S/T (or outside the protein)
and report per-kinase E-hard AUROC (M1 metric, A24 fold-normalised) per group for selected runs, via a subset-label
ksp.eval.evaluate (removed positives stay masked through KP(k)). M2 is not computed here (prediction rows are
restricted to the priming kinases). Output: results/strata/cv_priming_diag.tsv
Usage: python -I scripts/strata_priming_diag.py
"""
import dataclasses
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from ksp import eval as E  # noqa: E402
from ksp import strata as ST  # noqa: E402
import strata_cv_report as R  # noqa: E402

TAGS = ["R2", "R3", "B1i", "B1prim", "B7"]


def main():
    data = E.prepare("cv", fold=None)
    S = data.sites
    sid = S.site_id.values
    observed = set(sid[S.in_U_ms.values.astype(bool)])
    win = dict(zip(S.site_id, S.window15))
    cur = {k: v for k, v in ST.priming_curated(data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY_")]).items()
           if v and k in data.positives}
    pm = data.pos_meta[data.pos_meta.kinase_acc.isin(cur)].reset_index(drop=True)
    grp = []
    for k, si in zip(pm.kinase_acc, pm.site_idx):
        s = sid[si]
        g = "offset_not_ST"
        for off in cur[k]:
            o = ST.offset_site_ids([s], [win[s]], off)[0]
            if o is not None and o in observed:
                g = "primed"
                break
            if o is not None:
                g = "offset_ST_unobs"
        grp.append(g)
    pm["group"] = grp
    gene = dict(zip(data.kinases.uniprot_acc, data.kinases.gene))
    rows = []
    for tag in TAGS:
        run_id, _, _ = R.RUNSPEC[tag]
        Pn = R.load_normalised(run_id, data, set(cur))
        for g in ["primed", "offset_ST_unobs", "offset_not_ST"]:
            x = pm[pm.group.eq(g)].reset_index(drop=True)
            pos = {k: np.sort(v.site_idx.values) for k, v in x.groupby("kinase_acc")}
            d = dataclasses.replace(data, positives=pos, pos_meta=x.iloc[:0])
            pk = E.evaluate(Pn, d, missing="worst", fold_normalise=False, keep_items=False).per_kinase
            for r in pk.itertuples():
                rows.append({"tag": tag, "group": g, "kinase_acc": r.kinase_acc, "gene": gene[r.kinase_acc],
                             "n_pos": r.n_pos, "auroc_hard": r.auroc_hard, "top1_recall": r.top1_recall})
        print(tag, "done", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(R.OUT / "cv_priming_diag.tsv", sep="\t", index=False, float_format="%.4f")
    P = D.pivot_table(index=["gene", "group"], columns="tag", values="auroc_hard")
    N = D[D.tag.eq("R2")].set_index(["gene", "group"]).n_pos
    P.insert(0, "n_pos", N)
    print(P.round(3).to_string())


if __name__ == "__main__":
    main()
