#!/usr/bin/env python
"""F6b inhibitor-response features (not label-derived), DESIGN §7 F6b.

Outputs (kinase x site float16, via ksp.store; name suffix = --variant):
  cache/features/phospho/f6b_ridge_beta{_v}.npy   ridge coefficient (scaled log2FC per unit occupancy; <0 = decrease)
  cache/features/phospho/f6b_ridge_z{_v}.npy      beta / posterior se
  cache/features/phospho/f6b_contrast{_v}.npy     mean y (conditions inhibiting k >= 50 %) - mean y (|occ_k| < 10 %)
  cache/features/phospho/f6b_kinase_maxcos{_v}.npy  kinase vector: max cosine of k's inhibition profile with any other
results/features/phospho/f6b_conditions{_v}.tsv  every condition: parsed drugs, doses, source, used/dropped reason
results/features/phospho/f6b_kinases{_v}.tsv, f6b_nonseparable_pairs{_v}.tsv, f6b_report{_v}.json

Usage: CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_f6b.py
       [--variant clean_bench --exclude-pmids 31959955,... --exclude-drugs AC220,...] [--lam 1.0] [--max-minutes 120]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import store  # noqa: E402
from ksp.kinases import load_resolver  # noqa: E402
from ksp.features import phospho as ph  # noqa: E402
from ksp.features import phospho_inhib as pi  # noqa: E402


def run(variant="", exclude_pmids=(), exclude_drugs=(), lam=1000.0, max_minutes=120.0, device="cuda"):
    t0 = time.time()
    sfx = f"_{variant}" if variant else ""
    res = load_resolver()
    kin = ph.kinase_index()
    sites = ph.site_index()
    A, tab, rep, keep = pi.build_inhibition_matrix(res, kin, max_minutes=max_minutes, exclude_pmids=exclude_pmids,
                                                   exclude_drugs=exclude_drugs, verbose=False)
    tab.to_csv(ph.RESULTS / f"f6b_conditions{sfx}.tsv", sep="\t", index=False)
    Au = A[keep]
    cond_ids = tab.conditionID[keep].tolist()
    expo = np.abs(Au).sum(0)
    kcols = np.nonzero(expo >= 1.0)[0]
    rows, Y = pi.assemble_Y(cond_ids, len(sites), exclude_pmids=exclude_pmids)
    nobs = np.isfinite(Y).sum(1)
    ok = nobs >= 3
    rows, Y = rows[ok], Y[ok]
    print(f"F6b: {len(cond_ids)} conditions, {len(kcols)} exposed kinases, {len(rows)} sites (>=3 obs) "
          f"[{time.time()-t0:.0f}s]", flush=True)
    beta, z, con, s0 = pi.ridge_deconvolution(Y, Au[:, kcols], lam=lam, device=device)
    out = {}
    for name, M in (("f6b_ridge_beta", beta), ("f6b_ridge_z", z), ("f6b_contrast", con)):
        F = np.full((len(kin), len(sites)), np.nan, np.float16)
        F[np.ix_(kcols, rows)] = M.astype(np.float16)
        store.save_matrix(f"phospho/{name}{sfx}", F, "kinase", "site",
                          meta=dict(variant=variant, exclude_pmids=list(map(str, exclude_pmids)),
                                    exclude_drugs=list(exclude_drugs), lam=lam, max_minutes=max_minutes,
                                    label_derived=False))
        out[name] = dict(n_finite=int(np.isfinite(F).sum()), n_kinases=int(np.isfinite(F).any(1).sum()),
                         n_sites=int(np.isfinite(F).any(0).sum()))
    pairs, ktab = pi.nonseparable_pairs(Au, kin)
    pairs.to_csv(ph.RESULTS / f"f6b_nonseparable_pairs{sfx}.tsv", sep="\t", index=False)
    ktab.to_csv(ph.RESULTS / f"f6b_kinases{sfx}.tsv", sep="\t", index=False)
    store.save_matrix(f"phospho/f6b_kinase_maxcos{sfx}", ktab.max_cosine_other.to_numpy(np.float32), "kinase", "none",
                      meta=dict(desc="max cosine similarity of the kinase's inhibition profile to any other kinase"))
    rep.update(variant=variant, lam=lam, s0_sq=s0, n_exposed_kinases=int(len(kcols)), n_sites=int(len(rows)),
               n_nonseparable_pairs=int(len(pairs)), n_nonseparable_kinases=int(ktab.nonseparable.sum()),
               features=out, seconds=round(time.time() - t0, 1))
    (ph.RESULTS / f"f6b_report{sfx}.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({k: v for k, v in rep.items() if k not in ("affinity_sources",
                                                               "hijazi_kinhib_rows_unresolved")}, default=str))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="")
    ap.add_argument("--exclude-pmids", default="")
    ap.add_argument("--exclude-drugs", default="")
    ap.add_argument("--lam", type=float, default=1000.0)
    ap.add_argument("--max-minutes", type=float, default=120.0)
    a = ap.parse_args()
    run(a.variant, [int(x) for x in a.exclude_pmids.split(",") if x],
        [x for x in a.exclude_drugs.split(",") if x], a.lam, a.max_minutes)


if __name__ == "__main__":
    main()
