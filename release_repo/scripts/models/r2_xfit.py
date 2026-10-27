"""Build R2 dense cross-fitted label-derived features (ksp.models.r2.xfit) for CV folds / final.
Usage: r2_xfit.py --family f2a f3nb f6a --folds 0 1 2 3 4 [--toy]"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from ksp.models.r2 import xfit
from ksp import crossfit, index, canary

ap = argparse.ArgumentParser()
ap.add_argument("--family", nargs="+", default=["f2a", "f3nb", "f6a"])
ap.add_argument("--folds", nargs="+", default=["0", "1", "2", "3", "4"])
ap.add_argument("--toy", action="store_true", help="leak test with crossfit.LeakyLabelFlag through DenseExpand")
ap.add_argument("--overwrite", action="store_true")
a = ap.parse_args()

if a.toy:
    xfit.FAMILIES["toy"] = (lambda: crossfit.LeakyLabelFlag, "site", ["is_labelled"])
    xfit.CHUNK["toy"] = 6_000_000
    p = xfit.build("toy", 0, tag="toy", overwrite=True)
    X, meta = xfit.load("toy", 0, tag="toy")
    tr, _, _ = crossfit.outer_split(0, with_canaries=True)
    K = canary.kinases_with_canaries(); kp = {k: i for i, k in enumerate(K.uniprot_acc)}
    sp = {s: j for j, s in enumerate(index.load_sites().site_id)}
    ki = tr.kinase_acc.map(kp).values; si = tr.site_id.map(sp).values
    v = X[0][ki, si]
    print("toy leak test: training positives with is_labelled=1:", int(np.nansum(v)), "of", len(v),
          "| NaN:", int(np.isnan(v).sum()), "| blocks:", meta["blocks"])
    assert np.nansum(v) == 0, "LEAK: training pair saw its own label"
    sys.exit(0)
for f in a.folds:
    fold = None if f == "final" else int(f)
    for fam in a.family:
        t = time.time()
        p = xfit.build(fam, fold, overwrite=a.overwrite)
        print(f"{fam} fold {f}: {p} {time.time() - t:.0f}s", flush=True)
