"""Final-model (all DEV labels) dense cross-fitted features for deployment-style scoring of U_obs:
training rows = every U_obs site on DEV proteins (inner-fold fits by substrate cluster), eval rows = every U_obs site
on TEST proteins (A1: randomly assigned inner-fold fit). No TEST label is read (ksp.labels.load_dev only).
Usage: r2_final_xfit.py --family f2a_sep f3nb f6a f6a_bench"""
import argparse, time
import numpy as np
from ksp import index, splits
from ksp.models.r2 import xfit
ap = argparse.ArgumentParser(); ap.add_argument("--family", nargs="+"); a = ap.parse_args()
S = index.load_sites()
dev = S.acc.isin(set(splits.load_cv_folds().index)).values
test = S.acc.isin(splits.load_test_proteins()).values
tr, ev = np.where(dev)[0], np.where(test & ~dev)[0]
print(f"train rows {len(tr)}, eval (TEST-protein) rows {len(ev)}, other {len(S) - len(tr) - len(ev)}", flush=True)
for fam in a.family:
    t = time.time()
    p = xfit.build(fam, None, tag="final", train_site_idx=tr, eval_site_idx=ev)
    print(f"{fam} final: {p} {time.time() - t:.0f}s", flush=True)
