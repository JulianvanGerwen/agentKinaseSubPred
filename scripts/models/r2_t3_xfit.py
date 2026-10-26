"""Pre-build T3 cross-fitted stores (labels of LKO fold j removed; rows = DEV U_obs sites, eval rows = TEST-protein
sites). Same calls as r2_final.py --target t3 (which skips existing stores). DEV labels only.
Usage: r2_t3_xfit.py --family f2a_sep f3nb f6a --lko 0 1 2 3 4"""
import argparse, time
import numpy as np
from ksp import index, splits
from ksp.models.r2 import xfit
ap = argparse.ArgumentParser(); ap.add_argument("--family", nargs="+"); ap.add_argument("--lko", nargs="+", type=int); a = ap.parse_args()
S = index.load_sites()
dev = S.acc.isin(set(splits.load_cv_folds().index)).values
test = S.acc.isin(splits.load_test_proteins()).values
tr, ev = np.where(dev)[0], np.where(test & ~dev)[0]
kf = splits.load_kinase_folds()
for j in a.lko:
    held = set(kf.kinase_acc[kf.lko_fold == j])
    for fam in a.family:
        t = time.time()
        p = xfit.build(fam, None, tag=f"t3lko{j}", train_site_idx=tr, eval_site_idx=ev, exclude_kinases=held)
        print(f"{fam} t3lko{j}: {p} {time.time() - t:.0f}s", flush=True)
