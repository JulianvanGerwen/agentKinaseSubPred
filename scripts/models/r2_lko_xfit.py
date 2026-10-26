"""Dense cross-fitted features for leave-kinase-out T2: all DEV labels minus the held-out kinase fold's labels;
rows = all DEV U_ms sites (inner-fold fits). Output tag lko<j>/all. Usage: r2_lko_xfit.py --family f2a_sep --lko 0 1 2 3 4"""
import argparse, time
from ksp import splits
from ksp.models.r2 import xfit
ap = argparse.ArgumentParser(); ap.add_argument("--family", nargs="+"); ap.add_argument("--lko", nargs="+", type=int)
a = ap.parse_args()
kf = splits.load_kinase_folds()
for j in a.lko:
    held = set(kf.kinase_acc[kf.lko_fold == j])
    for fam in a.family:
        t = time.time()
        p = xfit.build(fam, None, tag=f"lko{j}", exclude_kinases=held)
        print(f"{fam} lko{j}: {p} {time.time() - t:.0f}s ({len(held)} held-out kinases)", flush=True)
