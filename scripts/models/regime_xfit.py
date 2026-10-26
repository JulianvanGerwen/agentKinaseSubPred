"""Refit the cross-fitted label-derived features of the frozen R2 config (F2a f2a_sep, F3-label f3nb, F6a f6a) from a
training-label VARIANT, through ksp.models.r2.xfit.build -> ksp.crossfit.crossfit (read-only use of R2 code).
Stores: cache/regime/xfit/<variant>/fold<f>/<family>.npy (+ .json with n_train_labels, variant, runtime).
Usage: regime_xfit.py --variant L1 --folds 0 1 --family f2a_sep f3nb f6a [--threads 3] [--compare]
--compare: compare the new store with R2's cache/r2/xfit/cv store (sanity for variant L2: must be identical)."""
import argparse, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--variant", required=True)
ap.add_argument("--folds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
ap.add_argument("--family", nargs="+", default=["f2a_sep", "f3nb", "f6a"])
ap.add_argument("--threads", type=int, default=3)
ap.add_argument("--compare", action="store_true")
ap.add_argument("--overwrite", action="store_true")
a = ap.parse_args()
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = str(a.threads)
import numpy as np
try:
    import torch
    torch.set_num_threads(a.threads)
except Exception:
    pass
from ksp.models.regime import variants as V
V.patch()
from ksp.models.r2 import xfit

tag = V.xfit_tag(a.variant)
for f in a.folds:
    for fam in a.family:
        t = time.time()
        p = xfit.build(fam, f, seed=0, variant=a.variant, tag=tag, overwrite=a.overwrite)
        meta = json.loads(p.with_suffix(".json").read_text())
        msg = f"{a.variant} {fam} fold {f}: {p} n_train_labels={meta['n_train_labels']} {time.time() - t:.0f}s"
        if a.compare:
            new, _ = xfit.load(fam, f, tag=tag)
            old, mo = xfit.load(fam, f)
            S = np.load(p.with_suffix(".block.npy")) != -2
            dif = 0.0
            nn = 0
            for c in range(new.shape[0]):
                x, y = np.asarray(new[c][:, S]), np.asarray(old[c][:, S])
                both = ~np.isnan(x) & ~np.isnan(y)
                nn += int((np.isnan(x) != np.isnan(y)).sum())
                if both.any():
                    dif = max(dif, float(np.abs(x[both] - y[both]).max()))
            msg += f" | vs R2 store: max|diff|={dif:.3g} nan-mismatch={nn} (R2 n_train_labels={mo['n_train_labels']})"
        print(msg, flush=True)
