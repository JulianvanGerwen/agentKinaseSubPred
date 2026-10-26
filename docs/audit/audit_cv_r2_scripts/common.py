import json, sys
sys.path.insert(0, "/workspace/src")
import numpy as np, pandas as pd
from ksp import splits, index
C = "/workspace/cache/runs/"
R = "/workspace/results/runs/"
ENS = "r2_final_ens5__20261026-191213__s-1"
FUS = "r2_fusion_r3__20261026-191901__s0"
LKO = "r2_final_lko__20261026-202728__s0"
R3F = json.load(open("/workspace/results/models/r3/FINAL_OOF.json"))["cv_fold_runs"]

def site_fold():
    S = index.load_sites()
    return S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values

def load(path, cols=("kinase_idx", "site_idx", "score")):
    p = pd.read_parquet(path, columns=list(cols))
    return p.astype({"kinase_idx": np.int64, "site_idx": np.int64})

def blend_fold(f, sf, ens=None):
    """R2 ensemble + R3 rows of fold f merged, with per-fold z-scores (as scripts/models/r2_fusion.py)."""
    if ens is None:
        ens = load(C + ENS + "/predictions.parquet")
    A = ens[sf[ens.site_idx.values] == f]
    B = load(C + R3F[f] + "/predictions.parquet")
    M = A.merge(B, on=["kinase_idx", "site_idx"], how="inner", suffixes=("_r2", "_r3"))
    info = {"n_r2": len(A), "n_r3": len(B), "n_common": len(M)}
    for c in ["score_r2", "score_r3"]:
        z = M[c].values.astype(np.float64)
        M[c + "_z"] = (z - z.mean()) / z.std()
    return M, info
