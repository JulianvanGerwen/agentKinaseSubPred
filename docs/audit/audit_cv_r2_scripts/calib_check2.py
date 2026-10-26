import sys, json
sys.path.insert(0, "/workspace/src")
import numpy as np, pandas as pd
from ksp import index, splits, canary
from scipy.optimize import minimize
A = pd.read_parquet("/workspace/cache/r2/assign_r2_fusion_r3__20261026-191901__s0__calibrated.parquet")
S = index.load_sites(); K = canary.kinases_with_canaries()
A["fold"] = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values[A.site_idx.values]
A["res"] = np.where(S.res.values[A.site_idx.values] == "Y", "Y", "ST")
lg = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / np.clip(1 - p, 1e-6, None))
def platt_fit(x, y):
    f = lambda t: np.sum(np.logaddexp(0, t[0] * x + t[1]) - y * (t[0] * x + t[1]))
    return minimize(f, [1.0, 0.0], method="L-BFGS-B").x
bins = [(0, .2), (.2, .5), (.5, .8), (.8, 1.01)]
def rel(D, col):
    L = D[D.labelled]; out = []
    for lo, hi in bins:
        g = L[(L[col] >= lo) & (L[col] < hi)]
        out.append((f"[{lo},{hi})", len(g), round(g[col].mean(), 3) if len(g) else None, round(g.top_is_true.mean(), 3) if len(g) else None))
    return out
res = {"raw_by_res": {r: rel(A[A.res == r], "p_top") for r in ["ST", "Y"]}}
A["p_rs"] = np.nan
for f in range(5):
    for r in ["ST", "Y"]:
        tr = A[(A.fold != f) & A.labelled & (A.res == r)]
        a, b = platt_fit(lg(tr.p_top.values), tr.top_is_true.values.astype(float))
        m = (A.fold.values == f) & (A.res.values == r)
        A.loc[m, "p_rs"] = 1 / (1 + np.exp(-(a * lg(A.p_top.values[m]) + b)))
        res.setdefault("params", []).append((f, r, round(a, 3), round(b, 3)))
res["platt_by_res_rel"] = rel(A, "p_rs")
res["platt_by_res_rel_ST"] = rel(A[A.res == "ST"], "p_rs"); res["platt_by_res_rel_Y"] = rel(A[A.res == "Y"], "p_rs")
res["coverage_platt_by_res"] = {t: int((A.p_rs >= t).sum()) for t in (0.5, 0.8)} | {f"{t}_unlab": int(((A.p_rs >= t) & ~A.labelled).sum()) for t in (0.5, 0.8)}
res["coverage_raw_by_res"] = {r: {t: int(((A.p_top >= t) & (A.res == r)).sum()) for t in (0.5, 0.8)} for r in ["ST", "Y"]}
json.dump(res, open("out/calib_check2.json", "w"), indent=1, default=float); print(json.dumps(res, default=float, indent=0))
