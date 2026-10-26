import sys, json
sys.path.insert(0, "/workspace/src")
import numpy as np, pandas as pd
from ksp import index, splits
from sklearn.isotonic import IsotonicRegression
A = pd.read_parquet("/workspace/cache/r2/assign_r2_fusion_r3__20261026-191901__s0__calibrated.parquet")
S = index.load_sites()
A["fold"] = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values[A.site_idx.values]
A["res"] = np.where(S.res.values[A.site_idx.values] == "Y", "Y", "ST")
A["p_iso_rs"] = np.nan
for f in range(5):
    for r in ["ST", "Y"]:
        tr = A[(A.fold != f) & A.labelled & (A.res == r)]
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(tr.p_top.values, tr.top_is_true.values.astype(float))
        m = (A.fold.values == f) & (A.res.values == r)
        A.loc[m, "p_iso_rs"] = iso.predict(A.p_top.values[m])
out = {}
for r in ["ST", "Y", "all"]:
    D = A if r == "all" else A[A.res == r]; L = D[D.labelled]; rows = []
    for lo, hi in [(0, .2), (.2, .5), (.5, .8), (.8, 1.01)]:
        g = L[(L.p_iso_rs >= lo) & (L.p_iso_rs < hi)]
        rows.append((f"[{lo},{hi})", len(g), round(g.p_iso_rs.mean(), 3) if len(g) else None, round(g.top_is_true.mean(), 3) if len(g) else None))
    out[r] = rows
out["coverage"] = {t: int((A.p_iso_rs >= t).sum()) for t in (0.5, 0.8)} | {f"{t}_unlab": int(((A.p_iso_rs >= t) & ~A.labelled).sum()) for t in (0.5, 0.8)}
out["coverage_raw_by_res_unlab"] = {r: {t: int(((A.p_top >= t) & (A.res == r) & ~A.labelled).sum()) for t in (0.5, 0.8)} for r in ["ST", "Y"]}
print(json.dumps(out, indent=0)); json.dump(out, open("out/calib_check3.json", "w"), indent=1)
