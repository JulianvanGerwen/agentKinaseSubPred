import sys, json
sys.path.insert(0, "/workspace/src")
import numpy as np, pandas as pd
from ksp import labels, index, splits, canary
from scipy.optimize import minimize
from sklearn.isotonic import IsotonicRegression
A = pd.read_parquet("/workspace/cache/r2/assign_r2_fusion_r3__20261026-191901__s0__calibrated.parquet")
S = index.load_sites(); K = canary.kinases_with_canaries()
fold = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values
A["fold"] = fold[A.site_idx.values]
A["top_acc"] = K.uniprot_acc.values[A.top_kinase_idx.values]
A["fam"] = K.manning_family.values[A.top_kinase_idx.values]
A["gene"] = K.gene.values[A.top_kinase_idx.values]
dev = labels.load_dev("all")
sid = S.site_id.values
lab_any = dev.groupby("site_id").kinase_acc.agg(set)
par = labels.paralog_sets()
L = A[A.labelled].copy()
L["any_src"] = [ta in lab_any.get(sid[s], set()) for s, ta in zip(L.site_idx, L.top_acc)]
L["par_of_true"] = [any(ta in par.get(k, set()) for k in lab_any.get(sid[s], set())) for s, ta in zip(L.site_idx, L.top_acc)]
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(c - h, 3), round(c + h, 3)
out = {}
bins = [(0, .2), (.2, .5), (.5, .8), (.8, 1.01)]
rows = []
for lo, hi in bins:
    g = L[(L.p_top >= lo) & (L.p_top < hi)]
    rows.append({"bin": f"[{lo},{hi})", "n": len(g), "mean_p": g.p_top.mean(), "obs_true": g.top_is_true.mean(),
                 "wilson95": wilson(g.top_is_true.sum(), len(g)), "obs_any_source_label": g.any_src.mean(),
                 "obs_true_or_paralog": (g.top_is_true | g.par_of_true).mean()})
out["reliability"] = rows
top = L[L.p_top >= 0.8]
out["top_tier_genes"] = top.gene.value_counts().head(12).to_dict()
out["top_tier_true_by_gene"] = top.groupby("gene").top_is_true.agg(["size", "mean"]).sort_values("size", ascending=False).head(12).round(2).to_dict("index")
# cross-fitted recalibration of p_top on labelled sites (fit folds != f, apply fold f); applied to ALL sites for coverage
lg = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / np.clip(1 - p, 1e-6, None))
def platt_fit(x, y):
    f = lambda t: np.sum(np.logaddexp(0, t[0] * x + t[1]) - y * (t[0] * x + t[1]))
    return minimize(f, [1.0, 0.0], method="L-BFGS-B").x
A["p_platt"] = np.nan; A["p_iso"] = np.nan
for f in range(5):
    tr = L[L.fold != f]
    a, b = platt_fit(lg(tr.p_top.values), tr.top_is_true.values.astype(float))
    m = A.fold.values == f
    A.loc[m, "p_platt"] = 1 / (1 + np.exp(-(a * lg(A.p_top.values[m]) + b)))
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(tr.p_top.values, tr.top_is_true.values.astype(float))
    A.loc[m, "p_iso"] = iso.predict(A.p_top.values[m])
    out.setdefault("platt_params", {})[f] = [round(a, 3), round(b, 3)]
Lr = A[A.labelled]
for col in ["p_top", "p_platt", "p_iso"]:
    rr = []
    for lo, hi in bins:
        g = Lr[(Lr[col] >= lo) & (Lr[col] < hi)]
        rr.append({"bin": f"[{lo},{hi})", "n": len(g), "mean_p": round(g[col].mean(), 3) if len(g) else None,
                   "obs": round(g.top_is_true.mean(), 3) if len(g) else None})
    out[f"rel_{col}"] = rr
    out[f"coverage_{col}"] = {t: int((A[col] >= t).sum()) for t in (0.5, 0.8)} | \
        {f"{t}_unlabelled": int(((A[col] >= t) & ~A.labelled).sum()) for t in (0.5, 0.8)}
json.dump(out, open("out/calib_check.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
