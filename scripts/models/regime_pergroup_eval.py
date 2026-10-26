"""Combine per-Manning-group predictions (own model where the group had one) with the global R2 OOF fallback and
evaluate pooled 5-fold CV (registered) - score_site = same score (one LightGBM scale per group; M2 also with qmap)."""
import json
import numpy as np, pandas as pd
from ksp import splits, index
from ksp.index import ROOT
from ksp.models.regime import perkinase as PK, run as R

base = json.loads((ROOT / "results" / "runs" / R.BASE_RUN / "manifest.json").read_text())
G = pd.read_parquet(ROOT / base["outputs"]["predictions"]["path"])
S = index.load_sites()
G["fold"] = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values[G.site_idx.values]
Pg = pd.concat([pd.read_parquet(ROOT / "cache" / "regime" / "pergroup" / f"fold{f}.parquet") for f in range(5)])
M = G.merge(Pg, on=["kinase_idx", "site_idx"], how="left")
has = M.score_grp.notna().values
sc = np.where(has, M.score_grp.values, M.score.values)
ss = sc.copy()
hidx = np.flatnonzero(has)
for (k, f), ix in M[has].groupby(["kinase_idx", "fold"]).indices.items():
    rows = hidx[ix]
    ss[rows] = PK.qmap(sc[rows], M.score.values[rows])
for m2, s_site in (("raw", sc), ("qmap", ss)):
    P = M[["kinase_idx", "site_idx"]].assign(score=sc.astype(np.float64), score_site=s_site)
    R.pool(f"regime_pergroup_{m2}_cv5", [], {"regime": "per-Manning-group LightGBM (final R2 params)", "min_pos": 200,
                                            "fallback": "global " + R.BASE_RUN, "m2_score": m2, "seed": 0},
           notes=f"per-group models, M2 score {m2}", preds=R.to_eval(P), register_row=(m2 == "raw"))
