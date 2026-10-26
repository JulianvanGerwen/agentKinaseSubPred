"""Partial pooling, form 2 (DESIGN §8.3 'global model + per-family residual boosting with shrinkage'): the frozen
global R2 logit (seed 0) + w_k * g_fam(k)(x), w_k = n_k / (n_k + lambda), where g_fam is a small LightGBM residual
booster per Manning family fitted with init_score = OUT-OF-SAMPLE R2 offset on the family's training rows (inner
3-fold cluster split, R2 refitted at the outer run's round count - identical offsets to regime_partial_pool.py).
Booster (pre-declared): 15 leaves, lr 0.05, min_data_in_leaf 50, feature_fraction 0.5, bagging 0.8/1, lambda_l2 1,
early stopping 50 (max 500) on an inner 20% cluster split of the family rows; families with < 20 training positives
get no residual. Same lambda grid {0,3,10,30,100,300,inf}; selection on fold-0 M1.
Writes cache/regime/resboost/fold<f>.parquet [kinase_idx, site_idx, score, score_lam<..>].
Usage: regime_resboost.py --fold 0 [--threads 6]"""
import argparse, json, os, time
ap = argparse.ArgumentParser()
ap.add_argument("--fold", type=int, required=True)
ap.add_argument("--threads", type=int, default=6)
a = ap.parse_args()
os.environ["R2_THREADS"] = str(a.threads)
import lightgbm as lgb
import numpy as np, pandas as pd
from ksp import splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train
from ksp.models.regime import perkinase as PK, run as R

LAMS = [0.0, 3.0, 10.0, 30.0, 100.0, 300.0, float("inf")]
BP = dict(objective="binary", metric="None", learning_rate=0.05, num_leaves=15, min_data_in_leaf=50,
          feature_fraction=0.5, bagging_fraction=0.8, bagging_freq=1, lambda_l2=1.0, deterministic=True,
          force_col_wise=True, verbose=-1)
OUT = ROOT / "cache" / "regime" / "resboost"
OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
m = json.loads((ROOT / "results" / "runs" / R.BASE_RUN / "manifest.json").read_text())
cfg = m["config"]
f = a.fold
fm = json.loads((ROOT / "results" / "runs" / m["fold_runs"][f] / "manifest.json").read_text())
n_rounds = list(fm["fold_logs"].values())[0]["best_iteration"]
P_main = pd.read_parquet(ROOT / fm["outputs"]["predictions"]["path"])
fd = sample.FoldData(f, cfg["variant"])
A = assemble.Assembler(cfg["groups"], f, clean=cfg["clean"], comp_know=cfg["comp_know"],
                       f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
d, X, cols, stats = train.build_train(fd, A, cfg)
cat = [c for c in A.cat_cols if c in cols]
inner = splits.inner_folds(fd.S.site_id.values[d.site_idx.values], 3, cfg["seed"] + 11)
off = np.empty(len(d))
params = train.lgb_params(cfg)
params["num_threads"] = a.threads
for i in range(3):
    tr = inner != i
    b = lgb.train(params, lgb.Dataset(X[tr], d.y.values[tr], weight=d.w.values[tr], feature_name=cols,
                                      categorical_feature=cat), num_boost_round=n_rounds)
    off[~tr] = b.predict(X[~tr], raw_score=True)
t_off = time.time() - t0
fam = A.K.manning_family.fillna("(none)").values
fam_row = fam[d.kinase_idx.values]
va = splits.train_valid_split_inner(fd.S.site_id.values[d.site_idx.values], 0.2, cfg["seed"] + 3)
boosters, logs = {}, []
for g in np.unique(fam_row):
    ix = np.flatnonzero(fam_row == g)
    y = d.y.values[ix].astype(int)
    if y.sum() < 20:
        continue
    w = d.w.values[ix]
    v = va[ix]
    p = {**BP, "seed": cfg["seed"], "num_threads": a.threads}
    best = 100
    if y[v].sum() >= 3 and y[~v].sum() >= 3:
        dtr = lgb.Dataset(X[ix][~v], y[~v], weight=w[~v], init_score=off[ix][~v], feature_name=cols,
                          categorical_feature=cat, free_raw_data=False)
        dva = lgb.Dataset(X[ix][v], y[v], weight=w[v], init_score=off[ix][v], reference=dtr, feature_name=cols,
                          categorical_feature=cat, free_raw_data=False)
        yv, ov = y[v], off[ix][v]
        bb = lgb.train(p, dtr, 500, valid_sets=[dva], feval=lambda pr, ds: [("auc", PK.auc(yv, pr), True)],
                       callbacks=[lgb.early_stopping(50, verbose=False)])
        best = bb.best_iteration or 100
    if best < 1:
        continue
    boosters[g] = lgb.train(p, lgb.Dataset(X[ix], y, weight=w, init_score=off[ix], feature_name=cols,
                                           categorical_feature=cat), best)
    logs.append({"family": g, "n_pos": int(y.sum()), "rounds": int(best)})
Pf = P_main.reset_index(drop=True)
res_s = np.zeros(len(Pf), np.float32)
for k, ix in Pf.groupby("kinase_idx").indices.items():
    g = fam[k]
    if g in boosters:
        res_s[ix] = boosters[g].predict(A.kinase_block(int(k), Pf.site_idx.values[ix]), raw_score=True)
ntr = stats.n_pos.to_dict()
nk = np.array([ntr.get(k, 0) for k in Pf.kinase_idx.values], float)
out = Pf[["kinase_idx", "site_idx", "score"]].copy()
for lam in LAMS:
    with np.errstate(invalid="ignore"):
        wk = np.where(nk > 0, nk / (nk + lam), 0.0) if np.isfinite(lam) else np.zeros_like(nk)   # n_k = 0 -> w = 0
    out[f"score_lam{lam:g}"] = (Pf.score.values + wk * res_s).astype(np.float32)
out.to_parquet(OUT / f"fold{f}.parquet")
pd.DataFrame(logs).to_csv(OUT / f"fold{f}_models.tsv", sep="\t", index=False)
meta = {"fold": f, "n_rounds_global": n_rounds, "n_family_boosters": len(boosters), "t_offsets_s": round(t_off, 1),
        "t_s": round(time.time() - t0, 1), "threads": a.threads}
(OUT / f"fold{f}.json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta), flush=True)
