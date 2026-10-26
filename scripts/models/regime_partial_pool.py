"""Partial pooling (DESIGN §8.3, ALG-3) on top of the frozen final R2 config, seed 0 (base = R2's own 5-fold run;
not retrained). Reuses ksp.models.r2.pool (per-kinase ridge-logistic residual on 8 fixed features, shrunk towards the
Manning-family residual, mu = 10) and the logic of scripts/models/r2_partial_pool.py (reviewed: sound - training-row
offsets are out-of-sample via an inner 3-fold cluster split with R2 refitted at the outer run's round count; eval rows
use the outer OOF logit), split into per-fold jobs.
final logit = R2 logit + w_k * r_k(x), w_k = n_k / (n_k + lambda). Pre-declared grid lambda in {0,3,10,30,100,300,inf}.
Writes cache/regime/ppool/fold<f>.parquet [kinase_idx, site_idx, score, score_lam<..>].
Usage: regime_partial_pool.py --fold 0 [--threads 6]"""
import argparse, json, os, time
ap = argparse.ArgumentParser()
ap.add_argument("--fold", type=int, required=True)
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--mu", type=float, default=10.0)
a = ap.parse_args()
os.environ["R2_THREADS"] = str(a.threads)
import lightgbm as lgb
import numpy as np, pandas as pd
from ksp import splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, pool, sample, train
from ksp.models.regime import run as R

LAMS = [0.0, 3.0, 10.0, 30.0, 100.0, 300.0, float("inf")]
OUT = ROOT / "cache" / "regime" / "ppool"
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
    ds = lgb.Dataset(X[tr], d.y.values[tr], weight=d.w.values[tr], feature_name=cols, categorical_feature=cat)
    b = lgb.train(params, ds, num_boost_round=n_rounds)
    off[~tr] = b.predict(X[~tr], raw_score=True)
pidx = [cols.index(c) for c in pool.POOL_FEATURES]
trd = d.assign(offset=off)
fam = dict(enumerate(A.K.manning_family.fillna("(none)").values))
coefs, std, _ = pool.fit_residuals(trd, X[:, pidx], fam, a.mu)
ntr = stats.n_pos.to_dict()
Pf = P_main[np.isin(P_main.site_idx.values, fd.eval_sites)].reset_index(drop=True)
assert len(Pf) == len(P_main), "fold predictions contain other folds' sites"
Xe = np.full((len(Pf), len(pidx)), np.nan, np.float32)
for k, ix in Pf.groupby("kinase_idx").indices.items():
    if k in coefs:
        Xe[ix] = A.kinase_block(int(k), Pf.site_idx.values[ix])[:, pidx]
out = Pf[["kinase_idx", "site_idx", "score"]].copy()
for lam in LAMS:
    out[f"score_lam{lam:g}"] = pool.apply_residuals(Pf.kinase_idx.values, Xe, Pf.score.values, coefs, std, ntr,
                                                    lam).astype(np.float32)
out.to_parquet(OUT / f"fold{f}.parquet")
meta = {"fold": f, "n_rounds": n_rounds, "mu": a.mu, "lams": LAMS, "n_kinase_models": len(coefs),
        "features": pool.POOL_FEATURES, "t_s": round(time.time() - t0, 1), "threads": a.threads}
(OUT / f"fold{f}.json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta), flush=True)
