"""R2 partial pooling sweep over lambda (DESIGN §8.3) on top of a finished 5-fold R2 config.
For each CV fold: rebuild the same training rows (same config/seed), inner 3-fold out-of-sample R2 offsets on training
rows (R2 refitted with the outer run's n_rounds), fit per-kinase family-shrunk ridge residuals, apply to the outer OOF
logits with w_k = n_k/(n_k+lambda). Evaluate pooled CV per lambda; report M1/M3 by training-count bin.
Usage: r2_partial_pool.py --run <pooled run_id> [--lams 0 3 10 30 100 300 inf]"""
import argparse, json, time
import lightgbm as lgb
import numpy as np, pandas as pd
from ksp import eval as kev, manifest, splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train, pool

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--lams", nargs="+", default=["0", "3", "10", "30", "100", "300", "inf"])
ap.add_argument("--mu", type=float, default=10.0)
a = ap.parse_args()
lams = [float(x) for x in a.lams]
m = json.load(open(ROOT / "results" / "runs" / a.run / "manifest.json"))
cfg = m["config"]
fold_runs = m["fold_runs"]
P_main = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
outs = {lam: [] for lam in lams}
t0 = time.time()
with manifest.RunManifest("r2_partial_pool", config={"base_run": a.run, "lams": lams, "mu": a.mu, "features": pool.POOL_FEATURES},
                          seed=cfg["seed"], split="cv", notes="partial pooling sweep (§8.3)") as run:
    for f in range(5):
        fm = json.load(open(ROOT / "results" / "runs" / fold_runs[f] / "manifest.json"))
        n_rounds = list(fm["fold_logs"].values())[0]["best_iteration"]
        fd = sample.FoldData(f, cfg["variant"])
        A = assemble.Assembler(cfg["groups"], f, clean=cfg["clean"], comp_know=cfg["comp_know"],
                               f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
        d, X, cols, stats = train.build_train(fd, A, cfg)
        cat = [c for c in A.cat_cols if c in cols]
        inner = splits.inner_folds(fd.S.site_id.values[d.site_idx.values], 3, cfg["seed"] + 11)
        off = np.empty(len(d))
        params = train.lgb_params(cfg)
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
        # outer OOF rows of this fold
        Pf = P_main[np.isin(P_main.site_idx.values, fd.eval_sites)].reset_index(drop=True)
        Xe = np.full((len(Pf), len(pidx)), np.nan, np.float32)
        for k, ix in Pf.groupby("kinase_idx").indices.items():
            if k in coefs:
                Xe[ix] = A.kinase_block(int(k), Pf.site_idx.values[ix])[:, pidx]
        for lam in lams:
            sc = pool.apply_residuals(Pf.kinase_idx.values, Xe, Pf.score.values, coefs, std, ntr, lam)
            outs[lam].append(Pf.assign(score=sc.astype(np.float32)))
        print(f"fold {f} done {time.time() - t0:.0f}s", flush=True)
    data = kev.prepare("cv", fold=None)
    rows = []
    for lam in lams:
        P = pd.concat(outs[lam], ignore_index=True)
        res = kev.evaluate(P, data, missing="error", keep_items=False)
        pk = res.per_kinase
        pk = pk[pk.in_K_primary]
        r = {"lambda": lam, **{k: res.endpoints[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany"]}}
        for bname, g in pk.groupby("train_bin"):
            r[f"M1_bin_{bname}"] = g.auroc_hard.mean()
            r[f"M3_bin_{bname}"] = g.top1_recall.mean()
            r[f"n_bin_{bname}"] = len(g)
        rows.append(r)
        print(r, flush=True)
    R = pd.DataFrame(rows)
    run.dir.mkdir(parents=True, exist_ok=True)
    R.to_csv(run.dir / "lambda_sweep.tsv", sep="\t", index=False)
    R.to_csv(ROOT / "results" / "models" / "r2" / "partial_pool_lambda_sweep.tsv", sep="\t", index=False)
    run.save_metrics({"lambda_sweep": rows})
    best = R.loc[R.M1_auroc_hard.idxmax()]
    train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)", "config_name": "r2_partial_pool_sweep",
                    "config_json": json.dumps({"base_run": a.run, "lams": lams, "mu": a.mu}),
                    "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(float(best.M1_auroc_hard), 5),
                    "notes": f"lambda sweep ({len(lams)} values counted as one registry row); best lambda={best['lambda']}"})
