"""Per-Manning-group LightGBM models (frozen final R2 hyper-parameters), same training rows / features as the global
model (seed 0). Groups with < 200 training positives and canaries (own group each, A25) fall back to the global OOF
score. alpha = 0.5 weights are renormalised within group (FoldData.sample weights are already n_k^-alpha; mean 1 per
group). Writes cache/regime/pergroup/fold<f>.parquet [kinase_idx, site_idx, score_grp].
Usage: regime_pergroup.py --folds 0 [--threads 6]"""
import argparse, json, os, time
ap = argparse.ArgumentParser()
ap.add_argument("--folds", nargs="+", type=int, default=[0])
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--min-pos", type=int, default=200)
a = ap.parse_args()
os.environ["R2_THREADS"] = str(a.threads)
import numpy as np, pandas as pd
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train
from ksp.models.regime import run as R

OUT = ROOT / "cache" / "regime" / "pergroup"
OUT.mkdir(parents=True, exist_ok=True)
cfg = R.final_config(0)
for f in a.folds:
    t0 = time.time()
    fd = sample.FoldData(f, "L2")
    A = assemble.Assembler(cfg["groups"], f, clean=cfg["clean"], comp_know=cfg["comp_know"],
                           f2a_family=cfg["f2a_family"], desc_v=cfg["desc_v"])
    grp = A.K.manning_group.fillna("(none)").values
    real = np.arange(len(A.K)) < assemble.N_REAL
    d_all = fd.sample(cfg["ratio"], cfg["hard"], cfg["alpha"], cfg["seed"])
    X_all = A.rows(d_all.kinase_idx.values, d_all.site_idx.values)
    g_row = grp[d_all.kinase_idx.values]
    outs, logs = [], []
    for g in sorted(set(grp[real])):
        m = (g_row == g) & real[d_all.kinase_idx.values]
        npos = int(d_all.y.values[m].sum())
        if npos < a.min_pos:
            logs.append({"group": g, "n_pos": npos, "model": "fallback_global"})
            continue

        class Sub:            # FoldData-like view restricted to the group's rows for train.train_fold
            pass
        d = d_all[m].reset_index(drop=True)
        d["w"] = d.w * len(d) / d.w.sum()
        Xg = X_all[m]

        def build(fd_, A_, cfg_, d=d, Xg=Xg):
            return d, Xg, list(A.cols), None
        orig = train.build_train
        train.build_train = build
        lg = {}
        try:
            bst, cols, _ = train.train_fold(fd, A, cfg, lg)
        finally:
            train.build_train = orig
        ks = [k for k in range(assemble.N_REAL) if grp[k] == g]
        P = train.predict_fold(fd, A, bst, cfg, kinases=ks)
        outs.append(P.rename(columns={"score": "score_grp"}))
        logs.append({"group": g, "n_pos": npos, "model": "own", "n_kinases": len(ks), **lg})
        print(f"fold {f} group {g}: {lg}", flush=True)
    pd.concat(outs, ignore_index=True).to_parquet(OUT / f"fold{f}.parquet")
    pd.DataFrame(logs).to_csv(OUT / f"fold{f}_models.tsv", sep="\t", index=False)
    (OUT / f"fold{f}.json").write_text(json.dumps({"fold": f, "t_s": round(time.time() - t0, 1), "threads": a.threads,
                                                   "min_pos": a.min_pos}))
    print(f"fold {f} done {time.time() - t0:.0f}s", flush=True)
