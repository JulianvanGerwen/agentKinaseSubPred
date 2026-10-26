"""Per-kinase models (LightGBM and L2-logistic) for kinases with >= N training positives in the outer-training fold;
same training rows / features as the frozen global R2 config (seed 0). Writes per-fold predictions for the kinases
that got their own model: cache/regime/perkinase/fold<f>.parquet [kinase_idx, site_idx, score_lgb, score_lr].
Usage: regime_perkinase.py --folds 0 [--threads 4]"""
import argparse, json, os, time
ap = argparse.ArgumentParser()
ap.add_argument("--folds", nargs="+", type=int, default=[0])
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--n-min", type=int, default=20)
a = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", str(a.threads))
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from ksp import splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample
from ksp.models.regime import perkinase as PK, run as R

OUT = ROOT / "cache" / "regime" / "perkinase"
OUT.mkdir(parents=True, exist_ok=True)
cfg = R.final_config(0)
for f in a.folds:
    t0 = time.time()
    fd = sample.FoldData(f, "L2")
    A = assemble.Assembler(cfg["groups"], f, clean=cfg["clean"], comp_know=cfg["comp_know"],
                           f2a_family=cfg["f2a_family"], desc_v=cfg["desc_v"])
    d = fd.sample(cfg["ratio"], cfg["hard"], cfg["alpha"], cfg["seed"])
    X = A.rows(d.kinase_idx.values, d.site_idx.values)
    t_gather = time.time() - t0
    cols = list(A.cols)
    drop = set(PK.const_cols(cols))
    keep = [j for j in range(len(cols)) if j not in drop]
    kcols = [cols[j] for j in keep]
    cat_idx = [i for i, c in enumerate(kcols) if c in A.cat_cols]
    va_all = splits.train_valid_split_inner(fd.S.site_id.values[d.site_idx.values], 0.2, cfg["seed"])
    npos = d[d.y == 1].groupby("kinase_idx").size()
    ks = [int(k) for k in npos.index if npos[k] >= a.n_min]
    res = fd.S.res.values[fd.eval_sites]
    site_sets = {c: fd.eval_sites[sample.residue_ok(res, c)] for c in ("ST", "Y", "STY")}
    groups = d.groupby("kinase_idx").indices
    logs = []

    def one(k):
        ix = groups[k]
        Xk, yk, vk = X[ix][:, keep], d.y.values[ix].astype(int), va_all[ix]
        b, it = PK.fit_lgb_kinase(Xk, yk, vk, cfg["seed"], 1, cat_idx)
        lr, C = PK.fit_lr_kinase(Xk, yk, vk, cfg["seed"])
        s = site_sets[A.K.residue_class.values[k]]
        Xe = A.kinase_block(k, s)[:, keep]
        return k, s, b.predict(Xe, raw_score=True), lr(Xe), {"kinase_idx": k, "n_pos": int(yk.sum()), "n_rows": len(ix),
                                                             "lgb_rounds": it, "lr_C": C}
    outs = []
    with ThreadPoolExecutor(a.threads) as ex:
        for k, s, sl, sr, lg in ex.map(one, ks):
            outs.append(pd.DataFrame({"kinase_idx": np.int32(k), "site_idx": s.astype(np.int32),
                                      "score_lgb": sl.astype(np.float32), "score_lr": sr.astype(np.float32)}))
            logs.append(lg)
    P = pd.concat(outs, ignore_index=True)
    P.to_parquet(OUT / f"fold{f}.parquet")
    L = pd.DataFrame(logs)
    L.to_csv(OUT / f"fold{f}_models.tsv", sep="\t", index=False)
    meta = {"fold": f, "n_min": a.n_min, "n_models": len(ks), "n_real": int(sum(k < assemble.N_REAL for k in ks)),
            "t_gather_s": round(t_gather, 1), "t_total_s": round(time.time() - t0, 1), "threads": a.threads}
    (OUT / f"fold{f}.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta), flush=True)
