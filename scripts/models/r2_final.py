"""R2 final model(s) trained on ALL DEV labels (DESIGN §8.6) and exhaustive scoring.

Modes
  --target test      : TEST candidate set = U_ms sites on TEST proteins x every REAL kinase of matching class
                       (= ksp.eval_test.test_candidate_pairs()); writes cache/runs/<run_id>/test_predictions.parquet
                       (kinase_idx, site_idx, score[, score_site]) - mean logit over --seeds.
  --target t3        : per LKO kinase fold j: labels of fold-j kinases removed from training AND from the label-derived
                       features (xfit tag t3lko<j>, rows = DEV U_obs sites, eval rows = TEST sites, A1 inner-fold
                       fits); scores fold-j kinases on their TEST candidates -> test_t3_predictions.parquet.
  --target uobs      : every U_obs site (DEV + TEST proteins) x every real kinase of matching class (deployment-style
                       library, e.g. the benchmarKIN clean variant with --set clean='"bench"') -> uobs_predictions.parquet.
Label-derived features for TEST-protein rows come from the inner-fold fits (A1), exactly as CV eval rows; DEV rows are
training rows. No TEST label is read anywhere (ksp.labels.load_dev only).
The manifest records train_wall_s and inference_wall_s (G4 runbook).
Usage: r2_final.py --name r2_final_test --target test --seeds 0 1 2 3 4 --set num_leaves=127 ...
"""
import argparse, json, time
import numpy as np, pandas as pd
from ksp import index, manifest, splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train, xfit

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True)
ap.add_argument("--target", choices=["test", "t3", "uobs"], required=True)
ap.add_argument("--seeds", nargs="+", type=int, default=[0])
ap.add_argument("--lko", nargs="+", type=int, default=[0, 1, 2, 3, 4])
ap.add_argument("--set", nargs="*", default=[])
a = ap.parse_args()
cfg = train.default_config()
for kv in a.set:
    k, v = kv.split("=", 1)
    try:
        cfg[k] = json.loads(v)
    except json.JSONDecodeError:
        cfg[k] = v
S = index.load_sites()
dev = S.acc.isin(set(splits.load_cv_folds().index)).values
test = S.acc.isin(splits.load_test_proteins()).values
FAMS = ["f2a_sep" if cfg.get("f2a_family") == "f2a_sep" else "f2a", "f3nb",
        "f6a_bench" if cfg["clean"] == "bench" else ("f6a_clean" if cfg["clean"] else "f6a")]
REAL = assemble.N_REAL


def score_sites(target):
    if target in ("test", "t3"):
        return np.where(test & S.in_U_ms.values)[0]
    return np.arange(len(S))                                  # uobs


def fit_predict(fd, A, seed, sites, kinases):
    c = dict(cfg, seed=seed)
    lg = {}
    t0 = time.time()
    bst, cols, stats = train.train_fold(fd, A, c, lg)
    t_train = time.time() - t0
    fd.eval_sites = sites
    t1 = time.time()
    P = train.predict_fold(fd, A, bst, c, kinases=kinases)
    return P, t_train, time.time() - t1, lg


with manifest.RunManifest(a.name, config=cfg | {"target": a.target, "seeds": a.seeds}, seed=a.seeds[0],
                          split="final_dev", notes=f"final DEV model, target={a.target}, seeds={a.seeds}") as run:
    t_all = time.time()
    sites = score_sites(a.target)
    tw, iw, logs = 0.0, 0.0, {}
    if a.target in ("test", "uobs"):
        fd = sample.FoldData(None, cfg["variant"])
        A = assemble.Assembler(cfg["groups"], None, xfit_tag="final", clean=cfg["clean"], comp_know=cfg["comp_know"],
                               f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
        acc = None
        for sd in a.seeds:
            P, t_tr, t_inf, lg = fit_predict(fd, A, sd, sites, list(range(REAL)))
            tw += t_tr; iw += t_inf; logs[sd] = lg
            if acc is None:
                keys = P[["kinase_idx", "site_idx"]].copy(); acc = P.score.values.astype(np.float64)
            else:
                if not (np.array_equal(keys.kinase_idx.values, P.kinase_idx.values) and np.array_equal(keys.site_idx.values, P.site_idx.values)):
                    raise RuntimeError("seed predictions not aligned")
                acc += P.score.values
            print(f"seed {sd}: {lg} train {t_tr:.0f}s infer {t_inf:.0f}s rows {len(P)}", flush=True)
        out = keys.assign(score=(acc / len(a.seeds)).astype(np.float32))
    else:
        kf = splits.load_kinase_folds()
        kpos = {x: i for i, x in enumerate(index.load_kinases().uniprot_acc)}
        tr_rows, ev_rows = np.where(dev)[0], np.where(test & ~dev)[0]
        parts = []
        for j in a.lko:
            held = set(kf.kinase_acc[kf.lko_fold == j])
            for fam in FAMS:
                xfit.build(fam, None, tag=f"t3lko{j}", train_site_idx=tr_rows, eval_site_idx=ev_rows, exclude_kinases=held)
            fd = sample.FoldData(None, cfg["variant"], exclude_kinases=held)
            A = assemble.Assembler(cfg["groups"], None, xfit_tag=f"t3lko{j}", clean=cfg["clean"], comp_know=cfg["comp_know"],
                                   f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
            hk = sorted(kpos[x] for x in held if x in kpos)
            acc = None
            for sd in a.seeds:
                P, t_tr, t_inf, lg = fit_predict(fd, A, sd, sites, hk)
                tw += t_tr; iw += t_inf; logs[f"{j}_{sd}"] = lg
                if acc is None:
                    keys = P[["kinase_idx", "site_idx"]].copy(); acc = P.score.values.astype(np.float64)
                else:
                    acc += P.score.values
            parts.append(keys.assign(score=(acc / len(a.seeds)).astype(np.float32)))
            print(f"T3 lko {j}: {len(hk)} kinases, {len(parts[-1])} rows", flush=True)
        out = pd.concat(parts, ignore_index=True)
    if out.kinase_idx.max() >= REAL or not np.isfinite(out.score).all():
        raise RuntimeError("canary rows or non-finite scores in final predictions")
    if out.duplicated(["kinase_idx", "site_idx"]).any():
        raise RuntimeError("duplicate rows")
    fname = {"test": "test_predictions", "t3": "test_t3_predictions", "uobs": "uobs_predictions"}[a.target]
    run.save_predictions(out, name=fname)
    run.set("train_wall_s", round(tw, 1)); run.set("inference_wall_s", round(iw, 1))
    run.set("total_wall_s", round(time.time() - t_all, 1)); run.set("seed_logs", logs)
    print("RUN", run.run_id, fname, len(out), f"train {tw:.0f}s infer {iw:.0f}s", flush=True)
