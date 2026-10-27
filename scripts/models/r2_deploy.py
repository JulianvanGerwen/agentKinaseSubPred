"""G5 deployment model for R2 (DESIGN §8.6): the FROZEN final configuration retrained on DEV+TEST labels
(ksp.deploy.load_deploy_labels; refuses until the TEST batch is COMPLETED), 5 seeds, scoring all U_obs x real kinases
of matching class. Canaries: the frozen config trains with canary pseudo-kinases (code path), so with_canaries=True is
used for training and xfit, and canary outputs are never written (only the 503 real kinases are scored).

Stages
  --stage xfit  --families f2a_sep f3nb f6a : cross-fitted stores tag 'deploy' (rows = every U_obs site; inner folds
                                              over substrate clusters of DEV+TEST proteins; no eval rows)
  --stage train                             : 5 seeds -> cache/release/<run_id>/r2_uobs_scores.parquet (written in
                                              kinase chunks) + manifest with train/inference wall-clock
--smoke: DEV labels (ksp.labels.load_dev) instead of deployment labels, toy kinase subset, existing 'final' stores;
         for testing the code path before GO (never touches TEST labels).
"""
import argparse, json, time
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from ksp import index, labels, manifest
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train, xfit

V3 = dict(num_leaves=127, learning_rate=0.03, f2a_family="f2a_sep", desc_v=2, es_real_only=True)
ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["xfit", "train"], required=True)
ap.add_argument("--families", nargs="*", default=["f2a_sep", "f3nb", "f6a"])
ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
ap.add_argument("--name", default="r2_deploy")
ap.add_argument("--smoke", action="store_true")
ap.add_argument("--smoke-kinases", type=int, default=6)
a = ap.parse_args()

if a.smoke:
    loader = lambda variant, with_canaries=False: labels.load_dev(variant, with_canaries=with_canaries)
    tag = "final"
else:
    from ksp.deploy import load_deploy_labels
    loader = lambda variant, with_canaries=False: load_deploy_labels(variant, with_canaries=with_canaries)
    tag = "deploy"
cfg = train.default_config() | V3
S = index.load_sites()

if a.stage == "xfit":
    if a.smoke:
        raise SystemExit("smoke uses the existing 'final' stores")
    tl = loader("L2", with_canaries=True)[["kinase_acc", "site_id"]].drop_duplicates().reset_index(drop=True)
    for fam in a.families:
        t = time.time()
        p = xfit.build(fam, None, tag="deploy", train_site_idx=np.arange(len(S)), eval_site_idx=np.zeros(0, np.int64),
                       train_labels_override=tl)
        print(f"{fam} deploy: {p} {time.time() - t:.0f}s ({len(tl)} label pairs)", flush=True)
    raise SystemExit(0)

with manifest.RunManifest(a.name + ("_smoke" if a.smoke else ""), config=cfg | {"seeds": a.seeds, "xfit_tag": tag},
                          seed=a.seeds[0], split="deploy" if not a.smoke else "smoke",
                          notes="G5 deployment R2 (DEV+TEST labels)" if not a.smoke else "SMOKE: DEV labels, toy size") as run:
    out_dir = ROOT / "cache" / "release" / run.run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    fd = sample.FoldData(None, cfg["variant"], deploy_loader=loader)
    run.set("n_label_rows_loaded", fd.n_label_rows_loaded); run.set("n_positive_pairs", int(len(fd.pos)))
    run.set("n_train_sites", int(len(fd.train_sites)))
    A = assemble.Assembler(cfg["groups"], None, xfit_tag=tag, clean=cfg["clean"], comp_know=cfg["comp_know"],
                           f2a_family=cfg["f2a_family"], desc_v=cfg["desc_v"])
    for k_, v_ in A.fallbacks.items():
        run.count_fallback(k_, v_)
    kin = list(range(assemble.N_REAL))
    if a.smoke:
        kin = kin[:a.smoke_kinases]
    tw = iw = 0.0
    acc = keys = None
    logs = {}
    for sd in a.seeds:
        c = dict(cfg, seed=sd)
        lg = {}
        t1 = time.time()
        bst, cols, stats = train.train_fold(fd, A, c, lg)
        tw += time.time() - t1
        fd.eval_sites = np.arange(len(S))                   # all U_obs sites
        t2 = time.time()
        P = train.predict_fold(fd, A, bst, c, kinases=kin)
        iw += time.time() - t2
        if acc is None:
            keys = P[["kinase_idx", "site_idx"]].copy(); acc = P.score.values.astype(np.float64)
        else:
            if not (np.array_equal(keys.kinase_idx.values, P.kinase_idx.values) and np.array_equal(keys.site_idx.values, P.site_idx.values)):
                raise RuntimeError("seed predictions not aligned")
            acc += P.score.values
        logs[sd] = lg
        print(f"seed {sd}: {lg} rows {len(P)}", flush=True)
        del P
    score = (acc / len(a.seeds)).astype(np.float32)
    del acc
    f = out_dir / "r2_uobs_scores.parquet"
    kk = keys.kinase_idx.values
    bounds = np.flatnonzero(np.r_[True, kk[1:] != kk[:-1], True])
    schema = pa.schema([("kinase_idx", pa.int32()), ("site_idx", pa.int32()), ("score", pa.float32())])
    with pq.ParquetWriter(f, schema) as w:                  # chunked write (kinase blocks of ~5M rows)
        i = 0
        while i < len(bounds) - 1:
            j = i
            while j < len(bounds) - 1 and bounds[j + 1] - bounds[i] < 5_000_000:
                j += 1
            j = max(j, i + 1)
            lo, hi = bounds[i], bounds[j]
            w.write_table(pa.table({"kinase_idx": kk[lo:hi].astype(np.int32),
                                    "site_idx": keys.site_idx.values[lo:hi].astype(np.int32), "score": score[lo:hi]},
                                   schema=schema))
            i = j
    run.m["outputs"]["r2_uobs_scores"] = {"path": str(f.relative_to(ROOT)), "sha256": manifest.sha256_file(f)}
    run.m["n_rows"]["r2_uobs_scores"] = int(len(score))
    run.set("train_wall_s", round(tw, 1)); run.set("inference_wall_s", round(iw, 1)); run.set("seed_logs", logs)
    print("RUN", run.run_id, f, len(score), f"train {tw:.0f}s infer {iw:.0f}s total {time.time() - t0:.0f}s", flush=True)
