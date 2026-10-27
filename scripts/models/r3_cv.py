"""R3 CV / LKO driver (DESIGN §8.4, §4.5, §6). Owner: R3-neural.

  CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --fold 0 --config w7_nce05
  CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --lko 0 --config w7_nce05
  python -I scripts/models/r3_cv.py pool --config w7_nce05 --runs <run_id> ...      (pooled 5-fold OOF evaluation)
  python -I scripts/models/r3_cv.py pool-lko --config w7_nce05 --runs <run_id> ...  (pooled LKO T2 evaluation)

Every run: ksp.manifest.RunManifest -> results/runs/<run_id>/{manifest,metrics}.json, predictions in
cache/runs/<run_id>/predictions.parquet (kinase_idx, site_idx, kinase_acc, site_id, score). Metrics only from ksp.eval.
Each CV training run appends a row to results/cv/REGISTRY.tsv (configs tried are counted there).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp import manifest  # noqa: E402




CONFIGS = {
    "w7_nce05": dict(flank=7, lambda_nce=0.5),
    "w10_nce05": dict(flank=10, lambda_nce=0.5),
    "w7_nce0": dict(flank=7, lambda_nce=0.0),
    # diagnostic ablation (not eligible for selection): no protein-context token, no phospho channel
    "w7_nce05_noctx": dict(flank=7, lambda_nce=0.5, use_prot=False, use_phos=False),
    # round 2 (PI): Cantley-atlas distillation pretraining (label-free targets, training-fold sites only)
    "w7_nce05_cpt": dict(flank=7, lambda_nce=0.5, pretrain_steps=3000),
    "w7_nce05_cpt_aux": dict(flank=7, lambda_nce=0.5, pretrain_steps=3000, distill_aux=0.2),
    # round 2 (PI): ESM-2 3B per-residue window embeddings (cache/r3/esm3b, +-7)
    "w7_nce05_esm3b": dict(flank=7, lambda_nce=0.5, esm="3b"),
    "w7_nce05_cptaux_esm3b": dict(flank=7, lambda_nce=0.5, pretrain_steps=3000, distill_aux=0.2, esm="3b"),
}
# Base hyper-parameters ("reg1"), chosen by the fold-0 INNER-validation dev sweep (scripts/models/r3_dev_sweep.sh;
# REGISTRY split=inner_dev); 5-seed ensemble (seeds [0..4]) per training run, logits averaged.
BASE: dict = dict(lr=3e-4, p_lab=0.25, dropout=0.2, esm_dropout=0.3, wd=0.05, eval_every=100, patience=10, n_seeds=5)
BASE_NAME = "reg1x5"
REGISTRY = ROOT / "results" / "cv" / "REGISTRY.tsv"
REG_COLS = ["timestamp", "model", "run_id", "split", "fold", "config_name", "config_json", "selection_metric",
            "value", "n_configs_tried_model", "notes"]
OUTDIR = ROOT / "results" / "models" / "r3"


def registry_append(row: dict):
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    new = not REGISTRY.exists()
    with open(REGISTRY, "a") as fh:
        if new:
            fh.write("\t".join(REG_COLS) + "\n")
        fh.write("\t".join(str(row.get(c, "")) for c in REG_COLS) + "\n")


def n_configs_tried() -> int:
    if not REGISTRY.exists():
        return 0
    r = pd.read_csv(REGISTRY, sep="\t")
    r = r[r.model.eq("R3") & r.split.eq("cv")]
    return int(r.config_name.nunique())


def inputs_list():
    return [ROOT / "cache/r3/win_embed_meta.json", ROOT / "cache/r3/win_rows.npy", ROOT / "cache/r3/win_sites.parquet",
            ROOT / "cache/features/esm2_650m/kinase_domain.npy", ROOT / "cache/labels/dev_pairs.parquet",
            ROOT / "cache/labels/canary_dev.parquet", ROOT / "results/splits/cv_folds.tsv",
            ROOT / "results/splits/kinase_folds.tsv", ROOT / "cache/features/motif/cantley_unprimed_score.json"]


def summarise(res) -> dict:
    E = res.endpoints
    keys = ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M5_auroc_hard_Kcv",
            "sec_auroc_all", "sec_mrr", "sec_hit1", "sec_hit10", "n_K_primary_evaluated", "n_K_any_evaluated",
            "n_missing_pairs"]
    from ksp import eval as kev
    bins = kev.strata(res, cols=("train_bin",), kinase_set="primary")
    return {k: E.get(k) for k in keys} | {"canary": res.canary,
                                          "M5_by_train_bin": bins.to_dict(orient="records")}


def cmd_train(a):
    import torch
    from ksp import eval as kev
    from ksp.models.r3 import train as r3train
    cfg = {**r3train.DEFAULT, **BASE, **CONFIGS[a.config], **json.loads(a.set), "seed": a.seed}
    if a.max_steps:
        cfg["max_steps"] = a.max_steps
    if a.inner_only:
        return cmd_inner(a, cfg)
    tag = f"fold{a.fold}" if a.fold is not None else f"lko{a.lko}"
    name = f"r3_{BASE_NAME}_{a.config}_{tag}"
    with manifest.RunManifest(name, config=cfg, seed=a.seed, inputs=inputs_list(),
                              split="cv" if a.fold is not None else "lko_T2",
                              notes="R3 neural kinase-site model (DESIGN §8.4)") as run:
        t0 = time.time()
        preds, info, _ = r3train.run(cfg, fold=a.fold, lko_fold=a.lko, log=lambda s: print(s, flush=True))
        info["gpu"] = torch.cuda.get_device_name(0)
        run.set("train_info", {k: v for k, v in info.items() if k != "history"})
        run.set("history", info["history"])
        run.set("runtime_train_and_score_s", round(time.time() - t0, 1))
        run.save_predictions(preds)
        print(f"[{name}] predictions {len(preds)} rows; train {info['t_train_s']}s score {info['t_score_s']}s",
              flush=True)
        t1 = time.time()
        if a.fold is not None:
            data = kev.prepare("cv", fold=a.fold)
            res = kev.evaluate(preds[["kinase_idx", "site_idx", "score"]], data)
        else:
            data = kev.prepare("cv", fold=None)
            held = sorted(set(preds.kinase_acc.astype(str)) & set(data.positives))
            res = kev.evaluate(preds[["kinase_idx", "site_idx", "score"]], data, kinases=held, missing="drop")
        run.save_metrics(res)
        res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
        summ = summarise(res)
        summ["t_eval_s"] = round(time.time() - t1, 1)
        summ["train_info"] = {k: info.get(k) for k in ["best_val_auroc_hard", "ensemble_val_auroc_hard", "best_step",
                                                        "distill", "t_train_s", "t_score_s",
                                                    "t_total_s", "t_targets_s", "t_gpu_load_s", "peak_gpu_mem_gb",
                                                    "n_trained_kinases", "n_pos_fit", "n_pred_rows"]}
        run.save_metrics(summ, name="summary")
        print(json.dumps(summ, indent=1, default=str), flush=True)
        if a.fold is not None:
            registry_append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R3", "run_id": run.run_id,
                             "split": "cv", "fold": a.fold, "config_name": f"{BASE_NAME}:{a.config}",
                             "config_json": json.dumps({**BASE, **CONFIGS[a.config]}), "selection_metric": "M1_auroc_hard(K_cv)",
                             "value": round(summ["M1_auroc_hard"], 5),
                             "n_configs_tried_model": n_configs_tried() + (0 if _seen(f"{BASE_NAME}:{a.config}") else 1),
                             "notes": a.note})


def cmd_inner(a, cfg):
    """Development run: train on fold-a.fold training rows and report ONLY the inner-validation metric (the CV eval
    fold is neither scored nor evaluated). Logged in REGISTRY as split=inner_dev."""
    from ksp.models.r3 import train as r3train
    t0 = time.time()
    _, info, _ = r3train.run(cfg, fold=a.fold, log=lambda s: print(s, flush=True), score=False)
    out = {"cfg": cfg, "best_val_auroc_hard": info["best_val_auroc_hard"], "best_step": info["best_step"],
           "ensemble_val_auroc_hard": info.get("ensemble_val_auroc_hard"),
           "t_train_s": info["t_train_s"], "history": info["history"]}
    d = OUTDIR / "inner_dev"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{a.config}_{a.set_name or 'base'}_fold{a.fold}.json").write_text(json.dumps(out, indent=1))
    registry_append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R3", "run_id": "",
                     "split": "inner_dev", "fold": a.fold, "config_name": f"dev:{a.config}+{a.set_name or 'base'}",
                     "config_json": json.dumps({**CONFIGS[a.config], **json.loads(a.set)}),
                     "selection_metric": "inner_val_auroc_hard(train rows of fold)", "value":
                     round(info["best_val_auroc_hard"], 5), "n_configs_tried_model": "",
                     "notes": f"dev, eval fold untouched; {time.time() - t0:.0f}s"})
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1))


def cmd_test(a):
    """Final-model TEST scoring component (DESIGN §8.6; G4 runbook): train on ALL DEV with the frozen config and score
    the TEST candidates (TEST-protein U_ms sites x real kinases of matching class). --lko J: T3 component (LKO fold J
    kinases removed from training; only they are scored). Reads no TEST labels."""
    import torch
    from ksp.models.r3 import train as r3train
    cfg = {**r3train.DEFAULT, **BASE, **CONFIGS[a.config], "seed": a.seed}
    tag = "final" if a.lko is None else f"t3lko{a.lko}"
    with manifest.RunManifest(f"r3_final_distill_test_{tag}", config=cfg, seed=a.seed, inputs=inputs_list(),
                              split="final_dev_to_test",
                              notes=f"R3 final ({BASE_NAME}:{a.config}) trained on all DEV; TEST component {tag}") as run:
        preds, info, _ = r3train.run(cfg, lko_fold=a.lko, test=True, log=lambda s: print(s, flush=True))
        info["gpu"] = torch.cuda.get_device_name(0)
        run.set("train_info", {k: v for k, v in info.items() if k != "history"})
        run.set("history", info["history"])
        run.set("train_wall_s", info["t_train_wall_s"])
        run.set("inference_wall_s", info["t_score_s"])
        run.save_predictions(preds[["kinase_idx", "site_idx", "score"]], name="component_predictions")
        print(json.dumps({"run_id": run.run_id, "n_rows": len(preds), "train_wall_s": info["t_train_wall_s"],
                          "inference_wall_s": info["t_score_s"], "n_kinases": int(preds.kinase_idx.nunique())}),
              flush=True)


def cmd_assemble_test(a):
    """Assemble the G4 deliverables into ONE run: test_predictions.parquet (final component) and
    test_t3_predictions.parquet (the 5 LKO components concatenated), with runtimes and the --check-predictions
    pre-flight result recorded in the manifest."""
    from ksp import eval_test
    comps = {"final": a.final, **{f"lko{j}": r for j, r in enumerate(a.lko_runs)}}
    files = {k: manifest.CACHE / r / "component_predictions.parquet" for k, r in comps.items()}
    mans = {k: json.loads((manifest.RESULTS / r / "manifest.json").read_text()) for k, r in comps.items()}
    for k, m in mans.items():
        if m.get("status") != "ok":
            raise RuntimeError(f"component {k} ({comps[k]}) status {m.get('status')}")
        want = None if k == "final" else int(k[3:])
        if m["train_info"]["lko_fold"] != want or not m["train_info"]["test"]:
            raise RuntimeError(f"component {k}: wrong lko_fold/test flag")
    with manifest.RunManifest("r3_final_distill_test", config={"variant": f"{BASE_NAME}:{a.config}",
                              "components": comps, "config": mans["final"]["config"]}, seed=a.seed,
                              inputs=list(files.values()), split="final_dev_to_test",
                              notes="R3 final TEST deliverables (G4 runbook items 1-3)") as run:
        P = pd.read_parquet(files["final"])
        run.save_predictions(P, name="test_predictions")
        T3 = pd.concat([pd.read_parquet(files[f"lko{j}"]) for j in range(5)], ignore_index=True)
        if T3.duplicated(["kinase_idx", "site_idx"]).any():
            raise RuntimeError("duplicate (kinase, site) rows across LKO components")
        run.save_predictions(T3, name="test_t3_predictions")
        tf = mans["final"]
        run.set("train_wall_s", tf["train_wall_s"])
        run.set("inference_wall_s", tf["inference_wall_s"])
        run.set("t3_train_wall_s_sum", round(sum(mans[f"lko{j}"]["train_wall_s"] for j in range(5)), 1))
        run.set("t3_inference_wall_s_sum", round(sum(mans[f"lko{j}"]["inference_wall_s"] for j in range(5)), 1))
        chk = {n: eval_test.check_predictions(manifest.CACHE / run.run_id / f"{n}.parquet")
               for n in ["test_predictions", "test_t3_predictions"]}
        run.set("check_predictions", chk)
        print(json.dumps({"run_id": run.run_id, "check": chk, "train_wall_s": tf["train_wall_s"],
                          "inference_wall_s": tf["inference_wall_s"]}, indent=1, default=str), flush=True)


def cmd_deploy(a):
    """G5 deployment (DESIGN §8.6): the frozen final config (BASE + w7_nce05_cpt_aux, 5 seeds) retrained on DEV+TEST
    labels (ksp.deploy.load_deploy_labels; canaries off) on U_ms rows, then scoring U_obs x real kinases of matching
    class (~131M pairs) from the U_obs window cache (identical ESM settings). Output (chunked parquet):
    cache/release/r3_deploy_<ts>/uobs_predictions.parquet [kinase_idx int32, site_idx int32, score float32 = raw
    5-seed mean logit]. --smoke: DEV labels (ksp.labels.load_dev), DEV rows, written under --out (never cache/release)."""
    import hashlib
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch
    from ksp import labels as klabels
    from ksp.models.r3 import train as r3train
    FROZEN = "w7_nce05_cpt_aux"                     # r3_final_distill (FREEZE)
    if not a.smoke and (a.max_steps or a.n_seeds or a.pretrain_steps is not None):
        raise SystemExit("real deployment uses the frozen config only (no overrides)")
    a.config = FROZEN
    cfg = {**r3train.DEFAULT, **BASE, **CONFIGS[a.config], "seed": a.seed, "esm": "650m_uobs"}
    if a.max_steps:
        cfg["max_steps"] = a.max_steps
    if a.n_seeds:
        cfg["n_seeds"] = a.n_seeds
    if a.pretrain_steps is not None:
        cfg["pretrain_steps"] = a.pretrain_steps
    ts = time.strftime("%Y%m%d-%H%M%S")
    if a.smoke:
        loader = lambda v: klabels.load_dev(v, with_canaries=False)  # noqa: E731
        rows_mode, out_dir = "dev", Path(a.out) / f"r3_deploy_smoke_{ts}"
    else:
        from ksp import deploy
        deploy.check_test_completed()               # refuse early (before any GPU work) if TEST is not COMPLETED
        loader = lambda v: deploy.load_deploy_labels(v)  # noqa: E731
        rows_mode, out_dir = "all", ROOT / "cache" / "release" / f"r3_deploy_{ts}"
    name = "r3_deploy_smoke" if a.smoke else "r3_deploy"
    with manifest.RunManifest(name, config=cfg, seed=a.seed, inputs=inputs_list() + [ROOT / "cache/r3/uobs/win_embed_meta.json"],
                              split="deploy_dev_test" if not a.smoke else "deploy_smoke_dev",
                              notes="R3 deployment model (frozen r3_final_distill config)") as run:
        um = json.loads((ROOT / "cache/r3/uobs/win_embed_meta.json").read_text())
        run.set("input_cache_note", {
            "what": "U_obs window-embedding cache cache/r3/uobs (rebuilt for deployment; PI-accepted 00:58): the U_ms "
                    "cache cache/r3 lacks the 3,382 non-U_ms U_obs sites. Identical ESM-2 650M settings "
                    "(scripts/models/r3_window_embed.py --variant 650m_uobs); input cache change only, config unchanged; "
                    "training rows remain U_ms; phospho channel = flanking U_ms sites only.",
            "n_sites": um["n_sites"], "n_residues": um["n_residues"],
            "centre_cos_vs_sty_cache": um["centre_cos_vs_sty_cache"],
            "window_residue_mismatches": um["window_residue_mismatches"], "gpu_minutes": um["gpu_minutes"]})
        res, info, _ = r3train.run(cfg, deploy=True, label_loader=loader, deploy_train_rows=rows_mode,
                                   log=lambda s: print(s, flush=True))
        t0 = time.time()
        S, sidx, ok = res["S"], res["site_idx"], res["ok"]
        if not np.isfinite(S[ok]).all():
            raise RuntimeError("non-finite deployment scores")
        out_dir.mkdir(parents=True, exist_ok=True)
        f = out_dir / "uobs_predictions.parquet"
        tmp = f.with_suffix(".tmp.parquet")
        schema = pa.schema([("kinase_idx", pa.int32()), ("site_idx", pa.int32()), ("score", pa.float32())])
        order = np.argsort(sidx)                   # rows ordered by site_idx, then kinase_idx
        n = 0
        with pq.ParquetWriter(tmp, schema, compression="zstd") as w:
            for c0 in range(0, len(order), 20000):
                rr = order[c0:c0 + 20000]
                ri, kk = np.nonzero(ok[rr])
                tb = pa.table({"kinase_idx": kk.astype(np.int32), "site_idx": sidx[rr][ri].astype(np.int32),
                               "score": S[rr][ri, kk].astype(np.float32)}, schema=schema)
                w.write_table(tb)
                n += len(kk)
        tmp.replace(f)
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 22), b""):
                h.update(b)
        t_write = time.time() - t0
        out = {"path": str(f.relative_to(ROOT)) if str(f).startswith(str(ROOT)) else str(f), "sha256": h.hexdigest(),
               "n_rows": n, "n_sites": int(len(sidx)), "n_kinases": int(ok.any(0).sum()),
               "train_wall_s": info["t_train_wall_s"], "inference_wall_s": round(info["t_score_s"] + t_write, 1),
               "score_s": info["t_score_s"], "write_s": round(t_write, 1), "peak_gpu_mem_gb": info["peak_gpu_mem_gb"],
               "gpu": torch.cuda.get_device_name(0), "train_rows": rows_mode, "targets": info["targets"],
               "ensemble_val_auroc_hard": info.get("ensemble_val_auroc_hard"), "best_step": info["best_step"]}
        run.set("train_info", {k: v for k, v in info.items() if k != "history"})
        run.set("history", info["history"])
        run.set("train_wall_s", out["train_wall_s"])
        run.set("inference_wall_s", out["inference_wall_s"])
        run.set("deploy_output", out)
        run.m["n_rows"]["uobs_predictions"] = n
        (out_dir / "README.json").write_text(json.dumps(out | {"run_id": run.run_id, "config": cfg}, indent=1,
                                                         default=str))
        print(json.dumps(out | {"run_id": run.run_id}, indent=1, default=str), flush=True)


def _seen(cfg_name) -> bool:
    if not REGISTRY.exists():
        return False
    r = pd.read_csv(REGISTRY, sep="\t")
    return bool((r.model.eq("R3") & r.config_name.eq(cfg_name)).any())


def _load_preds(run_ids):
    parts = []
    for rid in run_ids:
        p = pd.read_parquet(manifest.CACHE / rid / "predictions.parquet", columns=["kinase_idx", "site_idx", "score"])
        parts.append(p)
    return pd.concat(parts, ignore_index=True)


def cmd_pool(a, lko=False):
    from ksp import eval as kev
    t0 = time.time()
    P = _load_preds(a.runs)
    name = f"r3_{a.config}{'_' + a.tag if a.tag else ''}_{'lko_pooled' if lko else 'cv_pooled'}"
    with manifest.RunManifest(name, config={"config": CONFIGS[a.config], "source_runs": a.runs}, seed=a.seed,
                              inputs=[manifest.CACHE / r / "predictions.parquet" for r in a.runs],
                              split="lko_T2" if lko else "cv") as run:
        data = kev.prepare("cv", fold=None)
        if lko:
            from ksp import canary
            K = canary.kinases_with_canaries()
            held = sorted(set(K.uniprot_acc.values[np.unique(P.kinase_idx.values)]) & set(data.positives))
            res = kev.evaluate(P, data, kinases=held, missing="drop")
        else:
            res = kev.evaluate(P, data)
        run.save_metrics(res)
        res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
        kev.strata(res).to_csv(run.dir / "strata.tsv", sep="\t", index=False)
        summ = summarise(res)
        summ["n_pred_rows"] = len(P)
        if not lko:
            try:
                bs = kev.bootstrap(res, B=a.boot)
                summ["bootstrap"] = bs.to_dict(orient="records")
            except Exception as e:  # recorded, not hidden
                summ["bootstrap_error"] = repr(e)
                run.count_fallback("bootstrap_failed")
        summ["t_eval_s"] = round(time.time() - t0, 1)
        run.save_metrics(summ, name="summary")
        OUTDIR.mkdir(parents=True, exist_ok=True)
        (OUTDIR / f"{name}_summary.json").write_text(json.dumps(summ | {"run_id": run.run_id}, indent=1,
                                                               default=kev._json_default))
        res.per_kinase.to_csv(OUTDIR / f"{name}_per_kinase.tsv", sep="\t", index=False)
        print(json.dumps(summ, indent=1, default=kev._json_default), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "pool", "pool-lko", "test", "assemble-test", "deploy"])
    ap.add_argument("--smoke", action="store_true", help="deploy: DEV labels + DEV rows, output under --out")
    ap.add_argument("--out", default=str(ROOT / "cache" / "r3" / "deploy_smoke"))
    ap.add_argument("--n-seeds", type=int, default=0)
    ap.add_argument("--pretrain-steps", type=int, default=None)
    ap.add_argument("--final", help="assemble-test: final component run id")
    ap.add_argument("--lko-runs", nargs=5, default=[], help="assemble-test: T3 component run ids for LKO folds 0..4")
    ap.add_argument("--fold", type=int)
    ap.add_argument("--lko", type=int)
    ap.add_argument("--config", default="w7_nce05", choices=sorted(CONFIGS))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-steps", type=int, default=0)
    ap.add_argument("--runs", nargs="*", default=[])
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--note", default="")
    ap.add_argument("--tag", default="", help="pool: extra name tag (e.g. canaryv2)")
    ap.add_argument("--set", default="{}", help="JSON overrides of the training config")
    ap.add_argument("--set-name", default="")
    ap.add_argument("--inner-only", action="store_true")
    a = ap.parse_args()
    if a.cmd == "train":
        cmd_train(a)
    elif a.cmd == "test":
        cmd_test(a)
    elif a.cmd == "assemble-test":
        cmd_assemble_test(a)
    elif a.cmd == "deploy":
        cmd_deploy(a)
    else:
        cmd_pool(a, lko=a.cmd == "pool-lko")
