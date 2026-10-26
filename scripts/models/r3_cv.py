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

# Agents never run git (DESIGN §11.8 / AGENT_BRIEF): do not let the manifest call `git rev-parse`.
manifest.git_hash = lambda: "not_queried(agent no-git rule)"

CONFIGS = {
    "w7_nce05": dict(flank=7, lambda_nce=0.5),
    "w10_nce05": dict(flank=10, lambda_nce=0.5),
    "w7_nce0": dict(flank=7, lambda_nce=0.0),
}
BASE: dict = {}      # base hyper-parameters chosen on inner validation (filled after the dev sweep)
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
    r = r[r.model.eq("R3")]
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
    return {k: E.get(k) for k in keys} | {"canary": res.canary}


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
    name = f"r3_{a.config}_{tag}"
    with manifest.RunManifest(name, config=cfg, seed=a.seed, inputs=inputs_list(),
                              split="cv" if a.fold is not None else "lko_T2",
                              notes="R3 neural kinase-site model (DESIGN §8.4)") as run:
        t0 = time.time()
        preds, info, model = r3train.run(cfg, fold=a.fold, lko_fold=a.lko, log=lambda s: print(s, flush=True))
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
        summ["train_info"] = {k: info[k] for k in ["best_val_auroc_hard", "best_step", "t_train_s", "t_score_s",
                                                    "t_total_s", "t_targets_s", "t_gpu_load_s", "peak_gpu_mem_gb",
                                                    "n_trained_kinases", "n_pos_fit", "n_pred_rows"]}
        run.save_metrics(summ, name="summary")
        print(json.dumps(summ, indent=1, default=str), flush=True)
        if a.fold is not None:
            registry_append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R3", "run_id": run.run_id,
                             "split": "cv", "fold": a.fold, "config_name": a.config,
                             "config_json": json.dumps(CONFIGS[a.config]), "selection_metric": "M1_auroc_hard(K_cv)",
                             "value": round(summ["M1_auroc_hard"], 5),
                             "n_configs_tried_model": n_configs_tried() + (0 if _seen(a.config) else 1),
                             "notes": a.note})


def cmd_inner(a, cfg):
    """Development run: train on fold-a.fold training rows and report ONLY the inner-validation metric (the CV eval
    fold is neither scored nor evaluated). Logged in REGISTRY as split=inner_dev."""
    from ksp.models.r3 import train as r3train
    t0 = time.time()
    _, info, _ = r3train.run(cfg, fold=a.fold, log=lambda s: print(s, flush=True), score=False)
    out = {"cfg": cfg, "best_val_auroc_hard": info["best_val_auroc_hard"], "best_step": info["best_step"],
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
    name = f"r3_{a.config}_{'lko_pooled' if lko else 'cv_pooled'}"
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
                summ["bootstrap"] = bs if isinstance(bs, dict) else str(bs)
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
    ap.add_argument("cmd", choices=["train", "pool", "pool-lko"])
    ap.add_argument("--fold", type=int)
    ap.add_argument("--lko", type=int)
    ap.add_argument("--config", default="w7_nce05", choices=sorted(CONFIGS))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-steps", type=int, default=0)
    ap.add_argument("--runs", nargs="*", default=[])
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--note", default="")
    ap.add_argument("--set", default="{}", help="JSON overrides of the training config")
    ap.add_argument("--set-name", default="")
    ap.add_argument("--inner-only", action="store_true")
    a = ap.parse_args()
    if a.cmd == "train":
        cmd_train(a)
    else:
        cmd_pool(a, lko=a.cmd == "pool-lko")
