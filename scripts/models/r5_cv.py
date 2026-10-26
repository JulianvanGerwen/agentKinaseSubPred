"""R5 CV / LKO driver (same protocol and harness as scripts/models/r3_cv.py). Owner: R5 agent.

  CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_cv.py train --fold 0 --config li35_cpt_aux
  python -I scripts/models/r5_cv.py pool --config li35_cpt_aux --runs <5 fold run_ids>

Every run: ksp.manifest.RunManifest -> results/runs/<run_id>/{manifest,metrics}.json, predictions in
cache/runs/<run_id>/predictions.parquet. Metrics only from ksp.eval (frozen harness, DESIGN A28).
Each CV training run appends a row to results/cv/REGISTRY.tsv.
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
    # fold-0 selection configs (<= 3, per protocol)
    "li35_base":     dict(trunk="35m"),
    "li35_cpt_aux":  dict(trunk="35m", pretrain_steps=1500, distill_aux=0.2),
    "li150_cpt_aux": dict(trunk="150m", batch=256, pretrain_batch=128, pretrain_steps=1500, distill_aux=0.2,
                          score_chunk=512, emb_on_gpu=False),
    # diagnostic ablations (not eligible for selection)
    "li35_cpt_aux_noli":  dict(trunk="35m", pretrain_steps=1500, distill_aux=0.2, use_li=False),
    "li35_cpt_aux_lionly": dict(trunk="35m", pretrain_steps=1500, distill_aux=0.2, use_pair=False),
}
BASE_NAME = "r5li"
REGISTRY = ROOT / "results" / "cv" / "REGISTRY.tsv"
REG_COLS = ["timestamp", "model", "run_id", "split", "fold", "config_name", "config_json", "selection_metric",
            "value", "n_configs_tried_model", "notes"]
OUTDIR = ROOT / "results" / "models" / "r5"


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
    return int(r[r.model.eq("R5") & r.split.eq("cv")].config_name.nunique())


def _seen(cfg_name) -> bool:
    if not REGISTRY.exists():
        return False
    r = pd.read_csv(REGISTRY, sep="\t")
    return bool((r.model.eq("R5") & r.config_name.eq(cfg_name)).any())


def inputs_list():
    return [ROOT / "cache/r3/win_embed_meta.json", ROOT / "cache/r3/win_rows.npy", ROOT / "cache/r3/win_sites.parquet",
            ROOT / "cache/r5/kin_tok.npy", ROOT / "cache/r5/kin_meta.json", ROOT / "cache/r5/win31_tokens.json",
            ROOT / "cache/features/esm2_650m/kinase_domain.npy", ROOT / "cache/labels/dev_pairs.parquet",
            ROOT / "cache/labels/canary_dev.parquet", ROOT / "results/splits/cv_folds.tsv",
            ROOT / "results/splits/kinase_folds.tsv", ROOT / "cache/features/motif/cantley_unprimed_score.json"]


def summarise(res) -> dict:
    from ksp import eval as kev
    E = res.endpoints
    keys = ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M5_auroc_hard_Kcv",
            "sec_auroc_all", "sec_mrr", "sec_hit1", "sec_hit10", "n_K_primary_evaluated", "n_K_any_evaluated",
            "n_missing_pairs"]
    bins = kev.strata(res, cols=("train_bin",), kinase_set="primary")
    return {k: E.get(k) for k in keys} | {"canary": res.canary, "M5_by_train_bin": bins.to_dict(orient="records")}


def cmd_train(a):
    import torch
    from ksp import eval as kev
    from ksp.models.r5 import train as r5train
    cfg = {**r5train.DEFAULT, **CONFIGS[a.config], **json.loads(a.set), "seed": a.seed}
    if a.max_steps:
        cfg["max_steps"] = a.max_steps
    if a.n_seeds:
        cfg["n_seeds"] = a.n_seeds
    tag = f"fold{a.fold}" if a.fold is not None else f"lko{a.lko}"
    name = f"r5_{BASE_NAME}_{a.config}_{tag}"
    with manifest.RunManifest(name, config=cfg, seed=a.seed, inputs=inputs_list(),
                              split="cv" if a.fold is not None else "lko_T2",
                              notes="R5 fine-tuned ESM-2 late-interaction kinase-site model") as run:
        t0 = time.time()
        preds, info, _ = r5train.run(cfg, fold=a.fold, lko_fold=a.lko, log=lambda s: print(s, flush=True))
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
                                                       "distill", "t_train_s", "t_score_s", "t_total_s",
                                                       "t_targets_s", "t_gpu_load_s", "peak_gpu_mem_gb",
                                                       "n_trained_kinases", "n_pos_fit", "n_pred_rows"]}
        run.save_metrics(summ, name="summary")
        print(json.dumps(summ, indent=1, default=str), flush=True)
        if a.fold is not None:
            registry_append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R5", "run_id": run.run_id,
                             "split": "cv", "fold": a.fold, "config_name": f"{BASE_NAME}:{a.config}",
                             "config_json": json.dumps(CONFIGS[a.config]),
                             "selection_metric": "M1_auroc_hard(K_cv)", "value": round(summ["M1_auroc_hard"], 5),
                             "n_configs_tried_model": n_configs_tried() + (0 if _seen(f"{BASE_NAME}:{a.config}") else 1),
                             "notes": a.note})


def _load_preds(run_ids):
    return pd.concat([pd.read_parquet(manifest.CACHE / r / "predictions.parquet",
                                      columns=["kinase_idx", "site_idx", "score"]) for r in run_ids],
                     ignore_index=True)


def cmd_pool(a, lko=False):
    from ksp import eval as kev
    t0 = time.time()
    P = _load_preds(a.runs)
    name = f"r5_{a.config}{'_' + a.tag if a.tag else ''}_{'lko_pooled' if lko else 'cv_pooled'}"
    with manifest.RunManifest(name, config={"config": CONFIGS[a.config], "source_runs": a.runs}, seed=a.seed,
                              inputs=[manifest.CACHE / r / "predictions.parquet" for r in a.runs],
                              split="lko_T2" if lko else "cv") as run:
        run.save_predictions(P)          # the leaderboard (scripts/reeval_cv_runs.py) re-evaluates from these
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
        if not lko and a.boot:
            try:
                summ["bootstrap"] = kev.bootstrap(res, B=a.boot).to_dict(orient="records")
            except Exception as e:
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
    ap.add_argument("--config", default="li35_cpt_aux", choices=sorted(CONFIGS))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-steps", type=int, default=0)
    ap.add_argument("--n-seeds", type=int, default=0)
    ap.add_argument("--runs", nargs="*", default=[])
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--note", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--set", default="{}")
    a = ap.parse_args()
    if a.cmd == "train":
        cmd_train(a)
    else:
        cmd_pool(a, lko=a.cmd == "pool-lko")
