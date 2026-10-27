"""Phase B orchestrator for B3a (GPU 0, sequential).

1. Wait for cache/phosformer/pretrain/mlm_final.pt.
2. Fold 0 at neg ratio 1:16 and 1:50 (50% hard). Choose the ratio by INNER-validation macro AUROC(hard), not by the
   outer fold-0 OOF metric, so the outer fold stays untouched by selection. Fold-0 OOF M1 is logged for both, for
   information only.
3. Folds 1-4 at the chosen ratio.
4. Pooled OOF evaluation (ksp.eval.prepare('cv', fold=None)) of set-A predictions -> run b3a_phosformer_cv_pooled.
5. (--eall) re-score every fold with E-all (set B) from the saved best checkpoints, then re-pool.
Usage: python -I scripts/baselines/phosformer_cv_run.py [--eall] [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--init", default=str(ROOT / "cache/phosformer/pretrain/mlm_final.pt"))
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--ratios", default="16,50")
ap.add_argument("--eall", action="store_true")
ap.add_argument("--folds", default="0,1,2,3,4")
args = ap.parse_args()
from ksp.paths import PYTHON as PY  # noqa: E402  $KSP_PYTHON or the running interpreter
FT = str(REPO_ROOT / "scripts/baselines/phosformer_finetune.py")
state_f = ROOT / "results/baselines/phosformer/cv_state.json"
state = json.loads(state_f.read_text()) if state_f.exists() else {"runs": {}}


def log(m):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}", flush=True)


def save_state():
    state_f.write_text(json.dumps(state, indent=1, default=str))


def run_fold(fold: int, ratio: int, extra=()):
    key = f"f{fold}_r{ratio}" + ("_eall" if "--eall" in extra else "")
    if key in state["runs"]:
        log(f"skip {key} (done: {state['runs'][key]['run_id']})")
        return state["runs"][key]
    cmd = [PY, "-I", FT, "--fold", str(fold), "--neg-ratio", str(ratio), "--init", args.init,
           "--seed", str(args.seed), *extra]
    log("RUN " + " ".join(cmd))
    out = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
    logf = ROOT / f"cache/phosformer/finetune/{key}_s{args.seed}.log"
    logf.parent.mkdir(parents=True, exist_ok=True)
    logf.write_text(out.stdout + "\n--- stderr ---\n" + out.stderr)
    if out.returncode != 0:
        raise RuntimeError(f"{key} failed; see {logf}")
    res = json.loads(out.stdout.strip().splitlines()[-1])
    state["runs"][key] = res
    save_state()
    log(f"{key}: best inner-val {res['best_val']:.4f}; OOF M1 {res['endpoints'].get('M1_auroc_hard')}")
    return res


while not Path(args.init).exists():
    time.sleep(60)
log(f"init checkpoint present: {args.init}")
ratios = [int(r) for r in args.ratios.split(",")]
folds = [int(f) for f in args.folds.split(",")]
sel = {r: run_fold(0, r) for r in ratios}
chosen = max(ratios, key=lambda r: sel[r]["best_val"])
state["chosen_ratio"] = chosen
state["ratio_selection"] = {r: {"inner_val": sel[r]["best_val"], "fold0_oof_M1": sel[r]["endpoints"].get(
    "M1_auroc_hard")} for r in ratios}
save_state()
log(f"chosen ratio 1:{chosen} by inner validation")
for f in folds:
    if f != 0:
        run_fold(f, chosen)


def pool(tag: str, keys: list):
    from ksp import eval as E
    from ksp.manifest import RunManifest
    parts = []
    for k in keys:
        rid = state["runs"][k]["run_id"]
        parts.append(pd.read_parquet(ROOT / "cache/runs" / rid / "predictions.parquet"))
    preds = pd.concat(parts, ignore_index=True)
    data = E.prepare("cv", fold=None, with_canaries=False)
    ev = E.evaluate(preds, data, missing="drop")
    cfg = {"arm": "B3a", "chosen_ratio": chosen, "fold_runs": {k: state["runs"][k]["run_id"] for k in keys},
           "ratio_selection": state["ratio_selection"], "eall_complete": tag.endswith("eall")}
    with RunManifest(f"b3a_phosformer_cv_pooled{tag}", config=cfg, seed=args.seed,
                     inputs=[ROOT / "cache/runs" / state["runs"][k]["run_id"] / "predictions.parquet" for k in keys],
                     split="cv_pooled", notes="pooled OOF predictions of B3a fold runs; missing='drop'") as run:
        run.save_predictions(preds)
        run.save_metrics(ev)
        run.set("endpoints", ev.endpoints)
    state[f"pooled{tag}"] = {"run_id": run.run_id, "endpoints": ev.endpoints}
    save_state()
    log(f"pooled{tag}: {json.dumps(ev.endpoints, default=str)}")


pool("", [f"f{f}_r{chosen}" for f in folds])
if args.eall:
    for f in folds:
        run_fold(f, chosen, extra=("--eall", "--skip-train"))
    pool("_eall", [f"f{f}_r{chosen}_eall" for f in folds])
log("ALL DONE")
