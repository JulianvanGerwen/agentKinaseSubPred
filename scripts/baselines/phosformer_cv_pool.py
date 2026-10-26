"""B3a pooled 5-fold OOF: fill + pool.

Fold runs (phosformer_finetune.py --eall) score, per fold f: set A (fold positives U E-hard sites x all residue-matched
kinases) + E-all for kinases with >= 1 positive IN FOLD f. Pooled CV evaluation (ksp.eval.prepare('cv', fold=None))
needs E-all for every kinase with >= 1 pooled positive over every fold's universe. `fill` scores exactly the missing
(kinase, site) pairs of fold f with fold f's own model (OOF preserved: the model never saw fold f); `pool` concatenates
fold predictions + fills and evaluates.
Usage:
  CUDA_VISIBLE_DEVICES=g python -I scripts/baselines/phosformer_cv_pool.py fill --fold f
  python -I scripts/baselines/phosformer_cv_pool.py pool
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from ksp import eval as E
from ksp import index, splits

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=["fill", "pool"])
ap.add_argument("--fold", type=int)
ap.add_argument("--ratio", type=int, default=16)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
FILL = ROOT / "cache/phosformer/fill"
FILL.mkdir(parents=True, exist_ok=True)


def fold_run(f: int) -> dict:
    """Latest successful E-all-complete manifest of fold f."""
    best = None
    for m in sorted((ROOT / "results/runs").glob(f"b3a_phosformer_f{f}_r{args.ratio}__*/manifest.json")):
        d = json.loads(m.read_text())
        if d.get("status") == "ok" and d["config"].get("eall_complete") and d["seed"] == args.seed:
            best = d
    if best is None:
        raise FileNotFoundError(f"no ok E-all fold run for fold {f}")
    return best


S = index.load_sites()
K = index.load_kinases()
folds = splits.load_cv_folds()
data = E.prepare("cv", fold=None, with_canaries=False)
req_k = sorted(k for k in data.positives if not k.startswith("CANARY_"))
kcls = dict(zip(K.uniprot_acc, K.residue_class))
res = S.res.values


def cm(c):
    return np.ones(len(S), bool) if c == "STY" else (np.isin(res, ["S", "T"]) if c == "ST" else res == "Y")


if args.mode == "fill":
    import torch
    from ksp.baselines.phosformer.inputs import kinase_domains, window15_to_11mer
    from ksp.baselines.phosformer.model import PhosformerClassifier, phosformer_config
    from ksp.baselines.phosformer.score import score_pairs
    torch.set_num_threads(2)
    f = args.fold
    run = fold_run(f)
    preds = pd.read_parquet(ROOT / "cache/runs" / run["run_id"] / "predictions.parquet")
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    have = set(zip(preds.kinase_acc.values, spos[preds.site_id].values))
    funi = data.universe & S.acc.isin(set(folds.index[folds == f])).values
    rows = []
    for k in req_k:
        sel = np.where(funi & cm(kcls[k]))[0]
        miss = [s for s in sel if (k, s) not in have]
        if miss:
            rows.append(pd.DataFrame({"kinase_acc": k, "site_idx": np.array(miss, dtype=np.int64)}))
    pairs = pd.concat(rows, ignore_index=True)
    pep = np.array([window15_to_11mer(w) for w in S.window15.values], dtype=object)
    pairs["peptide"] = pep[pairs.site_idx.values]
    print(f"fold {f}: run {run['run_id']}; filling {len(pairs)} pairs over {pairs.kinase_acc.nunique()} kinases",
          flush=True)
    kd, _ = kinase_domains(K)
    dom = dict(zip(kd.uniprot_acc, kd.domain_seq))
    ck = ROOT / f"cache/phosformer/finetune/b3a_phosformer_f{f}_r{args.ratio}_s{args.seed}/best.pt"
    model = PhosformerClassifier(phosformer_config())
    model.load_state_dict(torch.load(ck, map_location="cpu", weights_only=True))
    model = model.cuda().eval()
    model.roberta = torch.compile(model.roberta, dynamic=True)
    t0 = time.time()
    pairs["score"] = score_pairs(model, dom, pairs, batch_size=1024, log_every=2000)
    out = pd.DataFrame({"kinase_acc": pairs.kinase_acc.values, "site_id": S.site_id.values[pairs.site_idx.values],
                        "score": pairs.score.values.astype(np.float64)})
    out.to_parquet(FILL / f"fill_f{f}_r{args.ratio}_s{args.seed}.parquet", index=False)
    (FILL / f"fill_f{f}_r{args.ratio}_s{args.seed}.json").write_text(json.dumps(
        {"fold": f, "fold_run": run["run_id"], "checkpoint": str(ck.relative_to(ROOT)), "n_pairs": len(out),
         "score_s": time.time() - t0}, indent=1))
    print(f"fold {f}: filled {len(out)} in {time.time() - t0:.0f}s", flush=True)
else:
    from ksp.manifest import RunManifest
    runs = {f: fold_run(f) for f in range(5)}
    parts, inputs = [], []
    for f, r in runs.items():
        p = ROOT / "cache/runs" / r["run_id"] / "predictions.parquet"
        q = FILL / f"fill_f{f}_r{args.ratio}_s{args.seed}.parquet"
        parts += [pd.read_parquet(p), pd.read_parquet(q)]
        inputs += [p, q]
    preds = pd.concat(parts, ignore_index=True)
    ev = E.evaluate(preds, data, missing="drop")
    fold_m = {f: {k: r["endpoints"].get(k) for k in ("M1_auroc_hard", "M2_hit5", "M3_top1_recall")}
              for f, r in runs.items()}
    cfg = {"arm": "B3a", "neg_ratio": args.ratio, "hard_frac": 0.5, "fold_runs": {f: r["run_id"] for f, r in runs.items()},
           "fill": "missing pooled E-all pairs scored by each fold's own model", "per_fold_endpoints": fold_m,
           "ratio_note": "1:16 fixed by PI (r2, 19:54) = paper value; 1:50 fold-0 run is an informational ablation"}
    with RunManifest("b3a_phosformer_cv_pooled", config=cfg, seed=args.seed, inputs=inputs, split="cv_pooled",
                     notes="B3a Phosformer reimplementation (own MLM on TEST-filtered corpus), pooled 5-fold OOF, "
                           "E-all complete for every kinase with a pooled positive") as run:
        run.save_predictions(preds)
        run.save_metrics(ev)
        run.set("endpoints", ev.endpoints)
        run.set("missing", {k: v for k, v in ev.missing.items() if k != "per_kinase"})
    print(json.dumps({"run_id": run.run_id, "endpoints": ev.endpoints}, default=str, indent=1))
