"""B3_retrained TEST deliverable (docs/G4_RUNBOOK.md): score every TEST candidate pair with the all-DEV B3a model.

Model: phosformer_finetune.py --fold -1 (all DEV, frozen config = f0_r16: 1:16, 50% hard, augmentations, focal
gamma 2, lr 2e-5, batch 64, early stopping on a cluster-grouped inner 20% DEV split), initialised from the
TEST-filtered MLM checkpoint cache/phosformer/pretrain/mlm_final.pt (A22). Never reads the sealed TEST labels.
Candidates: ksp.eval_test.test_candidate_pairs() (public, label-free). Sharded BY KINASE, so each kinase is scored
completely by one shard (all-or-nothing per kinase).
Usage:
  CUDA_VISIBLE_DEVICES=g python -I scripts/baselines/phosformer_test.py score --shard i --n-shards n
  python -I scripts/baselines/phosformer_test.py assemble --n-shards n
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

from ksp import eval_test, index

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=["score", "assemble"])
ap.add_argument("--shard", type=int, default=0)
ap.add_argument("--n-shards", type=int, default=4)
ap.add_argument("--ckdir", default=str(ROOT / "cache/phosformer/finetune/b3_retrained_devfull_s0"))
args = ap.parse_args()
OUT = ROOT / "cache/phosformer/test"
OUT.mkdir(parents=True, exist_ok=True)
ckdir = Path(args.ckdir)
K = index.load_kinases()


def shard_kinases(n: int) -> list[np.ndarray]:
    """Greedy balance of kinase_idx over shards by candidate count x domain length (deterministic)."""
    cand = eval_test.test_candidate_pairs()
    cnt = cand.groupby("kinase_idx").size()
    from ksp.baselines.phosformer.inputs import kinase_domains
    kd, _ = kinase_domains(K)
    L = kd.domain_seq.str.len().values
    cost = cnt * (L[cnt.index.values] + 15)
    load, out = np.zeros(n), [[] for _ in range(n)]
    for k in cost.sort_values(ascending=False).index:
        j = int(np.argmin(load))
        out[j].append(int(k))
        load[j] += cost[k]
    return [np.array(sorted(o)) for o in out], cand


if args.mode == "score":
    import torch
    from ksp.baselines.phosformer.inputs import kinase_domains, window15_to_11mer
    from ksp.baselines.phosformer.model import PhosformerClassifier, phosformer_config
    from ksp.baselines.phosformer.score import score_pairs
    torch.set_num_threads(2)
    shards, cand = shard_kinases(args.n_shards)
    mine = shards[args.shard]
    c = cand[cand.kinase_idx.isin(mine)].reset_index(drop=True)
    S = index.load_sites()
    pep = np.array([window15_to_11mer(w) for w in S.window15.values], dtype=object)
    kd, fb = kinase_domains(K)
    dom = dict(zip(kd.uniprot_acc, kd.domain_seq))
    pairs = pd.DataFrame({"kinase_acc": K.uniprot_acc.values[c.kinase_idx.values], "peptide": pep[c.site_idx.values]})
    print(f"shard {args.shard}/{args.n_shards}: {len(mine)} kinases, {len(pairs)} pairs", flush=True)
    model = PhosformerClassifier(phosformer_config())
    model.load_state_dict(torch.load(ckdir / "best.pt", map_location="cpu", weights_only=True))
    model = model.cuda().eval()
    model.roberta = torch.compile(model.roberta, dynamic=True)
    t0 = time.time()
    sc = score_pairs(model, dom, pairs, batch_size=1024, log_every=2000)
    wall = time.time() - t0
    out = pd.DataFrame({"kinase_idx": c.kinase_idx.values.astype(np.int64),
                        "site_idx": c.site_idx.values.astype(np.int64), "score": sc.astype(np.float64)})
    out.to_parquet(OUT / f"shard{args.shard}of{args.n_shards}.parquet", index=False)
    (OUT / f"shard{args.shard}of{args.n_shards}.json").write_text(json.dumps(
        {"shard": args.shard, "n_kinases": len(mine), "n_pairs": len(out), "inference_wall_s": wall,
         "gpu": __import__("os").environ.get("CUDA_VISIBLE_DEVICES"), "kinase_domain_fallbacks": len(fb)}, indent=1))
    print(f"shard {args.shard}: {len(out)} pairs in {wall:.0f}s ({len(out) / wall:,.0f}/s)", flush=True)
else:
    from ksp.manifest import RunManifest
    parts, meta = [], []
    for i in range(args.n_shards):
        parts.append(pd.read_parquet(OUT / f"shard{i}of{args.n_shards}.parquet"))
        meta.append(json.loads((OUT / f"shard{i}of{args.n_shards}.json").read_text()))
    preds = pd.concat(parts, ignore_index=True).sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
    cand = eval_test.test_candidate_pairs()
    if len(preds) != len(cand):
        raise SystemExit(f"row count {len(preds)} != candidate count {len(cand)}")
    if not np.isfinite(preds.score).all():
        raise SystemExit("non-finite scores")
    ts = json.loads((ckdir / "train_summary.json").read_text())
    inf_sum = float(sum(m["inference_wall_s"] for m in meta))
    cfg = {"baseline_id": "B3_retrained", "arm": "B3a (Phosformer reimplementation, retrained)",
           "train": ts["args"], "best_inner_val_macro_auroc_hard": ts["best_val"],
           "init": "cache/phosformer/pretrain/mlm_final.pt (MLM 50 ep on TEST-filtered corpus, A22)",
           "frozen_from": "CV config f0_r16 (neg 1:16, 50% hard); PI fixed ratio 19:54",
           "shards": meta, "inference_wall_s_note": "sum of per-shard GPU wall times (single-GPU equivalent)"}
    inputs = [ckdir / "best.pt", ROOT / "cache/phosformer/pretrain/mlm_final.pt",
              ROOT / "results/splits/test_proteins.tsv", ROOT / "cache/index/kinases.parquet",
              ROOT / "cache/index/sites.parquet"]
    with RunManifest("b3_retrained_test", config=cfg, seed=int(ts["args"]["seed"]), inputs=inputs, split="test",
                     notes="B3_retrained TEST predictions; exhaustive over ksp.eval_test.test_candidate_pairs()") as run:
        run.count_fallback("kinase_domain_full_seq", meta[0]["kinase_domain_fallbacks"])
        f = run.save_predictions(preds, name="test_predictions")
        run.set("train_wall_s", float(ts["train_s"]))
        run.set("inference_wall_s", inf_sum)
        run.set("runtime", {"train_wall_s": float(ts["train_s"]), "inference_wall_s": inf_sum})
    print(json.dumps({"run_id": run.run_id, "path": str(f.relative_to(ROOT)), "n": len(preds),
                      "train_wall_s": ts["train_s"], "inference_wall_s": inf_sum}, indent=1))
