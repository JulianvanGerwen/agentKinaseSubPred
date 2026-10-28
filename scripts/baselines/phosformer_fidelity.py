"""Training-pipeline fidelity check (PI 22:35): our TEST-filtered MLM checkpoint + our B3a fine-tuning objective
(focal gamma 2, AdamW lr 2e-5, batch 64, bf16, early stopping patience 5 on validation AUROC) trained on THEIR
published train split (hf_dataset_curated/splits/train.csv: pre-augmented by the authors: upsampling, boundary
shifts, 5% '.'-masking, 1:16 multilevel negatives) and evaluated on THEIR test split (test_pos_neg.csv; paper
Fig. 4A AUROC 0.885; released weights through our scorer 0.8855). Early stopping uses THEIR val_pos_neg.csv.
Rows are used as given (their sampling); --frac subsamples positives AND negatives uniformly (keeps 1:16).
Usage: CUDA_VISIBLE_DEVICES=g python -I scripts/baselines/phosformer_fidelity.py --frac 0.6
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from ksp.baselines.phosformer.model import PhosformerClassifier, focal_loss, phosformer_config
from ksp.baselines.phosformer.tokens import collate_pairs, encode_pair

ap = argparse.ArgumentParser()
ap.add_argument("--init", default=str(ROOT / "cache/phosformer/pretrain/mlm_final.pt"))
ap.add_argument("--frac", type=float, default=1.0)
ap.add_argument("--epochs", type=float, default=1.0)
ap.add_argument("--batch", type=int, default=64)
ap.add_argument("--lr", type=float, default=2e-5)
ap.add_argument("--eval-every", type=int, default=2000)
ap.add_argument("--patience", type=int, default=5)
ap.add_argument("--max-hours", type=float, default=1.3, help="hard wall-clock cap on training")
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
torch.set_num_threads(2)
torch.manual_seed(args.seed)
SPL = ROOT / "external/phosformer/hf_dataset_curated/splits"
OUT = ROOT / "results/baselines/phosformer"
CK = ROOT / "cache/phosformer/fidelity"
CK.mkdir(parents=True, exist_ok=True)
dev = "cuda"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


tr = pd.read_csv(SPL / "train.csv", usecols=["kinase", "substrate", "label"])
n_full = len(tr)
rng = np.random.default_rng(args.seed)
if args.frac < 1:
    keep = np.zeros(len(tr), bool)
    for y in (0, 1):
        idx = np.where(tr.label.values == y)[0]
        keep[rng.choice(idx, int(round(args.frac * len(idx))), replace=False)] = True
    tr = tr[keep].reset_index(drop=True)
va = pd.read_csv(SPL / "val_pos_neg.csv")
te = pd.read_csv(SPL / "test_pos_neg.csv")
neg = pd.read_csv(SPL / "test_neg_only.csv")
log(f"train rows {len(tr)} of {n_full} (frac {args.frac}; pos {int(tr.label.sum())}), val {len(va)}, test {len(te)}")
K = tr.kinase.values
P = tr.substrate.values
Y = tr.label.values.astype(np.int64)

cfg = phosformer_config()
ck = torch.load(args.init, map_location="cpu", weights_only=False)
model = PhosformerClassifier.from_mlm(ck["model"], cfg).to(dev)
opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.0)
total = int(args.epochs * int(np.ceil(len(tr) / args.batch)))
sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: max(0.0, 1.0 - s / total))


def batch(kin, pep):
    ids, att = collate_pairs([encode_pair(p, k) for k, p in zip(kin, pep)])
    return torch.from_numpy(ids).to(dev), torch.from_numpy(att).to(dev)


@torch.no_grad()
def predict(df):
    model.eval()
    out = np.empty(len(df), np.float32)
    for b0 in range(0, len(df), 256):
        ids, att = batch(df.kinase.values[b0:b0 + 256], df.substrate.values[b0:b0 + 256])
        with torch.autocast("cuda", dtype=torch.bfloat16):
            out[b0:b0 + 256] = model.prob(ids, att).float().cpu().numpy()
    model.train()
    return out


def per_kinase_mean(df, pr):
    a = []
    for _, g in df.assign(p=pr).groupby("uniprot"):
        if g.label.nunique() == 2:
            a.append(roc_auc_score(g.label, g.p))
    return float(np.mean(a)), len(a)


hist, best, bad, step = [], -1.0, 0, 0
t0 = time.time()
model.train()
stop = False
ep = 0
while step < total and not stop:
    perm = np.random.default_rng([args.seed, ep]).permutation(len(tr))
    for b0 in range(0, len(perm), args.batch):
        bi = perm[b0:b0 + args.batch]
        ids, att = batch(K[bi], P[bi])
        y = torch.from_numpy(Y[bi]).to(dev)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(ids, att)
        loss = focal_loss(logits, y, 2.0)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        sched.step()
        step += 1
        timeout = (time.time() - t0) / 3600 > args.max_hours
        if step % args.eval_every == 0 or step == total or timeout:
            v = float(roc_auc_score(va.label, predict(va)))
            hist.append({"step": step, "rows_seen": step * args.batch, "val_auroc": v, "loss": float(loss.item()),
                         "elapsed_s": time.time() - t0})
            log(f"step {step}/{total} ({step * args.batch / (time.time() - t0):,.0f} rows/s) val AUROC {v:.4f} "
                f"(best {best:.4f})")
            if v > best:
                best, bad = v, 0
                torch.save(model.state_dict(), CK / "best.pt")
            else:
                bad += 1
            if bad >= args.patience or timeout:
                stop = True
                log("early stop" if not timeout else "wall-clock cap reached")
                break
        if step >= total:
            break
    ep += 1
train_s = time.time() - t0
model.load_state_dict(torch.load(CK / "best.pt", map_location=dev, weights_only=True))
pt = predict(te)
pn = predict(neg)
pkm, nk = per_kinase_mean(te, pt)
res = {"init": str(Path(args.init).relative_to(ROOT)), "train_rows_used": len(tr), "train_rows_full": n_full,
       "frac": args.frac, "epochs_cap": args.epochs, "steps": step, "rows_seen": step * args.batch,
       "train_s": train_s, "best_val_auroc": best, "history": hist,
       "test_auroc": float(roc_auc_score(te.label, pt)), "test_auprc": float(average_precision_score(te.label, pt)),
       "test_per_kinase_mean_auroc": pkm, "test_n_kinases": nk,
       "test_neg_only_FPR@0.5": float((pn >= 0.5).mean()),
       "reference": {"paper_fig4A_auroc": 0.885, "paper_per_kinase_mean_auroc": 0.860,
                     "released_weights_our_scorer_auroc": 0.8855, "released_weights_auprc": 0.8770,
                     "released_weights_FPR@0.5": 0.0230}}
(OUT / "fidelity_check.json").write_text(json.dumps(res, indent=1))
log(json.dumps({k: v for k, v in res.items() if k != "history"}, indent=1))
