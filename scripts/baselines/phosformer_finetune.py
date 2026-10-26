"""B3a Phosformer fine-tuning on one DEV CV fold + out-of-fold scoring (paper §2.2, §2.3.2; recipe.md §3).

Per fold f (outer eval = DEV proteins of fold f; never used for training, early stopping or ratio choice):
  positives   : ksp.crossfit.outer_split(f, 'L2') training labels (U_ms sites, in class).
  inner valid : ksp.splits.train_valid_split_inner over ALL U_ms sites of training proteins (whole clusters, 20%);
                used for early stopping and for choosing the negative ratio.
  masks       : ksp.labels.NegativeMasks built from training-protein label rows only (load_dev('all') restricted to
                outer-training proteins) -> no held-out label can shape training negatives.
  negatives   : per sampled positive (k,s): m hard + n easy, residue class of k, minus excluded(k) =
                KP(k) U Hamming<=2 U Par(k)-labelled.
                hard = NegativeMasks.hard_candidates(k) (human-native PSP label to a kinase outside {k} U Par(k));
                easy = U_ms sites with no kinase label of any source (paper: "no evidence of being phosphorylated by
                any kinase"; our universe is observed sites, so 'easy' = unlabelled).
  augment     : upsampling (>=50 positives -> fold max, <50 -> 50), domain-boundary shift +-5 (each end, independent),
                5% kinase-token <mask>; resampled every epoch.
  train       : focal loss gamma 2, AdamW lr 2e-5, batch 64, bf16; validate every 1/4 epoch on inner valid
                (macro AUROC positives vs hard negatives over kinases with >=3 valid positives); patience 5; max 3 ep.
  OOF scoring : set A = (positive U E-hard sites of fold f) x all residue-matched kinases (M1/M2/M4/M5 complete);
                set B (--eall) = remaining E-all sites x kinases with >=1 fold-f positive.
Usage:
  CUDA_VISIBLE_DEVICES=0 python -I scripts/baselines/phosformer_finetune.py --fold 0 --neg-ratio 16 \
      --init cache/phosformer/pretrain/mlm_final.pt [--eall] [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
import torch

from ksp import eval as E
from ksp import index, labels, splits
from ksp.baselines.phosformer.inputs import canonical_seqs, kinase_domains, window15_to_11mer
from ksp.baselines.phosformer.model import PhosformerClassifier, focal_loss, phosformer_config
from ksp.baselines.phosformer.score import score_pairs
from ksp.baselines.phosformer.tokens import MASK, collate_pairs, encode_pair, encode_seq, BOS, EOS
from ksp.crossfit import outer_split
from ksp.manifest import RunManifest

ap = argparse.ArgumentParser()
ap.add_argument("--fold", type=int, required=True)
ap.add_argument("--neg-ratio", type=int, default=16, help="negatives per positive (total)")
ap.add_argument("--hard-frac", type=float, default=0.5, help="fraction of negatives that are hard")
ap.add_argument("--init", required=True, help="MLM checkpoint (.pt with 'model'), 'released' or 'random'")
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--lr", type=float, default=2e-5)
ap.add_argument("--gamma", type=float, default=2.0)
ap.add_argument("--batch", type=int, default=64)
ap.add_argument("--max-epochs", type=float, default=3.0)
ap.add_argument("--evals-per-epoch", type=int, default=4)
ap.add_argument("--patience", type=int, default=5)
ap.add_argument("--shift", type=int, default=5)
ap.add_argument("--mask-frac", type=float, default=0.05)
ap.add_argument("--no-aug", action="store_true", help="ablation: no upsampling/shift/mask")
ap.add_argument("--eall", action="store_true", help="also score E-all (set B)")
ap.add_argument("--max-train-rows", type=int, default=0, help="debug: cap rows per epoch")
ap.add_argument("--name", default=None)
ap.add_argument("--threads", type=int, default=6)
args = ap.parse_args()
torch.set_num_threads(args.threads)
dev = "cuda"
rng = np.random.default_rng([args.seed, args.fold, args.neg_ratio])
torch.manual_seed(args.seed)
name = args.name or f"b3a_phosformer_f{args.fold}_r{args.neg_ratio}" + ("_noaug" if args.no_aug else "")
init_tag = args.init if args.init in ("released", "random") else Path(args.init).name
t_start = time.time()


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ------------------------------------------------------------------------------------------------ inputs
K = index.load_kinases()
S = index.load_sites()
spos = pd.Series(np.arange(len(S)), index=S.site_id)
kd, kfb = kinase_domains(K)
full = canonical_seqs(K.uniprot_acc)
kacc = K.uniprot_acc.values
kidx = {a: i for i, a in enumerate(kacc)}
kcls = K.residue_class.values
dom_a = kd.dom_start.values.astype(int) - 1          # 0-based start
dom_b = kd.dom_end.values.astype(int)                # exclusive end
full_tok = [encode_seq(full[a]) for a in kacc]
dom_seq = dict(zip(kd.uniprot_acc, kd.domain_seq))
pep = np.array([window15_to_11mer(w) for w in S.window15.values], dtype=object)
pep_tok = np.stack([encode_seq(p) for p in pep])     # (nS, 11)
res = S.res.values
is_st = np.isin(res, ["S", "T"])
is_y = res == "Y"


def class_mask(cls):
    return np.ones(len(S), bool) if cls == "STY" else (is_st if cls == "ST" else is_y)


# ------------------------------------------------------------------------------------------------ fold data
tr_lab, train_prot, eval_prot = outer_split(args.fold, "L2")
tr_lab = tr_lab[tr_lab.kinase_acc.isin(kidx) & tr_lab.site_id.isin(spos.index)]
ums = S.in_U_ms.values.astype(bool)
in_train_prot = S.acc.isin(train_prot).values
tr_sites_all = S.site_id.values[in_train_prot & ums]
val_mask = splits.train_valid_split_inner(tr_sites_all, valid_frac=0.2, seed=args.seed)
val_site = np.zeros(len(S), bool)
val_site[spos[tr_sites_all[val_mask]].values] = True
fit_site = np.zeros(len(S), bool)
fit_site[spos[tr_sites_all[~val_mask]].values] = True
all_dev = labels.load_dev("all")
mask_pairs = all_dev[all_dev.acc.isin(train_prot)]   # training-protein label rows only
log(f"fold {args.fold}: train L2 pairs {len(tr_lab)}, fit sites {fit_site.sum()}, valid sites {val_site.sum()}")
masks = labels.NegativeMasks(mask_pairs, paralogs=labels.paralog_sets())
labelled_any = np.zeros(len(S), bool)
labelled_any[spos.reindex(mask_pairs.site_id.unique()).dropna().values.astype(int)] = True
tr_lab = tr_lab.assign(s=spos[tr_lab.site_id].values, k=tr_lab.kinase_acc.map(kidx).values)
fit_pos = tr_lab[fit_site[tr_lab.s.values]]
val_pos = tr_lab[val_site[tr_lab.s.values]]

pools = {}  # k -> (hard_fit, easy_fit)
fallbacks = {"kinase_domain_full_seq": len(kfb), "hard_pool_empty_used_easy": 0}
easy_base = fit_site & ~labelled_any
for k in fit_pos.k.unique():
    a = kacc[k]
    excl = np.zeros(len(S), bool)
    excl[masks.excluded(a)] = True
    cm = class_mask(kcls[k])
    hard = np.zeros(len(S), bool)
    hard[masks.hard_candidates(a)] = True
    pools[k] = (np.where(hard & fit_site & cm & ~excl)[0], np.where(easy_base & cm & ~excl)[0])

counts = fit_pos.groupby("k").size()
cmax = int(counts.max())
log(f"fit positives {len(fit_pos)} over {len(counts)} kinases (max {cmax}); valid positives {len(val_pos)}")


def epoch_rows(r: np.random.Generator):
    """(k, s, y) arrays for one epoch: upsampled positives + m hard + n easy negatives per positive."""
    ks, ss, ys = [], [], []
    n_neg = args.neg_ratio
    m_hard = int(round(n_neg * args.hard_frac))
    n_easy = n_neg - m_hard
    for k, g in fit_pos.groupby("k"):
        s = g.s.values
        if args.no_aug:
            tgt = len(s)
        else:
            tgt = cmax if len(s) >= 50 else 50
        reps = np.concatenate([r.permutation(s) for _ in range(math.ceil(tgt / len(s)))])[:tgt]
        hp, ep = pools[k]
        nh = m_hard * tgt
        ne = n_easy * tgt
        if len(hp) == 0:
            fallbacks["hard_pool_empty_used_easy"] += nh
            ne, nh = ne + nh, 0
        neg = np.concatenate([r.choice(hp, nh) if nh else np.zeros(0, int), r.choice(ep, ne)])
        ks.append(np.full(tgt + len(neg), k))
        ss.append(np.concatenate([reps, neg]))
        ys.append(np.r_[np.ones(tgt, int), np.zeros(len(neg), int)])
    k, s, y = np.concatenate(ks), np.concatenate(ss), np.concatenate(ys)
    p = r.permutation(len(k))
    if args.max_train_rows:
        p = p[: args.max_train_rows]
    return k[p], s[p], y[p]


def make_batch(k, s, r: np.random.Generator, augment: bool):
    seqs = []
    for ki, si in zip(k, s):
        a, b = dom_a[ki], dom_b[ki]
        ft = full_tok[ki]
        if augment:
            a = int(np.clip(a + r.integers(-args.shift, args.shift + 1), 0, len(ft) - 1))
            b = int(np.clip(b + r.integers(-args.shift, args.shift + 1), a + 1, len(ft)))
        kt = ft[a:b].copy()
        if augment and args.mask_frac > 0:
            kt[r.random(len(kt)) < args.mask_frac] = MASK
        seqs.append(np.concatenate([[BOS], pep_tok[si], [EOS, EOS], kt, [EOS]]))
    ids, att = collate_pairs(seqs)
    return torch.from_numpy(ids).to(dev, non_blocking=True), torch.from_numpy(att).to(dev, non_blocking=True)


# fixed inner-validation set: valid positives + up to 200 valid hard negatives per kinase (E-hard analogue)
vk, vs, vy = [], [], []
val_kin = []
for k, g in val_pos.groupby("k"):
    a = kacc[k]
    excl = np.zeros(len(S), bool)
    excl[masks.excluded(a)] = True
    hard = np.zeros(len(S), bool)
    hard[masks.hard_candidates(a)] = True
    hn = np.where(hard & val_site & class_mask(kcls[k]) & ~excl)[0]
    if len(g) < 3 or len(hn) == 0:
        continue
    hn = rng.choice(hn, min(200, len(hn)), replace=False)
    vk += [k] * (len(g) + len(hn))
    vs += list(g.s.values) + list(hn)
    vy += [1] * len(g) + [0] * len(hn)
    val_kin.append(k)
vk, vs, vy = np.array(vk), np.array(vs), np.array(vy)
log(f"inner-valid set: {len(vk)} pairs, {len(val_kin)} kinases")


@torch.no_grad()
def validate(model):
    model.eval()
    pr = np.empty(len(vk), np.float32)
    for b0 in range(0, len(vk), 512):
        ids, att = make_batch(vk[b0:b0 + 512], vs[b0:b0 + 512], None, augment=False)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            pr[b0:b0 + 512] = model.prob(ids, att).float().cpu().numpy()
    model.train()
    au = [E.auroc(pr[(vk == k) & (vy == 1)], pr[(vk == k) & (vy == 0)]) for k in val_kin]
    return float(np.nanmean(au))


# ------------------------------------------------------------------------------------------------ model
cfg = phosformer_config()
if args.init == "released":
    model = PhosformerClassifier.from_released(ROOT / "external/phosformer/hf_phosformer")
elif args.init == "random":
    model = PhosformerClassifier(cfg)
else:
    ck = torch.load(args.init, map_location="cpu", weights_only=False)
    model = PhosformerClassifier.from_mlm(ck["model"], cfg)
model = model.to(dev)
opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.0)
k0, s0, y0 = epoch_rows(rng)
steps_per_epoch = math.ceil(len(k0) / args.batch)
total_steps = int(args.max_epochs * steps_per_epoch)
eval_every = max(1, steps_per_epoch // args.evals_per_epoch)
sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda st: max(0.0, 1.0 - st / total_steps))
log(f"epoch rows {len(k0)} (pos {int(y0.sum())}), steps/epoch {steps_per_epoch}, eval every {eval_every}, "
    f"max steps {total_steps}")

ckdir = ROOT / "cache/phosformer/finetune" / f"{name}_s{args.seed}"
ckdir.mkdir(parents=True, exist_ok=True)
best, bad, step, hist = -1.0, 0, 0, []
v0 = validate(model)
hist.append({"step": 0, "val_macro_auroc_hard": v0})
log(f"init valid macro AUROC(hard) {v0:.4f}")
model.train()
t0 = time.time()
ep_rows, ep = (k0, s0, y0), 0
stop = False
while step < total_steps and not stop:
    kk, ss, yy = ep_rows
    for b0 in range(0, len(kk), args.batch):
        ids, att = make_batch(kk[b0:b0 + args.batch], ss[b0:b0 + args.batch], rng, augment=not args.no_aug)
        y = torch.from_numpy(yy[b0:b0 + args.batch]).to(dev)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(ids, att)
        loss = focal_loss(logits, y, args.gamma)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        sched.step()
        step += 1
        if step % 200 == 0:
            log(f"step {step}/{total_steps} ep {step / steps_per_epoch:.2f} loss {loss.item():.4f} "
                f"{step * args.batch / (time.time() - t0):,.0f} rows/s")
        if step % eval_every == 0 or step == total_steps:
            v = validate(model)
            hist.append({"step": step, "epoch": step / steps_per_epoch, "val_macro_auroc_hard": v,
                         "train_loss": float(loss.item())})
            log(f"VALID step {step}: macro AUROC(hard) {v:.4f} (best {best:.4f})")
            if v > best:
                best, bad = v, 0
                torch.save(model.state_dict(), ckdir / "best.pt")
            else:
                bad += 1
                if bad >= args.patience:
                    log("early stop")
                    stop = True
                    break
        if step >= total_steps:
            break
    ep += 1
    if not stop and step < total_steps:
        ep_rows = epoch_rows(rng)
train_s = time.time() - t0
model.load_state_dict(torch.load(ckdir / "best.pt", map_location=dev, weights_only=True))
model.eval()

# ------------------------------------------------------------------------------------------------ OOF scoring
data = E.prepare("cv", fold=args.fold, with_canaries=False)
uni = data.universe
hard_u = data.masks.hard_idx[uni[data.masks.hard_idx]]
pos_u = np.unique(np.concatenate([v for v in data.positives.values()])) if data.positives else np.zeros(0, int)
setA = np.zeros(len(S), bool)
setA[np.union1d(hard_u, pos_u)] = True
rows = []
for i, a in enumerate(kacc):
    sel = np.where(setA & class_mask(kcls[i]))[0]
    rows.append(pd.DataFrame({"kinase_acc": a, "site_idx": sel}))
if args.eall:
    for a in data.positives:
        i = kidx[a]
        sel = np.where(uni & ~setA & class_mask(kcls[i]))[0]
        rows.append(pd.DataFrame({"kinase_acc": a, "site_idx": sel}))
pairs = pd.concat(rows, ignore_index=True)
pairs["peptide"] = pep[pairs.site_idx.values]
log(f"OOF scoring {len(pairs)} pairs (setA {int(setA.sum())} sites; eall={args.eall})")
t1 = time.time()
model.roberta = torch.compile(model.roberta, dynamic=True)
pairs["score"] = score_pairs(model, dom_seq, pairs, batch_size=1024)
score_s = time.time() - t1
preds = pd.DataFrame({"kinase_acc": pairs.kinase_acc.values, "site_id": S.site_id.values[pairs.site_idx.values],
                      "score": pairs.score.values.astype(np.float64)})
res_ = E.evaluate(preds, data, missing="drop")
cfg_rec = {**vars(args), "init_tag": init_tag, "n_fit_pos": int(len(fit_pos)), "fold_max_count": cmax,
           "epoch_rows": int(len(k0)), "steps_per_epoch": steps_per_epoch, "best_val": best, "history": hist,
           "train_s": train_s, "score_s": score_s, "eall_complete": bool(args.eall),
           "pairs_per_s": len(pairs) / score_s}
inputs = [ROOT / "results/splits/cv_folds.tsv", ROOT / "cache/index/kinases.parquet",
          ROOT / "cache/index/sites.parquet"]
if args.init not in ("released", "random"):
    inputs.append(Path(args.init))
with RunManifest(name, config=cfg_rec, seed=args.seed, inputs=inputs, split=f"cv_fold{args.fold}",
                 notes="B3a Phosformer reimplementation; OOF predictions for CV fold; metrics with missing='drop' "
                       "(E-all metrics valid only if eall_complete)") as run:
    for kf, v in fallbacks.items():
        run.count_fallback(kf, v)
    run.save_predictions(preds)
    run.save_metrics(res_)
    run.set("endpoints", res_.endpoints)
log(f"done: run {run.run_id}; endpoints {json.dumps({k: v for k, v in res_.endpoints.items() if 'M' in k}, default=str)}")
print(json.dumps({"run_id": run.run_id, "best_val": best, "endpoints": res_.endpoints}, default=str))
