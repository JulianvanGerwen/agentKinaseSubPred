"""B3 Phosformer MLM pretraining (paper §2.3.1), resumable, single GPU.

Recipe (paper): RoBERTa-style encoder 6x12x768, MLM 15% masking, AdamW lr 5e-5, 50 epochs, batch 48/device x 8
devices = 384 sequences per optimiser step, corpus = 295,320 kinase domains + 5,338 substrate proteins.
Our choices where the paper is silent (documented in results/baselines/phosformer/recipe.md):
- 80/10/10 mask/random/keep replacement (BERT/RoBERTa + HF DataCollator default);
- linear LR decay to 0, no warmup, weight decay 0 (HF Trainer defaults, the most likely setup);
- sequences > 1022 residues are randomly cropped each epoch to 1022 (+2 specials = 1024 tokens);
- 1,000 sequences held out for MLM validation loss/accuracy monitoring;
- bf16 autocast, length-bucketed batches split into micro-batches by a token budget (gradient accumulation;
  loss = sum over masked tokens / masked tokens in the full 384-seq batch).

Usage:
  CUDA_VISIBLE_DEVICES=0 python -I scripts/baselines/phosformer_pretrain.py --run-dir cache/phosformer/pretrain \
      --epochs 50 [--pilot-steps 60] [--max-steps N]
Re-running the same command resumes from <run-dir>/ckpt_last.pt.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import torch
from transformers import RobertaForMaskedLM

from ksp.baselines.phosformer.corpus import load_corpus
from ksp.baselines.phosformer.model import phosformer_config
from ksp.baselines.phosformer.tokens import AA20_IDS, BOS, EOS, MASK, PAD

p = argparse.ArgumentParser()
p.add_argument("--run-dir", default=str(ROOT / "cache/phosformer/pretrain"))
p.add_argument("--epochs", type=float, default=50.0)
p.add_argument("--batch-seqs", type=int, default=384)
p.add_argument("--token-budget", type=int, default=49152, help="max padded tokens per micro-batch")
p.add_argument("--max-len", type=int, default=1022, help="max residues per sequence (random crop)")
p.add_argument("--lr", type=float, default=5e-5)
p.add_argument("--mask-prob", type=float, default=0.15)
p.add_argument("--n-val", type=int, default=1000)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--ckpt-every", type=int, default=500, help="optimizer steps between ckpt_last saves")
p.add_argument("--keep-every-epochs", type=int, default=5)
p.add_argument("--val-every", type=int, default=1000)
p.add_argument("--pilot-steps", type=int, default=0, help="run N steps, report throughput + ETA, exit (no ckpt)")
p.add_argument("--compile", action="store_true")
p.add_argument("--threads", type=int, default=4)
args = p.parse_args()

torch.set_num_threads(args.threads)
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
dev = torch.device("cuda")
run_dir = Path(args.run_dir)
run_dir.mkdir(parents=True, exist_ok=True)
log_f = open(run_dir / ("pilot.log" if args.pilot_steps else "train.log"), "a")


def log(msg):
    s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(s, flush=True)
    log_f.write(s + "\n")
    log_f.flush()


# ---------------- data ---------------------------------------------------------------------------------------
tokens, offsets, source = load_corpus()
n_all = len(offsets) - 1
lens_all = np.diff(offsets)
rng0 = np.random.default_rng(args.seed)
perm0 = rng0.permutation(n_all)
val_idx = np.sort(perm0[: args.n_val])
train_idx = np.sort(perm0[args.n_val:])
n_train = len(train_idx)
steps_per_epoch = math.ceil(n_train / args.batch_seqs)
total_steps = int(round(args.epochs * steps_per_epoch))
log(f"corpus: {n_all} seqs ({int((source == 0).sum())} kinase, {int((source == 1).sum())} substrate), "
    f"{int(lens_all.sum())} residues; train {n_train}, val {len(val_idx)}; steps/epoch {steps_per_epoch}; "
    f"epochs {args.epochs} -> total steps {total_steps}")


def epoch_batches(epoch: int):
    """Deterministic per-epoch list of global batches (arrays of seq indices) + crop offsets."""
    r = np.random.default_rng([args.seed, epoch])
    order = train_idx[r.permutation(n_train)]
    L = np.minimum(lens_all[order], args.max_len)
    chunk = args.batch_seqs * 64  # bucket within mega-chunks to keep randomness
    batches = []
    for c0 in range(0, n_train, chunk):
        sub = order[c0 : c0 + chunk]
        sub = sub[np.argsort(L[c0 : c0 + chunk], kind="stable")]
        for b0 in range(0, len(sub), args.batch_seqs):
            batches.append(sub[b0 : b0 + args.batch_seqs])
    bo = r.permutation(len(batches))
    batches = [batches[i] for i in bo]
    crop_seed = int(r.integers(1 << 62))
    return batches, crop_seed


def get_seq(i: int, r: np.random.Generator) -> np.ndarray:
    a, b = offsets[i], offsets[i + 1]
    n = b - a
    if n > args.max_len:
        s = int(r.integers(0, n - args.max_len + 1))
        a, b = a + s, a + s + args.max_len
    return np.asarray(tokens[a:b], dtype=np.int64)


def make_micro(seq_list, r: np.random.Generator):
    """seq_list sorted roughly by length; split into micro-batches under token budget; apply MLM masking."""
    micro, cur, cur_max = [], [], 0
    for s in sorted(seq_list, key=len):
        L = len(s) + 2
        if cur and max(cur_max, L) * (len(cur) + 1) > args.token_budget:
            micro.append(cur)
            cur, cur_max = [], 0
        cur.append(s)
        cur_max = max(cur_max, L)
    if cur:
        micro.append(cur)
    out = []
    for m in micro:
        n = max(len(s) for s in m) + 2
        ids = np.full((len(m), n), PAD, dtype=np.int64)
        for j, s in enumerate(m):
            ids[j, 0] = BOS
            ids[j, 1 : 1 + len(s)] = s
            ids[j, 1 + len(s)] = EOS
        att = (ids != PAD).astype(np.int64)
        cand = (ids != PAD) & (ids != BOS) & (ids != EOS)
        sel = cand & (r.random(ids.shape) < args.mask_prob)
        labels = np.where(sel, ids, -100)
        u = r.random(ids.shape)
        inp = ids.copy()
        inp[sel & (u < 0.8)] = MASK
        rnd = sel & (u >= 0.8) & (u < 0.9)
        inp[rnd] = AA20_IDS[r.integers(0, 20, size=int(rnd.sum()))]
        out.append((torch.from_numpy(inp), torch.from_numpy(att), torch.from_numpy(labels)))
    return out


# ---------------- model --------------------------------------------------------------------------------------
cfg = phosformer_config()
cfg._attn_implementation = "sdpa"
torch.manual_seed(args.seed)
model = RobertaForMaskedLM(cfg).to(dev)
n_params = sum(p.numel() for p in model.parameters())
opt = torch.optim.AdamW(model.parameters(), lr=args.lr, betas=(0.9, 0.999), eps=1e-8, weight_decay=0.0)
sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: max(0.0, 1.0 - s / max(1, total_steps)))
fwd = torch.compile(model, dynamic=True) if args.compile else model
log(f"model params {n_params / 1e6:.1f}M; config {cfg.num_hidden_layers}L/{cfg.num_attention_heads}H/"
    f"{cfg.hidden_size}d; lr {args.lr}; batch {args.batch_seqs} seqs; token budget {args.token_budget}")

step = 0
ckpt_last = run_dir / "ckpt_last.pt"
if ckpt_last.exists() and not args.pilot_steps:
    ck = torch.load(ckpt_last, map_location="cpu", weights_only=False)
    if ck["total_steps"] != total_steps:
        raise SystemExit(f"checkpoint total_steps {ck['total_steps']} != requested {total_steps}; refuse to "
                         f"resume with a different LR schedule")
    model.load_state_dict(ck["model"])
    opt.load_state_dict(ck["opt"])
    sched.load_state_dict(ck["sched"])
    step = ck["step"]
    log(f"resumed from step {step}")


def save(path: Path, extra=None):
    tmp = path.with_suffix(".tmp")
    d = {"model": model.state_dict(), "opt": opt.state_dict(), "sched": sched.state_dict(), "step": step,
         "total_steps": total_steps, "args": vars(args), "config": cfg.to_dict()}
    if extra:
        d.update(extra)
    torch.save(d, tmp)
    os.replace(tmp, path)


@torch.no_grad()
def validate():
    model.eval()
    r = np.random.default_rng(12345)
    tot_loss, tot_n, tot_corr = 0.0, 0, 0
    seqs = [get_seq(i, r) for i in val_idx]
    for inp, att, lab in make_micro(seqs, r):
        inp, att, lab = inp.to(dev), att.to(dev), lab.to(dev)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(input_ids=inp, attention_mask=att).logits
        m = lab != -100
        lg = logits[m].float()
        tot_loss += torch.nn.functional.cross_entropy(lg, lab[m], reduction="sum").item()
        tot_corr += (lg.argmax(-1) == lab[m]).sum().item()
        tot_n += int(m.sum())
    model.train()
    return tot_loss / tot_n, tot_corr / tot_n


# ---------------- train loop ---------------------------------------------------------------------------------
model.train()
t0 = time.time()
tok_seen, res_seen = 0, 0
target = total_steps if not args.pilot_steps else min(total_steps, step + args.pilot_steps)
cur_epoch, batches, crop_rng = -1, None, None
while step < target:
    ep, bi = divmod(step, steps_per_epoch)
    if ep != cur_epoch:
        batches, crop_seed = epoch_batches(ep)
        cur_epoch = ep
    if bi >= len(batches):  # should not happen (steps_per_epoch == len(batches))
        raise RuntimeError("batch index overflow")
    r = np.random.default_rng([args.seed, ep, bi, 7])
    seqs = [get_seq(int(i), r) for i in batches[bi]]
    micro = make_micro(seqs, r)
    n_masked = sum(int((lab != -100).sum()) for _, _, lab in micro)
    opt.zero_grad(set_to_none=True)
    loss_acc = 0.0
    for inp, att, lab in micro:
        inp, att, lab = inp.to(dev, non_blocking=True), att.to(dev, non_blocking=True), lab.to(dev, non_blocking=True)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = fwd(input_ids=inp, attention_mask=att).logits
        m = lab != -100
        loss = torch.nn.functional.cross_entropy(logits[m].float(), lab[m], reduction="sum") / n_masked
        loss.backward()
        loss_acc += loss.item()
        tok_seen += inp.numel()
        res_seen += int(att.sum())
    gn = torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    opt.step()
    sched.step()
    step += 1
    if step % 50 == 0 or args.pilot_steps:
        el = time.time() - t0
        sp = (step - (target - args.pilot_steps if args.pilot_steps else 0))
        msg = (f"step {step}/{total_steps} ep {step / steps_per_epoch:.2f} loss {loss_acc:.4f} gnorm {gn:.2f} "
               f"lr {sched.get_last_lr()[0]:.2e} | {res_seen / el:,.0f} real tok/s, {tok_seen / el:,.0f} padded tok/s "
               f"| {len(micro)} micro")
        if step % 50 == 0 or step == target:
            log(msg)
    if not args.pilot_steps:
        if step % args.val_every == 0 or step == total_steps:
            vl, va = validate()
            log(f"VAL step {step} ep {step / steps_per_epoch:.2f}: mlm loss {vl:.4f} acc {va:.4f}")
        if step % args.ckpt_every == 0 or step == total_steps:
            save(ckpt_last)
        if step % (steps_per_epoch * args.keep_every_epochs) == 0 or step == total_steps:
            save(run_dir / f"ckpt_ep{step // steps_per_epoch:03d}.pt")

el = time.time() - t0
if args.pilot_steps:
    n = args.pilot_steps
    warm = 0
    sec_per_step = el / n
    eta_h = sec_per_step * total_steps / 3600
    res = {"pilot_steps": n, "elapsed_s": el, "sec_per_step": sec_per_step, "real_tok_per_s": res_seen / el,
           "padded_tok_per_s": tok_seen / el, "total_steps_for_epochs": total_steps, "epochs": args.epochs,
           "eta_hours_full": eta_h, "compile": args.compile, "token_budget": args.token_budget,
           "max_mem_GB": torch.cuda.max_memory_allocated() / 1e9}
    log("PILOT " + json.dumps(res))
    (run_dir / "pilot.json").write_text(json.dumps(res, indent=1))
else:
    vl, va = validate()
    final = {"step": step, "total_steps": total_steps, "epochs": args.epochs, "val_mlm_loss": vl, "val_mlm_acc": va}
    (run_dir / "final.json").write_text(json.dumps(final, indent=1))
    torch.save({"model": model.state_dict(), "config": cfg.to_dict(), "step": step}, run_dir / "mlm_final.pt")
    log("DONE " + json.dumps(final))
