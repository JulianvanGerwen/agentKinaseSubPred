"""Batched (kinase domain, 11-mer) scoring for any PhosformerClassifier (released or ours)."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
import torch

from .tokens import collate_pairs, encode_pair


@torch.no_grad()
def score_pairs(model, kinase_seq: dict, pairs: pd.DataFrame, batch_size: int = 512, device="cuda",
                amp: bool = True, log_every: int = 0, max_tokens: int = 200_000) -> np.ndarray:
    """pairs: DataFrame[kinase_acc, peptide]; kinase_seq: acc -> domain sequence. Returns P(phospho) per row.
    Rows are processed grouped by kinase (equal lengths within a batch -> no padding waste)."""
    model.eval()
    out = np.empty(len(pairs), dtype=np.float32)
    order = np.argsort(pairs.kinase_acc.values, kind="stable")
    kin = pairs.kinase_acc.values[order]
    pep = pairs.peptide.values[order]
    t0, done = time.time(), 0
    klen = np.array([len(kinase_seq[k]) + 15 for k in kin])
    bounds, b0 = [], 0
    while b0 < len(order):  # batches of <= batch_size rows and <= max_tokens padded tokens
        b1, mx = b0, 0
        while b1 < len(order) and b1 - b0 < batch_size and max(mx, klen[b1]) * (b1 - b0 + 1) <= max_tokens:
            mx = max(mx, klen[b1])
            b1 += 1
        b1 = max(b1, b0 + 1)
        bounds.append((b0, b1))
        b0 = b1
    for bi, (b0, b1) in enumerate(bounds):
        sl = slice(b0, b1)
        ids, att = collate_pairs([encode_pair(p, kinase_seq[k]) for k, p in zip(kin[sl], pep[sl])])
        ids = torch.from_numpy(ids).to(device, non_blocking=True)
        att = torch.from_numpy(att).to(device, non_blocking=True)
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
            pr = model.prob(ids, att)
        out[order[sl]] = pr.float().cpu().numpy()
        done += len(ids)
        if log_every and bi % log_every == 0:
            print(f"  scored {done}/{len(order)} ({done / (time.time() - t0):,.0f} pairs/s)", flush=True)
    return out
