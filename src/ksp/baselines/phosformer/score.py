"""Batched (kinase domain, 11-mer) scoring for any PhosformerClassifier (released or ours)."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
import torch

from .tokens import BOS, EOS, _LUT, encode_seq


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
    # batches: contiguous runs of one kinase, each split into chunks of <= min(batch_size, max_tokens // len)
    chg = np.r_[0, np.where(kin[1:] != kin[:-1])[0] + 1, len(kin)]
    bounds = []
    for g0, g1 in zip(chg[:-1], chg[1:]):
        bs = max(1, min(batch_size, max_tokens // (len(kinase_seq[kin[g0]]) + 15)))
        bounds += [(c, min(c + bs, g1)) for c in range(g0, g1, bs)]
    ktok = {k: encode_seq(kinase_seq[k]) for k in np.unique(kin)}
    for bi, (b0, b1) in enumerate(bounds):
        sl = slice(b0, b1)
        kt = ktok[kin[b0]]  # one kinase per batch (bounds never mix kinases)
        n, L = b1 - b0, len(kt)
        pb = np.frombuffer("".join(pep[sl]).encode("ascii"), dtype=np.uint8)
        if len(pb) != 11 * n:
            raise ValueError("all peptides must be 11-mers")
        ids = np.empty((n, 15 + L), dtype=np.int64)
        ids[:, 0] = BOS
        ids[:, 1:12] = _LUT[pb].reshape(n, 11)
        ids[:, 12:14] = EOS
        ids[:, 14:14 + L] = kt
        ids[:, -1] = EOS
        att = np.ones_like(ids)
        ids = torch.from_numpy(ids).to(device, non_blocking=True)
        # all rows share one kinase -> no padding; attention_mask=None lets SDPA use flash kernels (identical math:
        # in-peptide <pad> tokens have attention_mask 1 anyway)
        att = None if att.all() else torch.from_numpy(att).to(device, non_blocking=True)
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
            pr = model.prob(ids, att)
        out[order[sl]] = pr.float().cpu().numpy()
        done += len(ids)
        if log_every and bi % log_every == 0:
            print(f"  scored {done}/{len(order)} ({done / (time.time() - t0):,.0f} pairs/s)", flush=True)
    return out
