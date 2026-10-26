"""R5 time pilot: measure cross-encoder vs late-interaction inference cost (DESIGN §8.4-like, R5 brief).

Decides whether full joint cross-encoding of ~21.3M eval pairs/fold is feasible (<= ~1.5 GPU-h/fold).

  CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_pilot.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import torch  # noqa: E402

N_EVAL_PAIRS = 21_300_000      # eval candidate pairs per CV fold (R3 docs §3)


def bench_cross(model_id: str, L_kin: int, L_sub: int = 31, batch: int = 64, steps: int = 12, train: bool = False):
    from transformers import AutoModel
    m = AutoModel.from_pretrained(model_id, torch_dtype=torch.bfloat16).cuda()
    m.train(train)
    L = L_kin + L_sub + 3
    ids = torch.randint(4, 29, (batch, L), device="cuda")
    att = torch.ones_like(ids)
    if train:
        m = AutoModel.from_pretrained(model_id).cuda()
        m.train()
    def one():
        if train:
            with torch.autocast("cuda", dtype=torch.bfloat16):
                out = m(input_ids=ids, attention_mask=att).last_hidden_state.float().mean()
            out.backward()
            m.zero_grad(set_to_none=True)
        else:
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                m(input_ids=ids, attention_mask=att)
    for _ in range(3):
        one()
    torch.cuda.synchronize()
    t0 = time.time()
    for _ in range(steps):
        one()
    torch.cuda.synchronize()
    dt = (time.time() - t0) / steps
    pps = batch / dt
    peak = torch.cuda.max_memory_allocated() / 2**30
    del m
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    return {"model": model_id, "mode": "train" if train else "infer", "seq_len": L, "batch": batch,
            "s_per_batch": round(dt, 4), "pairs_per_s": round(pps, 1),
            "gpu_h_per_fold": round(N_EVAL_PAIRS / pps / 3600, 2), "peak_gb": round(peak, 2)}


def bench_sub_encoder(model_id: str, L_sub: int = 31, batch: int = 1024, steps: int = 10, train: bool = True):
    """Late interaction: the ESM trunk sees only the +-15 substrate window (31 tokens); kinase tokens are cached."""
    from transformers import AutoModel
    m = AutoModel.from_pretrained(model_id).cuda()
    m.train(train)
    ids = torch.randint(4, 29, (batch, L_sub), device="cuda")
    att = torch.ones_like(ids)
    opt = torch.optim.AdamW(m.parameters(), lr=1e-5)
    def one():
        if train:
            with torch.autocast("cuda", dtype=torch.bfloat16):
                out = m(input_ids=ids, attention_mask=att).last_hidden_state.float().mean()
            out.backward()
            opt.step()
            opt.zero_grad(set_to_none=True)
        else:
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                m(input_ids=ids, attention_mask=att)
    for _ in range(3):
        one()
    torch.cuda.synchronize()
    t0 = time.time()
    for _ in range(steps):
        one()
    torch.cuda.synchronize()
    dt = (time.time() - t0) / steps
    peak = torch.cuda.max_memory_allocated() / 2**30
    del m, opt
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    return {"model": model_id, "mode": "train" if train else "infer", "seq_len": L_sub, "batch": batch,
            "s_per_batch": round(dt, 4), "sites_per_s": round(batch / dt, 1), "peak_gb": round(peak, 2)}


def bench_maxsim(n_sites: int = 4096, K: int = 516, L_sub: int = 31, m_kin: int = 16, d: int = 64, steps: int = 10):
    """Cost of the late-interaction head itself: (B, L_sub, d) x (K, m, d) -> MaxSim -> (B, K)."""
    S = torch.randn(n_sites, L_sub, d, device="cuda", dtype=torch.bfloat16)
    Q = torch.randn(K, m_kin, d, device="cuda", dtype=torch.bfloat16)
    def one():
        with torch.no_grad():
            sim = torch.einsum("bld,kmd->bkml", S, Q)
            sim.amax(-1).sum(-1)
    for _ in range(2):
        one()
    torch.cuda.synchronize()
    t0 = time.time()
    for _ in range(steps):
        one()
    torch.cuda.synchronize()
    dt = (time.time() - t0) / steps
    pps = n_sites * K / dt
    peak = torch.cuda.max_memory_allocated() / 2**30
    torch.cuda.empty_cache()
    return {"op": "maxsim", "sites_per_batch": n_sites, "K": K, "m_kin": m_kin, "d": d,
            "s_per_batch": round(dt, 4), "pairs_per_s": round(pps, 0),
            "gpu_h_per_fold": round(N_EVAL_PAIRS / pps / 3600, 4), "peak_gb": round(peak, 2)}


if __name__ == "__main__":
    out = {"gpu": torch.cuda.get_device_name(0), "n_eval_pairs_per_fold": N_EVAL_PAIRS, "results": []}
    for mid, Lk in [("facebook/esm2_t12_35M_UR50D", 300), ("facebook/esm2_t30_150M_UR50D", 300)]:
        for bs in (64, 256):
            r = bench_cross(mid, Lk, batch=bs)
            out["results"].append(r)
            print(json.dumps(r), flush=True)
    for mid, tb in [("facebook/esm2_t12_35M_UR50D", 1024), ("facebook/esm2_t30_150M_UR50D", 256)]:
        for r in (bench_sub_encoder(mid, batch=tb), bench_sub_encoder(mid, train=False, batch=2048)):
            out["results"].append(r)
            print(json.dumps(r), flush=True)
    for m_kin in (8, 16, 32):
        r = bench_maxsim(m_kin=m_kin)
        out["results"].append(r)
        print(json.dumps(r), flush=True)
    d = ROOT / "results" / "models" / "r5"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pilot_timing.json").write_text(json.dumps(out, indent=1))
    print("written", d / "pilot_timing.json")
