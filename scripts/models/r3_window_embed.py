"""R3 window embeddings: frozen ESM-2 650M per-residue embeddings (full-protein context) for every residue within
+-FLANK of a U_ms site (DESIGN §8.4). The struct+ESM cache (cache/esm2_650m) stores only S/T/Y residues, so
non-S/T/Y window positions are missing there; this script recomputes them with the SAME settings
(ksp.features.esm_embed.load_model / embed_sequences: HF facebook/esm2_t33_650M_UR50D, last_hidden_state, fp16,
1022-aa windows stride 511, overlaps averaged).

Outputs (cache/r3/):
  win_emb_<w>.npy            float16 (n_res_w, 1280) per worker
  win_res_<w>.parquet        acc, pos (1-based), aa, row   (row into win_emb_<w>.npy)
  (merge) win_emb.npy, win_res.parquet (acc,pos,aa,row), win_sites.parquet (site_id), win_rows.npy
          (n_sites, 2*FLANK+1) int64 row ids into win_emb.npy, -1 = beyond protein terminus; win_sites.site_ids.npy
  win_embed_meta.json        timings + centre-residue consistency check against cache/esm2_650m

Usage:
  CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_window_embed.py run --worker 0 --n-workers 2
  CUDA_VISIBLE_DEVICES=3 python -I scripts/models/r3_window_embed.py run --worker 1 --n-workers 2
  python -I scripts/models/r3_window_embed.py merge
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import index  # noqa: E402
from ksp.features import esm_embed  # noqa: E402

OUT = ROOT / "cache" / "r3"
FLANK = 10
MODEL = "facebook/esm2_t33_650M_UR50D"


def ums_sites() -> pd.DataFrame:
    s = index.load_sites()
    return s[s.in_U_ms][["site_id", "acc", "res", "pos"]].reset_index(drop=True)


def run(args):
    import torch
    OUT.mkdir(parents=True, exist_ok=True)
    s = ums_sites()
    prot = esm_embed.canonical_swissprot()
    prot = prot[prot.acc.isin(set(s.acc))].copy()
    missing = set(s.acc) - set(prot.acc)
    if missing:
        raise KeyError(f"{len(missing)} U_ms proteins missing from canonical Swiss-Prot, e.g. {sorted(missing)[:5]}")
    prot["L"] = prot.seq.str.len()
    prot = prot.sort_values(["L", "acc"], ascending=[False, True]).reset_index(drop=True)
    mine = prot.iloc[args.worker::args.n_workers].reset_index(drop=True)
    pos_by_acc = s.groupby("acc").pos.agg(lambda x: np.asarray(x, dtype=np.int64)).to_dict()
    tok, model = esm_embed.load_model(MODEL, "fp16")
    t0 = time.time()
    mats, rows = [], []
    nrow = 0
    off = np.arange(-FLANK, FLANK + 1)
    for g0 in range(0, len(mine), args.group):
        g = mine.iloc[g0:g0 + args.group]
        embs = esm_embed.embed_sequences(g.seq.tolist(), tok, model, stride=511, token_budget=16384)
        for (acc, seq), E in zip(g[["acc", "seq"]].itertuples(index=False), embs):
            p = pos_by_acc[acc]
            for pp in p:
                if seq[pp - 1] not in "STY":
                    raise ValueError(f"{acc}_{pp}: residue {seq[pp - 1]} is not S/T/Y")
            keep = np.unique((p[:, None] + off[None]).ravel())
            keep = keep[(keep >= 1) & (keep <= len(seq))]
            mats.append(E[keep - 1].astype(np.float16))
            rows.append(pd.DataFrame({"acc": acc, "pos": keep, "aa": [seq[k - 1] for k in keep],
                                      "row": np.arange(nrow, nrow + len(keep))}))
            nrow += len(keep)
        el = time.time() - t0
        print(f"[w{args.worker}] {g0 + len(g)}/{len(mine)} proteins, {nrow} residues kept, {el / 60:.1f} min",
              flush=True)
    np.save(OUT / f"win_emb_{args.worker}.npy", np.concatenate(mats))
    pd.concat(rows, ignore_index=True).to_parquet(OUT / f"win_res_{args.worker}.parquet", index=False)
    (OUT / f"win_timing_{args.worker}.json").write_text(json.dumps(
        dict(worker=args.worker, proteins=len(mine), residues_kept=nrow, seconds=time.time() - t0,
             gpu=torch.cuda.get_device_name(0)), indent=1))


def merge(args):
    s = ums_sites()
    parts, mats, off = [], [], 0
    for w in range(args.n_workers):
        r = pd.read_parquet(OUT / f"win_res_{w}.parquet")
        m = np.load(OUT / f"win_emb_{w}.npy", mmap_mode="r")
        assert len(r) == m.shape[0]
        parts.append(r.assign(row=r.row + off))
        mats.append(m)
        off += len(r)
    res = pd.concat(parts, ignore_index=True)
    assert not res.duplicated(["acc", "pos"]).any()
    E = np.concatenate(mats)
    order = np.lexsort((res.pos.values, res.acc.values))
    E = E[res.row.values[order]]
    res = res.iloc[order].reset_index(drop=True)
    res["row"] = np.arange(len(res))
    np.save(OUT / "win_emb.npy", E)
    res.to_parquet(OUT / "win_res.parquet", index=False)
    key = pd.Series(res.row.values, index=pd.MultiIndex.from_arrays([res.acc.values, res.pos.values]))
    offs = np.arange(-FLANK, FLANK + 1)
    W = np.full((len(s), len(offs)), -1, dtype=np.int64)
    for j, o in enumerate(offs):
        mi = pd.MultiIndex.from_arrays([s.acc.values, s.pos.values + o])
        W[:, j] = key.reindex(mi).fillna(-1).astype(np.int64).values
    assert (W[:, FLANK] >= 0).all()
    # residue identity check against index windows (window31 = +-15)
    s31 = index.load_sites().set_index("site_id").window31.reindex(s.site_id).values
    aa = res.aa.values
    bad = 0
    for j, o in enumerate(offs):
        ch = np.array([w[15 + o] for w in s31])
        ok = W[:, j] >= 0
        bad += int((aa[W[ok, j]] != ch[ok]).sum()) + int((ch[~ok] != "_").sum())
    if bad:
        raise RuntimeError(f"{bad} window residues disagree with cache/index window31")
    np.save(OUT / "win_rows.npy", W)
    s[["site_id"]].to_parquet(OUT / "win_sites.parquet", index=False)
    np.save(OUT / "win_sites.site_ids.npy", np.array(s.site_id.tolist(), dtype="U"))
    # consistency check vs struct+ESM S/T/Y cache (centre residues, random 5000 sites)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(s), size=5000, replace=False)
    idx, M = esm_embed.load_sty(ROOT / "cache" / "esm2_650m", site_ids=s.site_id.values[pick])
    assert (idx.site_id.values == s.site_id.values[pick]).all()
    A = E[W[pick, FLANK]].astype(np.float32)
    B = np.asarray(M, dtype=np.float32)
    cos = (A * B).sum(1) / np.linalg.norm(A, axis=1) / np.linalg.norm(B, axis=1)
    tim = [json.loads((OUT / f"win_timing_{w}.json").read_text()) for w in range(args.n_workers)]
    meta = dict(model=MODEL, flank=FLANK, n_sites=len(s), n_residues=len(res), shape=list(E.shape),
                centre_cos_vs_esm2_650m_cache=dict(min=float(cos.min()), mean=float(cos.mean()), n=5000),
                window_residue_mismatches=bad, timings=tim,
                gpu_minutes=sum(t["seconds"] for t in tim) / 60)
    (OUT / "win_embed_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "merge"])
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--n-workers", type=int, default=2)
    ap.add_argument("--group", type=int, default=400)
    a = ap.parse_args()
    run(a) if a.cmd == "run" else merge(a)
