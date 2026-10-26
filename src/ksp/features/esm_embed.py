"""ESM-2 per-residue embeddings for S/T/Y residues of all Swiss-Prot canonical human proteins (DESIGN §7 F0).

Model: HuggingFace `facebook/esm2_t33_650M_UR50D` (MIT) or `facebook/esm2_t36_3B_UR50D`.
Representation: last layer (`EsmModel.last_hidden_state`, i.e. after `emb_layer_norm_after`, identical to
fair-esm `repr_layers=[num_layers]`), BOS/EOS stripped.

Long proteins (> MAXRES = 1022 residues): sliding windows of 1022 residues with stride `--stride`
(default 511, i.e. ~50% overlap; the last window is right-aligned to the protein end); per-residue
embeddings of overlapping windows are averaged with uniform weight (each window gets its own BOS/EOS).

Outputs (in --out):
  sty_<worker>_<shard>.npy        float16 (n_sites, D) rows for S/T/Y residues, shard order
  sty_index_<worker>.parquet      site_id, acc, pos (1-based), res, file, row
  protmean_<worker>.npy           float32 (n_prot, D) mean over all residues of the (window-averaged) embedding
  protmean_index_<worker>.parquet acc, length, n_windows, row
  timing_<worker>.json            throughput
  (merge) sty_index.parquet, sty_<w>_<shard>.site_ids.npy (explicit row ids per shard), protein_mean.npy +
          protein_mean_index.parquet + protein_mean.accs.npy. Aligned exports for the feature store:
          scripts/export_structure_esm_features.py (src/ksp/store.py).
`merge` subcommand concatenates the per-worker index files into sty_index.parquet / protein_mean.{npy,parquet}.

Usage:
  CUDA_VISIBLE_DEVICES=2 python -I src/ksp/features/esm_embed.py run --worker 0 --n-workers 2 --out cache/esm2_650m
  python -I src/ksp/features/esm_embed.py merge --out cache/esm2_650m
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from ksp.uniprot import read_fasta  # noqa: E402

MAXRES = 1022
STY = set("STY")


def canonical_swissprot() -> pd.DataFrame:
    df = read_fasta(ROOT / "external/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz")
    c = df[(df.db == "sp") & (df.isoform == "")][["acc", "gene", "seq"]].reset_index(drop=True)
    assert c.acc.is_unique
    return c


def windows(L: int, win: int = MAXRES, stride: int = 511) -> list[tuple[int, int]]:
    """0-based half-open windows covering [0, L). Last window right-aligned."""
    if L <= win:
        return [(0, L)]
    starts = list(range(0, L - win + 1, stride))
    if starts[-1] + win < L:
        starts.append(L - win)
    return [(s, s + win) for s in starts]


def load_model(name: str, dtype: str = "fp16"):
    import torch
    from transformers import AutoTokenizer, EsmModel
    td = {"fp16": torch.float16, "bf16": torch.bfloat16, "fp32": torch.float32}[dtype]
    tok = AutoTokenizer.from_pretrained(name)
    try:
        model = EsmModel.from_pretrained(name, dtype=td, add_pooling_layer=False, attn_implementation="sdpa")
    except (ValueError, TypeError):
        model = EsmModel.from_pretrained(name, dtype=td, add_pooling_layer=False)
    model = model.cuda().eval()
    return tok, model


def embed_sequences(seqs: list[str], tok, model, stride: int = 511, token_budget: int = 16384,
                    max_batch: int = 64) -> list[np.ndarray]:
    """Return list of float32 arrays (L_i, D): last-layer per-residue embeddings, overlaps averaged."""
    import torch
    D = model.config.hidden_size
    acc = [np.zeros((len(s), D), dtype=np.float32) for s in seqs]
    cnt = [np.zeros(len(s), dtype=np.float32) for s in seqs]
    jobs = [(i, a, b) for i, s in enumerate(seqs) for (a, b) in windows(len(s), MAXRES, stride)]
    jobs.sort(key=lambda j: -(j[2] - j[1]))
    k = 0
    while k < len(jobs):
        Lmax = jobs[k][2] - jobs[k][1] + 2
        bs = max(1, min(max_batch, token_budget // Lmax))
        batch = jobs[k:k + bs]
        k += bs
        enc = tok([seqs[i][a:b] for i, a, b in batch], return_tensors="pt", padding=True, add_special_tokens=True)
        with torch.inference_mode():
            out = model(input_ids=enc["input_ids"].cuda(), attention_mask=enc["attention_mask"].cuda())
        h = out.last_hidden_state.float().cpu().numpy()
        for r, (i, a, b) in enumerate(batch):
            acc[i][a:b] += h[r, 1:1 + (b - a)]
            cnt[i][a:b] += 1
    return [x / c[:, None] for x, c in zip(acc, cnt)]


def run(args):
    import torch
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    prot = canonical_swissprot()
    prot["L"] = prot.seq.str.len()
    prot = prot.sort_values(["L", "acc"], ascending=[False, True]).reset_index(drop=True)
    if args.accs:
        keep = set(Path(args.accs).read_text().split())
        prot = prot[prot.acc.isin(keep)].reset_index(drop=True)
    if args.limit:
        rng = np.random.default_rng(0)
        prot = prot.iloc[np.sort(rng.choice(len(prot), size=min(args.limit, len(prot)), replace=False))]
        prot = prot.reset_index(drop=True)
    mine = prot.iloc[args.worker::args.n_workers].reset_index(drop=True)
    tok, model = load_model(args.model, args.dtype)
    D = model.config.hidden_size
    w = args.worker
    idx_rows, pm_rows, pm_vecs = [], [], []
    t0 = time.time()
    n_res = n_tok = 0
    torch.cuda.reset_peak_memory_stats()
    for sh, g0 in enumerate(range(0, len(mine), args.group)):
        g = mine.iloc[g0:g0 + args.group]
        embs = embed_sequences(g.seq.tolist(), tok, model, stride=args.stride, token_budget=args.token_budget)
        rows = []
        fname = f"sty_{w}_{sh:03d}.npy"
        for (acc_, seq), E in zip(g[["acc", "seq"]].itertuples(index=False), embs):
            pm_rows.append((acc_, len(seq), len(windows(len(seq), MAXRES, args.stride)), len(pm_rows)))
            pm_vecs.append(E.mean(0))
            n_res += len(seq)
            n_tok += sum(b - a + 2 for a, b in windows(len(seq), MAXRES, args.stride))
            pos = [p for p, ch in enumerate(seq) if ch in STY]
            if pos:
                rows.append(E[pos].astype(np.float16))
                for p in pos:
                    idx_rows.append((f"{acc_}_{seq[p]}{p + 1}", acc_, p + 1, seq[p], fname))
        if not args.no_sty:
            M = np.concatenate(rows) if rows else np.zeros((0, D), np.float16)
            np.save(out / fname, M)
        el = time.time() - t0
        print(f"[w{w}] shard {sh} proteins {g0 + len(g)}/{len(mine)} res {n_res} {n_res / el:.0f} res/s "
              f"elapsed {el / 60:.1f} min peakmem {torch.cuda.max_memory_allocated() / 2**30:.1f}GB", flush=True)
    el = time.time() - t0
    idx = pd.DataFrame(idx_rows, columns=["site_id", "acc", "pos", "res", "file"])
    idx["row"] = idx.groupby("file").cumcount()
    if not args.no_sty:
        idx.to_parquet(out / f"sty_index_{w}.parquet", index=False)
    np.save(out / f"protmean_{w}.npy", np.stack(pm_vecs).astype(np.float32))
    pd.DataFrame(pm_rows, columns=["acc", "length", "n_windows", "row"]).to_parquet(
        out / f"protmean_index_{w}.parquet", index=False)
    timing = dict(model=args.model, dtype=args.dtype, worker=w, n_workers=args.n_workers, proteins=len(mine),
                  residues=n_res, tokens=n_tok, sty_sites=len(idx), seconds=el, residues_per_s=n_res / el,
                  proteins_per_s=len(mine) / el, gpu=torch.cuda.get_device_name(0),
                  peak_mem_gb=torch.cuda.max_memory_allocated() / 2**30, stride=args.stride,
                  token_budget=args.token_budget, limit=args.limit)
    (out / f"timing_{w}.json").write_text(json.dumps(timing, indent=1))
    print(json.dumps(timing), flush=True)


def merge(args):
    out = Path(args.out)
    idx = pd.concat([pd.read_parquet(p) for p in sorted(out.glob("sty_index_*.parquet"))], ignore_index=True)
    assert idx.site_id.is_unique
    idx = idx.sort_values(["acc", "pos"]).reset_index(drop=True)
    idx.to_parquet(out / "sty_index.parquet", index=False)
    # explicit site_id array next to every shard (row r of sty_X.npy <-> site_ids[r]); guards against order hazards
    for f, g in idx.groupby("file"):
        g = g.sort_values("row")
        assert (g.row.values == np.arange(len(g))).all() and len(g) == np.load(out / f, mmap_mode="r").shape[0]
        np.save(out / f.replace(".npy", ".site_ids.npy"), g.site_id.values.astype(str))
    pis, vecs, off = [], [], 0
    for p in sorted(out.glob("protmean_index_*.parquet")):
        d = pd.read_parquet(p)
        v = np.load(out / p.name.replace("protmean_index_", "protmean_").replace(".parquet", ".npy"))
        pis.append(d.assign(row=d.row + off))
        vecs.append(v)
        off += len(d)
    pi = pd.concat(pis, ignore_index=True)
    V = np.concatenate(vecs)
    order = np.argsort(pi.acc.values, kind="stable")
    pi = pi.iloc[order].reset_index(drop=True)
    V = V[pi.row.values]
    pi["row"] = np.arange(len(pi))
    np.save(out / "protein_mean.npy", V)
    np.save(out / "protein_mean.accs.npy", pi.acc.values.astype(str))
    pi.to_parquet(out / "protein_mean_index.parquet", index=False)
    print(f"merged: {len(idx)} S/T/Y rows in {idx.file.nunique()} shards; {len(pi)} protein means", flush=True)


def load_sty(out: str | Path, site_ids=None) -> tuple[pd.DataFrame, np.ndarray]:
    """Convenience loader: returns (index, float16 matrix) for the requested site_ids (all if None)."""
    out = Path(out)
    idx = pd.read_parquet(out / "sty_index.parquet")
    if site_ids is not None:
        idx = idx.set_index("site_id").loc[list(site_ids)].reset_index()
    mats = []
    for f, g in idx.groupby("file", sort=False):
        m = np.load(out / f, mmap_mode="r")
        mats.append((g.index.values, np.asarray(m[g.row.values])))
    D = mats[0][1].shape[1]
    M = np.zeros((len(idx), D), np.float16)
    for ii, m in mats:
        M[ii] = m
    return idx, M


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    r.add_argument("--dtype", default="fp16", choices=["fp16", "bf16", "fp32"])
    r.add_argument("--out", required=True)
    r.add_argument("--worker", type=int, default=0)
    r.add_argument("--n-workers", type=int, default=1)
    r.add_argument("--group", type=int, default=500, help="proteins per output shard")
    r.add_argument("--stride", type=int, default=511)
    r.add_argument("--token-budget", type=int, default=16384)
    r.add_argument("--limit", type=int, default=0, help="random subset of N proteins (seed 0), for benchmarking")
    r.add_argument("--accs", default="", help="file with accessions to restrict to")
    r.add_argument("--no-sty", action="store_true", help="skip writing S/T/Y matrices (benchmark only)")
    m = sub.add_parser("merge")
    m.add_argument("--out", required=True)
    a = ap.parse_args()
    run(a) if a.cmd == "run" else merge(a)


if __name__ == "__main__":
    main()
