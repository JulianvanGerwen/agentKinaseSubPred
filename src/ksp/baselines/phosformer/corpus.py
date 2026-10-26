"""MLM pretraining corpus (paper §2.3.1): 295,320 NCBI-nr kinase domains + 5,338 full-length substrate proteins.

The authors released exactly this corpus (external/phosformer/hf_dataset_curated/corpus/{kinases,substrates}.fa,
counts verified: 295,320 / 5,338 sequences; 86,944,788 / 3,845,761 residues — identical to the paper), so we use it
unchanged rather than reconstructing it from UniProt/Pfam.

Cached form: cache/phosformer/corpus/{tokens.npy (uint8 ids), offsets.npy (int64, n+1), source.npy (0=kinase,
1=substrate)} + meta.json with input SHA-256s.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .tokens import encode_seq

ROOT = Path(__file__).resolve().parents[4]
CORPUS_SRC = ROOT / "external" / "phosformer" / "hf_dataset_curated" / "corpus"
CORPUS_DIR = ROOT / "cache" / "phosformer" / "corpus"


def read_fasta_seqs(path: Path) -> list[str]:
    out, cur = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    out.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        out.append("".join(cur))
    return out


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build_corpus(out_dir: Path = CORPUS_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = [CORPUS_SRC / "kinases.fa", CORPUS_SRC / "substrates.fa"]
    toks, offs, src = [], [0], []
    meta = {"inputs": {}, "counts": {}}
    for si, p in enumerate(files):
        seqs = read_fasta_seqs(p)
        meta["inputs"][p.name] = _sha(p)
        meta["counts"][p.name] = {"n_seq": len(seqs), "n_res": int(sum(map(len, seqs)))}
        for s in seqs:
            t = encode_seq(s).astype(np.uint8)
            toks.append(t)
            offs.append(offs[-1] + len(t))
            src.append(si)
    np.save(out_dir / "tokens.npy", np.concatenate(toks))
    np.save(out_dir / "offsets.npy", np.asarray(offs, dtype=np.int64))
    np.save(out_dir / "source.npy", np.asarray(src, dtype=np.uint8))
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def load_corpus(out_dir: Path = CORPUS_DIR):
    if not (out_dir / "offsets.npy").exists():
        build_corpus(out_dir)
    return (np.load(out_dir / "tokens.npy", mmap_mode="r"), np.load(out_dir / "offsets.npy"),
            np.load(out_dir / "source.npy"))
