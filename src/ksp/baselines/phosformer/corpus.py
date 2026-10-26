"""MLM pretraining corpus (paper §2.3.1): 295,320 NCBI-nr kinase domains + 5,338 full-length substrate proteins.

The authors released exactly this corpus (external/phosformer/hf_dataset_curated/corpus/{kinases,substrates}.fa,
counts verified: 295,320 / 5,338 sequences; 86,944,788 / 3,845,761 residues — identical to the paper), so we use it
unchanged rather than reconstructing it from UniProt/Pfam.

FILTERED corpus (used for B3a; PI AIFM2-8): cache/phosformer/corpus_filtered/ built by
scripts/baselines/phosformer_corpus_filter.py (TEST proteins + homologs removed).

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


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """[(header without '>', sequence)]"""
    out, hdr, cur = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(cur)))
                hdr, cur = line[1:].rstrip("\n"), []
            else:
                cur.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(cur)))
    return out


def read_fasta_seqs(path: Path) -> list[str]:
    return [s for _, s in read_fasta(path)]


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build_corpus(out_dir: Path = CORPUS_DIR, src_dir: Path = CORPUS_SRC) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = [src_dir / "kinases.fa", src_dir / "substrates.fa"]
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


def load_corpus(out_dir: Path = CORPUS_DIR, src_dir: Path = CORPUS_SRC):
    if not (out_dir / "offsets.npy").exists():
        build_corpus(out_dir, src_dir)
    return (np.load(out_dir / "tokens.npy", mmap_mode="r"), np.load(out_dir / "offsets.npy"),
            np.load(out_dir / "source.npy"))
