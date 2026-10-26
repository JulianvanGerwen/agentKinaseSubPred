"""Readers for the UniProt files downloaded by scripts/download_uniprot.py (external/uniprot/)."""
from __future__ import annotations

import gzip
import os
import re
from pathlib import Path

import pandas as pd

DEFAULT_DIR = Path(os.environ.get("KSP_UNIPROT", Path(__file__).resolve().parents[2] / "external" / "uniprot"))
_HDR = re.compile(r"^>(sp|tr)\|([^|]+)\|(\S+)\s?(.*)$")


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path=DEFAULT_DIR / "proteome" / "UP000005640_all_isoforms.fasta.gz") -> pd.DataFrame:
    """Return DataFrame: acc_full, acc (base), isoform ('' for canonical), db (sp/tr), entry_name, gene, sv, seq.
    In UniProt FASTA the canonical isoform carries the bare accession; other isoforms 'ACC-n'."""
    rows, hdr, seq = [], None, []

    def flush():
        if hdr is None:
            return
        m = _HDR.match(hdr)
        db, acc_full, entry, desc = m.groups()
        gn = re.search(r"\bGN=(\S+)", desc)
        sv = re.search(r"\bSV=(\d+)", desc)
        base, _, iso = acc_full.partition("-")
        rows.append((acc_full, base, iso, db, entry, gn.group(1) if gn else "", int(sv.group(1)) if sv else None,
                     "".join(seq)))

    with _open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                hdr, seq = line, []
            else:
                seq.append(line.strip())
        flush()
    return pd.DataFrame(rows, columns=["acc_full", "acc", "isoform", "db", "entry_name", "gene", "sv", "seq"])


def read_annotations(path=DEFAULT_DIR / "human_tsv" / "UP000005640_annotations.tsv.gz") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def parse_pkinfam(path=DEFAULT_DIR / "pkinfam" / "pkinfam.txt", species: str = "HUMAN") -> pd.DataFrame:
    """Parse UniProt pkinfam.txt (protein kinase family classification, Hanks/Manning based).
    Returns rows (entry_name, acc, group, family) for the given species suffix.
    'group' = the top-level '===' underlined section title; 'family' = most recent subfamily line."""
    lines = _open(path).read().splitlines()
    rows, group, family = [], None, None
    pat = re.compile(r"(\w+_" + species + r")\s+\((\w+)\s*\)")
    for i, ln in enumerate(lines):
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if nxt.startswith("====") and ln.strip() and not ln.startswith("="):
            group, family = ln.strip(), None
            continue
        if group is None:
            continue
        s = ln.strip()
        if s and not pat.search(ln) and not s.startswith("=") and not s.startswith("-") and "(" not in s:
            family = s
        for m in pat.finditer(ln):
            rows.append((m.group(1), m.group(2), group, family))
    return pd.DataFrame(rows, columns=["entry_name", "acc", "group", "family"]).drop_duplicates("acc")
