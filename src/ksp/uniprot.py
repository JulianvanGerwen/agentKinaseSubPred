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
    """Parse UniProt pkinfam.txt ("Human and mouse protein kinases: classification and index").
    Sections are titled between '====' lines (e.g. 'AGC Ser/Thr protein kinase family',
    'Atypical: PI3/PI4-kinase family'); rows are 'GENE  ENTRY_HUMAN (ACC) ENTRY_MOUSE (ACC)'.
    Returns columns gene, gene_uc, entry_name, acc, group (section title), group_short (e.g. 'AGC', 'Atypical: PI3/PI4')."""
    lines = _open(path).read().splitlines()
    rows, group = [], None
    row = re.compile(r"^(\S+)\s+(\w+_" + species + r")\s+\((\w+)\s*\)")
    for i, ln in enumerate(lines):
        if 0 < i < len(lines) - 1 and lines[i - 1].startswith("====") and lines[i + 1].startswith("===="):
            group = ln.strip()
            continue
        m = row.match(ln)
        if m and group:
            rows.append((m.group(1), m.group(1).upper(), m.group(2), m.group(3), group))
    df = pd.DataFrame(rows, columns=["gene", "gene_uc", "entry_name", "acc", "group"]).drop_duplicates("acc")
    df["group_short"] = df.group.map(lambda g: "Atypical-" + g.split(":", 1)[1].split()[0] if g.startswith("Atypical:")
                                     else g.split()[0])
    return df
