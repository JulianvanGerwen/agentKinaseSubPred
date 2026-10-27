"""Phosformer inputs from our canonical universe.

- Kinase domain: canonical UniProt sequence (external/uniprot, via ksp.uniprot) cut at the UniProt "Protein kinase"
  domain coordinates stored in cache/index/kinases.parquet (1-based, inclusive). Kinases with no annotated domain
  (5 of 503 at the time of writing: atypical FASTK, PAN3, PEAK3, FAM20C, FAM20A) fall back to the full canonical
  sequence; every fallback is returned in `fallbacks` so callers log it in their run manifest.
- Peptide: central 11-mer of window15 (positions 3..13 of the ±7 window), '_' padding -> '-' (encoded as <pad>).
"""
from __future__ import annotations

import pandas as pd

from ksp.index import load_kinases
from ksp.uniprot import read_fasta

MAX_KINASE_LEN = 1000  # only hit by full-length fallbacks; inputs are ~250-300 aa domains otherwise


def canonical_seqs(accs) -> dict:
    fa = read_fasta()
    fa = fa[(fa.isoform == "") & fa.acc.isin(set(accs))]
    return dict(zip(fa.acc, fa.seq))


def kinase_domains(kinases: pd.DataFrame | None = None) -> tuple[pd.DataFrame, list[dict]]:
    """Returns (df[uniprot_acc, gene, domain_seq, dom_start, dom_end, source], fallbacks)."""
    k = load_kinases() if kinases is None else kinases
    seqs = canonical_seqs(k.uniprot_acc)
    rows, fallbacks = [], []
    for r in k.itertuples(index=False):
        if r.uniprot_acc not in seqs:
            raise KeyError(f"no canonical UniProt sequence for kinase {r.uniprot_acc} ({r.gene})")
        full = seqs[r.uniprot_acc]
        if pd.notna(r.kinase_domain_start) and pd.notna(r.kinase_domain_end):
            a, b = int(r.kinase_domain_start), int(r.kinase_domain_end)
            if not (1 <= a < b <= len(full)):
                raise ValueError(f"bad domain coords {a}-{b} for {r.uniprot_acc} (len {len(full)})")
            dom, src = full[a - 1 : b], "uniprot_domain"
        else:
            a, b = 1, min(len(full), MAX_KINASE_LEN)
            dom, src = full[:b], "full_sequence_fallback"
            fallbacks.append({"kinase": r.uniprot_acc, "gene": r.gene, "reason": "no kinase-domain coords",
                              "used": f"canonical 1-{b} of {len(full)}"})
        rows.append((r.uniprot_acc, r.gene, dom, a, b, src))
    return pd.DataFrame(rows, columns=["uniprot_acc", "gene", "domain_seq", "dom_start", "dom_end", "source"]), \
        fallbacks


def window15_to_11mer(w15: str) -> str:
    if len(w15) != 15:
        raise ValueError(f"window15 must have length 15: {w15!r}")
    return w15[2:13].replace("_", "-")
