"""Canonical row/column orders for all feature stores (DESIGN §7).

Every kinase x site / kinase x protein feature matrix in cache/features/ MUST use
these orders. Build once with scripts/build_index.py; never re-sort elsewhere.
"""
from pathlib import Path
import pandas as pd

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
INDEX_DIR = ROOT / "cache" / "index"


def load_kinases() -> pd.DataFrame:
    """Included kinases (+ canary pseudo-kinases appended later by ksp.canary), row i = kinase index i."""
    return pd.read_parquet(INDEX_DIR / "kinases.parquet")


def load_sites() -> pd.DataFrame:
    """U_obs sites, row j = site index j (sorted by site_id)."""
    return pd.read_parquet(INDEX_DIR / "sites.parquet")


def load_proteins() -> pd.DataFrame:
    """Proteins carrying >=1 U_obs site, row p = protein index p (sorted by acc)."""
    return pd.read_parquet(INDEX_DIR / "proteins.parquet")


def kinase_pos() -> dict:
    return {a: i for i, a in enumerate(load_kinases().uniprot_acc)}


def site_pos() -> dict:
    return {s: j for j, s in enumerate(load_sites().site_id)}


def protein_pos() -> dict:
    return {a: p for p, a in enumerate(load_proteins().acc)}
