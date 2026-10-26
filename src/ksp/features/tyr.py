"""TK-specific features (F9): SH2-mediated docking/processivity, receptor/membrane proximity,
kinase adaptor-domain descriptors, ITAM/ITIM motif context.

None of this is label-derived: the only phospho information used is *other* sites' membership of
U_ms (mass-spec observation), exactly as the F1 `primed_soft` variant already does (DESIGN §7 F1).
The query site's own observation status is never used (A6 / DL2-7).

Sources
  SH2 specificity matrices : Yaron-Barir 2024 Nature Suppl. Table 5 (MOESM7), sheet
                             `sh2_norm_scaled_matrices` -- 76 SH2 domains x 6 positions
                             (-2,-1,+1,+2,+3,+4 relative to pY) x 20 aa, position-normalised and
                             scaled to random amino-acid frequency (so 1.0 = no preference).
  Adaptor domains          : UniProt 2026_03 `Domain [FT]` / `Keywords` (SH2, SH3, PID/PTB,
                             Transmembrane, Receptor).
  Localisation             : existing F4 stores `context/loc/vec_{hpa,hcm}_protein`.

Public API
  load_sh2_matrices()            -> (DataFrame log2 matrices indexed by SH2 domain name, positions)
  sh2_domain_to_kinase()         -> {sh2 name: kinase gene} for the 24 kinase SH2 domains
  kinase_adaptor_domains()       -> DataFrame (acc, n_sh2, n_sh3, n_ptb, is_tm, is_receptor)
  score_windows(mat, windows)    -> log2 SH2 score per 15-mer window
  itam_features(seqs, sites)     -> per-site ITAM/ITIM descriptors
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ksp.index import ROOT

SH2_XLSX = ROOT / "external/cantley_yaronbarir2024/41586_2024_7407_MOESM7_ESM.xlsx"
SH2_SHEET = "sh2_norm_scaled_matrices"
UNIPROT_ANN = ROOT / "external/uniprot/human_tsv/UP000005640_annotations.tsv.gz"

#: SH2 matrix name -> kinase HGNC symbol. Only domains belonging to a human protein kinase.
#: SYK/ZAP70/PLCG etc. have N- and C-terminal SH2s; both Syk-family tandem domains are kept.
SH2_TO_KINASE = {
    "ABL1": "ABL1", "ABL2": "ABL2", "BLK": "BLK", "BMX": "BMX", "BRK": "PTK6", "BTK": "BTK",
    "CSK": "CSK", "FER": "FER", "FES": "FES", "FGR": "FGR", "FRK": "FRK", "FYN": "FYN",
    "HCK": "HCK", "ITK": "ITK", "LCK": "LCK", "LYN": "LYN", "MATK": "MATK", "SRC": "SRC",
    "SRMS": "SRMS", "SYK_N": "SYK", "SYK_C": "SYK", "TEC": "TEC", "TXK": "TXK", "YES": "YES1",
    "ZAP70_N": "ZAP70",
}
#: SH2 offsets relative to the phosphotyrosine, as shipped by Yaron-Barir 2024.
SH2_OFFSETS = (-2, -1, 1, 2, 3, 4)
AAS = "PGACSTVILMFYWHKRQNDE"
FLANK = 7                      # window15 = +/-7


def load_sh2_matrices(path: Path = SH2_XLSX, sheet: str = SH2_SHEET) -> pd.DataFrame:
    """log2 of the position-normalised, random-frequency-scaled SH2 matrices.

    Returns a DataFrame indexed by SH2 domain name with MultiIndex columns (offset, aa).
    0 means 'no preference', because the shipped values are already scaled so that 1.0 = random.
    """
    d = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    d = d.rename(columns={d.columns[0]: "sh2"}).set_index("sh2")
    cols, keep = [], []
    for c in d.columns:
        m = re.match(r"^(-?\d+)([A-Z])$", str(c))
        if m and int(m.group(1)) in SH2_OFFSETS and m.group(2) in AAS:
            cols.append((int(m.group(1)), m.group(2)))
            keep.append(c)
    d = d[keep].apply(pd.to_numeric, errors="coerce")
    d.columns = pd.MultiIndex.from_tuples(cols, names=["offset", "aa"])
    out = np.log2(d.clip(lower=1e-3))
    return out.sort_index(axis=1)


def sh2_lookup(mats: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    """(domain names, array[n_domains, n_offsets, 256]) for fast byte-indexed scoring."""
    names = list(mats.index)
    tab = np.zeros((len(names), len(SH2_OFFSETS), 256), dtype=np.float32)
    for oi, off in enumerate(SH2_OFFSETS):
        for aa in AAS:
            if (off, aa) in mats.columns:
                tab[:, oi, ord(aa)] = mats[(off, aa)].to_numpy(dtype=np.float32)
    return names, tab


def score_windows(tab: np.ndarray, windows: np.ndarray) -> np.ndarray:
    """SH2 log2 score of each 15-mer (centre = index FLANK) for every domain.

    tab: (n_domains, n_offsets, 256) from `sh2_lookup`. windows: array of 15-char strings.
    Returns (n_domains, n_windows) float32. Out-of-range flanks ('_', 'X', 'U') contribute 0,
    as in the F1 Cantley engine.
    """
    w = np.frombuffer("".join(windows).encode("ascii", "replace"), dtype=np.uint8)
    w = w.reshape(len(windows), -1)
    out = np.zeros((tab.shape[0], len(windows)), dtype=np.float32)
    for oi, off in enumerate(SH2_OFFSETS):
        col = w[:, FLANK + off]
        out += tab[:, oi, :][:, col]
    return out


def sh2_domain_to_kinase() -> dict[str, str]:
    return dict(SH2_TO_KINASE)


def kinase_adaptor_domains(accs, path: Path = UNIPROT_ANN) -> pd.DataFrame:
    """Per accession: SH2 / SH3 / PID(PTB) domain counts and transmembrane / receptor keywords."""
    u = pd.read_csv(path, sep="\t", dtype=str,
                    usecols=["Entry", "Domain [FT]", "Keywords"]).set_index("Entry")
    dom, kw = u["Domain [FT]"].fillna(""), u["Keywords"].fillna("")
    rows = []
    for a in accs:
        s, k = dom.get(a, ""), kw.get(a, "")
        rows.append(dict(
            acc=a,
            n_sh2=len(re.findall(r'/note="SH2', s)),
            n_sh3=len(re.findall(r'/note="SH3', s)),
            n_ptb=len(re.findall(r'/note="(?:PID|PTB)', s)),
            is_tm=int("Transmembrane" in k),
            is_receptor=int("Receptor" in k)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ ITAM / ITIM context
def itam_features(seqs: dict[str, str], acc: np.ndarray, pos: np.ndarray) -> pd.DataFrame:
    """ITAM/ITIM-like context of each Y site (SYK/ZAP70/SRC-family recruitment).

    ITAM half-site  : Yxx[LIV]                       (hydrophobic at +3)
    ITAM full       : a half-site here AND a second Yxx[LIV] 6-12 residues up- or downstream
    ITIM            : [ILVS]xYxx[LV]                 (hydrophobic at -2 and +3)
    ITAM order      : 0 none, 1 half only, 2 full (this Y is in a tandem Yxx[LIV] pair)
    """
    hyd3, hyd3i = set("LIV"), set("LV")
    n2 = set("ILVS")
    rows = []
    for a, p in zip(acc, pos):
        s = seqs.get(a)
        if s is None or p < 1 or p > len(s) or s[p - 1] != "Y":
            rows.append((np.nan, np.nan, np.nan, np.nan))
            continue
        i = p - 1

        def half(j):                                   # Y at 0-based j with Yxx[LIV]
            return 0 <= j < len(s) and s[j] == "Y" and j + 3 < len(s) and s[j + 3] in hyd3

        h = half(i)
        partner = any(half(i + d) for d in list(range(6, 13)) + list(range(-12, -5)))
        itim = int(i - 2 >= 0 and s[i - 2] in n2 and i + 3 < len(s) and s[i + 3] in hyd3i)
        rows.append((int(h), int(h and partner), itim, 2 if (h and partner) else int(h)))
    return pd.DataFrame(rows, columns=["itam_half", "itam_full", "itim", "itam_order"])
