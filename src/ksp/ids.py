"""Canonical identifiers: Swiss-Prot canonical human sequence store, site windows, 15-mer -> canonical site mapping,
and site_id formatting/parsing.

Conventions (docs/AGENT_BRIEF.md "Canonical identifiers"):
  * protein = UniProt Swiss-Prot canonical accession (no isoform suffix), UniProt release in external/uniprot;
  * site_id = f"{acc}_{res}{pos}" (res in S/T/Y, pos 1-based on the canonical sequence), e.g. "P05198_S52";
  * windows are centred on the site, padded with '_' beyond the protein termini (±7 -> 15-mer, ±15 -> 31-mer).

Mapping a source site (accession, stated position, 15-mer) to a canonical site (`CanonicalStore.map_window`):
  1. candidate canonical accessions: the source base accession if it is a Swiss-Prot canonical entry; otherwise the
     Swiss-Prot canonical entries of the source gene symbol (given explicitly, or looked up from the UniProt entry of
     the source accession, e.g. a TrEMBL accession) — "the same gene's canonical";
  2. 'exact_pos'  : the window matches the canonical window at the stated position;
  3. 'shifted'    : otherwise, the window matches at exactly one other position over all candidate canonicals;
  4. otherwise unmapped with a reason: 'bad_window' (length != 15 or centre not the stated residue),
     'no_candidate' (no Swiss-Prot canonical for the accession/gene), 'no_match' (window not found),
     'multiple_matches' (window found at >1 position — e.g. repeats/paralogous duplicates within the gene).
  Window comparison: case-insensitive; '_' (padding) and 'X' in the SOURCE window are wildcards (a source window
  padded where the canonical has residues, e.g. an initiator Met removed, still matches). The canonical centre
  residue must equal the source centre residue.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import uniprot

SITE_RE = re.compile(r"^([A-Z0-9]+)_([STY])(\d+)$")
STY = frozenset("STY")
WILDCARDS = frozenset("_X")


# ------------------------------------------------------------------------------------------------- site ids
def make_site_id(acc: str, res: str, pos: int) -> str:
    if "-" in acc:
        raise ValueError(f"site_id needs a canonical accession without isoform suffix: {acc}")
    res = res.upper()
    if res not in STY:
        raise ValueError(f"residue must be S/T/Y, got {res!r}")
    pos = int(pos)
    if pos < 1:
        raise ValueError(f"position must be >= 1, got {pos}")
    return f"{acc}_{res}{pos}"


def parse_site_id(site_id: str) -> tuple[str, str, int]:
    m = SITE_RE.match(site_id)
    if not m:
        raise ValueError(f"not a site_id: {site_id!r}")
    return m.group(1), m.group(2), int(m.group(3))


def site_ids(acc: pd.Series, res: pd.Series, pos: pd.Series) -> pd.Series:
    """Vectorised make_site_id (no validation beyond string formatting)."""
    return acc.astype(str) + "_" + res.astype(str) + pos.astype("Int64").astype(str)


# ------------------------------------------------------------------------------------------------- windows
def window(seq: str, pos: int, flank: int = 7, pad: str = "_") -> str:
    """Window of length 2*flank+1 centred on 1-based `pos`, padded with `pad` beyond the termini."""
    i = int(pos) - 1
    if i < 0 or i >= len(seq):
        raise IndexError(f"position {pos} outside sequence of length {len(seq)}")
    left = seq[max(0, i - flank): i]
    right = seq[i + 1: i + 1 + flank]
    return pad * (flank - len(left)) + left + seq[i] + right + pad * (flank - len(right))


def windows_match(src: str, ref: str) -> bool:
    """True if every non-wildcard character of `src` equals `ref` at the same index (same length required)."""
    if len(src) != len(ref):
        return False
    return all(a == b or a in WILDCARDS for a, b in zip(src, ref))


def find_window(src: str, seq: str) -> list[int]:
    """All 1-based centre positions in `seq` where the 15-mer (or any odd-length window) `src` matches with
    wildcards ('_'/'X' in src match anything, including beyond-terminus padding)."""
    src = src.upper()
    L = len(src)
    flank = L // 2
    core = src.strip("_")
    if not core:
        return []
    lead = len(src) - len(src.lstrip("_"))
    hits = []
    if "X" not in core:
        start = seq.find(core)
        while start >= 0:
            centre = start + (flank - lead) + 1  # 1-based
            if 1 <= centre <= len(seq):
                hits.append(centre)
            start = seq.find(core, start + 1)
        return hits
    # slow path for internal 'X' wildcards
    c = src[flank]
    for i, ch in enumerate(seq):
        if ch == c and windows_match(src, window(seq, i + 1, flank)):
            hits.append(i + 1)
    return hits


@dataclass
class MapResult:
    acc: str | None
    pos: int | None
    res: str | None
    method: str  # 'exact_pos' | 'shifted' | 'unmapped'
    reason: str  # '' or failure reason
    n_hits: int = 0

    @property
    def site_id(self) -> str | None:
        if self.acc is None or self.res not in STY:
            return None
        return make_site_id(self.acc, self.res, self.pos)


# ------------------------------------------------------------------------------------------------- store
class CanonicalStore:
    """Swiss-Prot canonical human sequences (acc -> seq) plus gene-symbol and accession indices.

    seqs        : dict acc -> canonical sequence (Swiss-Prot, isoform suffix '')
    gene_to_accs: dict upper-case gene symbol -> list of Swiss-Prot canonical accessions (primary gene names from
                  the FASTA 'GN=' field, plus optional extra symbol -> acc pairs, e.g. UniProt/HGNC synonyms)
    acc_to_gene : dict any UniProt accession of the proteome (Swiss-Prot or TrEMBL, base accession) -> gene symbol
    """

    def __init__(self, seqs: dict, gene_to_accs: dict | None = None, acc_to_gene: dict | None = None):
        self.seqs = dict(seqs)
        self.gene_to_accs = {k.upper(): list(dict.fromkeys(v)) for k, v in (gene_to_accs or {}).items()}
        self.acc_to_gene = dict(acc_to_gene or {})

    # -- construction
    @classmethod
    def from_uniprot(cls, reviewed_fasta: Path | None = None, all_fasta: Path | None = None) -> "CanonicalStore":
        """Build from the UniProt FASTA files in external/uniprot (see scripts/download_uniprot.py).
        `all_fasta` (Swiss-Prot + TrEMBL) is only used for the accession -> gene index (TrEMBL accessions)."""
        rf = reviewed_fasta or uniprot.DEFAULT_DIR / "proteome" / "UP000005640_reviewed_isoforms.fasta.gz"
        sp = uniprot.read_fasta(rf)
        can = sp[(sp.isoform == "") & (sp.db == "sp")]
        seqs = dict(zip(can.acc, can.seq))
        g2a: dict[str, list] = {}
        for g, a in zip(can.gene.str.upper(), can.acc):
            if g:
                g2a.setdefault(g, []).append(a)
        a2g = dict(zip(can.acc, can.gene))
        if all_fasta is not False:
            af = all_fasta or uniprot.DEFAULT_DIR / "proteome" / "UP000005640_all_isoforms.fasta.gz"
            if Path(af).exists():
                al = uniprot.read_fasta(af)
                al = al[al.isoform == ""]
                for a, g in zip(al.acc, al.gene):
                    a2g.setdefault(a, g)
        st = cls(seqs, g2a, a2g)
        st.isoform_seqs = dict(zip(sp.acc_full, sp.seq))  # all Swiss-Prot isoforms (incl. canonical, bare acc)
        return st

    def add_gene_synonyms(self, pairs) -> int:
        """Add (symbol, acc) pairs to the gene index without overriding primary symbols. Returns # added.
        A synonym already used as a primary symbol of another entry is NOT added (ambiguous)."""
        n = 0
        primary = set(self.gene_to_accs)
        for sym, acc in pairs:
            s = str(sym).upper()
            if not s or acc not in self.seqs or s in primary:
                continue
            lst = self.gene_to_accs.setdefault(s, [])
            if acc not in lst:
                lst.append(acc)
                n += 1
        return n

    # -- accessors
    def __contains__(self, acc) -> bool:
        return acc in self.seqs

    def window(self, acc: str, pos: int, flank: int = 7) -> str:
        return window(self.seqs[acc], pos, flank)

    def residue(self, acc: str, pos: int) -> str | None:
        s = self.seqs.get(acc)
        if s is None or pos is None or pos < 1 or pos > len(s):
            return None
        return s[pos - 1]

    def candidates(self, acc: str | None, gene: str | None = None) -> list[str]:
        base = (acc or "").split("-", 1)[0]
        if base in self.seqs:
            return [base]
        g = (gene or "").upper() or (self.acc_to_gene.get(base, "") or "").upper()
        return list(self.gene_to_accs.get(g, [])) if g else []

    # -- mapping
    def map_window(self, acc: str | None, pos, win: str, res: str | None = None, gene: str | None = None,
                   candidates: list[str] | None = None) -> MapResult:
        """Map a source site to a canonical site; see the module docstring for the rules."""
        if not isinstance(win, str) or len(win) % 2 == 0:
            return MapResult(None, None, None, "unmapped", "bad_window")
        w = win.upper()
        flank = len(w) // 2
        c = w[flank]
        if res is not None and str(res).upper() != c:
            return MapResult(None, None, None, "unmapped", "bad_window")
        cands = candidates if candidates is not None else self.candidates(acc, gene)
        cands = [a for a in cands if a in self.seqs]
        if not cands:
            return MapResult(None, None, c, "unmapped", "no_candidate")
        try:
            p = int(pos)
        except (TypeError, ValueError):
            p = None
        if p is not None:
            for a in cands:
                s = self.seqs[a]
                if 1 <= p <= len(s) and s[p - 1] == c and windows_match(w, window(s, p, flank)):
                    return MapResult(a, p, c, "exact_pos", "", 1)
        hits = [(a, h) for a in cands for h in find_window(w, self.seqs[a]) if self.seqs[a][h - 1] == c]
        hits = list(dict.fromkeys(hits))
        if len(hits) == 1:
            a, h = hits[0]
            return MapResult(a, h, c, "shifted", "", 1)
        if len(hits) > 1:
            return MapResult(None, None, c, "unmapped", "multiple_matches", len(hits))
        return MapResult(None, None, c, "unmapped", "no_match", 0)

    def map_position(self, acc: str | None, pos, res: str, gene: str | None = None,
                     expected_length: int | None = None) -> MapResult:
        """For sources WITHOUT a sequence window (e.g. Ochoa 2020, OmniPath, UniProt MOD_RES): accept the stated
        position on the canonical of `acc` if the residue matches (and, if `expected_length` is given, the canonical
        length equals it — a proxy for 'same sequence version'). Methods: 'pos_residue' or unmapped with reason
        'no_candidate' / 'pos_out_of_range' / 'residue_mismatch' / 'length_changed'."""
        base = (acc or "").split("-", 1)[0]
        if "-" in (acc or "") and not (acc or "").endswith("-1"):
            return MapResult(None, None, res, "unmapped", "isoform_accession")
        cands = [base] if base in self.seqs else []
        if not cands:
            return MapResult(None, None, res, "unmapped", "no_candidate")
        s = self.seqs[base]
        try:
            p = int(pos)
        except (TypeError, ValueError):
            return MapResult(None, None, res, "unmapped", "pos_out_of_range")
        if p < 1 or p > len(s):
            return MapResult(None, None, res, "unmapped", "pos_out_of_range")
        if expected_length is not None and pd.notna(expected_length) and int(expected_length) != len(s):
            return MapResult(None, None, res, "unmapped", "length_changed")
        if s[p - 1] != str(res).upper():
            return MapResult(None, None, res, "unmapped", "residue_mismatch")
        return MapResult(base, p, s[p - 1], "pos_residue", "", 1)


def map_table(store: CanonicalStore, df: pd.DataFrame, acc_col="acc", pos_col="pos", win_col="window",
              res_col: str | None = None, gene_col: str | None = None) -> pd.DataFrame:
    """Apply CanonicalStore.map_window row-wise. Returns DataFrame (index of df) with columns
    map_acc, map_pos (Int64), map_res, map_method, map_reason, map_n_hits, site_id."""
    out = []
    accs = df[acc_col].tolist()
    poss = df[pos_col].tolist()
    wins = df[win_col].tolist()
    ress = df[res_col].tolist() if res_col else [None] * len(df)
    genes = df[gene_col].tolist() if gene_col else [None] * len(df)
    cache: dict = {}
    for a, p, w, r, g in zip(accs, poss, wins, ress, genes):
        key = (a, p, w, r, g)
        m = cache.get(key)
        if m is None:
            m = store.map_window(a, p, w, r, g)
            cache[key] = m
        out.append((m.acc, m.pos, m.res, m.method, m.reason, m.n_hits, m.site_id))
    o = pd.DataFrame(out, index=df.index,
                     columns=["map_acc", "map_pos", "map_res", "map_method", "map_reason", "map_n_hits", "site_id"])
    o["map_pos"] = o["map_pos"].astype("Int64")
    return o


def add_windows(store: CanonicalStore, df: pd.DataFrame, acc_col="acc", pos_col="pos", flanks=(7, 15)) -> pd.DataFrame:
    """Add canonical windows window15 / window31 (generally window{2f+1}) for mapped rows."""
    d = df.copy()
    for f in flanks:
        d[f"window{2 * f + 1}"] = [store.window(a, int(p), f) if isinstance(a, str) and pd.notna(p) else None
                                   for a, p in zip(d[acc_col], d[pos_col])]
    return d


def sty_positions(seq: str) -> np.ndarray:
    """1-based positions of S/T/Y in seq."""
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    return np.flatnonzero((a == ord("S")) | (a == ord("T")) | (a == ord("Y"))) + 1
