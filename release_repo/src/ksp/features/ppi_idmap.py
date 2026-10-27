"""Identifier mapping to canonical Swiss-Prot accessions for the context features (PPI, localisation, co-expression).

All context sources are mapped to the canonical UniProt Swiss-Prot accession (no isoform suffix) of release
2026_03. Mapping tables:
  * external/uniprot_topology/UP000005640_reviewed_topology.tsv.gz (reviewed set; Entrez GeneID, RefSeq NP, gene primary)
  * external/hgnc/hgnc_complete_set.txt (symbol, prev/alias symbols, Entrez, Ensembl gene, UniProt)
  * external/string/9606.protein.aliases.v12.0.txt.gz (ENSP -> UniProt_AC; only loaded on demand)

Rules (stated, not silent):
  * UniProt accessions: isoform suffix stripped; kept only if in the reviewed set.
  * Gene-level IDs (symbol, Entrez, ENSG) map to *all* reviewed accessions linked to that gene (HGNC first, then
    UniProt cross-references). A gene that maps to several reviewed accessions yields several nodes (rare).
  * Symbols: HGNC approved symbol > UniProt primary gene name > HGNC previous symbol (unique only) >
    HGNC alias symbol (unique only). Case-insensitive.
Every mapping call can record its success rate via `MapLog`.
"""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import pandas as pd

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
EXT = ROOT / "external"
TOPO = EXT / "uniprot_topology" / "UP000005640_reviewed_topology.tsv.gz"
HGNC = EXT / "hgnc" / "hgnc_complete_set.txt"
STRING_ALIASES = EXT / "string" / "9606.protein.aliases.v12.0.txt.gz"


def _split(v, sep="|"):
    if v is None or (isinstance(v, float)) or v == "":
        return []
    return [x.strip() for x in str(v).split(sep) if x.strip()]


class IdMapper:
    def __init__(self):
        u = pd.read_csv(TOPO, sep="\t", dtype=str, keep_default_na=False)
        self.reviewed = set(u["Entry"])
        self.uniprot_table = u
        h = pd.read_csv(HGNC, sep="\t", dtype=str, keep_default_na=False)
        h = h[h.status == "Approved"]
        self.sym, self.entrez, self.ensg = defaultdict(set), defaultdict(set), defaultdict(set)
        prev, alias = defaultdict(set), defaultdict(set)
        for r in h.itertuples(index=False):
            accs = {a for a in _split(r.uniprot_ids) if a in self.reviewed}
            if not accs:
                continue
            self.sym[r.symbol.upper()] |= accs
            for e in _split(r.entrez_id):
                self.entrez[e] |= accs
            for e in _split(r.ensembl_gene_id):
                self.ensg[e] |= accs
            for s in _split(r.prev_symbol):
                prev[s.upper()] |= accs
            for s in _split(r.alias_symbol):
                alias[s.upper()] |= accs
        # UniProt Entrez xrefs fill genes absent from HGNC's UniProt column
        uni_entrez = defaultdict(set)
        for acc, gid in zip(u.Entry, u.GeneID):
            for e in _split(gid, ";"):
                uni_entrez[e].add(acc)
        for e, accs in uni_entrez.items():
            if e not in self.entrez:
                self.entrez[e] = accs
        u2 = u.rename(columns={"Gene Names (primary)": "gene"})
        uprim = defaultdict(set)
        for acc, g in zip(u2.Entry, u2.gene):
            for gg in _split(g, ";"):
                uprim[gg.upper()].add(acc)
        self.uprim = uprim
        self.prev = {k: v for k, v in prev.items() if len(v) == 1}
        self.alias = {k: v for k, v in alias.items() if len(v) == 1}
        # RefSeq protein (NP_) -> acc, version-stripped
        self.refseq = defaultdict(set)
        for acc, rs in zip(u.Entry, u.RefSeq):
            for x in _split(rs, ";"):
                np_ = x.split()[0].split(".")[0]
                self.refseq[np_].add(acc)
        self._ensp = None

    # ---- single-ID mappers (return set of accessions) ----
    def from_uniprot(self, a) -> set:
        out = set()
        for x in _split(a, ";"):
            b = x.split("-")[0]
            if b in self.reviewed:
                out.add(b)
        return out

    def from_symbol(self, s) -> set:
        if not s or not isinstance(s, str):
            return set()
        s = s.strip().upper()
        for d in (self.sym, self.uprim, self.prev, self.alias):
            if s in d:
                return set(d[s])
        return set()

    def from_entrez(self, e) -> set:
        if e is None or (isinstance(e, float) and e != e):
            return set()
        e = str(e).split(".")[0]
        return set(self.entrez.get(e, ()))

    def from_ensg(self, g) -> set:
        return set(self.ensg.get(str(g).split(".")[0], ()))

    def from_refseq(self, np_) -> set:
        return set(self.refseq.get(str(np_).split(".")[0], ()))

    def from_ensp(self, p) -> set:
        if self._ensp is None:
            d = defaultdict(set)
            with gzip.open(STRING_ALIASES, "rt") as f:
                next(f)
                for line in f:
                    sid, alias, src = line.rstrip("\n").split("\t")
                    if src == "UniProt_AC" and alias in self.reviewed:
                        d[sid.split(".", 1)[1]].add(alias)
            self._ensp = d
        return set(self._ensp.get(str(p).replace("9606.", ""), ()))

    def first(self, *cands) -> set:
        """Return the first non-empty set among (fn, value) candidates."""
        for fn, v in cands:
            r = fn(v)
            if r:
                return r
        return set()


@lru_cache(maxsize=1)
def get_mapper() -> IdMapper:
    return IdMapper()


class MapLog:
    """Collects mapping rates per source; dumped to results/features/context/id_mapping_rates.json."""

    MIN_RATE = 0.80  # a source mapping below this raises (never silently drop large parts of a source)

    def __init__(self):
        self.rows = {}

    def add(self, source, what, n_total, n_mapped, note="", min_rate=None):
        rate = n_mapped / max(n_total, 1)
        self.rows[f"{source}:{what}"] = dict(n_total=int(n_total), n_mapped=int(n_mapped), n_unmapped=int(n_total - n_mapped),
                                              rate=round(rate, 4), note=note)
        if n_total and rate < (self.MIN_RATE if min_rate is None else min_rate):
            raise ValueError(f"ID mapping rate for {source}:{what} = {rate:.3f} below threshold; check the mapping")

    def dump(self, path):
        path = Path(path)
        old = json.loads(path.read_text()) if path.exists() else {}
        old.update(self.rows)
        path.write_text(json.dumps(old, indent=1, sort_keys=True))
