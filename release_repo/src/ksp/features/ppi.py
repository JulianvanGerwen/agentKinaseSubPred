"""F3 PPI features from systematic, hypothesis-free screens only (DESIGN §7 F3).

Sources (all physical interaction screens; no curated/literature records):
  bioplex_293t   BioPlex 3.0 293T network (Huttlin 2021; external/bioplex/BioPlex_293T_Network_10K_Dec_2019.tsv)
  bioplex_hct116 BioPlex 3.0 HCT116 network (external/bioplex/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv)
  huri           HuRI Y2H (Luck 2020; external/huri/HuRI.tsv). Lit-BM and HI-union are NOT used.
  opencell       OpenCell IP-MS interactions (Cho 2022; external/opencell/opencell-protein-interactions.csv)
  buljan         Buljan 2020 kinase AP-MS, Table S2 final high-confidence interactome (mmc3.xlsx)
  hcm            Human Cell Map BioID (Go 2021) bait-prey pairs with SAINT BFDR <= 0.01

Excluded by design: STRING (all channels), BioGRID, IntAct, HIPPIE, OmniPath, SIGNOR, BioPlex 'interactionList_v4a'
(older combined set), HuRI Lit-BM / HI-union.

All edges are undirected, unweighted, self-loops removed, mapped to Swiss-Prot canonical accessions.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .ppi_idmap import EXT, ROOT, MapLog, get_mapper

SOURCES = ["bioplex_293t", "bioplex_hct116", "huri", "opencell", "buljan", "hcm"]
CACHE = ROOT / "cache" / "features" / "context"
HCM_BFDR = 0.01


def _pairs_from(a_sets, b_sets, source, log: MapLog, n_rows):
    rows = []
    n_ok = 0
    for A, B in zip(a_sets, b_sets):
        if A and B:
            n_ok += 1
            for a in A:
                for b in B:
                    if a != b:
                        rows.append((a, b) if a < b else (b, a))
    log.add(source, "edges(rows)", n_rows, n_ok, "row mapped if both endpoints map to >=1 reviewed acc")
    df = pd.DataFrame(rows, columns=["a", "b"]).drop_duplicates()
    df["source"] = source
    return df


def load_bioplex(cell: str, log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    f = {"293t": "BioPlex_293T_Network_10K_Dec_2019.tsv", "hct116": "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"}[cell]
    d = pd.read_csv(EXT / "bioplex" / f, sep="\t", dtype=str)
    A = [m.first((m.from_entrez, g), (m.from_uniprot, u), (m.from_symbol, s)) for g, u, s in zip(d.GeneA, d.UniprotA, d.SymbolA)]
    B = [m.first((m.from_entrez, g), (m.from_uniprot, u), (m.from_symbol, s)) for g, u, s in zip(d.GeneB, d.UniprotB, d.SymbolB)]
    return _pairs_from(A, B, f"bioplex_{cell}", log, len(d))


def load_huri(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    d = pd.read_csv(EXT / "huri" / "HuRI.tsv", sep="\t", header=None, names=["a", "b"], dtype=str)
    A = [m.from_ensg(x) for x in d.a]
    B = [m.from_ensg(x) for x in d.b]
    ids = pd.unique(pd.concat([d.a, d.b]))
    log.add("huri", "ensg_nodes", len(ids), sum(bool(m.from_ensg(x)) for x in ids))
    return _pairs_from(A, B, "huri", log, len(d))


def load_opencell(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    d = pd.read_csv(EXT / "opencell" / "opencell-protein-interactions.csv", dtype=str)
    meta = pd.read_csv(EXT / "opencell" / "opencell-library-metadata.csv", dtype=str)
    t_uni = dict(zip(meta.ensg_id, meta.uniprot_id))
    A = [m.first((m.from_uniprot, t_uni.get(e, "")), (m.from_ensg, e), (m.from_symbol, s))
         for e, s in zip(d.target_ensg_id, d.target_gene_name)]
    B = [m.first((m.from_ensg, e), (m.from_uniprot, u), (m.from_symbol, s))
         for e, u, s in zip(d.interactor_ensg_id, d.interactor_uniprot_ids.fillna(""), d.interactor_gene_name)]
    return _pairs_from(A, B, "opencell", log, len(d))


def buljan_table() -> pd.DataFrame:
    out = CACHE / "raw" / "buljan" / "mmc3.xlsx"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        zipfile.ZipFile(EXT / "buljan2020" / "supp_PMC7427327.zip").extract("mmc3.xlsx", out.parent)
    d = pd.read_excel(out, header=None, skiprows=2, dtype=str)
    d.columns = ["Bait_id", "Protein_id", "wdscore", "gfpratio", "NUMSPECSTOT_SUM", "NUMSPECSUNIQ_SUM"]
    return d.dropna(subset=["Bait_id", "Protein_id"])


def load_buljan(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    d = buljan_table()
    A = [m.from_uniprot(x) for x in d.Bait_id]
    B = [m.from_uniprot(x) for x in d.Protein_id]
    return _pairs_from(A, B, "buljan", log, len(d))


def load_hcm(log: MapLog, bfdr: float = HCM_BFDR) -> pd.DataFrame:
    m = get_mapper()
    s = pd.read_csv(EXT / "humancellmap" / "saint-latest.txt", sep="\t", dtype=str)
    baits = pd.read_csv(EXT / "humancellmap" / "baits-latest.txt", sep="\t", dtype=str)
    s = s[s.BFDR.astype(float) <= bfdr]
    b_uni = dict(zip(baits.bait, baits.UniProt))
    b_ent = dict(zip(baits.bait, baits.Entrez))
    bait_map = {b: m.first((m.from_uniprot, b_uni.get(b, "")), (m.from_entrez, b_ent.get(b, "")), (m.from_symbol, b))
                for b in s.Bait.unique()}
    A = [bait_map[b] for b in s.Bait]
    B = [m.first((m.from_refseq, p), (m.from_symbol, g)) for p, g in zip(s.Prey, s.PreyGene)]
    log.add("hcm", "baits", len(bait_map), sum(bool(v) for v in bait_map.values()))
    return _pairs_from(A, B, "hcm", log, len(s))


def load_all_edges(log: MapLog | None = None) -> pd.DataFrame:
    log = log or MapLog()
    parts = [load_bioplex("293t", log), load_bioplex("hct116", log), load_huri(log), load_opencell(log),
             load_buljan(log), load_hcm(log)]
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------- graph + RWR
class Graph:
    def __init__(self, edges: pd.DataFrame):
        e = edges[["a", "b"]].drop_duplicates()
        self.nodes = np.array(sorted(set(e.a) | set(e.b)))
        self.pos = {a: i for i, a in enumerate(self.nodes)}
        i = e.a.map(self.pos).to_numpy()
        j = e.b.map(self.pos).to_numpy()
        n = len(self.nodes)
        A = sp.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)).tocsr()
        A = ((A + A.T) > 0).astype(np.float64)
        A.setdiag(0)
        A.eliminate_zeros()
        self.A = A.tocsr()
        self.deg = np.asarray(self.A.sum(0)).ravel()

    def rwr_dense(self, restart: float) -> np.ndarray:
        """Full RWR matrix R (N x N), column s = stationary visiting distribution of a walk restarting at s:
        R = r (I - (1-r) W)^-1, W = A D^-1 (column-stochastic). Columns sum to 1. Returned float32."""
        n = len(self.nodes)
        W = self.A.multiply(1.0 / self.deg[None, :]).toarray()
        M = np.eye(n) - (1.0 - restart) * W
        del W
        R = restart * np.linalg.inv(M)
        del M
        R[R < 0] = 0.0  # numerical noise
        return R.astype(np.float32)

    def rwr_power(self, seeds: np.ndarray, restart: float, tol: float = 1e-10, max_iter: int = 200) -> np.ndarray:
        """Power iteration for a set of seed node indices; returns N x len(seeds)."""
        n = len(self.nodes)
        W = self.A.multiply(1.0 / self.deg[None, :]).tocsr()
        E = np.zeros((n, len(seeds)))
        E[seeds, np.arange(len(seeds))] = 1.0
        P = E.copy()
        for _ in range(max_iter):
            Pn = restart * E + (1 - restart) * (W @ P)
            if np.abs(Pn - P).sum(0).max() < tol:
                P = Pn
                break
            P = Pn
        return P
