"""F5 co-expression / co-dependency (low priority; DESIGN §7 F5). Gene-level data mapped to Swiss-Prot accessions.

  hpa_tissue  HPA v25 RNA tissue consensus (external/hpa/rna_tissue_consensus.tsv.zip): log2(nTPM+1), genes x tissues;
              Spearman correlation across tissues (rank per gene, then Pearson).
  depmap      DepMap 26Q1 Chronos gene effect (external/depmap/DepMap_26Q1_chronos_gene_effect.csv; cell lines x genes):
              Pearson co-essentiality (NaN cells mean-imputed per gene); variant with the top 5 principal components
              removed (screen-quality / lineage confounders), since raw co-essentiality is dominated by a few PCs.
A gene mapping to several accessions gives each accession the same profile; several genes mapping to one accession:
the first (HGNC order) is kept.
"""
from __future__ import annotations

import re
import zipfile

import numpy as np
import pandas as pd

from .ppi_idmap import EXT, MapLog, get_mapper


def hpa_tissue_profiles(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    with zipfile.ZipFile(EXT / "hpa" / "rna_tissue_consensus.tsv.zip") as z:
        d = pd.read_csv(z.open("rna_tissue_consensus.tsv"), sep="\t")
    X = d.pivot_table(index=["Gene", "Gene name"], columns="Tissue", values="nTPM", aggfunc="mean")
    X = np.log2(X.fillna(0) + 1)
    rows, n_ok = {}, 0
    for (g, s), v in zip(X.index, X.to_numpy()):
        accs = m.first((m.from_ensg, g), (m.from_symbol, s))
        if accs:
            n_ok += 1
        for a in accs:
            rows.setdefault(a, v)
    log.add("hpa_rna_tissue", "genes", len(X), n_ok, min_rate=0.5)
    P = pd.DataFrame.from_dict(rows, orient="index", columns=list(X.columns))
    return P.rank(axis=1)  # Spearman via ranks


def depmap_profiles(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    d = pd.read_csv(EXT / "depmap" / "DepMap_26Q1_chronos_gene_effect.csv", index_col=0)
    rows, n_ok = {}, 0
    for col in d.columns:
        mm = re.match(r"^(.*) \((\d+)\)$", col)
        sym, ent = (mm.group(1), mm.group(2)) if mm else (col, None)
        accs = m.first((m.from_entrez, ent), (m.from_symbol, sym))
        if accs:
            n_ok += 1
        for a in accs:
            rows.setdefault(a, col)
    log.add("depmap_chronos", "genes", d.shape[1], n_ok)
    cols = list(dict.fromkeys(rows.values()))
    X = d[cols].to_numpy(float)
    X = np.where(np.isnan(X), np.nanmean(X, 0, keepdims=True), X)
    P = pd.DataFrame(X.T, index=cols, columns=d.index)
    acc_rows = {a: c for a, c in rows.items()}
    out = P.loc[list(acc_rows.values())]
    out.index = list(acc_rows.keys())
    return out


def remove_top_pcs(P: pd.DataFrame, k: int) -> pd.DataFrame:
    X = P.to_numpy(float)
    X = X - X.mean(1, keepdims=True)
    U, S, Vt = np.linalg.svd(X - X.mean(0, keepdims=True), full_matrices=False)
    Xr = X - (U[:, :k] * S[:k]) @ Vt[:k]
    return pd.DataFrame(Xr, index=P.index, columns=P.columns)
