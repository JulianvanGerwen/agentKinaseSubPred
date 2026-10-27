"""Kinase activity inference from a site x sample matrix and a binary kinase -> site library (own implementations).

For each sample j only observed sites O_j (non-NaN) are used; kinase k is scored iff n_kj = |T_k ∩ O_j| >= min_targets
(benchmarKIN default 5), else NaN. With x = values on O_j, N = |O_j|, s_kj = sum of x over T_k ∩ O_j:

  zscore  RoKAI z-score as re-implemented in benchmarKIN `run_zscore`: s / (S * sqrt(n)), S = sd(x) (ddof 1);
          no centering (the data are fold changes / median-centred values).
  ksea    KSEA z-score (Casado et al. 2013): (s/n - mean(x)) * sqrt(n) / S.
  wmean   decoupler `wmean` raw score with unit weights: s / n (the target mean; no permutation normalisation).
  ulm     decoupler `ulm`: t-value of the slope of x ~ 1 + w over O_j, w = 1 for targets and 0 otherwise
          (closed form: b = Sxy/Sxx, se = sqrt(RSS/(N-2)/Sxx)).
All four are vectorised over samples with sparse incidence products.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

METHODS = ("zscore", "ksea", "wmean", "ulm")


def incidence(lib: pd.DataFrame, site_index: pd.Index, exclude_autophos: bool = False):
    """Sparse kinase x site (rows: sorted kinase_acc present in lib with >=1 site in site_index)."""
    d = lib[["kinase_acc", "site_id"]]
    if exclude_autophos:
        d = d[d.site_id.str.split("_").str[0] != d.kinase_acc]
    pos = pd.Series(np.arange(len(site_index)), index=site_index)
    d = d[d.site_id.isin(pos.index)].drop_duplicates()
    kin = np.array(sorted(d.kinase_acc.unique()))
    kpos = pd.Series(np.arange(len(kin)), index=kin)
    W = sparse.csr_matrix((np.ones(len(d)), (kpos[d.kinase_acc].to_numpy(), pos[d.site_id].to_numpy())),
                          shape=(len(kin), len(site_index)))
    return W, kin


def infer(mat: pd.DataFrame, lib: pd.DataFrame, method: str = "zscore", min_targets: int = 5,
          exclude_autophos: bool = False) -> pd.DataFrame:
    """mat: site_id x sample. Returns kinase_acc x sample activity (NaN = not scored)."""
    if method not in METHODS:
        raise ValueError(f"unknown method {method}; use {METHODS}")
    W, kin = incidence(lib, mat.index, exclude_autophos)
    if len(kin) == 0:
        return pd.DataFrame(index=pd.Index([], name="kinase_acc"), columns=mat.columns, dtype=float)
    X = mat.to_numpy(dtype=np.float64)
    M = np.isfinite(X)
    X0 = np.where(M, X, 0.0)
    N = M.sum(0).astype(float)                                    # per sample
    mean = X0.sum(0) / np.maximum(N, 1)
    Syy = (np.where(M, X - mean, 0.0) ** 2).sum(0)
    S = np.sqrt(Syy / np.maximum(N - 1, 1))
    n = np.asarray(W @ M.astype(float))                           # kinase x sample
    s = np.asarray(W @ X0)
    with np.errstate(divide="ignore", invalid="ignore"):
        if method == "zscore":
            A = s / (S * np.sqrt(n))
        elif method == "ksea":
            A = (s / n - mean) * np.sqrt(n) / S
        elif method == "wmean":
            A = s / n
        else:  # ulm
            Sxx = n * (N - n) / N
            Sxy = s - n * mean
            b = Sxy / Sxx
            rss = np.maximum(Syy - b * Sxy, 0.0)
            se = np.sqrt(rss / (N - 2) / Sxx)
            A = b / se
    A[(n < min_targets) | ~np.isfinite(A)] = np.nan
    return pd.DataFrame(A, index=pd.Index(kin, name="kinase_acc"), columns=mat.columns)
