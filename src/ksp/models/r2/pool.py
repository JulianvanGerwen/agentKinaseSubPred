"""R2 partial pooling (DESIGN §8.3, ALG-3): final logit = R2 logit + w_k * r_k(x), w_k = n_k / (n_k + lambda).

r_k = per-kinase ridge logistic RESIDUAL on a small fixed feature vector, with the R2 logit as offset, coefficients
shrunk towards the Manning-family mean vector (penalty mu * ||beta_k - beta_fam||^2, beta_fam fitted first as a pooled
family residual model with ridge mu).
The offset on TRAINING rows must be out-of-sample: inner K-fold over training clusters (ksp.splits.inner_folds), R2
refitted on each inner complement with the outer config's n_rounds, predicting the held-out inner rows. Eval rows use
the outer R2 OOF logit. lambda = inf -> global R2 only; lambda = 0 -> full per-kinase residual.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

POOL_FEATURES = ["f1_unprimed_pct_universe", "f2a_pssm_ph", "f6a_comp_med", "f6b_ridge_z", "f3_rwr_dc",
                 "f4_loc_cos_combined", "f0_idr", "f0_plddt"]


def _prep(X: np.ndarray, mu=None, sd=None):
    """NaN -> column mean (0 after standardisation) + missing indicators are NOT added (keeps 8 dims)."""
    if mu is None:
        mu = np.nanmean(X, 0)
        sd = np.nanstd(X, 0) + 1e-6
    Z = (X - mu) / sd
    Z = np.where(np.isnan(Z), 0.0, Z)
    return Z, mu, sd


def fit_ridge_offset(Z, y, off, w, prior, mu_pen):
    """argmin sum w * logloss(y, off + Z b + c) + mu_pen * ||b - prior||^2 (intercept c unpenalised)."""
    d = Z.shape[1]

    def f(theta):
        b, c = theta[:d], theta[d]
        eta = off + Z @ b + c
        p = 1 / (1 + np.exp(-eta))
        ll = w * (np.logaddexp(0, eta) - y * eta)
        g = w * (p - y)
        pen = mu_pen * np.sum((b - prior) ** 2)
        grad = np.r_[Z.T @ g + 2 * mu_pen * (b - prior), g.sum()]
        return ll.sum() + pen, grad
    r = minimize(f, np.r_[prior, 0.0], jac=True, method="L-BFGS-B")
    return r.x[:d], r.x[d]


def fit_residuals(train: pd.DataFrame, Xtr: np.ndarray, family: dict, mu_pen: float = 10.0):
    """train: DataFrame[kinase_idx, y, w, offset]; Xtr: pool features. Returns dict k -> (b, c), standardisation."""
    Z, m, s = _prep(Xtr)
    y, w, off = train.y.values.astype(float), train.w.values, train.offset.values
    fam_of = train.kinase_idx.map(family).values
    coefs_f = {}
    for f in np.unique(fam_of):
        ix = fam_of == f
        coefs_f[f] = fit_ridge_offset(Z[ix], y[ix], off[ix], w[ix], np.zeros(Z.shape[1]), mu_pen)
    out = {}
    for k, ix in train.groupby("kinase_idx").indices.items():
        bf, _ = coefs_f[family[k]]
        out[k] = fit_ridge_offset(Z[ix], y[ix], off[ix], w[ix], bf, mu_pen)
    return out, (m, s), coefs_f


def apply_residuals(kinase_idx, X, off, coefs, std, n_train: dict, lam: float):
    Z, _, _ = _prep(X, *std)
    out = off.astype(float).copy()
    for k in np.unique(kinase_idx):
        if k not in coefs:
            continue
        ix = kinase_idx == k
        b, c = coefs[k]
        wk = n_train.get(k, 0) / (n_train.get(k, 0) + lam) if np.isfinite(lam) else 0.0
        out[ix] = off[ix] + wk * (Z[ix] @ b + c)
    return out
