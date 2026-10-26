"""Training regimes beyond the global R2 model (DESIGN §8.3 / brief: 'one model per kinase, or something in between').

All regimes reuse EXACTLY the global model's training rows (ksp.models.r2.sample.FoldData(f,'L2').sample(50, 0.5,
0.5, seed) of the frozen final config) and feature assembler (L2 cross-fitted stores), so only the pooling regime
differs. Pre-declared choices (fixed before any result was seen):
  per-kinase  N = 20 training positives (in the outer-training fold) -> own model, else fall back to the global R2
              OOF score. Applies identically to canary pseudo-kinases (so the canary gate tests the regime).
              LightGBM: binary, 15 leaves, lr 0.05, min_data_in_leaf 20, feature_fraction 0.5, bagging 0.8/1,
              lambda_l2 1, early stopping (50 rounds, max 1000, min 30) on an inner 20% substrate-cluster split of the
              kinase's training rows (macro AUROC = AUROC), then refit on all rows with the stopped round count.
              L2-logistic: median-impute + standardise, C chosen per kinase from {0.01, 0.1, 1} on the same inner
              split (training data only), refit on all rows.
              Kinase-constant columns (F7 descriptors, f6b_kinase_maxcos) are dropped (no within-kinase variance).
  per-group   one LightGBM per Manning group with the frozen final R2 hyper-parameters (127 leaves, lr .03, alpha .5
              weights renormalised within group, inner ES on real kinases), fallback to global for groups whose
              training positives < 200. Canaries (one Manning group each, A25) fall back to the global model.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

N_MIN = 20
LGB_PK = dict(objective="binary", metric="None", learning_rate=0.05, num_leaves=15, min_data_in_leaf=20,
              feature_fraction=0.5, bagging_fraction=0.8, bagging_freq=1, lambda_l2=1.0, max_bin=255,
              deterministic=True, force_col_wise=True, verbose=-1)
LR_C = (0.01, 0.1, 1.0)


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, s))


def const_cols(cols):
    return [j for j, c in enumerate(cols) if c.startswith("f7_") or c == "f6b_kinase_maxcos"]


def fit_lgb_kinase(X, y, va, seed, threads, cat_idx):
    import lightgbm as lgb
    p = {**LGB_PK, "seed": seed, "num_threads": threads}
    best = 150
    if y[va].sum() >= 3 and (~va).sum() and y[~va].sum() >= 3:
        dtr = lgb.Dataset(X[~va], y[~va], categorical_feature=cat_idx, free_raw_data=False)
        dva = lgb.Dataset(X[va], y[va], reference=dtr, categorical_feature=cat_idx, free_raw_data=False)
        yv = y[va]
        b = lgb.train(p, dtr, 1000, valid_sets=[dva],
                      feval=lambda pr, ds: [("auc", auc(yv, pr), True)],
                      callbacks=[lgb.early_stopping(50, verbose=False)])
        best = max(b.best_iteration or 150, 30)
    b = lgb.train(p, lgb.Dataset(X, y, categorical_feature=cat_idx), best)
    return b, best


def fit_lr_kinase(X, y, va, seed):
    from sklearn.linear_model import LogisticRegression
    med = np.nanmedian(X, 0)
    med = np.where(np.isnan(med), 0.0, med)

    def prep(Z, m, mu=None, sd=None):
        Z = np.where(np.isnan(Z), m, Z)
        if mu is None:
            mu, sd = Z.mean(0), Z.std(0) + 1e-6
        return (Z - mu) / sd, mu, sd
    bestC, bestA = LR_C[1], -1
    if y[va].sum() >= 3 and y[~va].sum() >= 3:
        Ztr, mu, sd = prep(X[~va], med)
        Zva, _, _ = prep(X[va], med, mu, sd)
        for C in LR_C:
            m = LogisticRegression(C=C, max_iter=300).fit(Ztr, y[~va])
            a = auc(y[va], m.decision_function(Zva))
            if a > bestA:
                bestC, bestA = C, a
    Z, mu, sd = prep(X, med)
    m = LogisticRegression(C=bestC, max_iter=500).fit(Z, y)
    return (lambda Xn: m.decision_function(prep(Xn, med, mu, sd)[0])), bestC


def qmap(score: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Rank-preserving map of `score` onto the empirical distribution of `ref` (same rows; label-free): keeps the
    per-kinase model's within-kinase order, borrows the global model's cross-kinase scale for M2."""
    r = pd.Series(score).rank(method="average").values
    q = (r - 0.5) / len(score)
    srt = np.sort(np.asarray(ref, float))          # == np.quantile(ref, q) (linear), without per-q partitions
    return np.interp(q * (len(srt) - 1), np.arange(len(srt)), srt)
