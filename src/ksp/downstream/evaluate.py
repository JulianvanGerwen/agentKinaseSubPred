"""benchmarKIN evaluation (re-implemented from the paper's Methods and the benchmarKIN / decoupler-py 1.6 sources,
read as reference text only; nothing downloaded is executed).

Perturbation benchmark (benchmarKIN run_perturbBench / run_rank / run_phit):
  * AUROC: activities (kinase x experiment) are scaled per experiment by R `scale(center=FALSE)` (divide by the
    root-mean-square sqrt(sum x^2/(n-1)) over scored kinases), multiplied by the experiment's sign, 0 -> NaN, and
    flattened over experiments; positives = (experiment, perturbed kinase), negatives = every other scored
    (experiment, kinase) cell. rm_bg=True drops experiments in which no perturbed kinase was scored (the paper's
    "only experiments where >= 5 targets of the perturbed kinase were identified"); rm_bg=False is the package
    default. decoupler `mcauroc`: n_iter=1000 draws of n_pos negatives WITH replacement; mean (and median) AUROC.
  * Scaled rank: per experiment the (sign-adjusted, unscaled) activities are ranked (1 = most consistent with the
    perturbation); scaled rank = rank / #scored kinases; averaged per perturbed kinase over its experiments
    (average=T), then mean (paper) and median (vignette) over kinases. PHit(k) = fraction of kinases whose averaged
    rank <= k.
Tumour benchmark (benchmarKIN benchmarkROC): per cohort keep kinases with >= min_num (30) scored samples, z-score
  each kinase across samples (R scale: centre, sd ddof 1); GS+ / GS- values pooled over cohorts for kinases with
  >= 1 GS+ and >= 1 GS- value; 1000 draws of floor(0.8 * n) GS+ and floor(0.8 * n) GS- values without replacement;
  mean AUROC. matched=True mode (the package's multi-method call) restricts every compared activity set to the
  common kinases -- see api.match_*.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    ps = np.sort(pos)
    lo = np.searchsorted(ps, neg, "left")
    hi = np.searchsorted(ps, neg, "right")
    return float(((len(ps) - hi) + 0.5 * (hi - lo)).sum() / (len(ps) * len(neg)))


def mc_auroc(pos: np.ndarray, neg: np.ndarray, n_iter: int = 1000, seed: int = 42, n_neg: int | None = None,
             n_pos: int | None = None, replace: bool = True) -> np.ndarray:
    """AUROC distribution over n_iter subsamples. Default (decoupler mcauroc): all positives, n_pos negatives drawn
    with replacement. With n_pos/n_neg and replace=False: benchmarKIN tumour .ROC_sampler."""
    rng = np.random.default_rng(seed)
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.array([np.nan])
    out = np.empty(n_iter)
    for i in range(n_iter):
        p = pos if n_pos is None else rng.choice(pos, n_pos, replace=replace)
        k = len(pos) if n_neg is None else n_neg
        q = rng.choice(neg, k, replace=replace)
        out[i] = auroc(p, q)
    return out


# ------------------------------------------------------------------------------------------------ perturbation
def _rms_scale(A: pd.DataFrame) -> pd.DataFrame:
    n = A.notna().sum(0)
    rms = np.sqrt((A ** 2).sum(0) / (n - 1).clip(lower=1))
    return A / rms.replace(0, np.nan)


def perturbation_ranks(act: pd.DataFrame, bench) -> pd.DataFrame:
    """One row per (experiment, perturbed kinase) with a score: rank, n_scored, scaled_rank."""
    rows = []
    sign = dict(zip(bench.experiments.experiment, bench.experiments.sign))
    for e, grp in bench.meta.groupby("experiment", sort=False):
        if e not in act.columns:
            continue
        a = act[e].dropna() * sign[e]
        if a.empty:
            continue
        order = a.sort_values(ascending=False, kind="stable")
        rk = pd.Series(np.arange(1, len(order) + 1), index=order.index)
        for k, g in zip(grp.kinase_acc, grp.gene):
            if k in rk.index:
                rows.append((e, k, g, int(rk[k]), len(order), rk[k] / len(order)))
    return pd.DataFrame(rows, columns=["experiment", "kinase_acc", "gene", "rank", "n_scored", "scaled_rank"])


def perturbation_eval(act: pd.DataFrame, bench, rm_bg: bool = True, n_iter: int = 1000, seed: int = 42,
                      ks=(5, 10, 20)) -> tuple[dict, pd.DataFrame]:
    exps = [e for e in bench.experiments.experiment if e in act.columns]
    A = act.reindex(columns=exps)
    sign = bench.experiments.set_index("experiment").sign.reindex(exps)
    S = _rms_scale(A) * sign.to_numpy()[None, :]
    S = S.mask(S == 0)
    G = pd.DataFrame(False, index=S.index, columns=S.columns)
    for e, k in zip(bench.meta.experiment, bench.meta.kinase_acc):
        if e in G.columns and k in G.index:
            G.loc[k, e] = True
    if rm_bg:
        has = (G & S.notna()).any(axis=0)
        S, G = S.loc[:, has], G.loc[:, has]
    v, g = S.to_numpy(dtype=float).ravel(), G.to_numpy(dtype=bool).ravel()
    ok = np.isfinite(v)
    pos, neg = v[ok & g], v[ok & ~g]
    dist = mc_auroc(pos, neg, n_iter=n_iter, seed=seed)
    rk = perturbation_ranks(act, bench)
    per_k = rk.groupby(["kinase_acc", "gene"]).agg(rank=("rank", "mean"), scaled_rank=("scaled_rank", "mean"),
                                                   n_exp=("experiment", "nunique")).reset_index()
    out = {"auroc_mean": float(np.nanmean(dist)), "auroc_median": float(np.nanmedian(dist)),
           "auroc_sd": float(np.nanstd(dist, ddof=1)) if len(dist) > 1 else np.nan,
           "auroc_full": auroc(pos, neg), "n_pos": int(len(pos)), "n_neg": int(len(neg)),
           "n_experiments_auroc": int(S.shape[1]), "n_experiments_ranked": int(rk.experiment.nunique()),
           "n_kinases_scored": int(S.notna().any(axis=1).sum()),
           "n_perturbed_kinases_evaluable": int(len(per_k)),
           "scaled_rank_mean": float(per_k.scaled_rank.mean()) if len(per_k) else np.nan,
           "scaled_rank_median": float(per_k.scaled_rank.median()) if len(per_k) else np.nan,
           "scaled_rank_mean_pairs": float(rk.scaled_rank.mean()) if len(rk) else np.nan}
    for k in ks:
        out[f"phit{k}"] = float((per_k["rank"] <= k).mean()) if len(per_k) else np.nan
    return out, per_k


# ------------------------------------------------------------------------------------------------ tumour
def tumour_values(acts: dict, gs: pd.DataFrame, min_num: int = 30, kinases: dict | None = None) -> pd.DataFrame:
    """acts: cohort -> kinase x sample activity. Returns rows (cohort, kinase_acc, gene, sample, label, z)."""
    rows = []
    for c, A in acts.items():
        A = A[A.notna().sum(1) >= min_num]
        if kinases is not None:
            A = A[A.index.isin(kinases.get(c, set()))]
        Z = A.sub(A.mean(1), axis=0).div(A.std(1, ddof=1), axis=0)
        g = gs[gs.cohort.eq(c) & gs.kinase_acc.isin(Z.index) & gs["sample"].isin(Z.columns)]
        if g.empty:
            continue
        z = Z.to_numpy()[Z.index.get_indexer(g.kinase_acc), Z.columns.get_indexer(g["sample"])]
        rows.append(g.assign(z=z))
    if not rows:
        return pd.DataFrame(columns=["cohort", "kinase_acc", "gene", "sample", "label", "z"])
    d = pd.concat(rows, ignore_index=True).dropna(subset=["z"])
    both = d.groupby(["cohort", "kinase_acc"]).label.agg(lambda s: s.nunique() == 2)
    keep = set(both[both].index)
    return d[[(c, k) in keep for c, k in zip(d.cohort, d.kinase_acc)]].reset_index(drop=True)


def tumour_eval(acts: dict, gs: pd.DataFrame, min_num: int = 30, sub_frac: float = 0.8, n_iter: int = 1000,
                seed: int = 1, kinases: dict | None = None) -> tuple[dict, pd.DataFrame]:
    d = tumour_values(acts, gs, min_num, kinases)
    pos, neg = d.z[d.label.eq(1)].to_numpy(), d.z[d.label.eq(0)].to_numpy()
    dist = mc_auroc(pos, neg, n_iter=n_iter, seed=seed, n_pos=int(np.floor(sub_frac * len(pos))),
                    n_neg=int(np.floor(sub_frac * len(neg))), replace=False) if len(pos) and len(neg) else \
        np.array([np.nan])
    per_k = d.groupby(["kinase_acc", "gene"]).agg(n_cohorts=("cohort", "nunique"), n_pos=("label", "sum"),
                                                  n_pairs=("label", "size")).reset_index()
    out = {"auroc_mean": float(np.nanmean(dist)), "auroc_median": float(np.nanmedian(dist)),
           "auroc_sd": float(np.nanstd(dist, ddof=1)) if len(dist) > 1 else np.nan, "auroc_full": auroc(pos, neg),
           "n_pos": int(len(pos)), "n_neg": int(len(neg)), "n_kinases": int(d.kinase_acc.nunique()),
           "n_kinase_cohorts": int(d.groupby(["cohort", "kinase_acc"]).ngroups)}
    return out, per_k
