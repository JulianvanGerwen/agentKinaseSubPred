"""Paired comparisons of two libraries on matched kinase sets, with paired bootstrap CIs.

Matched set: perturbation = the (experiment, kinase) cells scored by BOTH libraries (so each experiment's ranking
and negatives are over the kinases evaluable by both); tumour = per cohort the kinases with >= 30 scored tumours
in BOTH, and (cohort, kinase, tumour) GS values present in both.
Metrics on the matched set (same preprocessing as evaluate.py):
  perturbation AUROC (full; equals the expectation of the paper's balanced-subsampled mcauroc) -- bootstrap unit
  = experiment (clusters of cells sharing one perturbation); kinase-averaged scaled rank -- bootstrap unit =
  perturbed kinase; tumour AUROC -- bootstrap unit = kinase (all its GS pairs across cohorts).
B = 1000 resamples, identical draws for both libraries (paired); percentile 95% CI of the difference b - a and
the bootstrap fraction of draws with difference <= 0 (one-sided, descriptive; no multiplicity correction --
downstream endpoints are secondary, DESIGN §6.7).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import evaluate


def _prep_auc(scores: np.ndarray, labels: np.ndarray):
    o = np.argsort(scores, kind="mergesort")
    s, y = scores[o], labels[o]
    starts = np.r_[0, np.flatnonzero(np.diff(s)) + 1]
    return o, y, starts


def _wauc(prep, w: np.ndarray) -> float:
    o, y, starts = prep
    ws = w[o]
    wp = np.add.reduceat(ws * y, starts)
    wn = np.add.reduceat(ws * (1 - y), starts)
    P, N = wp.sum(), wn.sum()
    if P <= 0 or N <= 0:
        return np.nan
    below = np.cumsum(wn) - wn
    return float((wp * (below + 0.5 * wn)).sum() / (P * N))


def _summ(a: np.ndarray, b: np.ndarray, pa: float, pb: float, prefix: str) -> dict:
    d = b - a
    ok = np.isfinite(d)
    lo, hi = (np.percentile(d[ok], [2.5, 97.5]) if ok.any() else (np.nan, np.nan))
    return {f"{prefix}_a": pa, f"{prefix}_b": pb, f"{prefix}_diff": pb - pa, f"{prefix}_ci_lo": float(lo),
            f"{prefix}_ci_hi": float(hi), f"{prefix}_frac_le0": float((d[ok] <= 0).mean()) if ok.any() else np.nan}


def _signed_scaled(A: pd.DataFrame, bench) -> pd.DataFrame:
    sign = bench.experiments.set_index("experiment").sign.reindex(A.columns)
    S = evaluate._rms_scale(A) * sign.to_numpy()[None, :]
    return S.mask(S == 0)


def perturbation_paired(act_a: pd.DataFrame, act_b: pd.DataFrame, bench, B: int = 1000, seed: int = 0) -> dict:
    exps = [e for e in bench.experiments.experiment if e in act_a.columns and e in act_b.columns]
    idx = act_a.index.intersection(act_b.index)
    A, Bm = act_a.reindex(index=idx, columns=exps), act_b.reindex(index=idx, columns=exps)
    m = A.notna() & Bm.notna()
    A, Bm = A.where(m), Bm.where(m)
    G = pd.DataFrame(False, index=idx, columns=exps)
    for e, k in zip(bench.meta.experiment, bench.meta.kinase_acc):
        if e in G.columns and k in G.index:
            G.loc[k, e] = True
    has = (G & m).any(axis=0)                                   # rm_bg on the matched set
    exps = list(has[has].index)
    out = {"n_experiments": len(exps)}
    if not exps:
        return out
    A, Bm, G, m = A[exps], Bm[exps], G[exps], m[exps]
    Sa, Sb = _signed_scaled(A, bench), _signed_scaled(Bm, bench)
    cell = (Sa.notna() & Sb.notna()).to_numpy()
    ei = np.broadcast_to(np.arange(len(exps)), cell.shape)[cell]
    y = G.to_numpy()[cell].astype(float)
    sa, sb = Sa.to_numpy()[cell], Sb.to_numpy()[cell]
    out.update(n_pos=int(y.sum()), n_neg=int((1 - y).sum()), n_kinases_scored=int(m.any(axis=1).sum()))
    pa, pb = _prep_auc(sa, y), _prep_auc(sb, y)
    rng = np.random.default_rng(seed)
    one = np.ones(len(y))
    ba, bb = np.empty(B), np.empty(B)
    for i in range(B):
        c = np.bincount(rng.integers(0, len(exps), len(exps)), minlength=len(exps)).astype(float)
        w = c[ei]
        ba[i], bb[i] = _wauc(pa, w), _wauc(pb, w)
    out.update(_summ(ba, bb, _wauc(pa, one), _wauc(pb, one), "auroc"))
    # scaled rank on the matched set, kinase-averaged; bootstrap over perturbed kinases
    sub = type(bench)(bench.data, bench.meta[bench.meta.experiment.isin(exps)], bench.experiments[
        bench.experiments.experiment.isin(exps)], bench.mapping)
    ra = evaluate.perturbation_ranks(A, sub).groupby("kinase_acc").scaled_rank.mean()
    rb = evaluate.perturbation_ranks(Bm, sub).groupby("kinase_acc").scaled_rank.mean()
    kk = ra.index.intersection(rb.index)
    ra, rb = ra[kk].to_numpy(), rb[kk].to_numpy()
    out["n_perturbed_kinases"] = int(len(kk))
    if len(kk):
        rng = np.random.default_rng(seed + 1)
        draws = rng.integers(0, len(kk), (B, len(kk)))
        out.update(_summ(ra[draws].mean(1), rb[draws].mean(1), float(ra.mean()), float(rb.mean()), "scaled_rank"))
        out["perturbed_kinases"] = ";".join(sorted(kk))
    return out


def tumour_paired(acts_a: dict, acts_b: dict, gs: pd.DataFrame, B: int = 1000, seed: int = 0,
                  min_num: int = 30) -> dict:
    common = {c: set(acts_a[c].index[acts_a[c].notna().sum(1) >= min_num]) &
              set(acts_b[c].index[acts_b[c].notna().sum(1) >= min_num]) for c in acts_a if c in acts_b}
    va = evaluate.tumour_values(acts_a, gs, min_num, common)
    vb = evaluate.tumour_values(acts_b, gs, min_num, common)
    key = ["cohort", "kinase_acc", "sample", "label"]
    v = va[key + ["z"]].merge(vb[key + ["z"]], on=key, suffixes=("_a", "_b"))
    out = {"n_kinases": int(v.kinase_acc.nunique()), "n_kinase_cohorts": int(v.groupby(["cohort", "kinase_acc"]).ngroups),
           "n_pos": int(v.label.sum()), "n_neg": int((1 - v.label).sum())}
    if v.empty or v.label.nunique() < 2:
        return out
    kin = v.kinase_acc.astype("category").cat.codes.to_numpy()
    nk = kin.max() + 1
    y = v.label.to_numpy(float)
    pa, pb = _prep_auc(v.z_a.to_numpy(), y), _prep_auc(v.z_b.to_numpy(), y)
    rng = np.random.default_rng(seed)
    ba, bb = np.empty(B), np.empty(B)
    for i in range(B):
        c = np.bincount(rng.integers(0, nk, nk), minlength=nk).astype(float)
        w = c[kin]
        ba[i], bb[i] = _wauc(pa, w), _wauc(pb, w)
    one = np.ones(len(y))
    out.update(_summ(ba, bb, _wauc(pa, one), _wauc(pb, one), "auroc"))
    return out
