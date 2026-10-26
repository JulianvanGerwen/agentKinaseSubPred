"""R2 tiered calibration (DESIGN §8.5, A16) and site-level assignment probability (A17), on EXHAUSTIVE pooled OOF
predictions over DEV U_ms (no sampling -> no case-control offset needed).

Target: y(k,s) = 1 if (k,s) is a DEV evaluation positive (human-native PSP, in class, not autophos), 0 for E-all
negatives (masks applied, ksp.eval definitions). Under PU labels this calibrates P(curated positive | score), a
lower bound of P(true substrate | score) - stated with every number.
Tiers:   n_pos(k) >= 30 -> per-kinase Platt (a_k + b_k * score)
         else            -> family-level slope b_fam (pooled Platt over the family's kinases, per-kinase centring)
                            + per-kinase intercept fitted with b fixed (kinases with 0 positives: family-mean
                            intercept; family without positives: global slope/intercept).
Honest evaluation: calibrators are CROSS-FITTED over CV folds (fit on the OOF rows of folds != f, apply to fold f),
so reliability/ECE are out-of-sample.
Assignment (A17): p(k|s) = softmax_k(logit_cal(k,s) / T) over real kinases of matching class at site s;
T fitted on labelled DEV sites (OOF) by maximising mean log sum_{k in true(s)} p(k|s), cross-fitted the same way.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from sklearn.linear_model import LogisticRegression

from ksp import splits


def _logit_fit(x, y, w=None, C=1e4):
    if y.sum() == 0 or y.sum() == len(y):
        return np.nan, np.nan
    m = LogisticRegression(C=C, max_iter=200).fit(x.reshape(-1, 1), y, sample_weight=w)
    return float(m.intercept_[0]), float(m.coef_[0, 0])


def _intercept_given_slope(x, y, b):
    """MLE intercept with fixed slope (1-d Newton)."""
    a = np.log((y.sum() + 0.5) / (len(y) - y.sum() + 0.5))
    for _ in range(50):
        p = 1 / (1 + np.exp(-(a + b * x)))
        g = (y - p).sum()
        h = (p * (1 - p)).sum() + 1e-9
        a += g / h
        if abs(g / h) < 1e-6:
            break
    return float(a)


def table_from_eval(P: pd.DataFrame, data, kinases=None, neg_subsample: int | None = 20000, seed: int = 0):
    """Rows for calibration: per real kinase, all positives + E-all negatives (optionally subsampled with weights
    so that fits stay fast; weights restore the exhaustive negative mass)."""
    rng = np.random.default_rng(seed)
    nS = len(data.sites)
    res = data.sites.res.values
    kacc = data.kinases.uniprot_acc.values
    kcls = data.kinases.residue_class.values
    P = P.sort_values("kinase_idx", kind="stable")
    kb = np.searchsorted(P.kinase_idx.values, np.arange(len(kacc) + 1))
    rows = []
    from ksp.eval import _res_ok
    for i, k in enumerate(kacc):
        if k.startswith("CANARY_") or (kinases is not None and k not in kinases):
            continue
        sl = slice(kb[i], kb[i + 1])
        if kb[i + 1] == kb[i]:
            continue
        s = P.site_idx.values[sl]
        v = P.score.values[sl]
        cand = data.universe & _res_ok(res, kcls[i])
        excl = np.zeros(nS, bool)
        excl[data.masks.excluded(k)] = True
        pos = np.zeros(nS, bool)
        pos[data.positives.get(k, np.zeros(0, np.int64))] = True
        keep = cand[s] & (pos[s] | ~excl[s])
        s, v = s[keep], v[keep]
        y = pos[s].astype(np.int8)
        w = np.ones(len(s))
        neg = np.flatnonzero(y == 0)
        if neg_subsample and len(neg) > neg_subsample:
            sel = rng.choice(neg, neg_subsample, replace=False)
            keepi = np.r_[np.flatnonzero(y == 1), sel]
            w = np.r_[np.ones((y == 1).sum()), np.full(len(sel), len(neg) / len(sel))]
            s, v, y = s[keepi], v[keepi], y[keepi]
        rows.append(pd.DataFrame({"kinase_idx": i, "site_idx": s, "score": v, "y": y, "w": w}))
    T = pd.concat(rows, ignore_index=True)
    T["cluster"] = data.site_cluster[T.site_idx.values]
    return T


def fit_calibrator(T: pd.DataFrame, kinases: pd.DataFrame, min_pos: int = 30) -> pd.DataFrame:
    """Per-kinase (a_k, b_k, tier) from calibration rows T (kinase_idx, score, y, w)."""
    fam = kinases.manning_family.fillna("(none)").values
    out = {}
    g_a, g_b = _logit_fit(T.score.values, T.y.values, T.w.values)
    by_k = dict(tuple(T.groupby("kinase_idx")))
    npos = {k: int(d.y.sum()) for k, d in by_k.items()}
    # family slopes from per-kinase-centred scores
    fam_b = {}
    for f in np.unique(fam):
        ks = [k for k in by_k if fam[k] == f]
        if not ks:
            continue
        d = pd.concat([by_k[k] for k in ks])
        if d.y.sum() < 5:
            continue
        mu = d.groupby("kinase_idx").score.transform("mean")
        _, b = _logit_fit((d.score - mu).values, d.y.values, d.w.values)
        if np.isfinite(b) and b > 0:
            fam_b[f] = b
    for k in range(len(kinases)):
        if k not in by_k:
            continue
        d = by_k[k]
        n = npos[k]
        if n >= min_pos:
            a, b = _logit_fit(d.score.values, d.y.values, d.w.values)
            tier = "kinase"
            if not (np.isfinite(b) and b > 0):
                a, b, tier = np.nan, np.nan, None
        else:
            a = b = np.nan
            tier = None
        if tier is None:
            b = fam_b.get(fam[k], g_b)
            tier = "family_slope" if fam[k] in fam_b else "global_slope"
            if n > 0:
                a = _intercept_given_slope(d.score.values, d.y.values, b) if d.w.nunique() == 1 else \
                    _intercept_given_slope_w(d.score.values, d.y.values, d.w.values, b)
            else:
                a = np.nan
                tier += "+no_pos"
        out[k] = (a, b, tier, n)
    C = pd.DataFrame.from_dict(out, orient="index", columns=["a", "b", "tier", "n_pos"])
    # intercepts for kinases without positives: family mean intercept, else global
    fam_a = C.assign(f=[fam[k] for k in C.index]).dropna(subset=["a"]).groupby("f").a.median()
    for k in C.index[C.a.isna()]:
        C.loc[k, "a"] = fam_a.get(fam[k], g_a)
    C["family"] = [fam[k] for k in C.index]
    C.attrs["global"] = (g_a, g_b)
    return C


def _intercept_given_slope_w(x, y, w, b):
    a = 0.0
    for _ in range(100):
        p = 1 / (1 + np.exp(-(a + b * x)))
        g = (w * (y - p)).sum()
        h = (w * p * (1 - p)).sum() + 1e-9
        a += g / h
        if abs(g / h) < 1e-7:
            break
    return float(a)


def apply(C: pd.DataFrame, kinase_idx: np.ndarray, score: np.ndarray) -> np.ndarray:
    ga, gb = C.attrs.get("global", (np.nan, np.nan))
    a = C.a.reindex(kinase_idx).fillna(ga).values
    b = C.b.reindex(kinase_idx).fillna(gb).values
    return a + b * score


def reliability(p: np.ndarray, y: np.ndarray, w: np.ndarray, bins=(0, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.5, 1.0)):
    out = []
    b = np.digitize(p, bins[1:-1])
    for i in range(len(bins) - 1):
        m = b == i
        if not m.any():
            continue
        out.append({"bin": f"[{bins[i]},{bins[i + 1]})", "n": float(w[m].sum()), "mean_pred": float(np.average(p[m], weights=w[m])),
                    "obs_rate": float(np.average(y[m], weights=w[m])), "n_pos": int(y[m].sum())})
    R = pd.DataFrame(out)
    ece = float((R.n * (R.mean_pred - R.obs_rate).abs()).sum() / R.n.sum()) if len(R) else np.nan
    return R, ece


# ------------------------------------------------------------------------------------------------ assignment (A17)
def site_softmax(site_idx: np.ndarray, logit: np.ndarray, T: float):
    """p(k|s) over the given rows grouped by site (rows must be sorted by site_idx)."""
    b = np.flatnonzero(np.r_[True, site_idx[1:] != site_idx[:-1]])
    z = logit / T
    mx = np.maximum.reduceat(z, b)
    cnt = np.diff(np.r_[b, len(z)])
    e = np.exp(z - np.repeat(mx, cnt))
    den = np.add.reduceat(e, b)
    return e / np.repeat(den, cnt)


def fit_temperature(site_idx, logit, is_true, grid=None):
    """T maximising mean over labelled sites of log sum_{true k} p(k|s). Rows sorted by site; only sites with >=1
    true kinase are used."""
    b = np.flatnonzero(np.r_[True, site_idx[1:] != site_idx[:-1]])
    has = np.add.reduceat(is_true.astype(float), b) > 0
    cnt = np.diff(np.r_[b, len(site_idx)])
    keep = np.repeat(has, cnt)
    s, z, t = site_idx[keep], logit[keep], is_true[keep]

    def nll(logT):
        p = site_softmax(s, z, np.exp(logT))
        bb = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
        pt = np.add.reduceat(np.where(t, p, 0.0), bb)
        return -np.mean(np.log(np.clip(pt, 1e-300, None)))
    r = minimize_scalar(nll, bounds=(np.log(0.05), np.log(50)), method="bounded")
    return float(np.exp(r.x)), float(r.fun)
