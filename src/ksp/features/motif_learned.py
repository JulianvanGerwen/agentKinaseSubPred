"""F2a / B2: learned per-kinase position-specific log-odds motif (label-derived; must be cross-fitted).

Model (per kinase k and residue type rt in {ST, Y} within k's residue class):
  window  : ±7 around the site (15 positions incl. the centre, so S-vs-T preference is learned).
  alphabet: 20 aa + '_' (outside protein) + 'X' (other) [+ 's','t','y' if phospho=True: a FLANK S/T/Y
            that is itself a U_obs phosphosite; the query site's own status is never used].
  q_k[p,a]: empirical-Bayes shrunk frequencies (Dirichlet-multinomial, hierarchical, leave-k-out priors)
            q_group(-k)  = (C_group - C_k + gamma*bg) / (n_group - n_k + gamma)
            q_family(-k) = (C_family - C_k + beta*q_group(-k)) / (n_family - n_k + beta)
            q_k          = (C_k + alpha*q_prior) / (n_k + alpha)
            prior='family' -> q_prior = q_family(-k); 'group' -> q_group(-k); 'none' -> bg.
            Kinases with 0 training labels get q_k = q_prior (family -> group -> background fallback).
  bg[p,a] : frequencies over U_ms windows of the same residue type (label-free).
  score   : sum_p log2(q_k[p, a_p] / bg[p, a_p]).
  alpha, beta (gamma fixed) are tuned by a quick internal CV over training pairs grouped by substrate
  protein (or by the caller's cluster ids): maximise the macro-average (over kinases) held-out
  log-likelihood of positive windows under q_k.

Pure functions `fit_pssm(...)` / `score(model, site_idx)` plus the ksp.crossfit contract class
`LearnedPSSM` (name 'f2a_pssm'). No labels are read here: the caller passes training pairs.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ksp.index import load_kinases, load_sites

BASE = "ACDEFGHIKLMNPQRSTVWY"
ALPHA = BASE + "_X" + "sty"          # indices 0..19 aa, 20 '_', 21 'X', 22..24 phospho flank
A_UNPH, A_PH = 22, 25
RTS = ("ST", "Y")
RT_RES = {"ST": ("S", "T"), "Y": ("Y",)}
CLASS_RTS = {"ST": ("ST",), "Y": ("Y",), "STY": ("ST", "Y")}
GAMMA = 10.0
ALPHA_GRID = (0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
BETA_GRID = (1, 4, 16, 64, 256, 1024)

_cache: dict = {}


def _site_codes(sites: pd.DataFrame, phospho: bool) -> np.ndarray:
    """N x 15 int8 codes of window15 (cached per (len(sites), phospho))."""
    key = ("codes", len(sites), phospho, id(sites))
    if key in _cache:
        return _cache[key]
    lut = np.full(256, 21, dtype=np.int8)
    for i, a in enumerate(BASE):
        lut[ord(a)] = i
    lut[ord("_")] = 20
    w = np.asarray(sites.window15.to_numpy(), dtype="U15")
    b = np.frombuffer(np.char.encode(w, "ascii"), dtype=np.uint8).reshape(len(w), 15)
    codes = lut[b]
    if phospho:
        known = set(sites.site_id)
        acc, pos = sites.acc.to_numpy(), sites.pos.to_numpy()
        for j, off in enumerate(range(-7, 8)):
            if off == 0:
                continue
            col = b[:, j]
            for ch, c in ((ord("S"), 22), (ord("T"), 23), (ord("Y"), 24)):
                idx = np.where(col == ch)[0]
                hit = np.array([f"{acc[i]}_{chr(ch)}{pos[i] + off}" in known for i in idx], dtype=bool)
                codes[idx[hit], j] = c
    _cache[key] = codes
    return codes


def _counts(codes: np.ndarray, A: int) -> np.ndarray:
    """codes (n x 15) -> 15 x A count matrix."""
    C = np.zeros((15, A))
    for p in range(15):
        C[p] = np.bincount(codes[:, p], minlength=A)[:A]
    return C


@dataclass
class PSSMModel:
    kinase_acc: np.ndarray               # index order (rows of score output)
    phospho: bool
    prior: str
    alpha: float
    beta: float
    gamma: float
    logodds: dict                        # rt -> K x 15 x A float32 (NaN rows where rt not in class)
    n_train: dict                        # rt -> K int
    level: dict                          # rt -> K array of 'kinase'|'family'|'group'|'background'|''
    bg: dict                             # rt -> 15 x A
    cv: pd.DataFrame | None = None       # tuning grid results
    info: dict = field(default_factory=dict)


def _hier(Ck, nk, kin, bg, alpha, beta, gamma, prior):
    """Vectorised hierarchical shrinkage. Ck: K x 15 x A counts, nk: K. Returns q (K x 15 x A), level."""
    fam = kin.manning_family.fillna("").to_numpy()
    grp = kin.manning_group.fillna("").to_numpy()
    K = len(kin)
    Cg = np.zeros_like(Ck); ng = np.zeros(K)
    Cf = np.zeros_like(Ck); nf = np.zeros(K)
    for key, Cx, nx in ((grp, Cg, ng), (fam, Cf, nf)):
        s = pd.Series(np.arange(K)).groupby(key)
        for g, idx in s:
            idx = idx.to_numpy()
            if g == "":            # unknown family/group: no pooling
                Cx[idx] = Ck[idx]; nx[idx] = nk[idx]
                continue
            Cx[idx] = Ck[idx].sum(0)
            nx[idx] = nk[idx].sum()
    Cg -= Ck; ng -= nk; Cf -= Ck; nf -= nk             # leave-k-out priors
    q_g = (Cg + gamma * bg) / (ng + gamma)[:, None, None]
    if prior == "family":
        q_f = (Cf + beta * q_g) / (nf + beta)[:, None, None]
        q_prior = q_f
    elif prior == "group":
        q_prior = q_g
    elif prior == "none":
        q_prior = np.broadcast_to(bg, Ck.shape)
    else:
        raise ValueError(prior)
    q = (Ck + alpha * q_prior) / (nk + alpha)[:, None, None]
    level = np.where(nk > 0, "kinase",
                     np.where((nf > 0) & (prior == "family"), "family",
                              np.where((ng > 0) & (prior != "none"), "group", "background")))
    return q, level


def _prepare(train_pairs, sites, kinases, phospho):
    spos = pd.Index(sites.site_id)
    kpos = pd.Index(kinases.uniprot_acc)
    tp = train_pairs[["kinase_acc", "site_id"]].drop_duplicates()
    si = spos.get_indexer(tp.site_id)
    ki = kpos.get_indexer(tp.kinase_acc)
    ok = (si >= 0) & (ki >= 0)
    tp = tp[ok].assign(si=si[ok], ki=ki[ok])
    res = sites.res.to_numpy()
    tp["rt"] = np.where(res[tp.si.to_numpy()] == "Y", "Y", "ST")
    cls = kinases.residue_class.to_numpy()
    in_cls = np.array([rt in CLASS_RTS.get(c, ()) for rt, c in zip(tp.rt, cls[tp.ki.to_numpy()])], dtype=bool)
    info = {"pairs_in": len(train_pairs), "pairs_unknown_site_or_kinase": int((~ok).sum()),
            "pairs_out_of_class": int((~in_cls).sum())}
    tp = tp[in_cls]
    codes = _site_codes(sites, phospho)
    A = A_PH if phospho else A_UNPH
    ms = sites.in_U_ms.to_numpy(bool)
    bg = {}
    for rt in RTS:
        m = ms & np.isin(res, RT_RES[rt])
        C = _counts(codes[m], A) + 0.5               # tiny pseudocount: no zero background cells
        bg[rt] = C / C.sum(1, keepdims=True)
    return tp, codes, A, bg, info


def _kinase_counts(tp, codes, K, A, rt):
    Ck = np.zeros((K, 15, A))
    sub = tp[tp.rt == rt]
    ki = sub.ki.to_numpy()
    c = codes[sub.si.to_numpy()]
    for p in range(15):
        np.add.at(Ck, (ki, p, c[:, p]), 1.0)
    nk = np.bincount(ki, minlength=K).astype(float)
    return Ck, nk


def _tune(tp, codes, kinases, A, bg, prior, groups, n_folds=5, seed=0):
    """Grid search alpha x beta maximising macro (per-kinase) held-out log-likelihood."""
    rng = np.random.default_rng(seed)
    ug = np.unique(groups)
    fold_of = dict(zip(ug, rng.integers(0, n_folds, len(ug))))
    f = np.array([fold_of[g] for g in groups])
    K = len(kinases)
    betas = BETA_GRID if prior == "family" else (np.nan,)
    rows = []
    for rt in RTS:
        for i in range(n_folds):
            tr, te = tp[f != i], tp[f == i]
            te = te[te.rt == rt]
            if len(te) == 0:
                continue
            Ck, nk = _kinase_counts(tr, codes, K, A, rt)
            tki, tc = te.ki.to_numpy(), codes[te.si.to_numpy()]
            for a in ALPHA_GRID:
                for b in betas:
                    q, _ = _hier(Ck, nk, kinases, bg[rt], a, 1.0 if np.isnan(b) else b, GAMMA, prior)
                    ll = np.zeros(len(te))
                    for p in range(15):
                        ll += np.log2(q[tki, p, tc[:, p]])
                    per_k = pd.Series(ll).groupby(tki).mean()
                    rows.append(dict(rt=rt, fold=i, alpha=a, beta=b, macro_ll=per_k.mean(), n_k=len(per_k),
                                     micro_ll=ll.mean()))
    cv = pd.DataFrame(rows)
    if cv.empty:
        return 8.0, 4.0, cv
    agg = cv.groupby(["alpha", "beta"], dropna=False).apply(
        lambda d: np.average(d.macro_ll, weights=d.n_k), include_groups=False).rename("score").reset_index()
    best = agg.loc[agg.score.idxmax()]
    return float(best.alpha), (4.0 if np.isnan(best.beta) else float(best.beta)), cv


def fit_pssm(train_pairs: pd.DataFrame, sites: pd.DataFrame | None = None, kinases: pd.DataFrame | None = None,
             prior: str = "family", phospho: bool = False, alpha: float | None = None, beta: float | None = None,
             groups: pd.Series | None = None, seed: int = 0) -> PSSMModel:
    """Fit learned PSSMs for ALL kinases in `kinases` (index order) from train_pairs[kinase_acc, site_id].
    alpha/beta None -> tuned by internal 5-fold CV grouped by `groups` (Series site_id -> group id;
    default: substrate protein accession). Pure: no global label access."""
    sites = load_sites() if sites is None else sites
    kinases = load_kinases() if kinases is None else kinases
    tp, codes, A, bg, info = _prepare(train_pairs, sites, kinases, phospho)
    K = len(kinases)
    cv = None
    if alpha is None or beta is None:
        if groups is None:
            g = tp.site_id.str.rsplit("_", n=1).str[0].to_numpy()
        else:
            g = tp.site_id.map(groups).fillna(tp.site_id.str.rsplit("_", n=1).str[0]).astype(str).to_numpy()
        a_t, b_t, cv = _tune(tp, codes, kinases, A, bg, prior, g, seed=seed)
        alpha = a_t if alpha is None else alpha
        beta = b_t if beta is None else beta
    lo, n_train, level = {}, {}, {}
    cls = kinases.residue_class.to_numpy()
    for rt in RTS:
        Ck, nk = _kinase_counts(tp, codes, K, A, rt)
        q, lev = _hier(Ck, nk, kinases, bg[rt], alpha, beta, GAMMA, prior)
        L = np.log2(q / bg[rt][None]).astype(np.float32)
        has = np.array([rt in CLASS_RTS.get(c, ()) for c in cls])
        L[~has] = np.nan
        lo[rt] = L
        n_train[rt] = nk.astype(int)
        level[rt] = np.where(has, lev, "")
    info.update(n_pairs_used=len(tp), alpha=alpha, beta=beta, gamma=GAMMA)
    return PSSMModel(kinases.uniprot_acc.to_numpy(), phospho, prior, alpha, beta, GAMMA, lo, n_train, level, bg, cv, info)


def score(model: PSSMModel, site_idx=None, sites: pd.DataFrame | None = None, chunk: int = 50_000) -> np.ndarray:
    """K x len(site_idx) float32 log-odds scores (NaN where the site residue type is outside the kinase's class)."""
    sites = load_sites() if sites is None else sites
    site_idx = np.arange(len(sites)) if site_idx is None else np.asarray(site_idx)
    codes = _site_codes(sites, model.phospho)
    res = sites.res.to_numpy()[site_idx]
    K = len(model.kinase_acc)
    out = np.full((K, len(site_idx)), np.nan, np.float32)
    for rt in RTS:
        cols = np.where(np.isin(res, RT_RES[rt]))[0]
        L = model.logodds[rt]
        rows = np.where(np.isfinite(L[:, 0, 0]))[0]
        if len(rows) == 0 or len(cols) == 0:
            continue
        Lr = L[rows]
        for a in range(0, len(cols), chunk):
            cc = cols[a:a + chunk]
            c = codes[site_idx[cc]]
            acc = np.zeros((len(rows), len(cc)), np.float32)
            for p in range(15):
                acc += Lr[:, p, :][:, c[:, p]]
            out[np.ix_(rows, cc)] = acc
    return out


def score_pairs(model: PSSMModel, kinase_idx, site_idx, sites: pd.DataFrame | None = None) -> np.ndarray:
    """Scores for explicit (kinase_idx, site_idx) pairs (vectorised gather)."""
    sites = load_sites() if sites is None else sites
    kinase_idx, site_idx = np.asarray(kinase_idx), np.asarray(site_idx)
    codes = _site_codes(sites, model.phospho)[site_idx]
    isY = sites.res.to_numpy()[site_idx] == "Y"
    out = np.full(len(site_idx), np.nan, np.float32)
    for rt, m in (("ST", ~isY), ("Y", isY)):
        if m.any():
            L = model.logodds[rt]
            k, c = kinase_idx[m], codes[m]
            out[m] = sum(L[k, p, c[:, p]] for p in range(15))
    return out


def universe_percentile(model: PSSMModel, kinase_idx, values, site_idx, sites: pd.DataFrame | None = None) -> np.ndarray:
    """Percentile (0-100) of `values` within the kinase's scores over U_ms sites of the same residue type."""
    sites = load_sites() if sites is None else sites
    kinase_idx, values, site_idx = np.asarray(kinase_idx), np.asarray(values, float), np.asarray(site_idx)
    res = sites.res.to_numpy()
    ms = sites.in_U_ms.to_numpy(bool)
    rtq = np.where(res[site_idx] == "Y", "Y", "ST")
    out = np.full(len(values), np.nan)
    for rt in RTS:
        ref_idx = np.where(ms & np.isin(res, RT_RES[rt]))[0]
        m = rtq == rt
        ks = np.unique(kinase_idx[m])
        if len(ks) == 0:
            continue
        sub = PSSMModel(model.kinase_acc[ks], model.phospho, model.prior, model.alpha, model.beta, model.gamma,
                        {r: model.logodds[r][ks] for r in RTS}, {}, {}, model.bg)
        S = score(sub, ref_idx, sites)
        for j, k in enumerate(ks):
            r = np.sort(S[j][np.isfinite(S[j])])
            if len(r) == 0:
                continue
            mm = m & (kinase_idx == k)
            out[mm] = np.searchsorted(r, values[mm], side="right") / len(r) * 100
    return out


def top_positions(model: PSSMModel, kinase_acc: str, rt: str = "ST", n: int = 8) -> pd.DataFrame:
    """Top-weighted (position, residue, log2-odds) cells of a kinase's learned PSSM (positions -7..+7)."""
    k = int(np.where(model.kinase_acc == kinase_acc)[0][0])
    L = model.logodds[rt][k]
    A = L.shape[1]
    cells = [(p - 7, ALPHA[a], float(L[p, a])) for p in range(15) for a in range(A) if p != 7 and ALPHA[a] not in "_X"]
    d = pd.DataFrame(cells, columns=["pos", "aa", "log2odds"]).sort_values("log2odds", ascending=False)
    return d.head(n).reset_index(drop=True)


# ----------------------------------------------------------------------------------------------
# ksp.crossfit contract
# ----------------------------------------------------------------------------------------------
try:
    from ksp.crossfit import LabelDerivedFeature, register
except Exception:  # pragma: no cover - crossfit not available
    LabelDerivedFeature, register = object, (lambda c: c)


@register
class LearnedPSSM(LabelDerivedFeature):
    """F2a learned PSSM as a cross-fitted feature family.
    Columns: f2a_pssm (log-odds), f2a_pssm_pct (percentile within U_ms, same residue type),
    f2a_pssm_ph (+ phospho channel), f2a_pssm_ph_pct, f2a_pssm_level (0 kinase-specific, 1 family, 2 group,
    3 background prior only), f2a_pssm_log_ntrain (log1p of the kinase's training-label count)."""
    name = "f2a_pssm"
    prior = "family"

    def fit(self, train_pairs: pd.DataFrame, ctx=None) -> "LearnedPSSM":
        groups = None
        if ctx is not None and getattr(ctx, "clusters", None) is not None and len(ctx.clusters):
            acc = train_pairs.site_id.str.rsplit("_", n=1).str[0]
            groups = pd.Series(acc.map(ctx.clusters).fillna(acc).astype(str).to_numpy(), index=train_pairs.site_id)
            groups = groups[~groups.index.duplicated()]
        seed = getattr(ctx, "seed", 0) if ctx is not None else 0
        self.sites = load_sites()
        self.kinases = load_kinases()
        self.m0 = fit_pssm(train_pairs, self.sites, self.kinases, self.prior, False, groups=groups, seed=seed)
        self.m1 = fit_pssm(train_pairs, self.sites, self.kinases, self.prior, True, groups=groups, seed=seed)
        return self

    def transform(self, query: pd.DataFrame, ctx=None) -> pd.DataFrame:
        if "kinase_idx" in query and "site_idx" in query:
            ki, si = query.kinase_idx.to_numpy(int), query.site_idx.to_numpy(int)
        else:
            ki = pd.Index(self.kinases.uniprot_acc).get_indexer(query.kinase_acc)
            si = pd.Index(self.sites.site_id).get_indexer(query.site_id)
            if (ki < 0).any() or (si < 0).any():
                raise KeyError("query kinase/site not in canonical index")
        out = {}
        for tag, m in (("f2a_pssm", self.m0), ("f2a_pssm_ph", self.m1)):
            v = score_pairs(m, ki, si, self.sites)
            out[tag] = v
            out[tag + "_pct"] = universe_percentile(m, ki, v, si, self.sites)
        isY = self.sites.res.to_numpy()[si] == "Y"
        lev_code = {"kinase": 0, "family": 1, "group": 2, "background": 3, "": np.nan}
        lev = np.where(isY, self.m0.level["Y"][ki], self.m0.level["ST"][ki])
        out["f2a_pssm_level"] = np.array([lev_code[x] for x in lev], float)
        ntr = np.where(isY, self.m0.n_train["Y"][ki], self.m0.n_train["ST"][ki])
        out["f2a_pssm_log_ntrain"] = np.log1p(ntr).astype(float)
        return pd.DataFrame(out)
