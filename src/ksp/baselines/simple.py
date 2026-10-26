"""Simple baselines B0 / B1 / B2 / B7 on the 5-fold DEV CV (DESIGN §6.5, §8.1, §16 A19). Owner: baselines-simple.

All predictions are out-of-fold (OOF) for EVERY evaluation candidate ksp.eval needs (`candidate_pairs`):
  * for each kinase with >= 1 evaluable positive in the pooled DEV universe: all universe (U_ms, DEV) sites of its
    residue class  (= positives + E-all negatives; E-hard is a subset);
  * for M2: every real kinase x every labelled positive site of a matching residue class.
Rows a method cannot score are omitted (ksp.eval counts them; default policy 'worst').

Methods
  b0_prior   : score = fold-training L2 label count of the kinase (constant per kinase -> only M2 is informative).
  b0_pop     : score = site popularity  LT_LIT + MS_LIT + MS_CST + compendium condition count (cache/sites/universe)
               (constant per site -> M2 is ties). Diagnostic floor only; never a feature.
  b0_miss    : LightGBM on missingness indicators only (pair-, kinase- and site-level NaN patterns of the label-free
               phospho / context stores + site data-availability flags); trained per CV fold like any model.
  b1_*       : Cantley (F1) arms read from cache/features/motif via ksp.store: (i) as-published percentile,
               (ii) universe-recalibrated percentile, (iii) priming-aware hard / soft (universe percentile), and the
               nearest-atlas imputed variant (own matrix where available, imputed rows elsewhere; separate arm).
  b2_pssm    : F2a learned PSSM (ksp.features.motif_learned.fit_pssm) fit on fold-training L2 labels, scored as the
               within-U_ms percentile (monotone in the log-odds per kinase, comparable across kinases for M2).
  b7_lr      : logistic regression on [logit Cantley best arm, logit B2 percentile (cross-fitted through
               ksp.crossfit), PPI direct edge any screen, localisation combined cosine] + missing indicators, with
               per-kinase intercepts (L2-shrunk one-hot = partial pooling) and inverse-sampling-probability weights
               on the sampled negatives, so intercepts estimate the kinase's prevalence in the universe (needed for M2).

Training rows (b0_miss, b7): per kinase, fold-training L2 positives (U_ms, in class) + sampled negatives at 1:50 with
50% hard (A15: if the hard pool is short, fill with E-all-type), drawn from training-protein U_ms sites of the
kinase's class with the DESIGN §5 masks built from training-protein labels only (ksp.labels.NegativeMasks).
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from .. import canary as _canary
from .. import crossfit, index, labels, splits, store

N_REAL = None   # set lazily: number of real kinases in the index


# ================================================================================================ index helpers
def n_real() -> int:
    global N_REAL
    if N_REAL is None:
        N_REAL = len(index.load_kinases())
    return N_REAL


def site_protein_idx(sites: pd.DataFrame) -> np.ndarray:
    P = index.load_proteins()
    pidx = pd.Index(P.acc).get_indexer(sites.acc.values)
    if (pidx < 0).any():
        raise KeyError(f"{(pidx < 0).sum()} sites on proteins absent from the protein index")
    return pidx


def site_fold(sites: pd.DataFrame) -> np.ndarray:
    """CV fold per index site (-1 = not a DEV protein)."""
    f = splits.load_cv_folds()
    return sites.acc.map(f).fillna(-1).astype(np.int64).values


def candidate_pairs(data, kinases: list | None = None, include_m2: bool = True) -> pd.DataFrame:
    """All (kinase_idx, site_idx) pairs ksp.eval.evaluate needs for `kinases` (default: every kinase with >=1
    positive in data, real only) - E-all U positives over the universe - plus, for M2, every real kinase x every
    labelled positive site of a matching class. Sorted by (kinase_idx, site_idx), unique."""
    from ..eval import _res_ok
    K = data.kinases
    kpos = {a: i for i, a in enumerate(K.uniprot_acc.values)}
    res = data.sites.res.values
    todo = kinases if kinases is not None else [k for k in sorted(data.positives) if not k.startswith("CANARY_")]
    parts = []
    for k in todo:
        i = kpos[k]
        s = np.where(data.universe & _res_ok(res, K.residue_class.values[i]))[0]
        parts.append(np.c_[np.full(len(s), i), s])
    if include_m2:
        pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
        us = np.unique(pm.site_idx.values)
        for i in range(n_real()):
            ok = _res_ok(res[us], K.residue_class.values[i])
            parts.append(np.c_[np.full(ok.sum(), i), us[ok]])
    a = np.concatenate(parts).astype(np.int64)
    key = np.unique(a[:, 0] * np.int64(len(data.sites)) + a[:, 1])
    return pd.DataFrame({"kinase_idx": (key // len(data.sites)).astype(np.int64),
                         "site_idx": (key % len(data.sites)).astype(np.int64)})


# ================================================================================================ store gathers
def gather_ks(name: str, kidx: np.ndarray, sidx: np.ndarray) -> np.ndarray:
    """Values of a kinase x site store at (kidx, sidx); NaN for canary kinases (kidx >= n_real). kidx sorted."""
    arr, desc = store.load_matrix(name)
    if desc["row_axis"] != "kinase" or desc["col_axis"] != "site":
        raise ValueError(f"{name} is not kinase x site")
    out = np.full(len(kidx), np.nan, np.float32)
    b = np.searchsorted(kidx, np.arange(n_real() + 1))
    for k in range(n_real()):
        if b[k + 1] > b[k]:
            row = np.asarray(arr[k], dtype=np.float32)
            out[b[k]:b[k + 1]] = row[sidx[b[k]:b[k + 1]]]
    return out


def gather_kp(name: str, kidx: np.ndarray, pidx: np.ndarray, fill=np.nan) -> np.ndarray:
    """Values of a kinase x protein store at (kidx, pidx); `fill` for canary kinases."""
    arr, desc = store.load_matrix(name)
    if desc["row_axis"] != "kinase" or desc["col_axis"] != "protein":
        raise ValueError(f"{name} is not kinase x protein")
    A = np.asarray(arr, dtype=np.float32)
    out = np.full(len(kidx), fill, np.float32)
    m = kidx < n_real()
    out[m] = A[kidx[m], pidx[m]]
    return out


# ================================================================================================ fold labels / rows
def fold_train_labels(fold: int, with_canaries: bool = False) -> pd.DataFrame:
    tr, _, _ = crossfit.outer_split(fold, "L2", with_canaries=with_canaries)
    return tr


def fold_train_counts(fold: int, kinases: pd.DataFrame, with_canaries: bool = False) -> np.ndarray:
    tr = fold_train_labels(fold, with_canaries)
    n = tr.groupby("kinase_acc").site_id.nunique()
    return kinases.uniprot_acc.map(n).fillna(0).values.astype(float)


def sample_training_rows(fold: int, data, ratio: int = 50, hard_frac: float = 0.5, seed: int = 0,
                         with_canaries: bool = False, log: dict | None = None) -> pd.DataFrame:
    """Training rows (kinase_idx, site_idx, y, w_ipw) for outer CV fold `fold` (positives = fold-training L2 labels on
    U_ms sites in the kinase's class; negatives sampled as described in the module docstring)."""
    from ..eval import _res_ok
    S, K = data.sites, data.kinases
    spos = pd.Series(np.arange(len(S)), index=S.site_id.values)
    sf = site_fold(S)
    tr_site = (sf >= 0) & (sf != fold) & S.in_U_ms.values
    tr_lab = fold_train_labels(fold, with_canaries)
    all_dev = labels.load_dev("all")
    if with_canaries:
        all_dev = pd.concat([all_dev, _canary.canary_labels()], ignore_index=True)
    tr_prot = set(S.acc.values[tr_site])
    mask_pairs = all_dev[all_dev.acc.isin(tr_prot)]
    masks = labels.NegativeMasks(mask_pairs, paralogs=data.par)
    rng = np.random.default_rng(seed + 1000 * fold)
    res = S.res.values
    kpos = {a: i for i, a in enumerate(K.uniprot_acc.values)}
    rows = []
    n_short = 0
    for k, g in tr_lab.groupby("kinase_acc"):
        if k not in kpos:
            raise KeyError(f"training label kinase {k} not in kinase index")
        i = kpos[k]
        cls_ok = _res_ok(res, K.residue_class.values[i])
        p = spos.reindex(g.site_id.values).dropna().values.astype(np.int64)
        p = np.unique(p[tr_site[p] & cls_ok[p]])
        if len(p) == 0:
            continue
        elig = tr_site & cls_ok
        elig[masks.excluded(k)] = False
        elig[p] = False
        hard = np.zeros(len(S), bool)
        hard[masks.hard_candidates(k)] = True
        hard &= elig
        n_neg = ratio * len(p)
        h_idx, a_idx = np.where(hard)[0], np.where(elig & ~hard)[0]
        nh = min(int(round(hard_frac * n_neg)), len(h_idx))
        if nh < int(round(hard_frac * n_neg)):
            n_short += 1
        na = min(n_neg - nh, len(a_idx))
        hs = rng.choice(h_idx, nh, replace=False) if nh else np.zeros(0, np.int64)
        as_ = rng.choice(a_idx, na, replace=False) if na else np.zeros(0, np.int64)
        # inverse sampling probability per stratum (hard / other) -> weights recover the universe composition
        wh = len(h_idx) / max(nh, 1)
        wa = len(a_idx) / max(na, 1)
        rows.append(pd.DataFrame({"kinase_idx": i, "site_idx": np.r_[p, hs, as_],
                                  "y": np.r_[np.ones(len(p)), np.zeros(nh + na)],
                                  "w_ipw": np.r_[np.ones(len(p)), np.full(nh, wh), np.full(na, wa)]}))
    if log is not None:
        log["kinases_hard_pool_short"] = log.get("kinases_hard_pool_short", 0) + n_short
    out = pd.concat(rows, ignore_index=True)
    out["kinase_acc"] = K.uniprot_acc.values[out.kinase_idx.values]
    out["site_id"] = S.site_id.values[out.site_idx.values]
    return out


# ================================================================================================ B0
def b0_prior(cand: pd.DataFrame, data) -> dict:
    """score = the kinase's DEV L2 count (constant per kinase -> within-kinase metrics exactly chance);
    score_site (used by ksp.eval for M2) = the fold-training L2 count (OOF). NOTE: using the fold-specific count as
    the within-kinase score gives M1 ~0.40, not 0.5 - in pooled OOF a kinase's positives sit in folds where its
    training count is lower (pooled-OOF offset artefact; see summary notes)."""
    sf = site_fold(data.sites)[cand.site_idx.values]
    site = np.full(len(cand), np.nan)
    for f in range(splits.N_FOLDS):
        n = fold_train_counts(f, data.kinases)
        m = sf == f
        site[m] = n[cand.kinase_idx.values[m]]
    const = data.kinases.uniprot_acc.map(data.train_counts).fillna(0).values.astype(float)
    return {"score": const[cand.kinase_idx.values], "score_site": site}


def site_popularity(sites: pd.DataFrame) -> np.ndarray:
    u = pd.read_parquet(index.ROOT / "cache" / "sites" / "universe.parquet",
                        columns=["site_id", "psp_lt_lit", "psp_ms_lit", "psp_ms_cst", "comp_n_conditions"])
    u = u.set_index("site_id").reindex(sites.site_id.values)
    v = sum(pd.to_numeric(u[c], errors="coerce").fillna(0).values.astype(float)
            for c in ["psp_lt_lit", "psp_ms_lit", "psp_ms_cst", "comp_n_conditions"])
    return v


def b0_pop(cand: pd.DataFrame, data) -> np.ndarray:
    return site_popularity(data.sites)[cand.site_idx.values]


# ---------------------------------------------------------------- B0-miss feature set (missingness only)
MISS_KS = ["phospho/f6b_ridge_z", "phospho/f6b_contrast", "phospho/f6c_prot", "phospho/f6c_aloop"]
MISS_KP = ["context/ppi/ppi_rwr", "context/loc/loc_cos_hpa", "context/loc/loc_cos_opencell",
           "context/loc/loc_cos_hcm", "context/loc/loc_cos_hein", "context/loc/loc_cos_comp_exp",
           "context/loc/loc_hein_profile_corr", "context/coexpr/coexpr_hpa_tissue_spearman",
           "context/coexpr/coexpr_depmap_pearson"]


class MissFeatures:
    """Missingness indicators for arbitrary (kinase_idx, site_idx) rows. Kinase- and site-level indicators are
    'row has any non-NaN value' of the same stores; site-level data availability: in compendium / CPTAC,
    structure NaN, topology NaN (PSP/Ochoa presence deliberately excluded). Canary kinases: all kinase-dependent indicators = missing."""

    def __init__(self, sites: pd.DataFrame):
        self.sites = sites
        self.pidx = site_protein_idx(sites)
        self.site_level = {}
        self.kin_level = {}
        for n in MISS_KS:
            arr, _ = store.load_matrix(n)
            anyk = np.zeros(arr.shape[1], bool)
            kany = np.zeros(arr.shape[0], bool)
            for k in range(arr.shape[0]):
                ok = ~np.isnan(np.asarray(arr[k], np.float32))
                anyk |= ok
                kany[k] = ok.any()
            self.site_level[f"site_any_{n.split('/')[-1]}"] = anyk
            self.kin_level[f"kin_any_{n.split('/')[-1]}"] = kany
        for n in MISS_KP:
            arr, _ = store.load_matrix(n)
            A = ~np.isnan(np.asarray(arr, np.float32))
            self.kin_level[f"kin_any_{n.split('/')[-1]}"] = A.any(1)
            self.site_level[f"prot_any_{n.split('/')[-1]}"] = A.any(0)[self.pidx]
        kis, _ = store.load_matrix("context/ppi/ppi_kinase_in_source")
        pis, _ = store.load_matrix("context/ppi/ppi_protein_in_source")
        cols = ["bioplex_293t", "bioplex_hct116", "huri", "opencell", "buljan", "hcm"]
        for j, c in enumerate(cols):
            self.kin_level[f"kin_in_ppi_{c}"] = np.asarray(kis[:, j]).astype(bool)
            self.site_level[f"prot_in_ppi_{c}"] = np.asarray(pis[:, j]).astype(bool)[self.pidx]
        st, _ = store.load_matrix("structure/site_structure")
        self.site_level["site_structure_nan"] = np.isnan(np.asarray(st, np.float32)).any(1)
        lum, _ = store.load_matrix("context/loc/site_lumenal")
        self.site_level["site_lumenal_nan"] = np.isnan(np.asarray(lum, np.float32))
        u = pd.read_parquet(index.ROOT / "cache" / "sites" / "universe.parquet",
                            columns=["site_id", "src_compendium", "src_cptac"])
        u = u.set_index("site_id").reindex(sites.site_id.values)
        # only the data sources that phospho features are computed from (F6a/F6b compendium, F6a/F6c CPTAC);
        # PSP / Ochoa presence is popularity (B0-pop territory), not feature missingness -> excluded
        for c in ["src_compendium", "src_cptac"]:
            self.site_level[f"site_{c}"] = u[c].fillna(False).astype(bool).values

    def transform(self, kidx: np.ndarray, sidx: np.ndarray) -> pd.DataFrame:
        o = np.argsort(kidx, kind="stable")
        ks, ss = kidx[o], sidx[o]
        X = {}
        for n in MISS_KS:
            v = gather_ks(n, ks, ss)
            X[f"pair_has_{n.split('/')[-1]}"] = ~np.isnan(v)
        ps = self.pidx[ss]
        for n in MISS_KP:
            v = gather_kp(n, ks, ps)
            X[f"pair_has_{n.split('/')[-1]}"] = ~np.isnan(v)
        real = ks < n_real()
        for c, v in self.kin_level.items():
            x = np.zeros(len(ks), bool)
            x[real] = v[ks[real]]
            X[c] = x
        for c, v in self.site_level.items():
            X[c] = v[ss]
        D = pd.DataFrame({c: v.astype(np.uint8) for c, v in X.items()})
        inv = np.empty_like(o)
        inv[o] = np.arange(len(o))
        return D.iloc[inv].reset_index(drop=True)


def b0_miss(cand: pd.DataFrame, data, seed: int = 0, threads: int | None = None, log: dict | None = None,
            params: dict | None = None) -> np.ndarray:
    import lightgbm as lgb
    mf = MissFeatures(data.sites)
    sf = site_fold(data.sites)[cand.site_idx.values]
    out = np.full(len(cand), np.nan)
    prm = {"objective": "binary", "learning_rate": 0.05, "num_leaves": 31, "min_data_in_leaf": 100,
           "feature_fraction": 0.9, "bagging_fraction": 0.8, "bagging_freq": 1, "verbose": -1,
           "num_threads": threads or int(os.environ.get("OMP_NUM_THREADS", 8)), "seed": seed}
    prm.update(params or {})
    for f in range(splits.N_FOLDS):
        t0 = time.time()
        tr = sample_training_rows(f, data, seed=seed, log=log)
        Xtr = mf.transform(tr.kinase_idx.values, tr.site_idx.values)
        # early stopping on an inner validation split grouped by substrate cluster (never the outer fold)
        va = splits.train_valid_split_inner(tr.site_id.values, 0.2, seed)
        dtr = lgb.Dataset(Xtr[~va], tr.y.values[~va])
        dva = lgb.Dataset(Xtr[va], tr.y.values[va], reference=dtr)
        bst = lgb.train(prm, dtr, 1000, valid_sets=[dva], callbacks=[lgb.early_stopping(50, verbose=False)])
        m = sf == f
        Xev = mf.transform(cand.kinase_idx.values[m], cand.site_idx.values[m])
        out[m] = bst.predict(Xev, num_iteration=bst.best_iteration)
        if log is not None:
            log.setdefault("folds", {})[f] = {"n_train": int(len(tr)), "n_pos": int(tr.y.sum()),
                                              "best_iter": int(bst.best_iteration), "s": round(time.time() - t0, 1),
                                              "features": list(Xtr.columns),
                                              "gain": dict(zip(Xtr.columns,
                                                               bst.feature_importance("gain").round(1).tolist()))}
    return out


# ================================================================================================ B1 Cantley
B1_ARMS = {
    "b1_published": ("motif/cantley_unprimed_pct_published", None),
    "b1_universe": ("motif/cantley_unprimed_pct_universe", None),
    "b1_primed_hard": ("motif/cantley_primed_hard_pct_universe", None),
    "b1_primed_soft": ("motif/cantley_primed_soft_pct_universe", None),
    "b1_primed_hard_published": ("motif/cantley_primed_hard_pct_published", None),
    "b1_primed_soft_published": ("motif/cantley_primed_soft_pct_published", None),
}


def cantley(cand: pd.DataFrame, arm_store: str, imputed: bool = False) -> np.ndarray:
    v = gather_ks(arm_store, cand.kinase_idx.values, cand.site_idx.values)
    if imputed:
        imp = arm_store.replace("motif/cantley_", "motif/cantley_imputed_")
        vi = gather_ks(imp, cand.kinase_idx.values, cand.site_idx.values)
        both = ~np.isnan(v) & ~np.isnan(vi)
        if both.any():
            raise RuntimeError(f"{both.sum()} cells have both an own and an imputed Cantley value")
        v = np.where(np.isnan(v), vi, v)
    return v


# ================================================================================================ B2 learned PSSM
def pssm_kinase_table(with_canaries: bool) -> pd.DataFrame:
    """Kinase table for PSSM fits. Canary rows get EMPTY Manning group/family (= no hierarchical pooling, prior =
    background): ksp.canary gives all canaries group 'CANARY', which would pool the 12 negative canaries' (v2: real
    labelled-site) windows into every canary's prior - a cross-canary 'labelled-site-likeness' channel that is not
    part of any real kinase's biology."""
    K = (_canary.kinases_with_canaries() if with_canaries else index.load_kinases()).copy()
    c = K.uniprot_acc.astype(str).str.startswith("CANARY_")
    K.loc[c, "manning_group"] = ""
    K.loc[c, "manning_family"] = ""
    return K


def _fit_pssm_canary_safe(ml, pairs, sites, kinases, phospho, groups, seed, tune_without_canaries, info=None):
    """fit_pssm; if canary kinases are present and tune_without_canaries, alpha/beta are tuned on the REAL kinases'
    pairs only (canaries must never steer shared hyperparameters), then all kinases are fitted with them fixed."""
    has_can = pairs.kinase_acc.astype(str).str.startswith("CANARY_").any()
    if has_can and tune_without_canaries:
        real = pairs[~pairs.kinase_acc.astype(str).str.startswith("CANARY_")]
        m0 = ml.fit_pssm(real[["kinase_acc", "site_id"]], sites, index.load_kinases(), "family", phospho,
                         groups=groups, seed=seed)
        if info is not None:
            info["tuned_without_canaries"] = {"alpha": m0.alpha, "beta": m0.beta}
        return ml.fit_pssm(pairs[["kinase_acc", "site_id"]], sites, kinases, "family", phospho, alpha=m0.alpha,
                           beta=m0.beta, groups=groups, seed=seed)
    return ml.fit_pssm(pairs[["kinase_acc", "site_id"]], sites, kinases, "family", phospho, groups=groups, seed=seed)


def b2_pssm(cand: pd.DataFrame, data, phospho: bool = False, with_canaries: bool = False, seed: int = 0,
            log: dict | None = None, tune_without_canaries: bool = True) -> np.ndarray:
    """OOF learned-PSSM percentile: fit on fold-training L2 labels (all of them: B2 is a model, not a feature)."""
    from ..features import motif_learned as ml
    K = pssm_kinase_table(with_canaries)
    S = data.sites
    sf = site_fold(S)[cand.site_idx.values]
    out = np.full(len(cand), np.nan)
    cl = splits.load_clusters()
    for f in range(splits.N_FOLDS):
        t0 = time.time()
        tr = fold_train_labels(f, with_canaries)
        acc = tr.site_id.str.rsplit("_", n=1).str[0]
        groups = pd.Series(acc.map(cl).fillna(acc).astype(str).values, index=tr.site_id.values)
        groups = groups[~groups.index.duplicated()]
        info = {}
        m = _fit_pssm_canary_safe(ml, tr, S, K, phospho, groups, seed, tune_without_canaries, info)
        r = sf == f
        ki, si = cand.kinase_idx.values[r], cand.site_idx.values[r]
        if ki.max() >= len(K):
            raise KeyError("candidate kinase_idx beyond the kinase table used for the PSSM fit")
        v = ml.score_pairs(m, ki, si, S)
        out[r] = ml.universe_percentile(m, ki, v, si, S)
        if log is not None:
            log.setdefault("folds", {})[f] = {"alpha": m.alpha, "beta": m.beta, "n_pairs": m.info.get("n_pairs_used"),
                                              "s": round(time.time() - t0, 1), **info}
    return out


class LearnedPSSMAnyKinases(crossfit.LabelDerivedFeature):
    """F2a through the crossfit contract for an arbitrary kinase table (real kinases + canaries); same model as
    ksp.features.motif_learned.LearnedPSSM (family prior, no phospho channel), percentile output only."""
    name = "b7_f2a_pssm_pct"

    def __init__(self, kinases: pd.DataFrame, sites: pd.DataFrame, tune_without_canaries: bool = True):
        self.kinases, self.sites, self.twc = kinases, sites, tune_without_canaries

    def fit(self, train_pairs, ctx=None):
        from ..features import motif_learned as ml
        acc = train_pairs.site_id.str.rsplit("_", n=1).str[0]
        groups = pd.Series(acc.map(ctx.clusters).fillna(acc).astype(str).to_numpy(), index=train_pairs.site_id)
        groups = groups[~groups.index.duplicated()]
        self.m = _fit_pssm_canary_safe(ml, train_pairs, self.sites, self.kinases, False, groups, ctx.seed, self.twc)
        return self

    def transform(self, query, ctx=None):
        from ..features import motif_learned as ml
        ki, si = query.kinase_idx.to_numpy(np.int64), query.site_idx.to_numpy(np.int64)
        v = ml.score_pairs(self.m, ki, si, self.sites)
        return pd.DataFrame({"b2_pct": ml.universe_percentile(self.m, ki, v, si, self.sites)})


# ================================================================================================ B7
def _logit_pct(p: np.ndarray) -> np.ndarray:
    q = np.clip(np.asarray(p, float) / 100.0, 1e-4, 1 - 1e-4)
    return np.log(q / (1 - q))


def b7_features(kidx, sidx, pidx, cantley_store: str, b2_pct: np.ndarray, pos_partners: set | None,
                protein_acc: np.ndarray) -> pd.DataFrame:
    o = np.argsort(kidx, kind="stable")
    inv = np.empty_like(o)
    inv[o] = np.arange(len(o))
    c = gather_ks(cantley_store, kidx[o], sidx[o])[inv]
    ppi = gather_kp("context/ppi/ppi_edge_any", kidx, pidx, fill=0.0)
    if pos_partners is not None:            # planted canary PPI (ksp.canary): an ordinary edge for CANARY_POS
        pos_i = n_real() + _canary.CANARY_IDS.index(_canary.POS_ID)
        m = kidx == pos_i
        ppi[m] = np.isin(protein_acc[pidx[m]], list(pos_partners)).astype(np.float32)
    loc = gather_kp("context/loc/loc_cos_combined", kidx, pidx)
    X = pd.DataFrame({"cantley_logit": np.where(np.isnan(c), 0.0, _logit_pct(np.nan_to_num(c, nan=50.0))),
                      "cantley_missing": np.isnan(c).astype(float),
                      "b2_logit": np.where(np.isnan(b2_pct), 0.0, _logit_pct(np.nan_to_num(b2_pct, nan=50.0))),
                      "b2_missing": np.isnan(b2_pct).astype(float),
                      "ppi_edge_any": ppi.astype(float),
                      "loc_cos_combined": np.nan_to_num(loc, nan=0.0).astype(float),
                      "loc_missing": np.isnan(loc).astype(float)})
    return X


def b7_lr(cand: pd.DataFrame, data, cantley_store: str, with_canaries: bool = False, seed: int = 0, C: float = 1.0,
          log: dict | None = None, tune_without_canaries: bool = True) -> dict:
    """OOF B7 scores (logit scale). Per outer fold: training rows sampled from fold-training data; the B2 feature is
    cross-fitted with ksp.crossfit.crossfit (inner 5-fold over substrate clusters; eval rows by A1 inner-fold fits)."""
    import scipy.sparse as sp
    from sklearn.linear_model import LogisticRegression
    K = pssm_kinase_table(with_canaries)
    S = data.sites
    pidx_all = site_protein_idx(S)
    pacc = index.load_proteins().acc.values
    partners = set(_canary.canary_ppi_partners()) if with_canaries else None
    cl = splits.load_clusters()
    sf = site_fold(S)[cand.site_idx.values]
    out = np.full(len(cand), np.nan)
    out_site = np.full(len(cand), np.nan)
    nK = len(K)
    for f in range(splits.N_FOLDS):
        t0 = time.time()
        tr = sample_training_rows(f, data, seed=seed, with_canaries=with_canaries, log=log)
        lab = fold_train_labels(f, with_canaries)
        r = np.where(sf == f)[0]
        ki, si = cand.kinase_idx.values[r], cand.site_idx.values[r]
        ev = pd.DataFrame({"kinase_acc": data.kinases.uniprot_acc.values[ki], "site_id": S.site_id.values[si],
                           "kinase_idx": ki, "site_idx": si})
        ftr, fev = crossfit.crossfit(lambda: LearnedPSSMAnyKinases(K, S, tune_without_canaries), lab[["kinase_acc", "site_id"]],
                                     tr[["kinase_acc", "site_id", "kinase_idx", "site_idx"]], ev, seed=seed,
                                     outer_fold=f, clusters=cl)
        t1 = time.time()
        Xtr = b7_features(tr.kinase_idx.values, tr.site_idx.values, pidx_all[tr.site_idx.values], cantley_store,
                          ftr.b2_pct.values, partners, pacc)
        mu, sd = Xtr.mean(), Xtr.std().replace(0, 1.0)
        w = tr.w_ipw.values.astype(float)
        w = w / w.mean()
        Ztr = sp.hstack([sp.csr_matrix(((Xtr - mu) / sd).values),
                         sp.csr_matrix((np.ones(len(tr)), (np.arange(len(tr)), tr.kinase_idx.values)),
                                       shape=(len(tr), nK))]).tocsr()
        lr = LogisticRegression(C=C, max_iter=5000, tol=1e-6)
        lr.fit(Ztr, tr.y.values, sample_weight=w)
        Xev = b7_features(ki, si, pidx_all[si], cantley_store, fev.b2_pct.values, partners, pacc)
        nf = Xev.shape[1]
        coef = lr.coef_.ravel()
        lin = ((Xev - mu) / sd).values @ coef[:nf]
        out[r] = lin                                             # within-kinase score (M1/M3/M4)
        out_site[r] = lin + coef[nf:][ki] + lr.intercept_[0]     # with per-kinase intercept (M2, A20)
        if log is not None:
            log.setdefault("folds", {})[f] = {
                "n_train": int(len(tr)), "n_pos": int(tr.y.sum()), "crossfit_s": round(t1 - t0, 1),
                "s": round(time.time() - t0, 1), "intercept": float(lr.intercept_[0]),
                "coef_std": dict(zip(Xtr.columns, coef[:nf].round(4).tolist())),
                "kinase_intercept_sd": float(np.std(coef[nf:][np.unique(tr.kinase_idx.values)]))}
    return {"score": out, "score_site": out_site}
