"""SELPHI2.0 model on our splits: negatives, cross-fitted features, RF grid search, bagging, OOF prediction.

Their scheme (paper1 Methods; src/model_training/learn.py, do-learn*.sh):
  * positives = PSP kinase-site pairs; negatives = random kinase-site pairs (paper: 10x/20x/50x "no differences",
    50x presented; code: 1x when > 10,000 positives), re-drawn in each of 100 runs; final probability = mean of the
    100 runs' outputs.
  * RandomForestClassifier(bootstrap=True), GridSearchCV(cv=3, scoring='roc_auc') per run over
    max_depth {10,20,50,70,80,100} x min_samples_split {8,10,12} x n_estimators {150,500,1000,1500}.
Adaptations (documented in docs/baseline_selphi2.md):
  * outer 5-fold CV = our substrate-cluster folds (results/splits/cv_folds.tsv) instead of StratifiedKFold(10) over
    pairs; inner grid-search CV = GroupKFold(3) over substrate clusters;
  * negatives: uniform random (kinase, site) pairs over training-fold U_ms sites with Ochoa features, residue-class
    matched, minus KP(k) / Hamming / paralog masks (ksp.labels.NegativeMasks built from training-fold labels only);
  * label-derived columns cross-fitted by ksp.crossfit (inner folds for training rows);
  * bags / grid scaled to the time box (cfg).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

from ksp import crossfit, labels, splits
from ksp.index import load_kinases, load_proteins, load_sites

from . import features as F

GRID = {"max_depth": [10, 20, 50, 70, 80, 100], "min_samples_split": [8, 10, 12],
        "n_estimators": [150, 500, 1000, 1500]}


@dataclass
class Config:
    neg_ratio: int = 10                 # negatives per positive per bag
    n_bags: int = 10                    # their 100 runs, scaled
    grid: dict = field(default_factory=lambda: dict(GRID))
    grid_bags: tuple = (0,)             # bags in which the grid search is run (others reuse its winner)
    n_estimators_cap: int | None = None  # cap at prediction time (documented if used)
    n_jobs: int = 12
    seed: int = 0
    drop_cols: tuple = ()               # ablation arms (e.g. F.LEAKY_COLS)
    label_variant: str = "L2"


# ------------------------------------------------------------------------------------------------ universe / rows
def selphi_universe() -> np.ndarray:
    """Bool over index sites: U_ms sites with Ochoa site features (SELPHI's feature-available universe)."""
    S = load_sites()
    return S.site_id.isin(set(F.load_site_table().index)).to_numpy() & S.in_U_ms.to_numpy(bool)


def _res_ok(res, cls):
    if cls == "STY":
        return np.ones(len(res), bool)
    if cls == "ST":
        return (res == "S") | (res == "T")
    return res == "Y"


def candidate_pairs(site_mask: np.ndarray, kinases: pd.DataFrame | None = None) -> pd.DataFrame:
    """All (kinase, site) pairs with site in `site_mask` and residue in the kinase's class."""
    S = load_sites()
    K = load_kinases() if kinases is None else kinases
    res = S.res.to_numpy()
    allK = load_kinases().uniprot_acc
    kpos = pd.Series(np.arange(len(allK)), index=allK)
    parts = []
    for acc, cls in zip(K.uniprot_acc, K.residue_class):
        s = np.where(site_mask & _res_ok(res, cls))[0]
        parts.append(pd.DataFrame({"kinase_idx": np.full(len(s), kpos[acc], np.int32), "site_idx": s.astype(np.int32)}))
    d = pd.concat(parts, ignore_index=True)
    d["kinase_acc"] = allK.to_numpy()[d.kinase_idx.to_numpy()]
    d["site_id"] = S.site_id.to_numpy()[d.site_idx.to_numpy()]
    return d[["kinase_acc", "site_id", "kinase_idx", "site_idx"]]


def sample_negatives(pos: pd.DataFrame, site_mask: np.ndarray, masks: labels.NegativeMasks, kinases: list,
                     n_total: int, rng: np.random.Generator) -> pd.DataFrame:
    """Uniform random (kinase, site) pairs (SELPHI's random negatives) over valid pairs: site in site_mask and in the
    kinase's residue class, not masked for k (KP / Hamming / paralogs), not a positive."""
    S = load_sites()
    res = S.res.to_numpy()
    K = load_kinases().set_index("uniprot_acc")
    kpos = pd.Series(np.arange(len(K)), index=K.index)
    pos_by_k = pos.groupby("kinase_acc").site_idx.agg(lambda x: set(x))
    valid = {}
    for k in kinases:
        ok = site_mask & _res_ok(res, K.residue_class[k])
        ok[masks.excluded(k)] = False
        v = np.where(ok)[0]
        if k in pos_by_k.index:
            v = v[~np.isin(v, list(pos_by_k[k]))]
        valid[k] = v
    sizes = np.array([len(valid[k]) for k in kinases], float)
    alloc = rng.multinomial(n_total, sizes / sizes.sum())
    parts = []
    for k, n in zip(kinases, alloc):
        n = min(int(n), len(valid[k]))
        if n == 0:
            continue
        s = rng.choice(valid[k], n, replace=False)
        parts.append(pd.DataFrame({"kinase_acc": k, "kinase_idx": kpos[k], "site_idx": s}))
    d = pd.concat(parts, ignore_index=True)
    d["site_id"] = S.site_id.to_numpy()[d.site_idx.to_numpy()]
    return d[["kinase_acc", "site_id", "kinase_idx", "site_idx"]]


# ------------------------------------------------------------------------------------------------ feature assembly
class Assembler:
    """Builds the 44-column SELPHI matrix for (kinase_idx, site_idx) rows + label-derived columns."""

    def __init__(self):
        S = load_sites()
        st = F.load_site_table()
        self.site_feat = st.reindex(S.site_id).to_numpy(np.float32)       # NaN rows = no Ochoa features
        prot = load_proteins().acc
        self.site_prot = pd.Series(np.arange(len(prot)), index=prot).reindex(S.acc).to_numpy()
        if np.isnan(self.site_prot).any():
            raise KeyError("sites without protein index")
        self.site_prot = self.site_prot.astype(np.int64)
        self.isY = S.res.eq("Y").to_numpy()
        self.expr = F.load_expression()
        self.cst, self.cy = F.load_cantley_norm()

    def matrix(self, ki: np.ndarray, si: np.ndarray, lab: pd.DataFrame, cols=F.SELPHI_COLS) -> np.ndarray:
        ki, si = np.asarray(ki, np.int64), np.asarray(si, np.int64)
        n = len(ki)
        out = {}
        sf = self.site_feat[si]
        if np.isnan(sf).any():
            raise ValueError("rows outside SELPHI's universe (no Ochoa site features)")
        for j, c in enumerate(F.SITE_COLS):
            out[c] = sf[:, j]
        pi = self.site_prot[si]
        for c in ("gtex", "rna.tiss", "rna.cell"):
            out[c] = self.expr[c][ki, pi].astype(np.float32)
        out["sel.1"] = self.expr["sel.1"][ki]
        out["sel.2"] = self.expr["sel.2"][pi]
        st = self.cst[ki, si]
        y = self.cy[ki, si]
        out["pssm_ser_thr"] = np.where(self.isY[si], 0, st).astype(np.float32)
        out["pssm_tyr"] = np.where(self.isY[si], y, 0).astype(np.float32)
        for c in F.LABEL_COLS:
            out[c] = lab[c].to_numpy(np.float32)
        X = np.column_stack([out[c] for c in cols]).astype(np.float32)
        X[~np.isfinite(X)] = 0.0                       # SELPHI: every remaining NA -> 0 (merge_kin_sites.r)
        assert X.shape == (n, len(cols))
        return X


def label_features(train_labels, train_rows, eval_rows, fold, clusters, seed=0, k_inner=5, log=print):
    """Cross-fitted label-derived SELPHI columns for training rows (inner folds) and eval rows (all outer-training
    labels), via ksp.crossfit (raises LeakageError on any protein/cluster overlap)."""
    outs_tr, outs_ev = [], []
    for fam in (F.SelphiPSSM, F.SelphiKSEACoreg):
        t0 = time.time()
        tr, ev = crossfit.crossfit(fam, train_labels, train_rows, eval_rows, k_inner=k_inner, seed=seed,
                                   outer_fold=fold, clusters=clusters)
        outs_tr.append(tr)
        outs_ev.append(ev)
        log(f"  crossfit {fam.name}: {time.time() - t0:.0f}s")
    tr = pd.concat(outs_tr, axis=1)
    ev = pd.concat(outs_ev, axis=1) if eval_rows is not None else None
    return tr, ev


# ------------------------------------------------------------------------------------------------ RF grid / bagging
def _tree_probas(rf: RandomForestClassifier, X: np.ndarray) -> np.ndarray:
    """(n_trees, n) positive-class probabilities per tree (prefix means = forests with fewer trees)."""
    from joblib import Parallel, delayed
    def one(t):
        p = t.predict_proba(X)
        return p[:, 1] if p.shape[1] > 1 else np.zeros(len(X))
    return np.array(Parallel(n_jobs=rf.n_jobs, prefer="threads")(delayed(one)(t) for t in rf.estimators_))


def grid_search(X, y, groups, grid, n_jobs=12, seed=0, n_splits=3, log=print) -> tuple[dict, pd.DataFrame]:
    """GridSearchCV(cv=3, roc_auc) equivalent with GroupKFold(3) over substrate clusters. n_estimators is evaluated
    exactly by prefixes of one max-size forest (sklearn draws per-tree seeds sequentially from random_state, so the
    first n trees ARE the n-tree forest with the same seed); max_depth values >= the deepest realised tree
    of the unrestricted fit give identical forests and are not refitted (exact)."""
    from sklearn.model_selection import GroupKFold
    gkf = GroupKFold(n_splits=n_splits)
    nmax = max(grid["n_estimators"])
    rows = []
    for mss in grid["min_samples_split"]:
        for fi, (tr, va) in enumerate(gkf.split(X, y, groups)):
            realised = None
            for md in sorted(grid["max_depth"], reverse=True):
                if realised is not None and md >= realised:
                    rows.extend({**r, "max_depth": md, "reused": True} for r in last)
                    continue
                t0 = time.time()
                rf = RandomForestClassifier(n_estimators=nmax, max_depth=md, min_samples_split=mss, bootstrap=True,
                                            n_jobs=n_jobs, random_state=seed + fi)
                rf.fit(X[tr], y[tr])
                P = _tree_probas(rf, X[va])
                cm = np.cumsum(P, 0)
                last = []
                for ne in grid["n_estimators"]:
                    r = {"min_samples_split": mss, "fold": fi, "n_estimators": ne,
                         "auc": roc_auc_score(y[va], cm[ne - 1] / ne)}
                    last.append(r)
                    rows.append({**r, "max_depth": md, "reused": False})
                d = max(t.get_depth() for t in rf.estimators_)
                if md == max(grid["max_depth"]):
                    realised = d
                log(f"    grid mss={mss} fold={fi} depth={md} (realised {d}) {time.time() - t0:.0f}s "
                    f"auc@{nmax}={last[-1]['auc']:.4f}")
    res = pd.DataFrame(rows)
    agg = res.groupby(["max_depth", "min_samples_split", "n_estimators"]).auc.mean().reset_index()
    best = agg.loc[agg.auc.idxmax()]
    return {"max_depth": int(best.max_depth), "min_samples_split": int(best.min_samples_split),
            "n_estimators": int(best.n_estimators)}, agg


def predict_chunked(models, assembler, rows: pd.DataFrame, lab: pd.DataFrame, cols, chunk=500_000) -> np.ndarray:
    out = np.zeros(len(rows), np.float64)
    ki, si = rows.kinase_idx.to_numpy(), rows.site_idx.to_numpy()
    for a in range(0, len(rows), chunk):
        b = min(a + chunk, len(rows))
        X = assembler.matrix(ki[a:b], si[a:b], lab.iloc[a:b], cols)
        acc = np.zeros(b - a)
        for m in models:
            acc += m.predict_proba(X)[:, 1]
        out[a:b] = acc / len(models)
    return out
