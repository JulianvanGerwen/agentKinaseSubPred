"""Cross-fitting driver for ALL label-derived feature families (DESIGN §7 rule; PI r2 update). Owner: core-infra.

Contract (feature agents implement this; register with @register):

    class MyFamily(LabelDerivedFeature):
        name = "f2a_pssm"                           # unique, used in manifests
        def fit(self, train_pairs: pd.DataFrame, ctx: FitContext) -> "MyFamily":
            # train_pairs: DataFrame[kinase_acc, site_id, (+ any label columns, e.g. provenance, L1, autophos)]
            #   = the ONLY labels you may use. Never call ksp.labels.load_dev() inside fit/transform.
            return self
        def transform(self, query: pd.DataFrame, ctx: FitContext) -> pd.DataFrame:
            # query: DataFrame[kinase_acc, site_id] (+ kinase_idx, site_idx: positions in cache/index orders)
            # return: DataFrame with len(query) rows, SAME ROW ORDER, feature columns only (float; NaN allowed)

Driver (`crossfit`), for one outer split (CV fold f, or final DEV->TEST):
  * eval rows (validation fold f / TEST): DEFAULT (A1) transformed by a randomly assigned inner-fold fit (the same
    fits used for training rows); eval_mode='all' = one fit on ALL outer-training labels (PI decision only).
  * training rows: inner k-fold (default 5) over substrate clusters of the training rows (ksp.splits.inner_folds);
    for inner fold i, fit on training labels whose substrate cluster is NOT in inner fold i, transform rows of fold i.
  * enforcement: for every fit/transform block the driver asserts that no fitted label lies on a protein or substrate
    cluster of a query site in that block (raises LeakageError). Training labels on clusters that carry no training
    row are assigned to an inner fold too (they are used by all other inner folds; never by their own).

`outer_split(fold)` returns (train_labels, train_sites, eval_sites) for CV fold f from ksp.labels.load_dev; for the
final model use fold=None (all DEV labels train; eval rows are supplied by the caller, e.g. TEST U_ms sites).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import splits

REGISTRY: dict = {}


class LeakageError(RuntimeError):
    pass


@dataclass
class FitContext:
    """What a feature family may know besides its training pairs.
    clusters: acc -> substrate cluster id; block: 'eval' or 'inner{i}'; outer_fold: CV fold or None (final);
    seed: run seed; extra: free dict (e.g. paths to feature stores)."""
    clusters: pd.Series
    block: str = "eval"
    outer_fold: int | None = None
    seed: int = 0
    extra: dict = field(default_factory=dict)


class LabelDerivedFeature:
    name = "base"

    def fit(self, train_pairs: pd.DataFrame, ctx: FitContext) -> "LabelDerivedFeature":
        raise NotImplementedError

    def transform(self, query: pd.DataFrame, ctx: FitContext) -> pd.DataFrame:
        raise NotImplementedError


def register(cls):
    if cls.name in REGISTRY:
        raise ValueError(f"feature family {cls.name!r} registered twice")
    REGISTRY[cls.name] = cls
    return cls


def _check_disjoint(fit_pairs: pd.DataFrame, query: pd.DataFrame, clusters: pd.Series, block: str) -> None:
    qp = set(splits.site_protein(query.site_id))
    fp = set(splits.site_protein(fit_pairs.site_id))
    if qp & fp:
        raise LeakageError(f"[{block}] {len(qp & fp)} query proteins also carry fitted labels, e.g. {list(qp & fp)[:3]}")
    qc = set(splits.site_cluster(query.site_id, clusters))
    fc = set(splits.site_cluster(fit_pairs.site_id, clusters))
    if qc & fc:
        raise LeakageError(f"[{block}] {len(qc & fc)} query clusters also carry fitted labels")


def _run(family_factory, fit_pairs, query, ctx) -> pd.DataFrame:
    _check_disjoint(fit_pairs, query, ctx.clusters, ctx.block)
    fam = family_factory()
    fam.fit(fit_pairs.reset_index(drop=True), ctx)
    out = fam.transform(query.reset_index(drop=True), ctx)
    if len(out) != len(query):
        raise ValueError(f"{getattr(fam, 'name', fam)}.transform returned {len(out)} rows for {len(query)} queries")
    return pd.DataFrame(out).reset_index(drop=True)


def crossfit(family_factory, train_labels: pd.DataFrame, train_rows: pd.DataFrame, eval_rows: pd.DataFrame | None,
             k_inner: int = 5, seed: int = 0, outer_fold: int | None = None, clusters: pd.Series | None = None,
             extra: dict | None = None, eval_mode: str = "inner") -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Return (features for train_rows, features for eval_rows), each aligned to the input row order.
    family_factory: zero-arg callable returning a fresh LabelDerivedFeature (e.g. the class itself, or REGISTRY[n]).
    train_labels: outer-training label pairs [kinase_acc, site_id, ...]; train_rows / eval_rows: [kinase_acc,
    site_id(, kinase_idx, site_idx)] rows needing features (pairs; training rows usually = positives + negatives).
    eval_mode='inner' (default, DESIGN §16 A1): each eval row's substrate cluster is randomly assigned (seed) to one
      inner fold i and the row is transformed by the SAME fit as training rows of fold i (labels of the other inner
      folds, ~80% coverage) -> count-type features have the same distribution for training and eval rows.
    eval_mode='all': eval rows use one fit on ALL outer-training labels (only if the PI decides so, e.g. deployment).
    Base (non-label) matrices must be read through ksp.store.load_matrix inside the family (A2)."""
    if eval_mode not in ("inner", "all"):
        raise ValueError("eval_mode must be 'inner' or 'all'")
    cl = clusters if clusters is not None else splits.load_clusters()
    ex = extra or {}
    have_eval = eval_rows is not None and len(eval_rows)
    eval_feat = None
    if have_eval and eval_mode == "all":
        ctx = FitContext(cl, "eval", outer_fold, seed, ex)
        eval_feat = _run(family_factory, train_labels, eval_rows, ctx)
    # inner folds over the clusters of training rows AND training labels (so every label has an inner fold)
    all_sites = pd.concat([train_rows.site_id, train_labels.site_id], ignore_index=True)
    uniq = pd.Series(all_sites.unique())
    f_of = pd.Series(splits.inner_folds(uniq, k_inner, seed, cl), index=uniq.values)
    row_f = train_rows.site_id.map(f_of).values
    lab_f = train_labels.site_id.map(f_of).values
    if have_eval and eval_mode == "inner":
        ec = splits.site_cluster(eval_rows.site_id, cl)
        uc = np.unique(ec)
        rng = np.random.default_rng(seed + 7919)
        cf = dict(zip(uc, rng.integers(0, k_inner, len(uc))))
        ev_f = np.array([cf[c] for c in ec])
    out = ev_out = None
    for i in range(k_inner):
        qi = np.where(row_f == i)[0]
        qe = np.where(ev_f == i)[0] if (have_eval and eval_mode == "inner") else np.zeros(0, int)
        if len(qi) == 0 and len(qe) == 0:
            continue
        ctx = FitContext(cl, f"inner{i}", outer_fold, seed, ex)
        q = pd.concat([train_rows.iloc[qi][["kinase_acc", "site_id"] + [c for c in ("kinase_idx", "site_idx")
                                                                          if c in train_rows]],
                       eval_rows.iloc[qe][["kinase_acc", "site_id"] + [c for c in ("kinase_idx", "site_idx")
                                                                         if c in eval_rows]]
                       if len(qe) else None], ignore_index=True)
        fi = _run(family_factory, train_labels[lab_f != i], q, ctx)
        if out is None:
            out = pd.DataFrame(np.nan, index=range(len(train_rows)), columns=fi.columns)
            if have_eval and eval_mode == "inner":
                ev_out = pd.DataFrame(np.nan, index=range(len(eval_rows)), columns=fi.columns)
        out.iloc[qi] = fi.values[:len(qi)]
        if len(qe):
            ev_out.iloc[qe] = fi.values[len(qi):]
    if have_eval and eval_mode == "inner":
        eval_feat = ev_out
    return (out if out is not None else pd.DataFrame(index=range(len(train_rows)))), eval_feat


def outer_split(fold: int | None, variant: str = "L2", with_canaries: bool = False):
    """(train_labels, train_proteins, eval_proteins) for CV fold `fold` (eval = DEV proteins of that fold).
    fold=None -> all DEV labels, all DEV proteins train, eval proteins = empty (caller supplies TEST rows)."""
    from . import labels
    lab = labels.load_dev(variant, with_canaries=with_canaries)
    folds = splits.load_cv_folds()
    acc = splits.site_protein(lab.site_id)
    if fold is None:
        return lab, set(folds.index), set()
    pf = pd.Series(acc).map(folds).values
    tr = lab[pf != fold].reset_index(drop=True)
    return tr, set(folds.index[folds != fold]), set(folds.index[folds == fold])


@register
class LeakyLabelFlag(LabelDerivedFeature):
    """Toy family used by tests: feature = 1 if (kinase, site) is in the fitted training pairs, else 0. Without
    cross-fitting this leaks the label perfectly (AUROC 1 on training rows); under crossfit() it must be ~0.5."""
    name = "toy_is_labelled"

    def fit(self, train_pairs, ctx):
        self.keys = set(zip(train_pairs.kinase_acc, train_pairs.site_id))
        return self

    def transform(self, query, ctx):
        return pd.DataFrame({"is_labelled": [float((k, s) in self.keys)
                                             for k, s in zip(query.kinase_acc, query.site_id)]})
