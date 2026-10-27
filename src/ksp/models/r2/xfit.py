"""R2: dense cross-fitted label-derived features, built ONLY through ksp.crossfit.crossfit (DESIGN §7, §16 A1).

Why a wrapper: R2 needs label-derived features for (a) every sampled training pair and (b) EXHAUSTIVE evaluation
pairs (all kinases x all U_ms sites of the eval fold, ~20M pairs/fold). The registered families (F2a LearnedPSSM,
F3 PPINeighbourhood, F6a CoregF6a) compute per *site* for all kinases anyway, so we call `crossfit()` once per
(outer fold, family) with ONE query row per site:
    train_rows = every U_ms site on the outer-training proteins (kinase_acc = '__ALL__')
    eval_rows  = every U_ms site on the outer-eval proteins
The driver assigns inner folds by substrate cluster, fits on labels outside the block, checks protein/cluster
disjointness of fitted labels vs. the block's query sites (LeakageError), and - for eval rows - uses a randomly
assigned inner-fold fit (A1 default eval_mode='inner'). Our `DenseExpand` family then expands each query SITE to
all kinases (516 = 503 real + 13 canaries) and writes the family's output into a dense K x S buffer. Expanding over
kinases at fixed query sites adds no label information: the disjointness guarantee is about sites/proteins/clusters,
and it holds for every kinase at those sites. Sites outside the outer split's U_ms sites stay NaN (never used).

Canary support: the families are subclassed so that canary pseudo-kinases (index 503..515, ksp.canary order) are
fitted and scored like real kinases (LearnedPSSM / CoregF6a otherwise index only the 503 real kinases).
F6a runs on CPU (R2 has no GPU budget); measured ~0.6 ms/query site/block-set.

Output per (fold, family): cache/r2/xfit/<tag>/fold<f>/<family>.npy  float32 (n_cols, K, S) + .json (columns,
inner-fold id per site, timings, input label count). fold=None -> tag 'final' (all DEV labels; eval rows = caller).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ksp import canary, crossfit, index, labels, splits
from ksp.index import ROOT

XFIT_DIR = ROOT / "cache" / "r2" / "xfit"
K_ALL = None  # filled lazily: len(kinases_with_canaries())


def kinases_all() -> pd.DataFrame:
    return canary.kinases_with_canaries()


# ------------------------------------------------------------------------------------------- canary-aware families
def _f2a_class():
    from ksp.features import motif_learned as ml

    class LearnedPSSMCanary(ml.LearnedPSSM):
        name = "r2_f2a_pssm_canary"   # not registered (subclass of the registered f2a_pssm family)

        def fit(self, train_pairs, ctx=None):
            groups = None
            if ctx is not None and getattr(ctx, "clusters", None) is not None and len(ctx.clusters):
                acc = train_pairs.site_id.str.rsplit("_", n=1).str[0]
                groups = pd.Series(acc.map(ctx.clusters).fillna(acc).astype(str).to_numpy(), index=train_pairs.site_id)
                groups = groups[~groups.index.duplicated()]
            seed = getattr(ctx, "seed", 0) if ctx is not None else 0
            self.sites = index.load_sites()
            self.kinases = kinases_all()
            self.m0 = ml.fit_pssm(train_pairs, self.sites, self.kinases, self.prior, False, groups=groups, seed=seed)
            self.m1 = ml.fit_pssm(train_pairs, self.sites, self.kinases, self.prior, True, groups=groups, seed=seed)
            return self
    return LearnedPSSMCanary


def _f2a_sep_class():
    """F2a with real kinases and canaries fitted SEPARATELY (PI 16:58): alpha/beta are tuned on real kinases only, and
    canary PSSMs never pool with, nor tune on, real kinases (each canary has its own group, A25)."""
    from ksp.features import motif_learned as ml
    base = _f2a_class()

    class LearnedPSSMSep(base):
        name = "r2_f2a_pssm_sep"

        def fit(self, train_pairs, ctx=None):
            groups = None
            if ctx is not None and getattr(ctx, "clusters", None) is not None and len(ctx.clusters):
                acc = train_pairs.site_id.str.rsplit("_", n=1).str[0]
                groups = pd.Series(acc.map(ctx.clusters).fillna(acc).astype(str).to_numpy(), index=train_pairs.site_id)
                groups = groups[~groups.index.duplicated()]
            seed = getattr(ctx, "seed", 0) if ctx is not None else 0
            self.sites = index.load_sites()
            self.kinases = kinases_all()
            isc = canary.is_canary(train_pairs.kinase_acc)
            real, can = train_pairs[~isc].reset_index(drop=True), train_pairs[isc].reset_index(drop=True)
            nreal = int((~canary.is_canary(self.kinases.uniprot_acc)).sum())
            for attr, ph in (("m0", False), ("m1", True)):
                mr = ml.fit_pssm(real, self.sites, self.kinases, self.prior, ph, groups=groups, seed=seed)
                if len(can):
                    mc = ml.fit_pssm(can, self.sites, self.kinases, self.prior, ph, groups=groups, seed=seed)
                    for rt in mr.logodds:
                        mr.logodds[rt][nreal:] = mc.logodds[rt][nreal:]
                        mr.n_train[rt][nreal:] = mc.n_train[rt][nreal:]
                        mr.level[rt][nreal:] = mc.level[rt][nreal:]
                    mr.info["canary_alpha_beta"] = (mc.alpha, mc.beta)
                setattr(self, attr, mr)
            return self
    return LearnedPSSMSep


def bench_pmids() -> tuple:
    p = ROOT / "results" / "features" / "phospho" / "exclude_lists.json"
    v = json.loads(p.read_text())["clean_bench"]
    if v.get("exclude_cptac") is not True:
        raise ValueError("clean_bench is expected to drop all CPTAC")
    return tuple(int(x) for x in v["exclude_pmids"])


def _f6a_class(exclude_pmids=(), name="r2_f6a_coreg_canary", no_cptac=False):
    from ksp.features import phospho as ph

    class CoregF6aCanary(ph.CoregF6a):
        device = "cpu"

        def transform(self, query, ctx=None):
            sites = ph.site_index()
            K = kinases_all()
            kpos = {a: i for i, a in enumerate(K.uniprot_acc)}
            spos = {s: j for j, s in enumerate(sites.site_id)}
            tp = pd.DataFrame({"kinase_idx": self.train.kinase_acc.map(kpos), "site_idx": self.train.site_id.map(spos)})
            if tp.isna().any().any():
                raise KeyError("F6a training pairs with unknown kinase/site")
            tp = tp.astype(np.int64)
            sidx = query.site_idx.to_numpy(np.int64)
            kidx = query.kinase_idx.to_numpy(np.int64)
            uq, inv = np.unique(sidx, return_inverse=True)
            clus = ph.site_clusters_from_acc(getattr(ctx, "clusters", None), sites)
            out = {}
            for kind in self.kinds:
                tag = "comp" if kind == "compendium" else "cptac"
                if kind == "cptac" and no_cptac:
                    out[f"f6a_{tag}_top"] = np.full(len(query), np.nan, np.float32)
                    out[f"f6a_{tag}_med"] = np.full(len(query), np.nan, np.float32)
                    continue
                blocks = ph.cached_blocks(kind, self.exclude_pmids if kind == "compendium" else (), False, ())
                r = ph.coreg_features(tp, uq, blocks, cluster_of_site=clus, m=self.m, min_shared=self.min_shared,
                                      min_blocks=self.min_blocks, n_kinases=len(K), device="cpu")
                out[f"f6a_{tag}_top"] = r["top"][kidx, inv]
                out[f"f6a_{tag}_med"] = r["med"][kidx, inv]
            return pd.DataFrame(out)
    CoregF6aCanary.name = name
    CoregF6aCanary.exclude_pmids = tuple(exclude_pmids)
    return CoregF6aCanary


def _f3_class():
    from ksp.features import ppi_label

    class PPINbCanary(ppi_label.PPINeighbourhood):
        name = "r2_f3_ppi_canary"
    return PPINbCanary


def f6_clean_pmids() -> tuple:
    p = ROOT / "results" / "features" / "phospho" / "exclude_lists.json"
    d = json.loads(p.read_text())
    v = d["f6_clean"]
    if v.get("exclude_cptac"):
        raise ValueError("f6_clean unexpectedly excludes CPTAC; R2 wrapper only implements PMID exclusion")
    return tuple(int(x) for x in v["exclude_pmids"])


FAMILIES = {
    # name: (factory-of-class, level, columns); CHUNK[name] = max expanded pairs per inner transform call
    "f2a": (lambda: _f2a_class(), "site",
            ["f2a_pssm", "f2a_pssm_pct", "f2a_pssm_ph", "f2a_pssm_ph_pct", "f2a_pssm_level", "f2a_pssm_log_ntrain"]),
    "f2a_sep": (lambda: _f2a_sep_class(), "site",
                ["f2a_pssm", "f2a_pssm_pct", "f2a_pssm_ph", "f2a_pssm_ph_pct", "f2a_pssm_level", "f2a_pssm_log_ntrain"]),
    "f3nb": (lambda: _f3_class(), "protein",
             ["f3nb_nb_max", "f3nb_nb_mean", "f3nb_nb_top3", "f3nb_nb_n_adj", "f3nb_nb_frac_adj", "f3nb_nb_n_used"]),
    "f6a": (lambda: _f6a_class(), "site", ["f6a_comp_top", "f6a_comp_med", "f6a_cptac_top", "f6a_cptac_med"]),
    "f6a_bench": (lambda: _f6a_class(bench_pmids(), "r2_f6a_coreg_bench_canary", no_cptac=True), "site",
                  ["f6a_comp_top", "f6a_comp_med", "f6a_cptac_top", "f6a_cptac_med"]),
    "f6a_clean": (lambda: _f6a_class(f6_clean_pmids(), "r2_f6a_coreg_f6clean_canary"), "site",
                  ["f6a_comp_top", "f6a_comp_med", "f6a_cptac_top", "f6a_cptac_med"]),
}


CHUNK = {"f2a": 6_000_000, "f2a_sep": 6_000_000, "f3nb": 6_000_000, "f6a": 40_000_000, "f6a_clean": 40_000_000, "f6a_bench": 40_000_000}


def _knn_class():
    from ksp.features import knn

    class KnnSimCPU(knn.KnnSim):
        name = "r2_knn_sim_cpu"

        def __init__(self):
            super().__init__(device="cpu")
    return KnnSimCPU


def _knn_cols():
    from ksp.features import knn
    return knn.columns()


FAMILIES["knn"] = (lambda: _knn_class(), "site", None)       # columns resolved lazily (ksp.features.knn.columns)
CHUNK["knn"] = 40_000_000


class DenseExpand(crossfit.LabelDerivedFeature):
    """Wraps one family: fit() delegates; transform() expands the query SITES to all kinases, writes the family's
    features into `buf` (n_cols, K, S) and returns a single bookkeeping column (the inner-block id)."""
    name = "r2_dense_expand"

    def __init__(self, inner_cls, level: str, buf: np.ndarray, cols: list, log: list, chunk: int = 6_000_000):
        self.inner = inner_cls()
        self.level, self.buf, self.cols, self.log, self.chunk = level, buf, cols, log, chunk

    def fit(self, train_pairs, ctx):
        t = time.time()
        self.inner.fit(train_pairs, ctx)
        self._tfit = time.time() - t
        self._nfit = len(train_pairs)
        return self

    def transform(self, query, ctx):
        t = time.time()
        K = kinases_all()
        nK = len(K)
        sidx = query.site_idx.to_numpy(np.int64)
        if self.level == "protein":
            # protein-level family: one representative site per protein, broadcast back
            acc = splits.site_protein(query.site_id)
            rep = pd.Series(np.arange(len(sidx))).groupby(acc).first()
            rep_rows = rep.values
            q_s = sidx[rep_rows]
            q_sid = query.site_id.values[rep_rows]
        else:
            q_s, q_sid = sidx, query.site_id.values
        n = len(q_s)
        CH = max(1, self.chunk // max(n, 1))       # kinases per sub-call (bounded memory)
        res = np.full((len(self.cols), nK, n), np.nan, np.float32)
        for k0 in range(0, nK, CH):
            ks = np.arange(k0, min(nK, k0 + CH))
            q = pd.DataFrame({"kinase_acc": np.repeat(K.uniprot_acc.values[ks], n),
                              "site_id": np.tile(q_sid, len(ks)),
                              "kinase_idx": np.repeat(ks, n), "site_idx": np.tile(q_s, len(ks))})
            f = self.inner.transform(q, ctx)
            missing = [c for c in self.cols if c not in f.columns]
            if missing:
                raise KeyError(f"family output lacks columns {missing}; has {list(f.columns)}")
            for ci, c in enumerate(self.cols):
                res[ci, ks, :] = f[c].to_numpy(np.float32).reshape(len(ks), n)
        if self.level == "protein":
            pos_of = pd.Series(np.arange(n), index=pd.Index(splits.site_protein(q_sid)))
            back = pos_of.reindex(splits.site_protein(query.site_id)).to_numpy()
            self.buf[:, :, sidx] = res[:, :, back]
        else:
            self.buf[:, :, sidx] = res
        self.log.append({"block": ctx.block, "n_fit_pairs": self._nfit, "n_query_sites": int(len(sidx)),
                         "fit_s": round(self._tfit, 1), "transform_s": round(time.time() - t, 1)})
        return pd.DataFrame({"blk": np.full(len(query), float(ctx.block[-1]) if ctx.block[-1].isdigit() else -1.0)})


def split_sites(fold: int | None):
    """(train_site_idx, eval_site_idx) = U_ms sites on outer-training / outer-eval DEV proteins."""
    S = index.load_sites()
    folds = splits.load_cv_folds()
    f = S.acc.map(folds)
    ums = S.in_U_ms.values
    if fold is None:
        return np.where(ums & f.notna().values)[0], np.zeros(0, np.int64)
    return (np.where(ums & f.notna().values & (f.values != fold))[0],
            np.where(ums & (f.values == fold))[0])


def build(family: str, fold: int | None, seed: int = 0, variant: str = "L2", tag: str | None = None,
          eval_site_idx: np.ndarray | None = None, overwrite: bool = False,
          exclude_kinases: set | None = None, train_site_idx: np.ndarray | None = None,
          train_labels_override: pd.DataFrame | None = None) -> Path:
    """exclude_kinases: kinase accessions whose labels are removed from the fitted labels (leave-kinase-out T2/T3)."""
    tag = tag or ("cv" if fold is not None else "final")
    out = XFIT_DIR / tag / (f"fold{fold}" if fold is not None else "all") / f"{family}.npy"
    if out.exists() and not overwrite:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cls_factory, level, cols = FAMILIES[family]
    if cols is None and family == "knn":
        cols = _knn_cols()
    inner_cls = cls_factory()
    if train_labels_override is not None:          # G5 deployment: DEV+TEST labels from ksp.deploy (fold must be None)
        if fold is not None:
            raise ValueError("train_labels_override requires fold=None")
        train_labels = train_labels_override
    else:
        train_labels, _, _ = crossfit.outer_split(fold, variant=variant, with_canaries=True)
    train_labels = train_labels[["kinase_acc", "site_id"]].drop_duplicates().reset_index(drop=True)
    n_excl = 0
    if exclude_kinases:
        keep = ~train_labels.kinase_acc.isin(set(exclude_kinases))
        n_excl = int((~keep).sum())
        train_labels = train_labels[keep].reset_index(drop=True)
    S = index.load_sites()
    tr_idx, ev_idx = split_sites(fold)
    if eval_site_idx is not None:
        ev_idx = np.asarray(eval_site_idx, np.int64)
    if train_site_idx is not None:
        tr_idx = np.asarray(train_site_idx, np.int64)
    K = kinases_all()
    buf = np.full((len(cols), len(K), len(S)), np.nan, np.float32)      # in RAM (scattered writes; NFS-safe)
    log: list = []
    tr_rows = pd.DataFrame({"kinase_acc": "__ALL__", "site_id": S.site_id.values[tr_idx], "site_idx": tr_idx})
    ev_rows = pd.DataFrame({"kinase_acc": "__ALL__", "site_id": S.site_id.values[ev_idx], "site_idx": ev_idx}) \
        if len(ev_idx) else None
    ftr, fev = crossfit.crossfit(lambda: DenseExpand(inner_cls, level, buf, cols, log, CHUNK[family]), train_labels, tr_rows, ev_rows,
                                 k_inner=5, seed=seed, outer_fold=fold)
    blk = np.full(len(S), -2, np.int8)
    blk[tr_idx] = ftr.blk.values.astype(np.int8)
    if fev is not None:
        blk[ev_idx] = 10 + fev.blk.values.astype(np.int8)     # 10+i: eval row transformed by inner fit i
    np.save(out.with_suffix(".tmp.npy"), buf)
    del buf
    Path(out.with_suffix(".tmp.npy")).rename(out)
    np.save(out.with_suffix(".block.npy"), blk)
    meta = {"family": family, "fold": fold, "seed": seed, "variant": variant, "columns": cols, "level": level,
            "n_kinases": len(K), "n_sites": len(S), "n_train_sites": int(len(tr_idx)), "n_eval_sites": int(len(ev_idx)),
            "n_train_labels": int(len(train_labels)), "n_canary_labels": int(canary.is_canary(train_labels.kinase_acc).sum()),
            "blocks": log, "runtime_s": round(time.time() - t0, 1), "inner_class": inner_cls.__name__,
            "eval_mode": "inner (A1 default)", "n_excluded_kinases": len(exclude_kinases or ()),
            "n_excluded_label_pairs": n_excl}
    out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    return out


def load(family: str, fold: int | None, tag: str | None = None) -> tuple[np.ndarray, dict]:
    tag = tag or ("cv" if fold is not None else "final")
    p = XFIT_DIR / tag / (f"fold{fold}" if fold is not None else "all") / f"{family}.npy"
    meta = json.loads(p.with_suffix(".json").read_text())
    return np.load(p, mmap_mode="r"), meta
