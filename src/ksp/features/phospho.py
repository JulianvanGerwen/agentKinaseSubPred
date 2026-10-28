"""F6 phosphoproteomics features: data blocks, GPU masked correlation, co-regulation (F6a).

DESIGN §7 F6. Docs: docs/features_phospho.md. Owner: feat-phospho agent.

Data blocks
-----------
A *block* is a site x column float32 matrix (NaN = not observed) whose rows are site indices into the canonical
site index (``ksp.index.load_sites()``). Two kinds:

* ``compendium``: perturbation log2FC (data/01-03), columns = conditions, one block per study (PMID); the one very
  large study (decryptM, PMID 36926954, 790 contrasts) is split per drug series; studies with < ``min_block_cols``
  conditions are pooled into ``pooled_small``. Each condition column is divided by its robust SD (1.4826*MAD), see
  ``conditions.parquet`` column ``scale``.
* ``cptac``: tumour phospho (BCM v1.2), one block per cohort, residualised per site on host-protein abundance and the
  available confounders (see ``ksp.features.phospho_cptac``).

Blocks live in ``cache/features/phospho/blocks/<block_id>.{X,rows}.npy`` with a manifest
``cache/features/phospho/blocks/manifest.parquet`` (block_id, kind, pmids, cohort, n_rows, n_cols).

Clean variants (DESIGN §6.7): every loader/builder takes ``exclude_pmids`` (drop compendium columns of those PMIDs)
and ``exclude_cptac`` (True = drop all CPTAC blocks; or an iterable of cohort names).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
CACHE = ROOT / "cache" / "features" / "phospho"
BLOCK_DIR = CACHE / "blocks"
RESULTS = ROOT / "results" / "features" / "phospho"


# ----------------------------------------------------------------------------------------------- small helpers
_TIME_RX = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*(s|sec|min|m|h|hr|hrs|d|days?)?\s*$", re.I)
_TIME_IN_TEXT = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)\s*(sec|s|min|h|hr|d)(?![a-zA-Z])", re.I)
_UNIT_MIN = {"s": 1 / 60, "sec": 1 / 60, "min": 1.0, "m": 1.0, "h": 60.0, "hr": 60.0, "hrs": 60.0,
             "d": 1440.0, "day": 1440.0, "days": 1440.0}


def parse_minutes(timepoint, condition: str | None = None) -> tuple[float, str]:
    """Treatment time in minutes from the compendium ``Timepoint`` field (fallback: the Condition text).
    Returns (minutes or NaN, how). ASSUMPTION: a unit-less number (Ochoa-2016-sourced rows: '60', '20', '240') is
    minutes."""
    tp = "" if timepoint is None or (isinstance(timepoint, float) and math.isnan(timepoint)) else str(timepoint)
    m = _TIME_RX.match(tp)
    if m:
        unit = (m.group(2) or "min").lower()
        return float(m.group(1)) * _UNIT_MIN[unit], "timepoint" if m.group(2) else "timepoint_unitless_as_min"
    if condition:
        hits = _TIME_IN_TEXT.findall(str(condition))
        if len(hits) == 1:
            v, u = hits[0]
            return float(v) * _UNIT_MIN[u.lower()], "condition_text"
    return float("nan"), "unknown" if tp in ("", "NA", "nan") else f"unparsed:{tp}"


def site_index() -> pd.DataFrame:
    from ksp.index import load_sites
    return load_sites()


def kinase_index() -> pd.DataFrame:
    from ksp.index import load_kinases
    return load_kinases()


def protein_of_site(sites: pd.DataFrame | None = None) -> np.ndarray:
    """int32 protein code per canonical site (codes = order of ksp.index.load_proteins())."""
    from ksp.index import protein_pos
    sites = site_index() if sites is None else sites
    pp = protein_pos()
    return sites.acc.map(pp).to_numpy().astype(np.int32)


# ----------------------------------------------------------------------------------------------- blocks
@dataclass
class Block:
    block_id: str
    kind: str                      # 'compendium' | 'cptac'
    rows: np.ndarray               # int32 site indices (sorted, unique)
    X: np.ndarray                  # float32 (len(rows), n_cols), NaN = missing
    cols: list = field(default_factory=list)
    pmids: tuple = ()
    cohort: str | None = None

    @property
    def n_cols(self) -> int:
        return self.X.shape[1]

    def subset_cols(self, keep: np.ndarray) -> "Block":
        X = self.X[:, keep]
        ok = ~np.isnan(X).all(axis=1)
        return Block(self.block_id, self.kind, self.rows[ok], X[ok], [c for c, k in zip(self.cols, keep) if k],
                     self.pmids, self.cohort)


def save_block(b: Block, out_dir: Path = BLOCK_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.+-]", "_", b.block_id)
    np.save(out_dir / f"{safe}.X.npy", b.X.astype(np.float32))
    np.save(out_dir / f"{safe}.rows.npy", b.rows.astype(np.int32))
    pd.Series(b.cols, dtype=object).to_frame("col").to_parquet(out_dir / f"{safe}.cols.parquet")
    return dict(block_id=b.block_id, kind=b.kind, file=safe, pmids=";".join(map(str, b.pmids)),
                cohort=b.cohort or "", n_rows=int(len(b.rows)), n_cols=int(b.n_cols),
                n_values=int(np.isfinite(b.X).sum()))


def load_blocks(kinds: Sequence[str] = ("compendium", "cptac"), exclude_pmids: Iterable = (),
                exclude_cptac: bool | Iterable[str] = False, block_dir: Path = BLOCK_DIR,
                conditions: pd.DataFrame | None = None, exclude_drugs: Iterable[str] = (),
                exclude_condition_ids: Iterable[str] = ()) -> list[Block]:
    """Load cached blocks. ``exclude_pmids`` removes the compendium columns (conditions) of those PMIDs from every
    block (a block that loses all columns is dropped); ``exclude_cptac`` = True drops all CPTAC blocks, or give a list
    of cohort names to drop only those."""
    man = pd.read_parquet(block_dir / "manifest.parquet")
    excl = {int(p) for p in exclude_pmids}
    excl_c = set(exclude_condition_ids)
    if exclude_drugs:
        from .phospho_inhib import condition_ids_for_drugs
        excl_c |= condition_ids_for_drugs(exclude_drugs)
    if conditions is None and excl:
        conditions = pd.read_parquet(CACHE / "conditions.parquet")
    if exclude_cptac is True:
        drop_coh = set(man.cohort)
    elif exclude_cptac:
        drop_coh = {str(c).upper() for c in exclude_cptac}
    else:
        drop_coh = set()
    out = []
    for r in man.itertuples():
        if r.kind not in kinds:
            continue
        if r.kind == "cptac" and r.cohort.upper() in drop_coh:
            continue
        X = np.load(block_dir / f"{r.file}.X.npy", mmap_mode="r")
        rows = np.load(block_dir / f"{r.file}.rows.npy")
        cols = pd.read_parquet(block_dir / f"{r.file}.cols.parquet").col.tolist()
        b = Block(r.block_id, r.kind, rows, np.asarray(X), cols,
                  tuple(int(p) for p in str(r.pmids).split(";") if p), r.cohort or None)
        if r.kind == "compendium" and (excl or excl_c):
            if conditions is None:
                conditions = pd.read_parquet(CACHE / "conditions.parquet")
            pm = conditions.set_index("conditionID").PMID.reindex(cols).astype("Int64")
            keep = ~pm.isin(list(excl)).to_numpy() & ~np.isin(np.array(cols, dtype=object), list(excl_c))
            if not keep.any():
                continue
            if not keep.all():
                b = b.subset_cols(keep)
        out.append(b)
    return out


def load_conditions() -> pd.DataFrame:
    return pd.read_parquet(CACHE / "conditions.parquet")


# ----------------------------------------------------------------------------------------------- GPU kernels
def _torch():
    import torch
    return torch


def to_dev(X: np.ndarray, device):
    """float32 tensor with NaN, rows centred by their own nan-mean (reduces cancellation in the sum formulas)."""
    torch = _torch()
    t = torch.as_tensor(np.ascontiguousarray(X, dtype=np.float32), device=device)
    mu = torch.nanmean(t, dim=1, keepdim=True)
    mu = torch.nan_to_num(mu, nan=0.0)
    return t - mu


def masked_corr(A, B, min_n: int = 8, prepared: tuple | None = None):
    """Pairwise-complete Pearson correlation between rows of A (a x n) and rows of B (b x n), torch tensors with NaN.
    Returns (r, n) as (a x b) tensors; r is NaN where n < min_n or a variance is ~0.
    Six masked matmuls (fp32, TF32 disabled by the caller's torch defaults)."""
    torch = _torch()
    MA = (~torch.isnan(A)).float()
    A0 = torch.nan_to_num(A, nan=0.0)
    if prepared is None:
        MB = (~torch.isnan(B)).float()
        B0 = torch.nan_to_num(B, nan=0.0)
        prepared = (MB, B0, B0 * B0)
    MB, B0, B02 = prepared
    n = MA @ MB.T
    sa = A0 @ MB.T
    sb = MA @ B0.T
    saa = (A0 * A0) @ MB.T
    sbb = MA @ B02.T
    sab = A0 @ B0.T
    nn = n.clamp(min=1.0)
    cov = sab - sa * sb / nn
    va = saa - sa * sa / nn
    vb = sbb - sb * sb / nn
    den = torch.sqrt((va * vb).clamp(min=0.0))
    r = cov / den
    bad = (n < min_n) | (va <= 1e-6 * nn) | (vb <= 1e-6 * nn)
    r = torch.where(bad, torch.full_like(r, float("nan")), r.clamp(-1.0, 1.0))
    return r, n


def fisher_z(r, clip: float = 0.999):
    torch = _torch()
    return torch.atanh(r.clamp(-clip, clip))


# ----------------------------------------------------------------------------------------------- F6a
def _pairs_to_idx(train_pairs: pd.DataFrame, kin: pd.DataFrame, sites: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    if {"kinase_idx", "site_idx"} <= set(train_pairs.columns):
        return train_pairs.kinase_idx.to_numpy(np.int64), train_pairs.site_idx.to_numpy(np.int64)
    kp = {a: i for i, a in enumerate(kin.uniprot_acc)}
    sp = {s: j for j, s in enumerate(sites.site_id)}
    ki = train_pairs.kinase_acc.map(kp)
    si = train_pairs.site_id.map(sp)
    ok = ki.notna() & si.notna()
    return ki[ok].to_numpy(np.int64), si[ok].to_numpy(np.int64)


def coreg_features(train_pairs: pd.DataFrame, query_site_idx: np.ndarray, blocks: Sequence[Block], *,
                   cluster_of_site: np.ndarray | None = None, protein_of_site_arr: np.ndarray | None = None,
                   exclude_same_protein: bool = True, exclude_pmids: Iterable = (),
                   exclude_cptac: bool | Iterable[str] = False, m: int = 3, min_shared: int = 8,
                   min_blocks: int = 2, n_kinases: int | None = None, device: str = "cuda",
                   chunk_elems: int = 160_000_000, return_numpy: bool = True, verbose: bool = False) -> dict:
    """F6a co-regulation with training substrates (label-derived; the CALLER cross-fits).

    train_pairs   : DataFrame with (kinase_idx, site_idx) or (kinase_acc, site_id) — the training KS pairs only.
    query_site_idx: int array of canonical site indices to score.
    blocks        : data blocks (``load_blocks``); pass compendium and CPTAC blocks in separate calls to get separate
                    features (recommended; F6a_comp / F6a_cptac).
    cluster_of_site: int array (n_sites,) substrate-cluster id per canonical site (DESIGN §4.1); training substrates in
                    the query's cluster are excluded. Default: no cluster exclusion beyond protein.
    exclude_same_protein: exclude training substrates on the query's protein (always also excludes the query itself).
    exclude_pmids / exclude_cptac: clean-variant filters applied to the given blocks (block metadata).

    Per block: Pearson r (pairwise complete, >= min_shared shared columns) of the query with every training substrate
    of kinase k (after exclusions), Fisher z; aggregate per kinase = mean of the top-m z ('top') and the median z
    ('med'). Across blocks: weighted mean with weight = number of columns in which the query is observed in that block;
    NaN unless >= min_blocks blocks contribute.

    Returns {'top': (n_kinases, n_query) float32, 'med': ..., 'n_blocks': uint8} (numpy unless return_numpy=False).
    """
    torch = _torch()
    sites = site_index()
    n_sites = len(sites)
    if n_kinases is None:
        n_kinases = len(kinase_index())
    kin = kinase_index() if not {"kinase_idx", "site_idx"} <= set(train_pairs.columns) else None
    ki, si = _pairs_to_idx(train_pairs, kin, sites)
    prot = protein_of_site(sites) if protein_of_site_arr is None else np.asarray(protein_of_site_arr)
    clus = None if cluster_of_site is None else np.asarray(cluster_of_site)
    query = np.asarray(query_site_idx, dtype=np.int64)
    nq = len(query)
    qpos_of_site = np.full(n_sites, -1, dtype=np.int64)
    qpos_of_site[query] = np.arange(nq)

    # filter blocks (clean variant)
    excl = {int(p) for p in exclude_pmids}
    use_blocks = []
    for b in blocks:
        if b.kind == "cptac" and (exclude_cptac is True or (exclude_cptac and b.cohort in set(exclude_cptac))):
            continue
        if b.kind == "compendium" and excl and set(b.pmids) & excl:
            if set(b.pmids) <= excl:
                continue
            cond = load_conditions().set_index("conditionID").PMID
            keep = ~cond.reindex(b.cols).isin(list(excl)).to_numpy()
            b = b.subset_cols(keep)
        use_blocks.append(b)

    # training substrates: unique sites, per kinase lists of positions into T
    T = np.unique(si)
    tpos = {s: i for i, s in enumerate(T)}
    nT = len(T)
    per_k: dict[int, list[int]] = {}
    for k, s in zip(ki, si):
        per_k.setdefault(int(k), []).append(tpos[int(s)])
    per_k = {k: sorted(set(v)) for k, v in per_k.items()}
    # bucket kinases by substrate count (padding index nT -> NaN column)
    buckets = []
    sizes = sorted(per_k.items(), key=lambda kv: len(kv[1]))
    edges = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 1 << 20]
    lo = 0
    for e in edges:
        grp = [(k, v) for k, v in sizes if lo < len(v) <= e]
        lo = e
        if not grp:
            continue
        P = max(len(v) for _, v in grp)
        idx = np.full((len(grp), P), nT, dtype=np.int64)
        for i, (_, v) in enumerate(grp):
            idx[i, :len(v)] = v
        buckets.append((np.array([k for k, _ in grp], dtype=np.int64), torch.as_tensor(idx, device=device)))

    dev = torch.device(device)
    acc_top = torch.zeros((n_kinases, nq), dtype=torch.float32, device=dev)
    acc_med = torch.zeros_like(acc_top)
    acc_w = torch.zeros_like(acc_top)
    acc_n = torch.zeros((n_kinases, nq), dtype=torch.uint8, device=dev)
    T_t = torch.as_tensor(T, device=dev)
    prot_t = torch.as_tensor(prot.astype(np.int64), device=dev)
    clus_t = None if clus is None else torch.as_tensor(clus.astype(np.int64), device=dev)

    for b in use_blocks:
        if b.n_cols < min_shared:
            continue
        row_of_site = {int(s): i for i, s in enumerate(b.rows)}
        # training-substrate matrix (nT x n_cols), NaN where substrate not in block
        XT = np.full((nT, b.n_cols), np.nan, dtype=np.float32)
        tin = np.array([row_of_site.get(int(s), -1) for s in T])
        have = tin >= 0
        if have.sum() == 0:
            continue
        XT[have] = b.X[tin[have]]
        nobsT = np.isfinite(XT).sum(1)
        XT[nobsT < min_shared] = np.nan
        if (nobsT >= min_shared).sum() == 0:
            continue
        Bt = to_dev(XT, dev)
        prep = ((~torch.isnan(Bt)).float(), torch.nan_to_num(Bt, nan=0.0))
        prep = (prep[0], prep[1], prep[1] * prep[1])
        # queries in this block
        qmask = qpos_of_site[b.rows] >= 0
        qrows = np.nonzero(qmask)[0]
        nobs_q = np.isfinite(b.X[qrows]).sum(1)
        qrows = qrows[nobs_q >= min_shared]
        nobs_q = nobs_q[nobs_q >= min_shared]
        if len(qrows) == 0:
            continue
        chunk = max(256, int(chunk_elems // max(nT, 1) // 4))
        for c0 in range(0, len(qrows), chunk):
            rr = qrows[c0:c0 + chunk]
            qs = b.rows[rr].astype(np.int64)
            A = to_dev(b.X[rr], dev)
            r, _ = masked_corr(A, None, min_n=min_shared, prepared=prep)
            z = fisher_z(r)
            del r
            qs_t = torch.as_tensor(qs, device=dev)
            # exclusions: the query itself, same protein, same cluster
            ex = qs_t[:, None] == T_t[None, :]
            if exclude_same_protein:
                ex |= prot_t[qs_t][:, None] == prot_t[T_t][None, :]
            if clus_t is not None:
                ex |= clus_t[qs_t][:, None] == clus_t[T_t][None, :]
            z = torch.where(ex, torch.full_like(z, float("nan")), z)
            z = torch.cat([z, torch.full((z.shape[0], 1), float("nan"), device=dev)], dim=1)
            w = torch.as_tensor(nobs_q[c0:c0 + chunk].astype(np.float32), device=dev)
            qp = torch.as_tensor(qpos_of_site[qs], device=dev)
            for kids, idx in buckets:
                nk, P = idx.shape
                sub = max(1, int(chunk_elems // max(nk * P, 1)))
                for s0 in range(0, z.shape[0], sub):
                    zz = z[s0:s0 + sub][:, idx]  # (c, nk, P)
                    valid = ~torch.isnan(zz)
                    nval = valid.sum(-1)
                    mm = min(m, P)
                    top = torch.topk(torch.where(valid, zz, torch.full_like(zz, -float("inf"))), mm, dim=-1).values
                    topv = torch.isfinite(top)
                    top_mean = torch.where(topv, top, torch.zeros_like(top)).sum(-1) / topv.sum(-1).clamp(min=1)
                    med = torch.nanmedian(zz, dim=-1).values
                    ok = nval > 0
                    ww = w[s0:s0 + sub][:, None].expand_as(top_mean) * ok
                    kk = torch.as_tensor(kids, device=dev)
                    qq = qp[s0:s0 + sub]
                    ii = (kk[None, :].expand_as(top_mean), qq[:, None].expand_as(top_mean))
                    acc_top.index_put_(ii, torch.where(ok, top_mean * ww, torch.zeros_like(top_mean)), accumulate=True)
                    acc_med.index_put_(ii, torch.where(ok, med * ww, torch.zeros_like(med)), accumulate=True)
                    acc_w.index_put_(ii, ww, accumulate=True)
                    acc_n.index_put_(ii, ok.to(torch.uint8), accumulate=True)
            del z, A, ex
        if verbose:
            print(f"  block {b.block_id}: {len(qrows)} queries, cols {b.n_cols}", flush=True)
        del Bt, prep
    good = acc_n >= min_blocks
    nan = torch.tensor(float("nan"), device=dev)
    top = torch.where(good, acc_top / acc_w.clamp(min=1e-9), nan)
    med = torch.where(good, acc_med / acc_w.clamp(min=1e-9), nan)
    out = {"top": top, "med": med, "n_blocks": acc_n}
    if return_numpy:
        out = {k: v.cpu().numpy() for k, v in out.items()}
    del acc_top, acc_med, acc_w
    torch.cuda.empty_cache()
    return out


# ----------------------------------------------------------------------------------------------- feature store
def feature_path(name: str) -> Path:
    return CACHE / f"{name}.npy"


def save_feature(name: str, M: np.ndarray, meta: dict | None = None) -> Path:
    """kinase x site float16 matrix aligned to ksp.index (503 x 350,736)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    p = feature_path(name)
    np.save(p, M.astype(np.float16))
    if meta is not None:
        import json
        (CACHE / f"{name}.json").write_text(json.dumps(meta, indent=1, default=str))
    return p


def load_feature(name: str, mmap: bool = True) -> np.ndarray:
    return np.load(feature_path(name), mmap_mode="r" if mmap else None)


# ----------------------------------------------------------------------------------------------- crossfit contract
_BLOCK_CACHE: dict = {}


def cached_blocks(kind: str, exclude_pmids=(), exclude_cptac=False, exclude_drugs=()) -> list:
    key = (kind, tuple(sorted(int(p) for p in exclude_pmids)),
           exclude_cptac if isinstance(exclude_cptac, bool) else tuple(sorted(exclude_cptac)),
           tuple(sorted(exclude_drugs)))
    if key not in _BLOCK_CACHE:
        _BLOCK_CACHE[key] = load_blocks(kinds=(kind,), exclude_pmids=exclude_pmids, exclude_cptac=exclude_cptac,
                                        exclude_drugs=exclude_drugs)
    return _BLOCK_CACHE[key]


def site_clusters_from_acc(clusters: pd.Series | None, sites: pd.DataFrame | None = None) -> np.ndarray | None:
    """acc -> cluster Series (FitContext.clusters) to an int64 cluster id per canonical site (-1 = unknown, never
    matches another site because unknowns get unique negative ids)."""
    if clusters is None:
        return None
    sites = site_index() if sites is None else sites
    c = sites.acc.map(clusters)
    out = np.empty(len(sites), np.int64)
    miss = c.isna().to_numpy()
    out[~miss] = c[~miss].astype(np.int64).to_numpy()
    out[miss] = -1 - np.arange(miss.sum())
    return out


try:
    from ksp.crossfit import LabelDerivedFeature, register  # noqa: E402
except Exception:  # pragma: no cover - crossfit module not available
    LabelDerivedFeature, register = object, (lambda c: c)


@register
class CoregF6a(LabelDerivedFeature):
    """F6a co-regulation with training substrates (ksp.crossfit contract).

    fit(train_pairs, ctx): keeps the training KS pairs (kinase_acc, site_id).
    transform(query, ctx): for the unique query sites, runs ``coreg_features`` on the compendium blocks and on the
    CPTAC blocks separately (substrates on the query's protein and substrate cluster excluded; ctx.clusters used) and
    returns per query pair: f6a_comp_top, f6a_comp_med, f6a_cptac_top, f6a_cptac_med (Fisher z; NaN = <2 blocks).
    Class attributes (override by subclassing or ``make_f6a``): kinds, exclude_pmids, exclude_cptac, m, min_shared,
    min_blocks, device.
    """
    name = "f6a_coreg"
    kinds = ("compendium", "cptac")
    exclude_pmids: tuple = ()
    exclude_cptac = False
    exclude_drugs: tuple = ()
    m, min_shared, min_blocks = 3, 8, 2
    device = "cuda"

    def fit(self, train_pairs: pd.DataFrame, ctx=None):
        self.train = train_pairs[["kinase_acc", "site_id"]].drop_duplicates().reset_index(drop=True)
        return self

    def transform(self, query: pd.DataFrame, ctx=None) -> pd.DataFrame:
        sites = site_index()
        kin = kinase_index()
        if "site_idx" in query.columns:
            sidx = query.site_idx.to_numpy(np.int64)
        else:
            sidx = query.site_id.map({s: j for j, s in enumerate(sites.site_id)}).to_numpy()
        if "kinase_idx" in query.columns:
            kidx = query.kinase_idx.to_numpy(np.int64)
        else:
            kidx = query.kinase_acc.map({a: i for i, a in enumerate(kin.uniprot_acc)}).to_numpy()
        uq, inv = np.unique(sidx, return_inverse=True)
        clus = site_clusters_from_acc(getattr(ctx, "clusters", None), sites) if ctx is not None else None
        out = {}
        for kind in self.kinds:
            tag = "comp" if kind == "compendium" else "cptac"
            if kind == "cptac" and self.exclude_cptac is True:
                out[f"f6a_{tag}_top"] = np.full(len(query), np.nan, np.float32)
                out[f"f6a_{tag}_med"] = np.full(len(query), np.nan, np.float32)
                continue
            blocks = cached_blocks(kind, self.exclude_pmids if kind == "compendium" else (),
                                   self.exclude_cptac if kind == "cptac" else False,
                                   self.exclude_drugs if kind == "compendium" else ())
            r = coreg_features(self.train, uq, blocks, cluster_of_site=clus, m=self.m, min_shared=self.min_shared,
                               min_blocks=self.min_blocks, n_kinases=len(kin), device=self.device)
            out[f"f6a_{tag}_top"] = r["top"][kidx, inv]
            out[f"f6a_{tag}_med"] = r["med"][kidx, inv]
        return pd.DataFrame(out)


def make_f6a(name: str, **attrs):
    """Factory for clean variants, e.g. make_f6a('f6a_coreg_clean', exclude_pmids=(...), exclude_cptac=True)."""
    cls = type(f"CoregF6a_{name}", (CoregF6a,), dict(name=name, **attrs))
    return cls
