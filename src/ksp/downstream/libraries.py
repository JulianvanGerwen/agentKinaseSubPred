"""Kinase -> substrate libraries for activity inference (DESIGN §6.7, §16 A21).

A library is a DataFrame with unique rows (kinase_acc, site_id) [+ optional 'score'/'rank'], kinase_acc = UniProt
accession of a kinase-list kinase, site_id = canonical site id. Builders:

  psp_all()          all PhosphoSitePlus human-kinase x human-substrate pairs (rows with KIN_ORGANISM = SUB_ORGANISM =
                     human), mapped as in ksp.labels (kinase by accession/symbol via the KinaseResolver; substrate by
                     15-mer on the canonical). INCLUDES TEST-protein pairs: allowed for the downstream benchmark only
                     (PI instruction: the benchmark evaluates libraries); built in memory, never written to disk, never
                     used by any model. Residue class is not enforced (the paper does not).
  psp_dev()          the same definition restricted to DEV proteins (ksp.labels.load_dev('all'), provenance PSP,
                     native_human) -- the leakage-free PSP arm.
  topn_from_matrix() top-N sites per kinase from a ksp.store kinase x site matrix (ties broken by a second matrix).
  cantley_topn()     Cantley/Kinase-Library-only arm: top-N by the as-published percentile (motif/
                     cantley_unprimed_pct_published), ties broken by the raw log2 score.
  from_predictions() top-N per kinase from a predictions parquet (long: kinase_acc|kinase_idx, site_id|site_idx,
                     score) or a store matrix ('store:<name>').
  union(a, b), minus(a, b)   'ours ∪ PSP' and 'ours minus PSP' arms (A21); minus_ranked() = top-N after removing b.
  shuffled(lib, seed) degree-preserving control: a random permutation of the library's distinct sites is applied to
                     all edges (each kinase keeps its set size, each site its number of kinases).
N per kinase is an int or a Series kinase_acc -> N (e.g. psp_sizes(psp_all())). Universe: None (= all U_obs index
sites), 'U_ms', or an explicit iterable of site_ids (e.g. sites measured in the benchmark).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CANTLEY_PCT = "motif/cantley_unprimed_pct_published"
CANTLEY_SCORE = "motif/cantley_unprimed_score"


def _lib(df: pd.DataFrame, name: str) -> pd.DataFrame:
    d = df.drop_duplicates(["kinase_acc", "site_id"]).reset_index(drop=True)
    d.attrs["name"] = name
    return d


def sizes(lib: pd.DataFrame) -> pd.Series:
    return lib.groupby("kinase_acc").site_id.nunique()


psp_sizes = sizes


# ------------------------------------------------------------------------------------------------ PSP
@lru_cache(maxsize=1)
def _psp_all_cached() -> pd.DataFrame:
    from ksp import io_psp, kinases, labels
    from .bench_data import _store
    raw = io_psp.load_psp()
    d = raw[raw.is_sty & raw.kin_organism.eq("human") & raw.sub_organism.eq("human")].copy()
    res = kinases.load_resolver()
    d = labels.resolve_psp_kinases(d, res)
    m = labels.map_psp_substrates(d, _store(), {}, 11)
    m = m.dropna(subset=["kinase", "site_id"])
    acc = dict(zip(res.kl.gene, res.kl.uniprot_acc))
    inc = dict(zip(res.kl.gene, res.kl.include))
    m["kinase_acc"] = m.kinase.map(acc)
    m = m[m.kinase.map(inc).astype(bool)]
    out = m.groupby(["kinase_acc", "site_id"]).agg(in_vivo=("in_vivo", "any"), in_vitro=("in_vitro", "any")) \
        .reset_index()
    out["autophos"] = out.site_id.str.split("_").str[0].eq(out.kinase_acc)
    return out


def psp_all() -> pd.DataFrame:
    return _lib(_psp_all_cached().copy(), "PSP_all")


def psp_dev() -> pd.DataFrame:
    from ksp.labels import load_dev
    d = load_dev("all", ums_only=False)
    d = d[d.provenance.eq("PSP") & d.native_human.astype(bool) & d.kinase_included.astype(bool)]
    out = d[["kinase_acc", "site_id", "in_vivo_native", "in_vitro_native"]].rename(
        columns={"in_vivo_native": "in_vivo", "in_vitro_native": "in_vitro"})
    out["autophos"] = out.site_id.str.split("_").str[0].eq(out.kinase_acc)
    return _lib(out, "PSP_dev")


# ------------------------------------------------------------------------------------------------ top-N builders
def _universe_mask(universe, site_ids: np.ndarray) -> np.ndarray:
    if universe is None:
        return np.ones(len(site_ids), bool)
    if isinstance(universe, str):
        if universe == "U_ms":
            from ksp import index
            return index.load_sites().in_U_ms.to_numpy().astype(bool)
        raise ValueError(universe)
    return np.isin(site_ids, np.asarray(list(universe)))


def _n_vector(n, kinase_accs: np.ndarray) -> np.ndarray:
    if isinstance(n, (int, np.integer)):
        return np.full(len(kinase_accs), int(n))
    s = pd.Series(n)
    return s.reindex(kinase_accs).fillna(0).astype(int).to_numpy()


def topn_from_matrix(name: str, n, universe=None, tiebreak: str | None = None, kinases=None,
                     lib_name: str | None = None) -> pd.DataFrame:
    """Top-n sites per kinase (row) of the store matrix `name` among `universe`; NaN = not scored."""
    from ksp import store
    from ksp.index import load_kinases, load_sites
    X, info = store.load_matrix(name)
    T = store.load_matrix(tiebreak)[0] if tiebreak else None
    kacc = load_kinases().uniprot_acc.to_numpy()
    sids = load_sites().site_id.to_numpy()
    mask = _universe_mask(universe, sids)
    cols = np.flatnonzero(mask)
    nv = _n_vector(n, kacc)
    keep = set(kinases) if kinases is not None else None
    rows = []
    for i, k in enumerate(kacc):
        if nv[i] <= 0 or (keep is not None and k not in keep):
            continue
        x = np.asarray(X[i, cols], dtype=np.float64)
        ok = np.isfinite(x)
        if not ok.any():
            continue
        c, x = cols[ok], x[ok]
        t = np.asarray(T[i, c], dtype=np.float64) if T is not None else np.zeros(len(c))
        t = np.nan_to_num(t, nan=-np.inf)
        m = min(nv[i], len(c))
        o = np.lexsort((t, x))[::-1][:m]              # primary x desc, then tiebreak desc
        rows.append(pd.DataFrame({"kinase_acc": k, "site_id": sids[c[o]], "score": x[o], "rank": np.arange(1, m + 1)}))
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["kinase_acc", "site_id", "score", "rank"])
    return _lib(out, lib_name or f"top{n if isinstance(n, int) else 'N'}:{name}")


def cantley_topn(n, universe=None, kinases=None) -> pd.DataFrame:
    tag = n if isinstance(n, (int, np.integer)) else "PSPsize"
    u = universe if isinstance(universe, (str, type(None))) else "measured"
    return topn_from_matrix(CANTLEY_PCT, n, universe, tiebreak=CANTLEY_SCORE, kinases=kinases,
                            lib_name=f"Cantley_top{tag}_{u or 'U_obs'}")


def cantley_top_kinases_per_site(k: int = 15, universe=None, lib_name: str | None = None) -> pd.DataFrame:
    """The paper's Kinase Library arm (Johnson 2023 convention): for every site in `universe` the k kinases with the
    highest as-published percentile (ties broken by the raw score)."""
    from ksp import store
    from ksp.index import load_kinases, load_sites
    kacc = load_kinases().uniprot_acc.to_numpy()
    sids = load_sites().site_id.to_numpy()
    cols = np.flatnonzero(_universe_mask(universe, sids))
    X = np.array(store.load_matrix(CANTLEY_PCT)[0])[:, cols]          # to RAM (column gathers on a mmap are slow)
    T = np.array(store.load_matrix(CANTLEY_SCORE)[0])[:, cols]
    cols_all, cols = cols, np.arange(len(cols))
    rows = []
    for c0 in range(0, len(cols), 20000):
        cc = cols[c0:c0 + 20000]
        x = np.asarray(X[:, cc], dtype=np.float64)
        t = np.asarray(T[:, cc], dtype=np.float64)
        key = np.where(np.isfinite(x), x + 1e-3 * np.tanh(np.nan_to_num(t) / 10.0), -np.inf)
        top = np.argsort(-key, axis=0, kind="stable")[:k]                       # k x n_cols kinase indices
        ok = np.isfinite(np.take_along_axis(key, top, 0))
        ki, ci = top[ok], np.broadcast_to(np.arange(len(cc)), top.shape)[ok]
        rows.append(pd.DataFrame({"kinase_acc": kacc[ki], "site_id": sids[cols_all[cc[ci]]],
                                  "score": np.take_along_axis(x, top, 0)[ok]}))
    return _lib(pd.concat(rows, ignore_index=True), lib_name or f"KL{k}site")


def read_predictions(path_or_store: str) -> pd.DataFrame:
    """Long predictions table -> columns kinase_acc, site_id, score."""
    if str(path_or_store).startswith("store:"):
        raise ValueError("use topn_from_matrix for store matrices")
    d = pd.read_parquet(path_or_store)
    from ksp.index import load_kinases, load_sites
    if "kinase_acc" not in d:
        if "kinase_idx" in d:
            d["kinase_acc"] = load_kinases().uniprot_acc.to_numpy()[d.kinase_idx.to_numpy()]
        else:
            raise KeyError("predictions need kinase_acc or kinase_idx")
    if "site_id" not in d:
        if "site_idx" in d:
            d["site_id"] = load_sites().site_id.to_numpy()[d.site_idx.to_numpy()]
        else:
            raise KeyError("predictions need site_id or site_idx")
    sc = [c for c in ("score", "prob", "p", "pred") if c in d]
    if not sc:
        raise KeyError("predictions need a score column (score/prob/p/pred)")
    return d[["kinase_acc", "site_id", sc[0]]].rename(columns={sc[0]: "score"})


def topn_from_long(pred: pd.DataFrame, n, universe=None, exclude: pd.DataFrame | None = None,
                   lib_name: str = "pred") -> pd.DataFrame:
    p = pred
    if universe is not None:
        if isinstance(universe, str):
            from ksp.index import load_sites
            s = load_sites()
            p = p[p.site_id.isin(set(s.site_id[s.in_U_ms.astype(bool)]))]
        else:
            p = p[p.site_id.isin(set(universe))]
    if exclude is not None:
        key = set(zip(exclude.kinase_acc, exclude.site_id))
        p = p[[(k, s) not in key for k, s in zip(p.kinase_acc, p.site_id)]]
    p = p[np.isfinite(p.score.astype(float))].sort_values(["kinase_acc", "score", "site_id"],
                                                           ascending=[True, False, True])
    p = p.assign(rank=p.groupby("kinase_acc").cumcount() + 1)
    nv = pd.Series(_n_vector(n, p.kinase_acc.unique()), index=p.kinase_acc.unique()) if not isinstance(n, int) \
        else None
    p = p[p["rank"] <= (n if isinstance(n, int) else p.kinase_acc.map(nv))]
    return _lib(p.reset_index(drop=True), lib_name)


def from_predictions(source: str, n, universe=None, exclude: pd.DataFrame | None = None,
                     lib_name: str | None = None) -> pd.DataFrame:
    """Top-n per kinase from predictions over U_obs x kinases. `source` = parquet path (long) or 'store:<name>'.
    `exclude` (e.g. psp_all()) removes those pairs BEFORE ranking ('ours minus PSP', novel predictions only)."""
    tag = n if isinstance(n, (int, np.integer)) else "PSPsize"
    name = lib_name or f"pred_top{tag}"
    if str(source).startswith("store:"):
        if exclude is not None:
            raise NotImplementedError("exclude with store matrices: export to long parquet first")
        return topn_from_matrix(source[6:], n, universe, lib_name=name)
    return topn_from_long(read_predictions(source), n, universe, exclude, name)


# ------------------------------------------------------------------------------------------------ set operations
def union(a: pd.DataFrame, b: pd.DataFrame, name: str | None = None) -> pd.DataFrame:
    return _lib(pd.concat([a[["kinase_acc", "site_id"]], b[["kinase_acc", "site_id"]]]),
                name or f"{a.attrs.get('name', 'a')}+{b.attrs.get('name', 'b')}")


def minus(a: pd.DataFrame, b: pd.DataFrame, name: str | None = None) -> pd.DataFrame:
    key = set(zip(b.kinase_acc, b.site_id))
    keep = [(k, s) not in key for k, s in zip(a.kinase_acc, a.site_id)]
    return _lib(a[keep], name or f"{a.attrs.get('name', 'a')}-{b.attrs.get('name', 'b')}")


def restrict_kinases(lib: pd.DataFrame, kinases) -> pd.DataFrame:
    return _lib(lib[lib.kinase_acc.isin(set(kinases))], lib.attrs.get("name", "lib"))


def shuffled(lib: pd.DataFrame, seed: int = 0, name: str | None = None) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    s = np.array(sorted(lib.site_id.unique()))
    perm = dict(zip(s, rng.permutation(s)))
    out = pd.DataFrame({"kinase_acc": lib.kinase_acc.to_numpy(), "site_id": lib.site_id.map(perm).to_numpy()})
    return _lib(out, name or f"{lib.attrs.get('name', 'lib')}_shuffled_s{seed}")
