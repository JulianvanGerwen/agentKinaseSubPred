"""Feature store with mandatory ID alignment (fixes r2 AIFM2-1: row-order hazards).

Every feature matrix is saved together with its row/column ID arrays and the SHA-256 of
the canonical index it is aligned to. `load_matrix` re-checks alignment and raises on any
mismatch. Never index a cached array by position without going through this module.

Axes: 'kinase' (cache/index/kinases.parquet: uniprot_acc), 'site' (sites.parquet: site_id),
'protein' (proteins.parquet: acc), or 'none' (free axis; IDs still stored).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ksp.index import INDEX_DIR, ROOT

STORE = ROOT / "cache" / "features"
_AXIS_FILE = {"kinase": ("kinases.parquet", "uniprot_acc"),
              "site": ("sites.parquet", "site_id"),
              "protein": ("proteins.parquet", "acc")}


def _index_ids(axis: str) -> np.ndarray:
    fname, col = _AXIS_FILE[axis]
    return pd.read_parquet(INDEX_DIR / fname, columns=[col])[col].to_numpy()


def _index_hash(axis: str) -> str:
    fname, _ = _AXIS_FILE[axis]
    return hashlib.sha256((INDEX_DIR / fname).read_bytes()).hexdigest()


def reorder_to_index(arr: np.ndarray, ids, axis: str, dim: int, fill=np.nan) -> np.ndarray:
    """Reorder `arr` along `dim` from the given `ids` into canonical index order.
    IDs absent from `ids` are filled with `fill`; IDs not in the index raise."""
    ids = np.asarray(ids)
    canon = _index_ids(axis)
    pos = pd.Index(canon).get_indexer(ids)
    if (pos < 0).any():
        bad = ids[pos < 0][:5]
        raise KeyError(f"{(pos < 0).sum()} ids not in canonical {axis} index, e.g. {bad}")
    shape = list(arr.shape)
    shape[dim] = len(canon)
    out = np.full(shape, fill, dtype=arr.dtype if np.issubdtype(arr.dtype, np.floating) else np.float32)
    sl = [slice(None)] * arr.ndim
    sl[dim] = pos
    out[tuple(sl)] = arr
    return out


def save_matrix(name: str, arr: np.ndarray, row_axis: str, col_axis: str,
                meta: dict | None = None, row_ids=None, col_ids=None) -> Path:
    """Save arr (rows x cols) already in canonical order for its axes.
    For axis 'none', pass row_ids/col_ids explicitly."""
    path = STORE / name
    path.parent.mkdir(parents=True, exist_ok=True)
    rid = _index_ids(row_axis) if row_axis != "none" else np.asarray(row_ids)
    cid = _index_ids(col_axis) if col_axis != "none" else np.asarray(col_ids) if col_ids is not None else None
    if arr.shape[0] != len(rid):
        raise ValueError(f"{name}: {arr.shape[0]} rows but {len(rid)} {row_axis} ids")
    if cid is not None and arr.ndim > 1 and arr.shape[1] != len(cid):
        raise ValueError(f"{name}: {arr.shape[1]} cols but {len(cid)} {col_axis} ids")
    tmp = path.with_suffix(".tmp.npy")
    np.save(tmp, arr)
    os.replace(tmp, path.with_suffix(".npy"))  # atomic
    desc = {"name": name, "shape": list(arr.shape), "dtype": str(arr.dtype),
            "row_axis": row_axis, "col_axis": col_axis,
            "row_index_sha256": _index_hash(row_axis) if row_axis != "none" else None,
            "col_index_sha256": _index_hash(col_axis) if col_axis not in ("none", None) else None,
            "nan_fraction": float(np.isnan(arr).mean()) if np.issubdtype(arr.dtype, np.floating) else None,
            "meta": meta or {}}
    if row_axis == "none":
        np.save(path.with_suffix(".rows.npy"), rid.astype(str))
    if col_axis == "none" and cid is not None:
        np.save(path.with_suffix(".cols.npy"), cid.astype(str))
    path.with_suffix(".json").write_text(json.dumps(desc, indent=1))
    return path.with_suffix(".npy")


def load_matrix(name: str, mmap: bool = True) -> tuple[np.ndarray, dict]:
    """Load a stored matrix; raises if the canonical index changed since it was written."""
    path = STORE / name
    desc = json.loads(path.with_suffix(".json").read_text())
    for ax_key, h_key in (("row_axis", "row_index_sha256"), ("col_axis", "col_index_sha256")):
        ax = desc[ax_key]
        if ax in _AXIS_FILE and desc[h_key] != _index_hash(ax):
            raise RuntimeError(f"{name}: {ax} index changed since this matrix was written; rebuild it")
    arr = np.load(path.with_suffix(".npy"), mmap_mode="r" if mmap else None)
    if list(arr.shape) != desc["shape"]:
        raise RuntimeError(f"{name}: shape {arr.shape} != descriptor {desc['shape']}")
    return arr, desc
