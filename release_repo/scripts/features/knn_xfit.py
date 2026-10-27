"""Dense cross-fitted knn_sim features for R2, built ONLY through ksp.crossfit.crossfit (DESIGN §7, §16 A1).

Same protocol and on-disk layout as ksp.models.r2.xfit (R2's DenseExpand), but GPU-dense and without expanding
pairs into DataFrames:
  train_rows = every U_ms site on the outer-training proteins (one row per site, kinase_acc '__ALL__')
  eval_rows  = every U_ms site on the outer-eval proteins (CV fold f) / on TEST proteins ('final')
The driver assigns inner folds by substrate cluster, fits on labels outside the block, raises LeakageError if any
fitted label shares a protein/cluster with a query site of the block, and transforms eval rows with a randomly
assigned inner-fold fit (A1 default). The writer expands each query SITE to all kinases_with_canaries() (516) via
KnnSim.dense(); this adds no label information (the disjointness is about sites/proteins/clusters).
Labels: crossfit.outer_split(fold, with_canaries=True) (as R2).

Output: cache/features/knn/xfit/<cv|final>/<fold{f}|all>/knn_sim.npy  float16 (n_cols, 516, 350736), NaN = not
computed / undefined; .block.npy (int8 per site: inner block i for training rows, 10+i for eval rows, -2 unused);
.json (columns, tau, timings, label counts). Loader: R2's xfit.load() layout, or knn_xfit.load(fold).
Usage: CUDA_VISIBLE_DEVICES=2 python -I scripts/features/knn_xfit.py --folds 0 1 2 3 4 final
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "6")
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from ksp import canary, crossfit, index, splits  # noqa: E402
from ksp.features import knn  # noqa: E402

torch.set_num_threads(6)
OUT = ROOT / "cache" / "features" / "knn" / "xfit"
SITE_CHUNK = 60000


class KnnDenseWriter(crossfit.LabelDerivedFeature):
    """fit() delegates to KnnSim; transform() writes all-kinase features of the query sites into `buf`."""
    name = "knn_dense_writer"          # not registered (wrapper of the registered knn_sim family)

    def __init__(self, buf, cols, kin_accs, log, spaces):
        self.inner = knn.KnnSim(spaces=spaces)
        self.buf, self.cols, self.kin_accs, self.log = buf, cols, kin_accs, log

    def fit(self, train_pairs, ctx):
        t = time.time()
        self.inner.fit(train_pairs, ctx)
        self._tfit, self._nfit = time.time() - t, len(train_pairs)
        return self

    def transform(self, query, ctx):
        t = time.time()
        sidx = query.site_idx.to_numpy(np.int64)
        sid = query.site_id.to_numpy()
        for a in range(0, len(sid), SITE_CHUNK):
            D = self.inner.dense(sid[a:a + SITE_CHUNK], self.kin_accs, ctx)
            for ci, c in enumerate(self.cols):
                self.buf[ci][:, sidx[a:a + SITE_CHUNK]] = D[c].astype(np.float16)
        self.log.append({"block": ctx.block, "n_fit_pairs": self._nfit, "n_query_sites": int(len(sidx)),
                         "n_fit_pairs_without_rep": self.inner.n_pairs_no_rep,
                         "fit_s": round(self._tfit, 2), "transform_s": round(time.time() - t, 1)})
        b = ctx.block[-1]
        return pd.DataFrame({"blk": np.full(len(query), float(b) if b.isdigit() else -1.0)})


def split_sites(fold):
    S = index.load_sites()
    folds = splits.load_cv_folds()
    f = S.acc.map(folds)
    ums = S.in_U_ms.values
    if fold is None:
        return np.where(ums & f.notna().values)[0], np.where(ums & f.isna().values)[0]   # eval = TEST proteins
    return np.where(ums & f.notna().values & (f.values != fold))[0], np.where(ums & (f.values == fold))[0]


def build(fold, seed=0, spaces=knn.SPACES, overwrite=False):
    tag = "cv" if fold is not None else "final"
    out = OUT / tag / (f"fold{fold}" if fold is not None else "all") / "knn_sim.npy"
    if out.exists() and not overwrite:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cols = knn.columns(spaces)
    train_labels, _, _ = crossfit.outer_split(fold, with_canaries=True)
    train_labels = train_labels[["kinase_acc", "site_id"]].drop_duplicates().reset_index(drop=True)
    S = index.load_sites()
    K = canary.kinases_with_canaries()
    tr_idx, ev_idx = split_sites(fold)
    buf = np.full((len(cols), len(K), len(S)), np.nan, np.float16)
    log = []
    tr_rows = pd.DataFrame({"kinase_acc": "__ALL__", "site_id": S.site_id.values[tr_idx], "site_idx": tr_idx})
    ev_rows = pd.DataFrame({"kinase_acc": "__ALL__", "site_id": S.site_id.values[ev_idx], "site_idx": ev_idx})
    ftr, fev = crossfit.crossfit(lambda: KnnDenseWriter(buf, cols, list(K.uniprot_acc), log, spaces),
                                 train_labels, tr_rows, ev_rows if len(ev_idx) else None, k_inner=5, seed=seed,
                                 outer_fold=fold)
    blk = np.full(len(S), -2, np.int8)
    blk[tr_idx] = ftr.blk.values.astype(np.int8)
    if fev is not None:
        blk[ev_idx] = 10 + fev.blk.values.astype(np.int8)
    tmp = out.with_suffix(".tmp.npy")
    np.save(tmp, buf)
    tmp.rename(out)
    np.save(out.with_suffix(".block.npy"), blk)
    meta = {"family": "knn_sim", "fold": fold, "seed": seed, "variant": "L2", "columns": cols, "level": "site",
            "dtype": "float16", "n_kinases": len(K), "kinase_order": "ksp.canary.kinases_with_canaries()",
            "n_sites": len(S), "site_order": "cache/index/sites.parquet", "n_train_sites": int(len(tr_idx)),
            "n_eval_sites": int(len(ev_idx)), "n_train_labels": int(len(train_labels)),
            "n_canary_labels": int(canary.is_canary(train_labels.kinase_acc).sum()), "tau": knn.load_tau(),
            "spaces": list(spaces), "blocks": log, "runtime_s": round(time.time() - t0, 1),
            "eval_mode": "inner (A1 default)",
            "eval_rows": "U_ms sites of CV fold f" if fold is not None else "U_ms sites on TEST proteins"}
    out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    return out


def load(fold, tag=None):
    tag = tag or ("cv" if fold is not None else "final")
    p = OUT / tag / (f"fold{fold}" if fold is not None else "all") / "knn_sim.npy"
    return np.load(p, mmap_mode="r"), json.loads(p.with_suffix(".json").read_text())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="+", default=["0", "1", "2", "3", "4"])
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    for f in a.folds:
        t = time.time()
        p = build(None if f == "final" else int(f), overwrite=a.overwrite)
        print(f"knn_sim fold {f}: {p} {time.time() - t:.0f}s", flush=True)
