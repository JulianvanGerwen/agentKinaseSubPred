"""F3-label: label-derived, cross-fitted network neighbourhood feature (DESIGN §7 F3, last bullet).

substrate_neighbourhood(train_pairs, query_pairs, exclude_same_cluster=cluster_map)
  For each (kinase k, query protein q): proximity of q to k's TRAINING substrate proteins S_k on the screens-only
  PPI graph (ksp.features.ppi), excluding
    * q itself,
    * any s in S_k with cluster_map[s] == cluster_map[q] (substrate homology cluster; proteins absent from the map
      are their own singleton cluster),
    * s == k (autophosphorylation partner = the kinase itself).
  Proximity of q from seed s = degree-corrected RWR, log2((R[q,s]+eps)/(bg[q]+eps)), R = full RWR matrix
  (restart 0.5) precomputed by scripts/features/context_ppi.py (cache/features/context/ppi/rwr_full.npy),
  bg[q] = mean_s R[q,s].
  Outputs per pair: nb_max, nb_mean, nb_top3 (mean of the 3 largest), nb_n_adj (number of usable training substrates
  that are direct neighbours of q), nb_frac_adj (= nb_n_adj / nb_n_used), nb_n_used (usable training substrates in graph).
  Missing (NaN): q not in graph, or nb_n_used == 0.

The function is pure w.r.t. labels: it uses only `train_pairs`; graph data are label-free. The caller is responsible
for cross-fitting (train_pairs must never contain pairs whose substrate is in the evaluation fold).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
PPI = ROOT / "cache" / "features" / "context" / "ppi"
OUT_NAMES = ["nb_max", "nb_mean", "nb_top3", "nb_n_adj", "nb_frac_adj", "nb_n_used"]


class GraphData:
    def __init__(self, ppi_dir: Path = PPI):
        self.nodes = np.array((ppi_dir / "rwr_nodes.txt").read_text().split())
        self.pos = {a: i for i, a in enumerate(self.nodes)}
        R = np.load(ppi_dir / "rwr_full.npy")
        bg = np.load(ppi_dir / "rwr_background.npy").astype(np.float32)
        eps = np.float32(0.01 / len(self.nodes))
        # DC[q, s] = log2((R[q,s]+eps)/(bg[q]+eps)); computed once, in place (float32, N x N)
        R += eps
        R /= (bg + eps)[:, None]
        np.log2(R, out=R)
        self.DC = R
        e = pd.read_parquet(ppi_dir / "edges.parquet")
        i = e.a.map(self.pos).to_numpy()
        j = e.b.map(self.pos).to_numpy()
        n = len(self.nodes)
        A = sp.coo_matrix((np.ones(len(i), np.float32), (i, j)), shape=(n, n)).tocsr()
        self.A = ((A + A.T) > 0).astype(np.float32).tocsr()


@lru_cache(maxsize=1)
def load_graph_data() -> GraphData:
    return GraphData()


def _as_pairs(x) -> pd.DataFrame:
    if isinstance(x, pd.DataFrame):
        cols = [c for c in ("kinase_acc", "kinase") if c in x.columns][:1] + [c for c in ("acc", "protein") if c in x.columns][:1]
        d = x[cols].copy()
        d.columns = ["kinase_acc", "acc"]
        return d.astype(str)
    return pd.DataFrame(list(x), columns=["kinase_acc", "acc"]).astype(str)


def substrate_neighbourhood(train_pairs, query_pairs=None, exclude_same_cluster: dict | None = None,
                            graph: GraphData | None = None, kinases=None, proteins=None):
    """train_pairs: DataFrame[kinase_acc, acc] (or iterable of (kinase_acc, acc)).
    query_pairs: DataFrame[kinase_acc, acc] -> returns DataFrame with OUT_NAMES columns (row-aligned to query_pairs);
                 None -> dense mode over `kinases` x `proteins` (default: cache/index order) -> dict name -> (K x P) float32.
    exclude_same_cluster: dict acc -> cluster id; None = only the identical protein is excluded."""
    g = graph or load_graph_data()
    cmap = exclude_same_cluster or {}
    tr = _as_pairs(train_pairs).drop_duplicates()
    tr = tr[tr.kinase_acc != tr.acc]
    tr = tr[tr.acc.isin(g.pos)]
    train_by_k = {k: grp.acc.to_numpy() for k, grp in tr.groupby("kinase_acc")}

    # cluster ids as integers (singletons for unmapped proteins)
    cid = {}

    def cl(a):
        c = cmap.get(a)
        key = ("c", c) if c is not None else ("p", a)
        return cid.setdefault(key, len(cid))

    def block(k, q_accs):
        """features for kinase k and an array of query accs (all in graph) -> (len(q), 6)."""
        out = np.full((len(q_accs), 6), np.nan, np.float32)
        S = train_by_k.get(k)
        if S is None or len(S) == 0:
            out[:, 3:6] = 0
            out[:, 4] = np.nan
            return out
        s_idx = np.array([g.pos[a] for a in S])
        q_idx = np.array([g.pos[a] for a in q_accs])
        s_cl = np.array([cl(a) for a in S])
        q_cl = np.array([cl(a) for a in q_accs])
        D = g.DC[np.ix_(q_idx, s_idx)]  # q x s
        mask = (q_cl[:, None] != s_cl[None, :]) & (q_idx[:, None] != s_idx[None, :])
        n_used = mask.sum(1)
        Dm = np.where(mask, D, -np.inf)
        mx = Dm.max(1)
        sm = np.where(mask, D, 0).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = sm / n_used
            if Dm.shape[1] >= 3:
                top = -np.partition(-Dm, 2, axis=1)[:, :3]
            else:
                top = Dm
            topf = np.where(np.isfinite(top), top, np.nan)
            with np.errstate(all="ignore"):
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    top3 = np.nanmean(topf, 1) if topf.size else np.full(len(q_idx), np.nan)
        adj = g.A[q_idx][:, s_idx].toarray() > 0
        n_adj = (adj & mask).sum(1)
        ok = n_used > 0
        out[ok, 0] = mx[ok]
        out[ok, 1] = mean[ok]
        out[ok, 2] = top3[ok]
        out[:, 3] = n_adj
        out[ok, 4] = n_adj[ok] / n_used[ok]
        out[:, 5] = n_used
        return out

    if query_pairs is not None:
        q = _as_pairs(query_pairs).reset_index(drop=True)
        res = np.full((len(q), 6), np.nan, np.float32)
        ing = q.acc.isin(g.pos).to_numpy()
        for k, idx in q[ing].groupby("kinase_acc").groups.items():
            idx = np.asarray(idx)
            res[idx] = block(k, q.acc.to_numpy()[idx])
        return pd.DataFrame(res, columns=OUT_NAMES)

    if kinases is None or proteins is None:
        from ksp.index import load_kinases, load_proteins
        kinases = load_kinases().uniprot_acc.tolist() if kinases is None else kinases
        proteins = load_proteins().acc.tolist() if proteins is None else proteins
    proteins = np.asarray(proteins)
    ing = np.array([a in g.pos for a in proteins])
    q_accs = proteins[ing]
    out = {n: np.full((len(kinases), len(proteins)), np.nan, np.float32) for n in OUT_NAMES}
    for i, k in enumerate(kinases):
        b = block(k, q_accs)
        for c, n in enumerate(OUT_NAMES):
            out[n][i, ing] = b[:, c]
    return out


# ---------------------------------------------------------------- ksp.crossfit contract (PI r2)
from ksp.crossfit import LabelDerivedFeature, register  # noqa: E402


@register
class PPINeighbourhood(LabelDerivedFeature):
    """F3-label under the ksp.crossfit contract. fit() keeps the training (kinase, substrate protein) pairs;
    transform() returns OUT_NAMES (prefixed 'f3nb_') for query (kinase_acc, site_id) rows. Protein-level: every site
    of a protein gets the same value. ctx.clusters (acc -> cluster) is used as exclude_same_cluster (redundant with the
    driver's disjointness check, kept as a second guard)."""
    name = "f3_ppi_neighbourhood"

    def fit(self, train_pairs, ctx):
        acc = train_pairs.site_id.astype(str).str.rsplit("_", n=1).str[0] if "acc" not in train_pairs else train_pairs.acc
        self.train = pd.DataFrame({"kinase_acc": train_pairs.kinase_acc.astype(str).values, "acc": np.asarray(acc)})
        return self

    def transform(self, query, ctx):
        qacc = query.site_id.astype(str).str.rsplit("_", n=1).str[0].values
        qp = pd.DataFrame({"kinase_acc": query.kinase_acc.astype(str).values, "acc": qacc})
        uq = qp.drop_duplicates().reset_index(drop=True)
        cmap = ctx.clusters.to_dict() if getattr(ctx, "clusters", None) is not None else None
        f = substrate_neighbourhood(self.train, uq, exclude_same_cluster=cmap)
        f.columns = ["f3nb_" + c for c in f.columns]
        f.index = pd.MultiIndex.from_frame(uq)
        return f.reindex(pd.MultiIndex.from_frame(qp)).reset_index(drop=True)
