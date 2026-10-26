"""Leakage-safe splits (DESIGN §4): substrate clusters, TEST, DEV 5-fold CV, kinase folds (LKO), inner folds.

Frozen outputs (results/splits/, written by scripts/build_splits.py):
  substrate_clusters.tsv  acc, cluster_id (int), n_sites, n_label_sites        (every protein of cache/index/proteins)
  test_proteins.tsv       acc, cluster_id                                       (TEST proteins)
  cv_folds.tsv            acc, cluster_id, fold (0..4)                          (DEV proteins only)
  kinase_folds.tsv        kinase_acc, gene, kinase_cluster, manning_family, lko_fold (0..4)
  eval_kinase_sets.tsv    kinase_acc, gene, in_K_eval, in_K_any, in_K_cv, n_test_pos, n_test_clusters, ...
  split_report.md         numbers + constraint checks

Readers here never touch the sealed TEST labels (read only by ksp.eval_test).

Definitions
  * Substrate clusters = connected components of the protein graph whose edges are (a) mmseqs2 all-vs-all hits with
    fident >= 0.40 and min(qcov, tcov) >= 0.50, and (b) "identical-window" edges: two proteins share any identical
    U_obs 15-mer (window15, padding '_' included). If one component holds > 3% of label sites, its mmseqs edges are
    re-thresholded at 50% identity (identical-window edges are always kept) and this is reported.
  * Labels used for split balancing: L2 pairs of included kinases (in class), plus the human-native evaluable TEST
    positives definition (L1, not autophos, in class) for the per-kinase constraint.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[2]
SPLIT_DIR = ROOT / "results" / "splits"
CACHE = ROOT / "cache" / "splits"
MMSEQS = "/home/tools/mmseqs/bin/mmseqs"
N_FOLDS = 5


# ------------------------------------------------------------------------------------------------ sequences / mmseqs
def write_protein_fasta(path: Path, accs, seqs: dict) -> None:
    missing = [a for a in accs if a not in seqs]
    if missing:
        raise KeyError(f"{len(missing)} proteins without canonical sequence, e.g. {missing[:5]}")
    with open(path, "w") as fh:
        for a in accs:
            fh.write(f">{a}\n{seqs[a]}\n")


def run_mmseqs_allvsall(fasta: Path, out_m8: Path, threads: int = 12, sens: float = 7.5, max_seqs: int = 2000,
                        min_id: float = 0.40, cov: float = 0.50) -> None:
    """mmseqs2 easy-search all-vs-all. Output columns: query,target,fident,alnlen,qcov,tcov,evalue,bits.
    --cov-mode 0 (coverage of query AND target >= cov), --min-seq-id min_id, e-value 1e-3."""
    tmp = out_m8.parent / "mmseqs_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    cmd = [MMSEQS, "easy-search", str(fasta), str(fasta), str(out_m8), str(tmp), "--threads", str(threads),
           "-s", str(sens), "--max-seqs", str(max_seqs), "-e", "1e-3", "--min-seq-id", str(min_id),
           "-c", str(cov), "--cov-mode", "0",
           "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    shutil.rmtree(tmp, ignore_errors=True)


def read_m8(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None,
                       names=["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"])


# ------------------------------------------------------------------------------------------------ frozen-split readers
def load_clusters(path: Path | None = None) -> pd.Series:
    """acc -> substrate cluster_id (int) for every protein of cache/index/proteins.parquet."""
    p = path or SPLIT_DIR / "substrate_clusters.tsv"
    if not Path(p).exists():
        raise FileNotFoundError(f"{p} missing: splits not built yet (scripts/build_splits.py)")
    d = pd.read_csv(p, sep="\t", comment="#")
    return pd.Series(d.cluster_id.values, index=d.acc.values, name="cluster_id")


def load_test_proteins() -> set:
    d = pd.read_csv(SPLIT_DIR / "test_proteins.tsv", sep="\t", comment="#")
    return set(d.acc)


def load_cv_folds() -> pd.Series:
    """acc -> CV fold (0..4) for DEV proteins (TEST proteins absent)."""
    d = pd.read_csv(SPLIT_DIR / "cv_folds.tsv", sep="\t", comment="#")
    return pd.Series(d.fold.values, index=d.acc.values, name="fold")


def load_kinase_folds() -> pd.DataFrame:
    return pd.read_csv(SPLIT_DIR / "kinase_folds.tsv", sep="\t", comment="#")


def load_eval_kinase_sets() -> pd.DataFrame:
    return pd.read_csv(SPLIT_DIR / "eval_kinase_sets.tsv", sep="\t", comment="#")


def site_protein(site_ids_) -> np.ndarray:
    return pd.Series(list(site_ids_), dtype=object).str.rsplit("_", n=1).str[0].values


def site_cluster(site_ids_, clusters: pd.Series | None = None) -> np.ndarray:
    """Substrate cluster of each site (raises on proteins without a cluster)."""
    cl = clusters if clusters is not None else load_clusters()
    acc = site_protein(site_ids_)
    c = pd.Series(acc).map(cl)
    if c.isna().any():
        bad = pd.Series(acc)[c.isna()].unique()[:5]
        raise KeyError(f"{c.isna().sum()} sites on proteins without substrate cluster, e.g. {list(bad)}")
    return c.values.astype(np.int64)


def inner_folds(train_sites, k: int = 5, seed: int = 0, clusters: pd.Series | None = None) -> np.ndarray:
    """Inner fold (0..k-1) per training site for cross-fitting (DESIGN §7): whole substrate clusters are assigned to
    folds (sites of one protein / homologous proteins / identical 15-mers always share a fold). Clusters are shuffled
    with `seed`, then assigned greedily (largest first) to the fold with fewest sites -> balanced sizes.
    Returns an int array aligned to `train_sites`."""
    cid = site_cluster(train_sites, clusters)
    u, cnt = np.unique(cid, return_counts=True)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(u))
    u, cnt = u[perm], cnt[perm]
    order = np.argsort(-cnt, kind="stable")
    load = np.zeros(k)
    fold_of = {}
    for i in order:
        f = int(np.argmin(load + rng.random(k) * 1e-9))
        fold_of[u[i]] = f
        load[f] += cnt[i]
    return np.array([fold_of[c] for c in cid], dtype=np.int64)


# ------------------------------------------------------------------------------------------------ cluster construction
def substrate_clusters(m8: pd.DataFrame, proteins: pd.Series, sites: pd.DataFrame, label_sites: pd.Series,
                       min_id: float = 0.40, min_cov: float = 0.50, max_label_frac: float = 0.03,
                       rethreshold_id: float = 0.50) -> tuple[pd.DataFrame, dict]:
    """Connected components of (mmseqs edges fident>=min_id & min(qcov,tcov)>=min_cov) U (identical U_obs 15-mer
    edges). Components holding > max_label_frac of `label_sites` (one entry per label pair) get their mmseqs edges
    re-thresholded at `rethreshold_id`. Returns (DataFrame acc, cluster_id ordered by protein index; report)."""
    proteins = pd.Series(list(proteins))
    n = len(proteins)
    pos = pd.Series(np.arange(n), index=proteins.values)
    m = m8[(m8["query"] != m8["target"]) & (m8.fident >= min_id) & (np.minimum(m8.qcov, m8.tcov) >= min_cov)]
    m = m[m["query"].isin(pos.index) & m["target"].isin(pos.index)]
    qi, ti, fid = pos[m["query"]].values, pos[m["target"]].values, m.fident.values
    g = sites.groupby("window15").acc.unique()
    g = g[g.map(len) > 1]
    ei, ej = [], []
    for accs in g:
        a = pos[accs].values
        ei += [a[0]] * (len(a) - 1)
        ej += list(a[1:])
    ei, ej = np.array(ei, dtype=np.int64), np.array(ej, dtype=np.int64)

    def comps(keep):
        i = np.r_[qi[keep], ei]
        j = np.r_[ti[keep], ej]
        A = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
        return connected_components(A, directed=False)[1]

    keep = np.ones(len(qi), bool)
    c = comps(keep)
    lab_acc = label_sites.str.rsplit("_", n=1).str[0]
    lab_acc = lab_acc[lab_acc.isin(pos.index)]
    fr = pd.Series(c[pos[lab_acc].values]).value_counts() / max(len(lab_acc), 1)
    big = list(fr.index[fr > max_label_frac])
    rep = {"n_proteins": n, "mmseqs_edges": int(len(qi)), "identical_window_groups": int(len(g)),
           "identical_window_edges": int(len(ei)), "max_component_label_frac_before": float(fr.max()),
           "components_rethresholded": len(big)}
    if big:
        keep = (fid >= rethreshold_id) | ~np.isin(c[qi], big)
        c = comps(keep)
        fr = pd.Series(c[pos[lab_acc].values]).value_counts() / max(len(lab_acc), 1)
    # relabel clusters by first protein occurrence (deterministic)
    _, first = np.unique(c, return_index=True)
    order = np.argsort(first)
    remap = np.empty(c.max() + 1, dtype=np.int64)
    remap[np.unique(c)[order]] = np.arange(len(order))
    cid = remap[c]
    sizes = pd.Series(cid).value_counts()
    rep.update({"n_clusters": int(len(sizes)), "max_component_label_frac_after": float(fr.max()),
                "largest_cluster_proteins": int(sizes.max()), "singleton_clusters": int((sizes == 1).sum())})
    return pd.DataFrame({"acc": proteins.values, "cluster_id": cid}), rep


# ------------------------------------------------------------------------------------------------ annealing optimiser
class _Anneal:
    """Assign groups (clusters) to bins minimising soft penalties:
      per kinase k with total n_k >= kin_min, per penalised bin b: W_band * hinge^2 / n_k + W_tgt * (t/n - tgt_b)^2,
        hinge = distance of count t outside [floor((tgt_b - band) n), ceil((tgt_b + band) n)]
      per global quantity g, per penalised bin: W_hard * max(0, |frac - tgt| - tol_g)^2 + W_soft * (frac - tgt)^2."""

    def __init__(self, kin_counts: list, glob: np.ndarray, n_kin: int, target: np.ndarray, pen_bins: np.ndarray,
                 tol: np.ndarray, kin_min: int = 5, band: float = 0.05, W_band=1.0, W_tgt=0.5, W_hard=1e6,
                 W_soft=1e3, glob_w: np.ndarray | None = None):
        self.kc = kin_counts                      # per group: (kin_idx array, count array)
        self.G = glob.astype(float)               # groups x n_glob
        self.nb = len(target)
        self.tgt = target
        self.pen = pen_bins
        self.tol = tol
        self.band = band
        self.Wb, self.Wt, self.Wh, self.Ws = W_band, W_tgt, W_hard, W_soft
        self.gw = glob_w if glob_w is not None else np.ones(glob.shape[1])
        N = np.zeros(n_kin)
        for ki, cnt in kin_counts:
            np.add.at(N, ki, cnt)
        self.N = N
        self.active = N >= kin_min
        self.Gtot = self.G.sum(0)
        self.lo = np.floor((target[None, :] - band) * N[:, None])
        self.hi = np.ceil((target[None, :] + band) * N[:, None])

    def kpen(self, k, b, t):
        if not self.active[k] or not self.pen[b]:
            return 0.0
        n = self.N[k]
        h = max(0.0, self.lo[k, b] - t, t - self.hi[k, b])
        return self.Wb * h * h / n + self.Wt * (t / n - self.tgt[b]) ** 2

    def gpen(self, b, gb):
        if not self.pen[b]:
            return 0.0
        f = gb / np.maximum(self.Gtot, 1e-12)
        d = np.abs(f - self.tgt[b])
        return float((self.gw * (self.Wh * np.maximum(0, d - self.tol) ** 2 + self.Ws * d ** 2)).sum())

    def run(self, assign: np.ndarray, n_iter: int, seed: int, T0=1.0, T1=1e-4, movable: np.ndarray | None = None,
            p_labeled: float = 0.85):
        rng = np.random.default_rng(seed)
        a = assign.copy()
        T = np.zeros((len(self.N), self.nb))
        for g, (ki, cnt) in enumerate(self.kc):
            np.add.at(T, (ki, np.full(len(ki), a[g])), cnt)
        GB = np.zeros((self.nb, self.G.shape[1]))
        np.add.at(GB, a, self.G)
        lab = np.array([len(ki) > 0 for ki, _ in self.kc])
        mov = movable if movable is not None else np.ones(len(a), bool)
        lab_idx = np.where(lab & mov)[0]
        unl_idx = np.where(~lab & mov)[0]
        temps = T0 * (T1 / T0) ** (np.arange(n_iter) / max(n_iter - 1, 1))
        for it in range(n_iter):
            if len(unl_idx) == 0 or (len(lab_idx) and rng.random() < p_labeled):
                g = lab_idx[rng.integers(len(lab_idx))]
            else:
                g = unl_idx[rng.integers(len(unl_idx))]
            src = a[g]
            dst = (src + 1 + rng.integers(self.nb - 1)) % self.nb
            ki, cnt = self.kc[g]
            d = 0.0
            for k, c in zip(ki, cnt):
                d += (self.kpen(k, src, T[k, src] - c) + self.kpen(k, dst, T[k, dst] + c)
                      - self.kpen(k, src, T[k, src]) - self.kpen(k, dst, T[k, dst]))
            gs_new, gd_new = GB[src] - self.G[g], GB[dst] + self.G[g]
            d += (self.gpen(src, gs_new) + self.gpen(dst, gd_new) - self.gpen(src, GB[src]) - self.gpen(dst, GB[dst]))
            if d <= 0 or rng.random() < np.exp(-d / temps[it]):
                a[g] = dst
                for k, c in zip(ki, cnt):
                    T[k, src] -= c
                    T[k, dst] += c
                GB[src], GB[dst] = gs_new, gd_new
        return a

    def total(self, a):
        T = np.zeros((len(self.N), self.nb))
        for g, (ki, cnt) in enumerate(self.kc):
            np.add.at(T, (ki, np.full(len(ki), a[g])), cnt)
        GB = np.zeros((self.nb, self.G.shape[1]))
        np.add.at(GB, a, self.G)
        return (sum(self.kpen(k, b, T[k, b]) for k in range(len(self.N)) for b in range(self.nb)) +
                sum(self.gpen(b, GB[b]) for b in range(self.nb)))


def _group_kin_counts(site_cluster_: np.ndarray, kin_idx: np.ndarray, n_groups: int) -> list:
    df = pd.DataFrame({"g": site_cluster_, "k": kin_idx}).groupby(["g", "k"]).size().reset_index(name="n")
    out = [(np.zeros(0, np.int64), np.zeros(0)) for _ in range(n_groups)]
    for g, sub in df.groupby("g"):
        out[g] = (sub.k.values.astype(np.int64), sub.n.values.astype(float))
    return out


def optimise_assignment(cluster_ids: np.ndarray, kin_site_cluster: np.ndarray, kin_idx: np.ndarray, n_kin: int,
                        glob: np.ndarray, target: np.ndarray, pen_bins: np.ndarray, tol: np.ndarray, seed: int = 0,
                        n_iter: int = 400_000, kin_min: int = 5, init_p: np.ndarray | None = None,
                        glob_w: np.ndarray | None = None) -> np.ndarray:
    """Return bin per cluster (aligned to cluster_ids = 0..n-1). kin_* arrays: one entry per (kinase, positive)
    pair: its cluster and kinase index. glob: clusters x quantities (e.g. labels, U_ms sites, residue labels)."""
    n = len(cluster_ids)
    kc = _group_kin_counts(kin_site_cluster, kin_idx, n)
    rng = np.random.default_rng(seed)
    p = init_p if init_p is not None else target
    a0 = rng.choice(len(target), size=n, p=p / p.sum())
    an = _Anneal(kc, glob, n_kin, target, pen_bins, tol, kin_min=kin_min, glob_w=glob_w)
    return an.run(a0, n_iter, seed)


# ------------------------------------------------------------------------------------------------ kinase folds (LKO)
def kinase_clusters(kinases: pd.DataFrame, identity: pd.DataFrame, min_id: float = 0.50) -> np.ndarray:
    """Connected components of (KinCore domain identity >= min_id) U (same Manning family; '*-Unique' never joins)."""
    from .labels import _real_family
    acc = kinases.uniprot_acc.values
    pos = pd.Series(np.arange(len(acc)), index=acc)
    e = identity[(identity.identity >= min_id) & identity.acc_a.isin(pos.index) & identity.acc_b.isin(pos.index)]
    i, j = list(pos[e.acc_a].values), list(pos[e.acc_b].values)
    fam = pd.Series([_real_family(f) for f in kinases.manning_family])
    for f, idx in fam.dropna().groupby(fam.dropna()).groups.items():
        idx = list(idx)
        i += [idx[0]] * (len(idx) - 1)
        j += idx[1:]
    n = len(acc)
    A = coo_matrix((np.ones(len(i)), (np.array(i, dtype=np.int64), np.array(j, dtype=np.int64))), shape=(n, n))
    return connected_components(A, directed=False)[1]


def kinase_folds(kclust: np.ndarray, n_labels: np.ndarray, k: int = N_FOLDS, seed: int = 0) -> np.ndarray:
    """Greedy balanced assignment of kinase clusters to k folds (largest label count first; ties by kinase count),
    balancing labels and kinases. No cluster (hence no family) is split."""
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({"c": kclust, "n": n_labels, "one": 1}).groupby("c").agg(n=("n", "sum"), m=("one", "sum"))
    df = df.iloc[rng.permutation(len(df))]
    df = df.sort_values(["n", "m"], ascending=False, kind="stable")
    lab, kin = np.zeros(k), np.zeros(k)
    tl, tk = max(df.n.sum(), 1), max(df.m.sum(), 1)
    fold_of = {}
    for c, r in df.iterrows():
        f = int(np.argmin(lab / tl + kin / tk))
        fold_of[c] = f
        lab[f] += r.n
        kin[f] += r.m
    return np.array([fold_of[c] for c in kclust], dtype=np.int64)


# ------------------------------------------------------------------------------------------------ training helpers
def train_valid_split_inner(train_sites, valid_frac: float = 0.2, seed: int = 0,
                            clusters: pd.Series | None = None) -> np.ndarray:
    """Boolean mask (True = early-stopping/validation row) INSIDE an outer training fold, grouped by substrate
    cluster (AIFM2-7): use for early stopping / tuning; never stop or tune on the outer evaluation fold."""
    k = max(2, int(round(1 / valid_frac)))
    return inner_folds(train_sites, k=k, seed=seed, clusters=clusters) == 0


def cv_eval_rows(fold: int):
    """Human-native CV evaluation positives (identical for every label variant, AIFM2-7) for CV fold `fold`:
    DEV 'eval' labels on U_ms sites of the fold's proteins."""
    from . import index, labels
    ev = labels.load_dev("eval")
    ums = set(index.load_sites().query("in_U_ms").site_id)
    ev = ev[ev.site_id.isin(ums)]
    folds = load_cv_folds()
    return ev[pd.Series(site_protein(ev.site_id)).map(folds).values == fold].reset_index(drop=True)
