"""knn_sim: similarity of a query site to kinase k's TRAINING substrate sites (label-derived; MUST be cross-fitted).

DESIGN §7 cross-fitting rule / §16 A1: registered `ksp.crossfit` family. Labels enter ONLY through fit(train_pairs).
Docs: docs/features_knn.md. Owner: feat-knn.

Representation spaces (all label-free; built by scripts/features/knn_build_reps.py, cached in cache/features/knn/):
  esmc    ESM-2 650M centre-residue embedding (full-protein context), centred per residue type (S/T vs Y mean over
          U_ms sites), L2-normalised -> cosine similarity.
  esmw    ESM-2 650M embedding mean-pooled over the +-7 window (in-protein positions), centred + L2-normalised ->
          cosine.
  esmcat  concatenated +-7 window: per-position ESM-2 embedding projected on a global residue PCA-128, L2-normalised,
          scaled by sqrt(w_p / sum w) (w_p as for blosum; padding -> 0) -> dot product = centre-weighted mean of
          per-position cosines. DEFAULT window space (beat esmw on fold-0 inner data); esmw is optional.
  blosum  BLOSUM62 15-mer similarity: sum_p w_p B62(q_p, t_p) / sum_p w_p, centre-weighted w_p = 1 / (1 + |p|/4)
          (p = -7..+7, centre included; padding '_' / non-standard residues contribute 0). Computed as a matmul of
          one-hot windows through a block-diagonal metric.
  (R3's site-encoder output was not saved by R3 - only predictions - so that optional space is skipped.)

For query site s (residue type r = S/T or Y) and kinase k, T = k's training sites of residue type r MINUS every
training site on the query's own substrate cluster (and protein; done inside transform, so transform is robust even
outside the driver - the driver additionally raises LeakageError if such labels are present at all). Per space:
  knn_<sp>_max     max_t sim(s, t)
  knn_<sp>_top3    mean of the top-3 similarities (fewer if |T| < 3)
  knn_<sp>_kde     log( sum_t exp(sim/tau) / (|T| + c) ),  c = 1 (shrinks kinases with few labels)
  knn_<sp>_kde_bg  log( mean_b exp(sim/tau) ) over a fixed label-free background: n_bg = 2000 random U_ms sites of
                   the same residue type (seed 0), minus sites on the query's own cluster. Site-only (kinase-free).
  knn_<sp>_lr      kde - kde_bg (log density ratio: kinase's substrates vs random sites of the same residue type)
  knn_log_ntrain   log1p(#training sites of k with residue type r) (before exclusion)
NaN if k has no remaining training site of the query's residue type or the site has no representation (non-U_ms).
tau per space is read from results/features/knn/tau.json (chosen ONCE on CV fold 0 inner data by
scripts/features/knn_tune_tau.py; documented in docs/features_knn.md).

Speed: query sites x unique training sites similarity by chunked fp16 matmul on GPU (fp32 on CPU), then grouped
gathers over kinases bucketed by training-set size; all 516 kinases x 347k U_ms sites in well under a minute on one
RTX 4090 (measured numbers in docs/features_knn.md).
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ksp import splits
from ksp.crossfit import LabelDerivedFeature, register
from ksp.index import ROOT, load_sites

KNN_DIR = ROOT / "cache" / "features" / "knn"
TAU_FILE = ROOT / "results" / "features" / "knn" / "tau.json"
ALL_SPACES = ("esmc", "esmw", "esmcat", "blosum")
SPACES = ("esmc", "esmcat", "blosum")    # default set, chosen on fold-0 inner data (results/features/knn/tau.json)
STATS = ("max", "top3", "kde", "kde_bg", "lr")
AA20 = "ACDEFGHIKLMNPQRSTVWY"
FLANK = 7


def columns(spaces=SPACES) -> list:
    return [f"knn_{sp}_{st}" for sp in spaces for st in STATS] + ["knn_log_ntrain"]


def blosum_weights(flank: int = FLANK) -> np.ndarray:
    p = np.arange(-flank, flank + 1)
    return 1.0 / (1.0 + np.abs(p) / 4.0)


def blosum_metric(flank: int = FLANK) -> np.ndarray:
    """(L*20, L*20) block-diagonal metric: sim(q, t) = onehot(q) @ M @ onehot(t)^T."""
    from Bio.Align import substitution_matrices
    B = substitution_matrices.load("BLOSUM62")
    B20 = np.array([[B[a][b] for b in AA20] for a in AA20], np.float32)
    w = blosum_weights(flank)
    L = len(w)
    M = np.zeros((L * 20, L * 20), np.float32)
    for i in range(L):
        M[i * 20:(i + 1) * 20, i * 20:(i + 1) * 20] = B20 * (w[i] / w.sum())
    return M


def onehot_windows(windows) -> np.ndarray:
    """(n, L*20) float16 one-hot of the windows (non-standard residues / padding -> zero block)."""
    w = np.asarray(windows).astype("U")
    L = len(w[0])
    lut = np.full(256, -1, np.int64)
    for i, a in enumerate(AA20):
        lut[ord(a)] = i
    codes = lut[np.frombuffer("".join(w).encode("ascii"), np.uint8).reshape(len(w), L)]
    out = np.zeros((len(w), L * 20), np.float16)
    r, p = np.nonzero(codes >= 0)
    out[r, p * 20 + codes[r, p]] = 1
    return out


def _centre_normalise(X: np.ndarray, rt: np.ndarray) -> np.ndarray:
    """Subtract the per-residue-type mean (S/T vs Y), L2-normalise rows; float16."""
    out = np.empty(X.shape, np.float16)
    for t in (0, 1):
        idx = np.where(rt == t)[0]
        mu = np.zeros(X.shape[1], np.float64)
        for a in range(0, len(idx), 50000):
            mu += X[idx[a:a + 50000]].astype(np.float64).sum(0)
        mu /= max(len(idx), 1)
        for a in range(0, len(idx), 50000):
            Z = X[idx[a:a + 50000]].astype(np.float32) - mu.astype(np.float32)
            Z /= np.maximum(np.linalg.norm(Z, axis=1, keepdims=True), 1e-6)
            out[idx[a:a + 50000]] = Z.astype(np.float16)
    return out


class SiteReps:
    """Site representations, ID-indexed (never by position: rows matched via site_id).
    mats[space]: (n, d) query-side matrix; metric[space]: None (cosine of unit vectors) or (d, d) metric applied
    to the training side."""

    def __init__(self, site_ids, rt, mats: dict, metric: dict | None = None):
        self.site_ids = np.asarray(site_ids).astype(str)
        self.rt = np.asarray(rt, np.int8)
        self.mats = mats
        self.metric = metric or {}
        self.pos = pd.Index(self.site_ids)
        if not self.pos.is_unique:
            raise ValueError("duplicate site ids in SiteReps")
        self.acc = splits.site_protein(self.site_ids)
        self._dev: dict = {}
        self._codes: dict = {}

    @property
    def spaces(self):
        return tuple(self.mats)

    @classmethod
    def from_cache(cls, spaces=SPACES) -> "SiteReps":
        sids = np.load(KNN_DIR / "reps_site_ids.npy", allow_pickle=True).astype(str)
        S = load_sites()
        spos = pd.Index(S.site_id).get_indexer(sids)
        if (spos < 0).any():
            raise KeyError("knn reps contain site ids absent from cache/index/sites.parquet")
        rt = (S.res.to_numpy()[spos] == "Y").astype(np.int8)
        mats, metric = {}, {}
        for sp in spaces:
            if sp == "esmc":
                mats[sp] = _centre_normalise(np.load(KNN_DIR / "reps_esm_centre.npy", mmap_mode="r"), rt)
            elif sp == "esmw":
                mats[sp] = _centre_normalise(np.load(KNN_DIR / "reps_esm_win7.npy", mmap_mode="r"), rt)
            elif sp == "esmcat":
                mats[sp] = np.load(KNN_DIR / "reps_esm_cat7.npy")
            elif sp == "blosum":
                w15 = S.window15.to_numpy()[spos]
                mats[sp] = onehot_windows(w15)
                metric[sp] = blosum_metric(FLANK)
            else:
                raise ValueError(f"unknown space {sp!r}")
        return cls(sids, rt, mats, metric)

    def tensors(self, device: str):
        import torch
        if device not in self._dev:
            dt = torch.float16 if device.startswith("cuda") else torch.float32
            self._dev[device] = {
                "mats": {sp: torch.as_tensor(np.asarray(m)).to(device=device, dtype=dt) for sp, m in self.mats.items()},
                "metric": {sp: torch.as_tensor(m).to(device=device, dtype=dt) for sp, m in self.metric.items()}}
        return self._dev[device]

    def cluster_codes(self, clusters: pd.Series | None) -> np.ndarray:
        """int code per row: substrate cluster (ctx.clusters) or, for proteins without a cluster, the protein."""
        key = id(clusters) if clusters is not None else None
        if key not in self._codes:
            acc = pd.Series(self.acc)
            if clusters is None:
                k = "p:" + acc
            else:
                c = acc.map(clusters)
                k = np.where(c.isna(), "p:" + acc, "c:" + c.astype(str))
            self._codes = {key: pd.factorize(pd.Series(k))[0].astype(np.int64)}   # keep one entry (memory)
        return self._codes[key]


_REPS: dict = {}


def get_reps(spaces=SPACES) -> SiteReps:
    key = tuple(spaces)
    if key not in _REPS:
        _REPS.clear()
        _REPS[key] = SiteReps.from_cache(spaces)
    return _REPS[key]


def load_tau(path: Path = TAU_FILE) -> dict:
    if not Path(path).exists():
        raise FileNotFoundError(f"{path} missing: run scripts/features/knn_tune_tau.py (CV fold 0 inner data)")
    return {k: float(v) for k, v in json.loads(Path(path).read_text())["tau"].items()}


def _buckets(sizes: np.ndarray, ratio: float = 1.25) -> list:
    """Group slot indices (sorted by size) so that padded length <= ratio * smallest member (+ 4)."""
    order = np.argsort(sizes, kind="stable")
    out, cur, lo = [], [], None
    for i in order:
        if cur and sizes[i] > ratio * lo + 4:
            out.append(np.array(cur))
            cur = []
        if not cur:
            lo = sizes[i]
        cur.append(i)
    if cur:
        out.append(np.array(cur))
    return out


@register
class KnnSim(LabelDerivedFeature):
    """Cross-fitted nearest-training-substrate similarity family (see module docstring).
    exclude_same_cluster=False disables the in-transform own-cluster exclusion: TEST/DIAGNOSTIC USE ONLY (the leaky
    variant used by tests/test_knn.py); the driver still raises LeakageError on overlapping labels."""
    name = "knn_sim"

    def __init__(self, spaces=SPACES, tau: dict | float | None = None, c: float = 1.0, n_bg: int = 2000,
                 bg_seed: int = 0, reps: SiteReps | None = None, device: str | None = None,
                 exclude_same_cluster: bool = True, chunk_elems: int = 2 ** 27):
        self.spaces = tuple(spaces)
        self._tau_arg, self.c, self.n_bg, self.bg_seed = tau, float(c), int(n_bg), int(bg_seed)
        self._reps_arg, self.exclude_same_cluster, self.chunk_elems = reps, exclude_same_cluster, int(chunk_elems)
        if device is None:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = device
        self._bg_memo: dict = {}

    # -------------------------------------------------------------------------------------------------- fit
    def _reps(self) -> SiteReps:
        return self._reps_arg if self._reps_arg is not None else get_reps(self.spaces)

    def fit(self, train_pairs: pd.DataFrame, ctx=None) -> "KnnSim":
        reps = self._reps()
        if self._tau_arg is None:
            t = load_tau()
            self.tau = {sp: t[sp] for sp in self.spaces}
        elif isinstance(self._tau_arg, dict):
            self.tau = {sp: float(self._tau_arg[sp]) for sp in self.spaces}
        else:
            self.tau = {sp: float(self._tau_arg) for sp in self.spaces}
        p = train_pairs[["kinase_acc", "site_id"]].drop_duplicates()
        r = reps.pos.get_indexer(p.site_id)
        self.n_pairs_in, self.n_pairs_no_rep = len(p), int((r < 0).sum())
        p = pd.DataFrame({"k": p.kinase_acc.to_numpy()[r >= 0], "r": r[r >= 0]})
        p["rt"] = reps.rt[p.r.to_numpy()]
        self.members = {key: np.unique(g.r.to_numpy()) for key, g in p.groupby(["k", "rt"], sort=False)}
        self.bg_rows = {}
        rng = np.random.default_rng(self.bg_seed)
        for t in (0, 1):
            pool = np.where(reps.rt == t)[0]
            self.bg_rows[t] = np.sort(rng.choice(pool, min(self.n_bg, len(pool)), replace=False)) if len(pool) \
                else pool
        return self

    # -------------------------------------------------------------------------------------------- transform
    def transform(self, query: pd.DataFrame, ctx=None) -> pd.DataFrame:
        reps = self._reps()
        qr = reps.pos.get_indexer(query.site_id)
        kcode, kuniq = pd.factorize(query.kinase_acc)
        valid = qr >= 0
        us, sinv = np.unique(qr[valid], return_inverse=True)
        D = self._compute(us, list(kuniq), ctx)
        out = {}
        for col, M in D.items():
            v = np.full(len(query), np.nan, np.float32)
            v[valid] = M[kcode[valid], sinv]
            out[col] = v
        return pd.DataFrame(out)

    def dense(self, site_ids, kinase_accs, ctx=None) -> dict:
        """{column: (len(kinase_accs), len(site_ids)) float32}; NaN for sites without a representation."""
        reps = self._reps()
        qr = reps.pos.get_indexer(pd.Index(site_ids))
        valid = np.where(qr >= 0)[0]
        D = self._compute(qr[valid], list(kinase_accs), ctx)
        out = {}
        for col, M in D.items():
            full = np.full((len(kinase_accs), len(qr)), np.nan, np.float32)
            full[:, valid] = M
            out[col] = full
        return out

    def _compute(self, q_rows: np.ndarray, kin_list: list, ctx) -> dict:
        import torch
        reps = self._reps()
        dev = self.device
        T = reps.tensors(dev)
        codes = reps.cluster_codes(getattr(ctx, "clusters", None) if ctx is not None else None)
        nK, nq = len(kin_list), len(q_rows)
        out = {col: np.full((nK, nq), np.nan, np.float32) for col in columns(self.spaces)}
        q_rt = reps.rt[q_rows]
        ninf = float("-inf")
        for rt in (0, 1):
            qsel = np.where(q_rt == rt)[0]
            if not len(qsel):
                continue
            mem = [self.members.get((k, rt)) for k in kin_list]
            nk = np.array([0 if m is None else len(m) for m in mem])
            out["knn_log_ntrain"][:, qsel] = np.log1p(nk)[:, None]
            slots = np.where(nk > 0)[0]
            if len(slots):
                trows = np.unique(np.concatenate([mem[i] for i in slots]))
                m = len(trows)
                bks = _buckets(nk[slots])
                bidx = []
                for b in bks:
                    sl = slots[b]
                    L = int(nk[sl].max())
                    I = np.full((len(sl), L), m, np.int64)          # pad -> dummy column m (-inf)
                    for j, s in enumerate(sl):
                        I[j, :nk[s]] = np.searchsorted(trows, mem[s])
                    bidx.append((sl, torch.as_tensor(I, device=dev)))
                padded = sum(I.numel() for _, I in bidx)
            else:
                trows, m, bidx, padded = np.zeros(0, np.int64), 0, [], 0
            brows = self.bg_rows[rt]
            C = int(max(64, min(len(qsel), self.chunk_elems // max(m + len(brows) + padded, 1))))
            ct = torch.as_tensor(codes[trows], device=dev)
            cb = torch.as_tensor(codes[brows], device=dev)
            Tt = {sp: (T["mats"][sp][torch.as_tensor(trows, device=dev)] @ T["metric"][sp]
                       if sp in T["metric"] else T["mats"][sp][torch.as_tensor(trows, device=dev)])
                  for sp in self.spaces}
            Bt = {sp: (T["mats"][sp][torch.as_tensor(brows, device=dev)] @ T["metric"][sp]
                       if sp in T["metric"] else T["mats"][sp][torch.as_tensor(brows, device=dev)])
                  for sp in self.spaces}
            bg_key = (rt, hashlib.blake2b(q_rows[qsel].tobytes(), digest_size=16).hexdigest(),
                      id(getattr(ctx, "clusters", None)))
            bg_done = self._bg_memo.get(bg_key)
            bg_new = {sp: np.full(len(qsel), np.nan, np.float32) for sp in self.spaces}
            for a in range(0, len(qsel), C):
                ch = qsel[a:a + C]
                qr_t = torch.as_tensor(q_rows[ch], device=dev)
                cq = torch.as_tensor(codes[q_rows[ch]], device=dev)
                same_t = (cq[:, None] == ct[None, :]) if self.exclude_same_cluster else None
                same_b = cq[:, None] == cb[None, :]
                for sp in self.spaces:
                    tau = self.tau[sp]
                    Q = T["mats"][sp][qr_t]
                    if bg_done is None:
                        Sb = (Q @ Bt[sp].T).float().masked_fill_(same_b, ninf)
                        nb = (~same_b).sum(1).float()
                        kb = torch.logsumexp(Sb / tau, 1) - torch.log(nb)
                        bg_new[sp][a:a + C] = kb.cpu().numpy()
                    kde_bg = bg_done[sp][a:a + C] if bg_done is not None else bg_new[sp][a:a + C]
                    out[f"knn_{sp}_kde_bg"][:, ch] = kde_bg[None, :]
                    if not len(bidx):
                        continue
                    S = Q @ Tt[sp].T
                    if same_t is not None:
                        S.masked_fill_(same_t, ninf)
                    S = torch.cat([S, torch.full((len(ch), 1), ninf, dtype=S.dtype, device=dev)], 1)
                    for sl, I in bidx:
                        V = S[:, I].float()                                      # (c, G, L)
                        fin = V > ninf
                        ne = fin.sum(-1)
                        none = ne == 0
                        mx = V.amax(-1)
                        top = V.topk(min(3, V.shape[-1]), -1).values
                        tf = top > ninf
                        t3 = torch.where(tf, top, torch.zeros_like(top)).sum(-1) / tf.sum(-1).clamp(min=1)
                        kde = torch.logsumexp(V / tau, -1) - torch.log(ne.float() + self.c)
                        nan = torch.tensor(float("nan"), device=dev)
                        mx, t3, kde = (torch.where(none, nan, x).T.cpu().numpy() for x in (mx, t3, kde))
                        out[f"knn_{sp}_max"][np.ix_(sl, ch)] = mx
                        out[f"knn_{sp}_top3"][np.ix_(sl, ch)] = t3
                        out[f"knn_{sp}_kde"][np.ix_(sl, ch)] = kde
                        out[f"knn_{sp}_lr"][np.ix_(sl, ch)] = kde - kde_bg[None, :]
            if bg_done is None:
                self._bg_memo = {bg_key: bg_new}       # keep only the latest site set (R2 DenseExpand reuses it)
        # kde_bg is site-level but reported only where the kinase has training sites (consistent NaN pattern)
        for sp in self.spaces:
            out[f"knn_{sp}_kde_bg"][np.isnan(out[f"knn_{sp}_kde"])] = np.nan
        return out
