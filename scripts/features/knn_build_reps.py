"""Build the label-free site representations used by the knn_sim family (ksp.features.knn).

Source: R3's ESM-2 650M per-residue window embeddings from FULL-protein context (cache/r3/win_emb.npy, rows via
cache/r3/win_rows.npy, +-10 window, -1 = beyond the protein end; site order = cache/r3/win_sites.site_ids.npy;
docs/model_r3.md §3). Written once; no labels involved.

Outputs (cache/features/knn/):
  reps_site_ids.npy        (n,) site_id of every row (U_ms sites; matched to cache/index/sites.parquet BY ID)
  reps_esm_centre.npy      (n, 1280) float16  centre-residue embedding (raw; centring/L2 norm done at load)
  reps_esm_win7.npy        (n, 1280) float16  mean over the +-7 window positions inside the protein (raw)
  reps_esm_cat7.npy        (n, 15*128) float16  (--cat) concatenated +-7 window: each position's embedding is
                           projected on a global residue PCA-128 (fit on 200k random window residues, label-free),
                           L2-normalised and scaled by sqrt(w_p / sum w), w_p = 1/(1+|p|/4); padding -> 0. Dot product
                           = centre-weighted mean of per-position cosines (PCA space).
  reps_meta.json           provenance + timings
Usage: python -I scripts/features/knn_build_reps.py [--cat]
"""
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "6")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np  # noqa: E402

from ksp import index  # noqa: E402

R3 = ROOT / "cache" / "r3"
OUT = ROOT / "cache" / "features" / "knn"
FLANK_SRC, FLANK = 10, 7


def build_cat(n_pca=128, n_fit=200_000, seed=0):
    import torch
    t0 = time.time()
    rows = np.load(R3 / "win_rows.npy")
    emb = np.load(R3 / "win_emb.npy", mmap_mode="r")
    sel = rows[:, FLANK_SRC - FLANK: FLANK_SRC + FLANK + 1]
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rng = np.random.default_rng(seed)
    allr = np.unique(sel[sel >= 0])
    fit_r = np.sort(rng.choice(allr, min(n_fit, len(allr)), replace=False))
    X = torch.as_tensor(np.asarray(emb[fit_r], dtype=np.float32), device=dev)
    mu = X.mean(0)
    Xc = X - mu
    cov = Xc.T @ Xc / (len(Xc) - 1)
    ev, V = torch.linalg.eigh(cov.double())
    V = V[:, torch.argsort(ev, descending=True)[:n_pca]].float()
    expl = float(ev.sort(descending=True).values[:n_pca].sum() / ev.sum())
    del X, Xc
    p = np.arange(-FLANK, FLANK + 1)
    w = 1.0 / (1.0 + np.abs(p) / 4.0)
    sw = torch.as_tensor(np.sqrt(w / w.sum()), device=dev, dtype=torch.float32)
    n = len(sel)
    out = np.zeros((n, len(p) * n_pca), np.float16)
    CH = 20000
    for a in range(0, n, CH):
        b = min(n, a + CH)
        r = sel[a:b]
        u, inv = np.unique(r[r >= 0], return_inverse=True)
        E = torch.as_tensor(np.asarray(emb[u], dtype=np.float32), device=dev)
        Z = (E - mu) @ V
        Z = Z / Z.norm(dim=1, keepdim=True).clamp(min=1e-6)
        idx = np.full(r.shape, -1, np.int64)
        idx[r >= 0] = inv
        it = torch.as_tensor(np.where(idx >= 0, idx, 0), device=dev)
        G = Z[it] * torch.as_tensor(idx >= 0, device=dev)[..., None] * sw[None, :, None]   # (m, 15, n_pca)
        out[a:b] = G.reshape(b - a, -1).half().cpu().numpy()
        print(f"cat {b}/{n} {time.time() - t0:.0f}s", flush=True)
    np.save(OUT / "reps_esm_cat7.npy", out)
    mp = OUT / "reps_meta.json"
    meta = json.loads(mp.read_text()) if mp.exists() else {}
    meta["cat7"] = {"n_pca": n_pca, "pca_fit_residues": int(len(fit_r)), "pca_explained_var": round(expl, 4),
                    "weights": "sqrt(w_p/sum w), w_p = 1/(1+|p|/4)", "runtime_s": round(time.time() - t0, 1)}
    mp.write_text(json.dumps(meta, indent=1))
    print(meta["cat7"])


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = np.load(R3 / "win_rows.npy")
    sids = np.load(R3 / "win_sites.site_ids.npy", allow_pickle=True).astype(str)
    emb = np.load(R3 / "win_emb.npy", mmap_mode="r")
    assert rows.shape == (len(sids), 2 * FLANK_SRC + 1)
    S = index.load_sites()
    known = set(S.site_id)
    assert all(s in known for s in sids), "R3 window sites not in canonical index"
    n, d = len(sids), emb.shape[1]
    cen = np.zeros((n, d), np.float16)
    win = np.zeros((n, d), np.float16)
    sel = rows[:, FLANK_SRC - FLANK: FLANK_SRC + FLANK + 1]          # (n, 15)
    assert (sel[:, FLANK] >= 0).all()
    CH = 20000
    for a in range(0, n, CH):
        b = min(n, a + CH)
        r = sel[a:b]
        u, inv = np.unique(r[r >= 0], return_inverse=True)
        E = np.asarray(emb[u], dtype=np.float32)                     # contiguous-ish rows (sites sorted by protein)
        idx = np.full(r.shape, -1, np.int64)
        idx[r >= 0] = inv
        valid = idx >= 0
        Eg = E[np.where(valid, idx, 0)] * valid[..., None]           # (m, 15, d)
        win[a:b] = (Eg.sum(1) / valid.sum(1, keepdims=True)).astype(np.float16)
        cen[a:b] = E[idx[:, FLANK]].astype(np.float16)
        print(f"{b}/{n} {time.time() - t0:.0f}s", flush=True)
    np.save(OUT / "reps_site_ids.npy", sids)
    np.save(OUT / "reps_esm_centre.npy", cen)
    np.save(OUT / "reps_esm_win7.npy", win)
    meta = {"source": "cache/r3/win_emb.npy (ESM-2 650M facebook/esm2_t33_650M_UR50D, full-protein context)",
            "n_sites": int(n), "dim": int(d), "window_flank": FLANK, "pooling": "mean over in-protein positions",
            "runtime_s": round(time.time() - t0, 1)}
    (OUT / "reps_meta.json").write_text(json.dumps(meta, indent=1))
    print(meta)


if __name__ == "__main__":
    if "--cat" in sys.argv:
        build_cat()
    else:
        main()
