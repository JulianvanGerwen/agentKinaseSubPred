"""Build the label-free site representations used by the knn_sim family (ksp.features.knn).

Source: R3's ESM-2 650M per-residue window embeddings from FULL-protein context (cache/r3/win_emb.npy, rows via
cache/r3/win_rows.npy, +-10 window, -1 = beyond the protein end; site order = cache/r3/win_sites.site_ids.npy;
docs/model_r3.md §3). Written once; no labels involved.

Outputs (cache/features/knn/):
  reps_site_ids.npy        (n,) site_id of every row (U_ms sites; matched to cache/index/sites.parquet BY ID)
  reps_esm_centre.npy      (n, 1280) float16  centre-residue embedding (raw; centring/L2 norm done at load)
  reps_esm_win7.npy        (n, 1280) float16  mean over the +-7 window positions inside the protein (raw)
  reps_meta.json           provenance + timings
Usage: python -I scripts/features/knn_build_reps.py
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
    main()
