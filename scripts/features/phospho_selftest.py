#!/usr/bin/env python
"""Self-tests for the F6 kernels (synthetic data only; no labels).
  1. masked_corr == pandas pairwise-complete Pearson (min_periods) on random data with NaN.
  2. coreg_features excludes (a) the query itself, (b) substrates on the same protein, (c) substrates in the same
     cluster (dummy cluster map), and requires >= min_blocks blocks.
  3. crossfit contract: CoregF6a.fit/transform returns one row per query pair in order.
Usage: CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_selftest.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp.features import phospho as ph  # noqa: E402


def test_masked_corr(dev):
    import torch
    rng = np.random.default_rng(0)
    A = rng.normal(size=(30, 40)).astype(np.float32)
    B = (A[:20] * 0.5 + rng.normal(size=(20, 40))).astype(np.float32)
    A[rng.random(A.shape) < 0.3] = np.nan
    B[rng.random(B.shape) < 0.3] = np.nan
    r, n = ph.masked_corr(ph.to_dev(A, dev), ph.to_dev(B, dev), min_n=8)
    r = r.cpu().numpy()
    ref = np.full_like(r, np.nan)
    for i in range(A.shape[0]):
        for j in range(B.shape[0]):
            ref[i, j] = pd.Series(A[i]).corr(pd.Series(B[j]), min_periods=8)
    ok = np.isfinite(ref)
    assert np.array_equal(ok, np.isfinite(r)), "NaN pattern differs"
    assert np.nanmax(np.abs(r - ref)) < 1e-4, np.nanmax(np.abs(r - ref))
    print("masked_corr OK, max diff", float(np.nanmax(np.abs(r - ref))))


def test_exclusions(dev):
    rng = np.random.default_rng(1)
    n_sites = len(ph.site_index())
    base = rng.normal(size=50).astype(np.float32)
    # sites 0 (query), 1 (same protein as 0), 2 (same cluster, other protein), 3 (unrelated protein/cluster),
    # 4 (noise substrate)
    X = np.stack([base, base + 0.01 * rng.normal(size=50), base + 0.01 * rng.normal(size=50),
                  base * 0.8 + 0.6 * rng.normal(size=50), rng.normal(size=50)]).astype(np.float32)
    rows = np.array([0, 1, 2, 3, 4], np.int32)
    blocks = [ph.Block(f"b{i}", "compendium", rows, X[:, i * 25:(i + 1) * 25].copy(), [f"c{j}" for j in range(25)],
                       (i,)) for i in range(2)]
    prot = np.arange(n_sites, dtype=np.int32)
    prot[1] = prot[0]
    clus = np.arange(n_sites, dtype=np.int64)
    clus[2] = clus[0]
    clus[1] = clus[0]
    tp = pd.DataFrame(dict(kinase_idx=[0, 0, 0, 0, 0, 1], site_idx=[0, 1, 2, 3, 4, 0]))
    r_nocl = ph.coreg_features(tp, np.array([0]), blocks, protein_of_site_arr=prot, n_kinases=3, device=dev, m=1)
    r_cl = ph.coreg_features(tp, np.array([0]), blocks, protein_of_site_arr=prot, cluster_of_site=clus,
                             n_kinases=3, device=dev, m=1)
    # without cluster map: best substrate = site 2 (same cluster, r ~ 1); with: site 3 (r ~ 0.8)
    top_nocl, top_cl = r_nocl["top"][0, 0], r_cl["top"][0, 0]
    assert np.tanh(top_nocl) > 0.99, top_nocl
    assert 0.3 < np.tanh(top_cl) < 0.95, top_cl
    # kinase 1 has only the query itself as substrate -> NaN; kinase 2 has none -> NaN
    assert np.isnan(r_cl["top"][1, 0]) and np.isnan(r_cl["top"][2, 0])
    # min_blocks: with one block only -> NaN
    r1 = ph.coreg_features(tp, np.array([0]), blocks[:1], protein_of_site_arr=prot, cluster_of_site=clus,
                           n_kinases=3, device=dev)
    assert np.isnan(r1["top"][0, 0])
    # exclude_pmids drops block 0 -> only one block -> NaN
    r2 = ph.coreg_features(tp, np.array([0]), blocks, protein_of_site_arr=prot, cluster_of_site=clus,
                           n_kinases=3, device=dev, exclude_pmids=[0])
    assert np.isnan(r2["top"][0, 0])
    print(f"exclusions OK: r(no cluster map)={np.tanh(top_nocl):.3f}, r(cluster excluded)={np.tanh(top_cl):.3f}")


def test_contract(dev):
    sites, kin = ph.site_index(), ph.kinase_index()
    blocks = ph.cached_blocks("compendium")
    big = max(blocks, key=lambda b: b.n_cols)
    s = sites.site_id.to_numpy()[big.rows[:40]]
    train = pd.DataFrame(dict(kinase_acc=[kin.uniprot_acc[0]] * 20 + [kin.uniprot_acc[1]] * 20, site_id=s))
    q = pd.DataFrame(dict(kinase_acc=[kin.uniprot_acc[1], kin.uniprot_acc[0], kin.uniprot_acc[0]],
                          site_id=sites.site_id.to_numpy()[big.rows[[100, 100, 101]]]))
    fam = ph.make_f6a("f6a_selftest", kinds=("compendium",), device=dev)().fit(train)
    out = fam.transform(q)
    assert len(out) == 3 and list(out.columns) == ["f6a_comp_top", "f6a_comp_med"]
    print("contract OK\n", out)


if __name__ == "__main__":
    dev = "cuda"
    test_masked_corr(dev)
    test_exclusions(dev)
    test_contract(dev)
    print("ALL OK")
