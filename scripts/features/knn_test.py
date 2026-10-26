"""Unit tests for the knn_sim crossfit family (ksp.features.knn). Run: python -m pytest -q scripts/features/knn_test.py
(kept under scripts/features/ = feat-knn's paths; pytest collects *_test.py).

Planted toy: proteins come in homologous pairs (one substrate cluster per pair); site j of the two proteins of a
pair has (near-)identical embeddings and the same label, so the label "leaks" through (a) the site itself and (b)
its homologue copy. Across clusters there is no signal.
  * leaky variant (fit on ALL labels, no in-transform exclusion, no driver)  -> AUROC ~ 1
  * transform-robust variant (fit on ALL labels, own-cluster exclusion ON)   -> AUROC ~ 0.5
  * through the ksp.crossfit driver                                          -> AUROC ~ 0.5
Positive control: a planted direction shared by positives across clusters is recovered under the driver.
Correctness: GPU/CPU engine == brute-force numpy for max / top3 / kde / kde_bg / lr with exclusions.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import crossfit  # noqa: E402
from ksp.eval import auroc  # noqa: E402
from ksp.features import knn  # noqa: E402

TAU = 0.1


def _toy(n_pairs=150, sites_per=4, frac=0.3, d=32, signal=0.0, seed=0):
    rng = np.random.default_rng(seed)
    prots = [f"P{i:05d}" for i in range(2 * n_pairs)]
    clusters = pd.Series([i // 2 for i in range(2 * n_pairs)], index=prots)
    base = rng.standard_normal((n_pairs, sites_per, d))
    lab_pair = rng.random((n_pairs, sites_per)) < frac
    motif = rng.standard_normal(d)
    sids, X, y = [], [], []
    for i, p in enumerate(prots):
        for j in range(sites_per):
            v = base[i // 2, j] + 0.05 * rng.standard_normal(d)       # homologue copies ~ identical
            if lab_pair[i // 2, j]:
                v = v + signal * motif
            sids.append(f"{p}_S{j + 1}")
            X.append(v / np.linalg.norm(v))
            y.append(int(lab_pair[i // 2, j]))
    X = np.asarray(X, np.float32)
    reps = knn.SiteReps(sids, np.zeros(len(sids), np.int8), {"esmc": X})
    rows = pd.DataFrame({"kinase_acc": "K1", "site_id": sids, "y": y})
    return rows, rows[rows.y == 1][["kinase_acc", "site_id"]].reset_index(drop=True), clusters, reps


def _fam(reps, **kw):
    return lambda: knn.KnnSim(spaces=("esmc",), tau=TAU, reps=reps, device="cpu", n_bg=200, **kw)


def _auc(f, rows, col="knn_esmc_max"):
    v = f[col].to_numpy()
    return auroc(v[rows.y.values == 1], v[rows.y.values == 0])


def test_leaky_variant_is_perfect():
    rows, labs, cl, reps = _toy()
    ctx = crossfit.FitContext(cl)
    fam = _fam(reps, exclude_same_cluster=False)().fit(labs, ctx)
    assert _auc(fam.transform(rows[["kinase_acc", "site_id"]], ctx), rows) > 0.99


def test_transform_excludes_own_cluster_even_without_driver():
    rows, labs, cl, reps = _toy()
    ctx = crossfit.FitContext(cl)
    fam = _fam(reps)().fit(labs, ctx)
    f = fam.transform(rows[["kinase_acc", "site_id"]], ctx)
    for col in ("knn_esmc_max", "knn_esmc_top3", "knn_esmc_kde", "knn_esmc_lr"):
        assert abs(_auc(f, rows, col) - 0.5) < 0.07, col


def test_driver_gives_chance():
    rows, labs, cl, reps = _toy(seed=1)
    tr, _ = crossfit.crossfit(_fam(reps), labs, rows[["kinase_acc", "site_id"]], None, clusters=cl)
    for col in ("knn_esmc_max", "knn_esmc_top3", "knn_esmc_kde", "knn_esmc_lr"):
        assert abs(_auc(tr, rows, col) - 0.5) < 0.07, col


def test_driver_raises_on_overlap():
    rows, labs, cl, reps = _toy()
    with pytest.raises(crossfit.LeakageError):
        crossfit.crossfit(_fam(reps), labs, rows, rows, clusters=cl, eval_mode="all")


def test_positive_control_signal_survives_driver():
    rows, labs, cl, reps = _toy(signal=0.8, seed=2)
    tr, _ = crossfit.crossfit(_fam(reps), labs, rows[["kinase_acc", "site_id"]], None, clusters=cl)
    assert _auc(tr, rows, "knn_esmc_kde") > 0.8
    assert _auc(tr, rows, "knn_esmc_lr") > 0.8


def _brute(reps, labs, q_sids, cl, kin, tau, c, bg_rows):
    X = reps.mats["esmc"].astype(np.float64)
    pos = reps.pos
    clus = pd.Series(reps.acc).map(cl).to_numpy()
    T = pos.get_indexer(labs[labs.kinase_acc == kin].site_id.unique())
    res = {k: [] for k in ("max", "top3", "kde", "kde_bg", "lr")}
    for s in q_sids:
        q = pos.get_loc(s)
        t = T[clus[T] != clus[q]]
        b = bg_rows[clus[bg_rows] != clus[q]]
        sb = X[b] @ X[q]
        kb = np.log(np.exp(sb / tau).mean())
        if len(t) == 0:
            for k in res:
                res[k].append(np.nan)
            continue
        sim = np.sort(X[t] @ X[q])[::-1]
        kde = np.log(np.exp(sim / tau).sum() / (len(t) + c))
        for k, v in zip(res, (sim[0], sim[:3].mean(), kde, kb, kde - kb)):
            res[k].append(v)
    return {k: np.array(v) for k, v in res.items()}


@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_engine_matches_bruteforce(device):
    import torch
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("no GPU")
    rows, labs, cl, reps = _toy(n_pairs=60, seed=3)
    labs2 = pd.concat([labs, pd.DataFrame({"kinase_acc": "K2", "site_id": labs.site_id.iloc[:2]})])
    ctx = crossfit.FitContext(cl)
    fam = knn.KnnSim(spaces=("esmc",), tau=TAU, reps=reps, device=device, n_bg=50, chunk_elems=3000).fit(labs2, ctx)
    q = rows.site_id.iloc[::3].to_numpy()
    D = fam.dense(q, ["K1", "K2", "K_none"], ctx)
    tol = 2e-3 if device == "cpu" else 3e-2           # fp16 matmul on GPU
    for ki, kin in enumerate(["K1", "K2"]):
        B = _brute(reps, labs2, q, cl, kin, TAU, 1.0, fam.bg_rows[0])
        for st in ("max", "top3", "kde", "kde_bg", "lr"):
            got = D[f"knn_esmc_{st}"][ki]
            assert np.array_equal(np.isnan(got), np.isnan(B[st])), (kin, st)
            m = ~np.isnan(got)
            assert np.allclose(got[m], B[st][m], atol=tol * (1 / TAU if st != "max" and st != "top3" else 1)), (kin, st)
    assert np.isnan(D["knn_esmc_max"][2]).all() and (D["knn_log_ntrain"][2] == 0).all()
    # pair transform == dense
    qq = pd.DataFrame({"kinase_acc": np.repeat(["K1", "K2"], len(q)), "site_id": np.tile(q, 2)})
    f = fam.transform(qq, ctx)
    assert np.allclose(f.knn_esmc_kde.to_numpy(), np.concatenate([D["knn_esmc_kde"][0], D["knn_esmc_kde"][1]]),
                       equal_nan=True)


def test_registered():
    assert crossfit.REGISTRY["knn_sim"] is knn.KnnSim


def test_blosum_metric_matches_direct_sum():
    from Bio.Align import substitution_matrices
    B = substitution_matrices.load("BLOSUM62")
    w = knn.blosum_weights()
    a, b = "RRKSPSLSGEYAX__", "RRASVSDAGEYQKLP"
    O = knn.onehot_windows([a, b]).astype(np.float32)
    got = O[0] @ knn.blosum_metric() @ O[1]
    std = set(knn.AA20)
    exp = sum(w[i] * B[x][y] for i, (x, y) in enumerate(zip(a, b)) if x in std and y in std) / w.sum()
    assert abs(got - exp) < 1e-5
