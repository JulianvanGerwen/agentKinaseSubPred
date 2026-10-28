"""Toy-data tests for src/ksp/features/{structure,esm_embed}.py.
Run: python tests/test_features_structure_esm.py"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp.features import structure as S  # noqa: E402
from ksp.features.esm_embed import windows  # noqa: E402


def ref_neighbors(ca, cb, c, n, pae, max_dist, max_angle):
    """Literal re-statement of structuremap 0.0.10 get_neighbors/get_angle/get_gly_vector loops."""
    L = len(ca)
    res = []
    for i in range(L):
        if np.isnan(cb[i, 0]):
            un = (n[i] - ca[i]) / np.linalg.norm(n[i] - ca[i])
            uc = (c[i] - ca[i]) / np.linalg.norm(c[i] - ca[i])
            u = S._rotation(uc[None], -120.0)[0] @ un
        else:
            u = (cb[i] - ca[i]) / np.linalg.norm(cb[i] - ca[i])
        k = 0
        for j in range(L):
            if i == j or pae[i, j] > max_dist:
                continue
            d = np.linalg.norm(ca[i] - ca[j])
            if d + pae[i, j] <= max_dist:
                ang = np.degrees(np.arccos(np.dot(u, (ca[j] - ca[i]) / d)))
                if ang <= max_angle:
                    k += 1
        res.append(k)
    return np.array(res)


def test_pse_matches_reference():
    rng = np.random.default_rng(0)
    L = 60
    ca = np.cumsum(rng.normal(0, 2.2, (L, 3)), 0)
    n = ca + rng.normal(0, 1, (L, 3))
    c = ca + rng.normal(0, 1, (L, 3))
    cb = ca + rng.normal(0, 1, (L, 3))
    cb[::7] = np.nan  # glycines
    pae = rng.integers(0, 15, (L, L)).astype(float)
    u = S.side_chain_vectors(ca, cb, c, n)
    for p in (None, pae):
        got = S.pse_counts(ca, u, p)
        for (r, th), v in got.items():
            exp = ref_neighbors(ca, cb, c, n, np.zeros((L, L)) if p is None else p, r, th)
            assert (v == exp).all(), (r, th, p is None)


def test_smooth():
    x = np.array([1.0, 2, 3, 4, 5])
    assert np.allclose(S.smooth(x, 1), [1.5, 2, 3, 4, 4.5])
    assert np.allclose(S.smooth(x, 10), [3] * 5)


def test_sasa_isolated_atom():
    a = S.sasa_atoms(np.zeros((1, 3)), np.array([1.7]))
    assert np.isclose(a[0], 4 * np.pi * (1.7 + 1.4) ** 2)
    two = S.sasa_atoms(np.array([[0, 0, 0], [100.0, 0, 0]]), np.array([1.7, 1.7]))
    assert np.allclose(two, a[0])
    buried = S.sasa_atoms(np.array([[0, 0, 0], [0.1, 0, 0]]), np.array([1.7, 1.7]))
    assert (buried < a[0]).all()


def test_windows():
    assert windows(100) == [(0, 100)]
    w = windows(2000, 1022, 511)
    assert w[0] == (0, 1022) and w[-1][1] == 2000 and all(b - a == 1022 for a, b in w)
    cov = np.zeros(2000, int)
    for a, b in w:
        cov[a:b] += 1
    assert (cov >= 1).all()


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f()
            print("ok", k)
