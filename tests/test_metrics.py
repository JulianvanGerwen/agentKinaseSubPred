"""Metric functions vs hand-computed toy examples (DESIGN §11.4) and sklearn cross-checks."""
import numpy as np
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from ksp import eval as E


def test_auroc_hand():
    # pos {3, 1}, neg {2, 0}: pairs (3>2,3>0,1<2,1>0) -> 3/4
    assert E.auroc([3, 1], [2, 0]) == pytest.approx(0.75)
    # tie counts 1/2: pos {1}, neg {1, 0} -> (0.5 + 1)/2
    assert E.auroc([1], [1, 0]) == pytest.approx(0.75)
    assert np.isnan(E.auroc([], [1]))
    assert E.auroc([5, 6], [1, 2, 3]) == 1.0


def test_pauroc_hand_and_sklearn():
    # perfect -> 1, reversed -> standardized minimum
    assert E.pauroc([2, 3], [0, 1], 0.5) == pytest.approx(1.0)
    rng = np.random.default_rng(1)
    for _ in range(5):
        p, n = rng.normal(1, 1, 40), rng.normal(0, 1, 400)
        y = np.r_[np.ones(40), np.zeros(400)]
        s = np.r_[p, n]
        for f in (0.01, 0.1, 0.3):
            assert E.pauroc(p, n, f) == pytest.approx(roc_auc_score(y, s, max_fpr=f), abs=1e-9)
        assert E.auroc(p, n) == pytest.approx(roc_auc_score(y, s))


def test_average_precision_hand_and_sklearn():
    # ranking: P N P N -> AP = (1/1 + 2/3)/2
    assert E.average_precision([4, 2], [3, 1]) == pytest.approx((1 + 2 / 3) / 2)
    rng = np.random.default_rng(2)
    p, n = np.round(rng.normal(1, 1, 30), 1), np.round(rng.normal(0, 1, 300), 1)   # rounding -> ties
    y = np.r_[np.ones(30), np.zeros(300)]
    assert E.average_precision(p, n) == pytest.approx(average_precision_score(y, np.r_[p, n]))


def test_topfrac_recall_hand():
    # N = 200 -> top 2. pos scores 10, 5, 0.5 ; negatives 199 with one at 7 -> top2 = {10, 7}
    neg = np.r_[7, np.zeros(196)]
    assert E.topfrac_recall([10, 5, 0.5], neg, 0.01) == pytest.approx(1 / 3)
    # tie at the cut-off: N=100 -> top 1; pos 1 tied with 1 neg at the max -> credit 1/2
    assert E.topfrac_recall([1.0], np.r_[1.0, np.zeros(98)], 0.01) == pytest.approx(0.5)


def test_expected_hit_at():
    assert E.expected_hit_at(0, 0, 1) == 1
    assert E.expected_hit_at(4, 0, 5) == 1
    assert E.expected_hit_at(5, 0, 5) == 0
    # 3 greater, tied with 3 others -> rank uniform on {4,5,6,7}: P(rank<=5) = 2/4
    assert E.expected_hit_at(3, 3, 5) == pytest.approx(0.5)


def test_weighted_auroc_matrix_matches_unweighted_and_replication():
    rng = np.random.default_rng(3)
    p, n = np.round(rng.normal(1, 1, 20), 1), np.round(rng.normal(0, 1, 50), 1)
    Wp, Wn = np.ones((1, 20)), np.ones((1, 50))
    assert E.weighted_auroc_matrix(p, n, Wp, Wn)[0] == pytest.approx(E.auroc(p, n))
    # integer weights == replicating items
    wp, wn = rng.integers(0, 3, 20), rng.integers(0, 3, 50)
    rep = E.auroc(np.repeat(p, wp), np.repeat(n, wn))
    assert E.weighted_auroc_matrix(p, n, wp[None].astype(float), wn[None].astype(float))[0] == pytest.approx(rep)


def test_holm():
    adj = E.holm([0.01, 0.04, 0.03])
    assert adj == pytest.approx([0.03, 0.06, 0.06])
