"""Frozen split invariants (DESIGN §4; §11.4): no protein / cluster / identical-window overlap between TEST and DEV or
across CV folds; constraints met. Reads only public split tables + DEV labels."""
import json

import numpy as np
import pandas as pd
import pytest

from ksp import index, labels, splits

pytestmark = pytest.mark.skipif(not (splits.SPLIT_DIR / "SEAL.json").exists(), reason="needs project data that is not part of this repository (frozen splits under KSP_DATA_DIR); see DATA_SOURCES.md")


@pytest.fixture(scope="module")
def tabs():
    cl = pd.read_csv(splits.SPLIT_DIR / "substrate_clusters.tsv", sep="\t")
    test = splits.load_test_proteins()
    folds = splits.load_cv_folds()
    return cl, test, folds


def test_every_protein_assigned_once(tabs):
    cl, test, folds = tabs
    P = set(index.load_proteins().acc)
    assert set(cl.acc) == P and cl.acc.is_unique
    assert not (test & set(folds.index))
    assert test | set(folds.index) == P
    assert set(folds.unique()) == {0, 1, 2, 3, 4}


def test_no_cluster_overlap(tabs):
    cl, test, folds = tabs
    c = cl.set_index("acc").cluster_id
    tc = set(c[list(test)])
    dc = set(c[list(folds.index)])
    assert not tc & dc
    fc = pd.DataFrame({"c": c[folds.index].values, "f": folds.values}).groupby("c").f.nunique()
    assert (fc == 1).all(), "a substrate cluster is split across CV folds"


def test_no_identical_window_overlap(tabs):
    cl, test, folds = tabs
    s = index.load_sites()
    grp = np.where(s.acc.isin(test), "T", s.acc.map(folds).astype("Int64").astype(str))
    g = pd.DataFrame({"w": s.window15.values, "g": grp}).groupby("w").g.nunique()
    assert (g == 1).all(), f"{(g > 1).sum()} identical 15-mers shared between TEST/DEV or CV folds"


def test_mmseqs_edges_within_clusters(tabs):
    cl, _, _ = tabs
    m = splits.read_m8(splits.CACHE / "allvsall.m8")
    m = m[(m.fident >= 0.4) & (np.minimum(m.qcov, m.tcov) >= 0.5)]
    c = cl.set_index("acc").cluster_id
    assert (c[m["query"]].values == c[m["target"]].values).all()


def test_split_constraints():
    R = json.loads((splits.SPLIT_DIR / "split_numbers.json").read_text())
    assert 0.19 <= R["frac_test:L2_labels"] <= 0.21
    assert 0.19 <= R["frac_test:U_ms_sites"] <= 0.21
    for r in "STY":
        assert 0.18 <= R[f"frac_test:{r}_labels"] <= 0.22
    for f in range(5):
        assert 0.19 <= R[f"frac_cv{f}:L2_labels"] <= 0.21
    assert 0.15 <= R["clusters_frac_test"] <= 0.25
    # soft per-kinase band (PI r2): >= 75% of kinases with >=5 eval positives in [15%, 25%], >= 95% in count band
    assert R["kinases_ge5_in_15_25pct"] >= 0.75 * R["kinases_ge5_eval_pos"]
    assert R["kinases_ge5_in_count_band(floor15%,ceil25%)"] >= 0.95 * R["kinases_ge5_eval_pos"]
    assert R["clusters:max_component_label_frac_after"] <= 0.03


def test_dev_table_excludes_test(tabs):
    _, test, _ = tabs
    d = labels.load_dev("all")
    assert not set(d.acc) & test
    c = labels.load_dev("L2", with_canaries=True)
    assert not set(splits.site_protein(c.site_id)) & test


def test_kinase_folds_no_family_split():
    kf = splits.load_kinase_folds()
    assert set(kf.lko_fold) == {0, 1, 2, 3, 4}
    fam = kf.dropna(subset=["manning_family"])
    fam = fam[~fam.manning_family.str.endswith("Unique")]
    assert (fam.groupby("manning_family").lko_fold.nunique() == 1).all()
    assert (kf.groupby("kinase_cluster").lko_fold.nunique() == 1).all()


def test_eval_kinase_sets_definitions():
    es = splits.load_eval_kinase_sets()
    k = es[es.in_K_eval]
    assert (k.n_test_pos >= 3).all() and (k.n_test_clusters >= 2).all()
    assert (es[es.in_K_any].n_test_pos >= 1).all()
    k = es[es.in_K_cv]
    assert (k.n_dev_pos >= 5).all() and (k.n_dev_clusters >= 2).all()
    assert es.in_K_eval.sum() <= es.in_K_any.sum()


def test_inner_folds_group_by_cluster():
    d = labels.load_dev("L2")
    f = splits.inner_folds(d.site_id, k=5, seed=0)
    c = splits.site_cluster(d.site_id)
    assert (pd.DataFrame({"c": c, "f": f}).groupby("c").f.nunique() == 1).all()
    share = np.bincount(f, minlength=5) / len(f)
    assert share.min() > 0.12
    assert np.array_equal(f, splits.inner_folds(d.site_id, k=5, seed=0))   # deterministic
    v = splits.train_valid_split_inner(d.site_id, 0.2, seed=1)
    assert 0.1 < v.mean() < 0.3
