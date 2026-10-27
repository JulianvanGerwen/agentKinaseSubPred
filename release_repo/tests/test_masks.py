"""Negative masks (DESIGN §5): Hamming search exactness, KP sets, paralog sets, NegativeMasks semantics."""
import numpy as np
import pandas as pd
import pytest

from ksp import index, labels


def _brute(q, ref, d):
    R = np.array([list(x) for x in ref])
    return [np.where((R != np.array(list(w))).sum(1) <= d)[0] for w in q]


def test_hamming_neighbors_exact_vs_bruteforce():
    rng = np.random.default_rng(0)
    alpha = list("ACDEFS")
    ref = ["".join(rng.choice(alpha, 15)) for _ in range(800)]
    # plant near-copies
    q = []
    for i in range(60):
        w = list(ref[rng.integers(len(ref))])
        for p in rng.choice(15, rng.integers(0, 4), replace=False):
            w[p] = "Y"
        q.append("".join(w))
    got = labels.hamming_neighbors(q, ref, 2)
    exp = _brute(q, ref, 2)
    for g, e in zip(got, exp):
        assert np.array_equal(np.sort(g), np.sort(e))


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet", "cache/labels/dev_pairs.parquet")
def test_paralog_examples():
    k = index.load_kinases()
    a = dict(zip(k.gene, k.uniprot_acc))
    g = dict(zip(k.uniprot_acc, k.gene))
    par = labels.paralog_sets()
    src = {g[x] for x in par[a["SRC"]]}
    assert {"YES1", "FYN", "FGR", "LCK", "HCK", "BLK", "LYN"} <= src          # DESIGN §5 example
    assert {g[x] for x in par[a["MAPK1"]]} >= {"MAPK3"}
    assert "PDK2" in {g[x] for x in par[a["PDK1"]]}
    for x, ps in par.items():                                              # symmetric, no self
        assert x not in ps
        for y in ps:
            assert x in par[y]


@pytest.fixture(scope="module")
def masks():
    d = labels.load_dev("all")
    return d, labels.NegativeMasks(d)


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet", "cache/labels/dev_pairs.parquet")
def test_kp_and_masks_semantics(masks):
    d, m = masks
    s = index.load_sites()
    spos = pd.Series(np.arange(len(s)), index=s.site_id)
    k = index.load_kinases()
    acc = dict(zip(k.gene, k.uniprot_acc))["CSNK2A1"]
    kp = labels.known_positive_sets(d)[acc]
    ex = set(m.excluded(acc))
    in_idx = [spos[x] for x in kp if x in spos.index]
    assert set(in_idx) <= ex                                                   # KP never negative
    # every site within Hamming <= 2 of a KP site window is excluded
    W = s.window15.values
    kp_w = labels.site_windows(sorted(kp))
    nb = labels.hamming_neighbors(list(kp_w[:50]), list(W), 2)
    for arr in nb:
        assert set(arr) <= ex
    # sites labelled to a paralog are excluded
    for j in m.par[acc]:
        for x in m.kp.get(j, ()):
            if x in spos.index:
                assert spos[x] in ex
    # E-hard candidates carry a human-native label to a non-paralog kinase
    hc = m.hard_candidates(acc)
    own = {acc} | m.par[acc]
    lab = m.hard_labelers
    sid = s.site_id.values[hc[:200]]
    for x in sid:
        assert len(lab[x] - own) > 0


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet", "cache/labels/dev_pairs.parquet")
def test_identical_window_expansion():
    toy = pd.DataFrame({"kinase_acc": ["K"], "site_id": [index.load_sites().site_id.iloc[0]]})
    s = index.load_sites()
    w = s.window15.iloc[0]
    kp = labels.known_positive_sets(toy)["K"]
    assert set(s.site_id[s.window15 == w]) <= kp


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet", "cache/labels/dev_pairs.parquet")
def test_load_dev_variants():
    d = labels.load_dev("L2")
    assert d.in_U_ms.all() and d.L2.all()
    assert not d.duplicated(["kinase_acc", "site_id"]).any()
    e = labels.load_dev("eval")
    assert (e.provenance == "PSP").all() and e.L1.all() and not e.autophos.any()
    with pytest.raises(ValueError):
        labels.load_dev("L9")
