"""Canaries (DESIGN §11.3, §16 A5)."""
import numpy as np
import pandas as pd
import pytest

from ksp import canary, index, labels, splits
from ksp import eval as E


def test_canary_kinases_and_descriptors():
    k = canary.kinases_with_canaries()
    assert len(k) == len(index.load_kinases()) + 13
    assert list(k.uniprot_acc[-13:]) == canary.CANARY_IDS
    d1, d2 = canary.canary_descriptors(), canary.canary_descriptors()
    assert np.array_equal(d1.values, d2.values)


def test_canary_labels_structure():
    c = canary.canary_labels()
    assert set(c.kinase_acc) == set(canary.CANARY_IDS)
    s = index.load_sites().set_index("site_id")
    assert s.loc[c.site_id, "in_U_ms"].all()
    neg = c[c.kinase_acc.str.startswith("CANARY_NEG")]
    per_prot = neg.groupby(["kinase_acc", "acc"]).size()
    assert per_prot.max() <= 4 and (per_prot >= 2).mean() > 0.3        # clustered by protein
    cls = dict(zip(canary.CANARY_IDS, canary.NEG_CLASSES + ["ST"]))
    for k, g in c.groupby("kinase_acc"):
        assert E._res_ok(g.res.values, cls[k]).all()
    pos = c[c.kinase_acc.eq(canary.POS_ID)]
    w = s.loc[pos.site_id, "window15"]
    assert all(x[4] == "R" and x[8] in "FILMV" for x in w)
    assert set(pos.acc) <= set(canary.canary_ppi_partners())
    assert not set(c.acc) & splits.load_test_proteins()


def test_load_dev_with_canaries_flag():
    a = labels.load_dev("L2")
    b = labels.load_dev("L2", with_canaries=True)
    assert len(b) > len(a) and b.is_canary.sum() == len(b) - len(a)
    assert not a.is_canary.any()


def test_canary_check_thresholds():
    pk = pd.DataFrame({"kinase_acc": canary.NEG_IDS + [canary.POS_ID],
                       "auroc_hard": [0.5] * 12 + [0.95], "auroc_all": [0.52] * 12 + [0.93]})
    r = E.canary_check(pk)
    assert r["neg_pass"] and r["pos_pass"]
    pk.loc[:11, "auroc_all"] = 0.6
    assert not E.canary_check(pk)["neg_pass"]
    pk.loc[12, "auroc_hard"] = 0.8
    assert not E.canary_check(pk)["pos_pass"]
