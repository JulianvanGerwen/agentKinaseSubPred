"""A11 (DESIGN §16): a feature source whose (kinase_acc, site_id) pairs overlap DEV PSP KS pairs at > 3x chance
fails. Feature agents call `check_feature_overlap(pairs)` from their own tests:
    from tests.test_feature_overlap import check_feature_overlap   (or copy the call)
Chance = |F| * |PSP| / |candidate pairs| (candidates: kinases x U_ms DEV sites of matching residue class)."""
import numpy as np
import pandas as pd
import pytest

from ksp import index, labels, splits

MAX_RATIO = 3.0


def overlap_ratio(pairs: pd.DataFrame) -> dict:
    psp = labels.load_dev("L2")
    psp_keys = set(zip(psp.kinase_acc, psp.site_id))
    test = splits.load_test_proteins()
    s = index.load_sites()
    s = s[s.in_U_ms & ~s.acc.isin(test)]
    k = index.load_kinases()
    n_st, n_y = int(s.res.isin(["S", "T"]).sum()), int((s.res == "Y").sum())
    n_cand = sum(n_st if c == "ST" else n_y if c == "Y" else n_st + n_y for c in k.residue_class)
    f = pairs.drop_duplicates(["kinase_acc", "site_id"])
    f = f[~splits.pd.Series(splits.site_protein(f.site_id)).isin(test).values]
    obs = sum((a, b) in psp_keys for a, b in zip(f.kinase_acc, f.site_id))
    exp = len(f) * len(psp_keys) / max(n_cand, 1)
    return {"n_feature_pairs": len(f), "observed_overlap": obs, "expected_overlap": exp,
            "ratio": obs / exp if exp > 0 else np.inf if obs else 0.0}


def check_feature_overlap(pairs: pd.DataFrame, max_ratio: float = MAX_RATIO) -> dict:
    r = overlap_ratio(pairs)
    if r["ratio"] > max_ratio:
        raise AssertionError(f"feature source overlaps PSP KS pairs at {r['ratio']:.1f}x chance (> {max_ratio}x): {r}")
    return r


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet", "cache/labels/dev_pairs.parquet", "results/splits/SEAL.json")
def test_helper_flags_label_copy_and_passes_random():
    d = labels.load_dev("L2")
    with pytest.raises(AssertionError):
        check_feature_overlap(d[["kinase_acc", "site_id"]].sample(500, random_state=0))
    s = index.load_sites()
    s = s[s.in_U_ms & ~s.acc.isin(splits.load_test_proteins()) & s.res.isin(["S", "T"])]
    k = index.load_kinases()
    k = k[k.residue_class.eq("ST")]
    rng = np.random.default_rng(0)
    rnd = pd.DataFrame({"kinase_acc": rng.choice(k.uniprot_acc.values, 20000),
                        "site_id": rng.choice(s.site_id.values, 20000)})
    assert check_feature_overlap(rnd)["ratio"] < MAX_RATIO
