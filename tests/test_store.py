import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import store
from ksp.index import load_kinases


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet")
def test_roundtrip_and_reorder(tmp_path, monkeypatch):
    monkeypatch.setattr(store, "STORE", tmp_path)
    k = load_kinases().uniprot_acc.to_numpy()
    perm = np.random.default_rng(0).permutation(len(k))
    vals = np.arange(len(k), dtype=np.float32)
    shuffled_ids, shuffled_vals = k[perm], vals[perm]
    canon = store.reorder_to_index(shuffled_vals[:, None], shuffled_ids, "kinase", 0)
    assert np.array_equal(canon[:, 0], vals)          # reorder recovers canonical order
    store.save_matrix("t/x", canon, "kinase", "none", col_ids=["f"])
    arr, desc = store.load_matrix("t/x")
    assert np.array_equal(np.asarray(arr)[:, 0], vals)


@pytest.mark.requires_data("cache/index/kinases.parquet", "cache/index/sites.parquet")
def test_unknown_id_raises():
    import pytest
    with pytest.raises(KeyError):
        store.reorder_to_index(np.zeros((1,), np.float32), ["NOT_A_KINASE"], "kinase", 0)
