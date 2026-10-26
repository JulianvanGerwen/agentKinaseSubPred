"""G4 batch path (ksp.eval_test) on a toy 'TEST' = DEV CV fold 0, with every path monkeypatched to tmp (the real
FREEZE.json / LEDGER.md / sealed/ are never touched). Known answers: perfect scorer -> 1, reversed -> 0,
random -> ~0.5; Holm over exactly the 9 primary tests; single-batch and seal guards hold."""
import json

import numpy as np
import pandas as pd
import pytest

from ksp import eval as E
from ksp import eval_test as T

pytestmark = pytest.mark.skipif(not (E.splits.SPLIT_DIR / "SEAL.json").exists(), reason="split not frozen")
RT = {"train_wall_s": 1.0, "inference_wall_s": 1.0}


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    w = tmp_path_factory.mktemp("g4")
    build = T.toy_data_builder(0)
    data = build()
    # 8 ST kinases of K_eval(toy) with most positives
    K = data.kinases.set_index("uniprot_acc")
    ks = [k for k in sorted(data.k_primary, key=lambda k: -len(data.positives[k])) if K.residue_class[k] == "ST"][:8]
    S = data.sites
    uni = np.where(data.universe & S.res.isin(["S", "T"]).values)[0]
    rng = np.random.default_rng(0)
    files = {}
    for name, f in [("perfect", lambda y, r: y * 2 + r), ("random", lambda y, r: r), ("reversed", lambda y, r: -2 * y + r)]:
        rows = []
        for k in ks:
            y = np.isin(uni, data.positives[k]).astype(float)
            rows.append(pd.DataFrame({"kinase_acc": k, "site_id": S.site_id.values[uni],
                                      "score": f(y, rng.random(len(uni)))}))
        p = w / f"{name}.parquet"
        pd.concat(rows, ignore_index=True).to_parquet(p, index=False)
        files[name] = str(p)
    cfg = {"version": 1, "harness_commit": "test", "B": 2000,   # min p = 2/2001 -> Holm x9 < 0.05
           "configs": [{"name": n, "role": "final" if n == "perfect" else "baseline", "predictions": files[n],
                        "runtime": RT} for n in files],
           "comparisons": [{"a": "perfect", "b": "random", "endpoints": ["M1", "M2", "M3"], "primary": True},
                           {"a": "perfect", "b": "reversed", "endpoints": ["M1", "M2", "M3"], "primary": True},
                           {"a": "random", "b": "reversed", "endpoints": ["M1", "M2", "M3"], "primary": True},
                           {"a": "random", "b": "perfect", "endpoints": ["M1"], "primary": False}]}
    fz = w / "FREEZE.json"
    fz.write_text(json.dumps(cfg))
    T.write_toy_sealed(w / "sealed", w / "SEAL.json", 0)
    summ = T.main(fz, w / "SEAL.json", w / "sealed", w / "LEDGER.md", out_root=w, require_git=False,
                  data_builder=build, kinases=ks)
    return w, cfg, summ, ks


def test_known_answers(toy):
    _, _, s, _ = toy
    assert s["perfect"]["M1"] == pytest.approx(1.0) and s["perfect"]["M3"] == pytest.approx(1.0)
    assert s["perfect"]["M2"] == pytest.approx(1.0)
    assert s["reversed"]["M1"] == pytest.approx(0.0) and s["reversed"]["M3"] == pytest.approx(0.0)
    assert s["reversed"]["M2"] < 0.1     # filtering can leave < 5 scored competitors at a few sites
    assert abs(s["random"]["M1"] - 0.5) < 0.1


def test_superiority_holm_and_outputs(toy):
    w, cfg, s, _ = toy
    C = pd.DataFrame(s["_comparisons"])
    prim = C[C.primary]
    assert len(prim) == T.N_PRIMARY
    assert np.allclose(prim.p_holm.values, E.holm(prim.p_boot.values))
    assert np.allclose(C[~C.primary].p_holm.values, E.holm(C[~C.primary].p_boot.values))
    pr = prim.set_index(["a", "b", "endpoint"])
    assert pr.loc[("perfect", "reversed", "M1"), "diff"] == pytest.approx(1.0)
    assert bool(pr.loc[("perfect", "random", "M1"), "beats"])
    outd = next(p for p in w.iterdir() if p.is_dir() and (p / "summary.md").exists())
    for n in ["perfect", "random", "reversed"]:
        st = pd.read_csv(outd / f"{n}_strata.tsv", sep="\t")
        assert {"manning_group", "residue_class", "train_bin", "priming", "pikk", "fam20c", "evidence", "temporal",
                "hub", "primed", "external_curated"} <= set(st.stratum)
    st = pd.read_csv(outd / "perfect_strata.tsv", sep="\t")
    assert (st[st.stratum != "external_curated"].M1.dropna() == 1.0).all()   # ext pairs are not eval positives
    assert (outd / "superiority.tsv").exists() and (outd / "secondary.tsv").exists()
    assert "FREEZE.json sha256" in (w / "LEDGER.md").read_text()


def test_guards(toy, tmp_path):
    w, cfg, _, ks = toy
    build = T.toy_data_builder(0)
    with pytest.raises(T.FreezeError, match="already scored"):           # single batch
        T.main(w / "FREEZE.json", w / "SEAL.json", w / "sealed", w / "LEDGER.md", out_root=w, require_git=False,
               data_builder=build, kinases=ks)
    with pytest.raises(T.FreezeError, match="not committed"):            # real path requires git-tracked FREEZE
        T.check_freeze(w / "FREEZE.json", require_git=True)
    bad = json.loads(json.dumps(cfg))
    bad["comparisons"] = bad["comparisons"][:2]                          # 6 primary tests != 9
    with pytest.raises(T.FreezeError, match="exactly 9"):
        T.validate_freeze(bad)
    bad = json.loads(json.dumps(cfg))
    del bad["configs"][0]["runtime"]
    with pytest.raises(T.FreezeError, match="runtime"):
        T.validate_freeze(bad)
    sd = tmp_path / "sealed"
    sd.mkdir()
    (sd / "test_labels.parquet").write_bytes((w / "sealed" / "test_labels.parquet").read_bytes() + b"x")
    with pytest.raises(T.FreezeError, match="SHA-256"):                  # tampered seal
        T.check_seal(w / "SEAL.json", sd)


def test_real_paths_untouched(toy):
    assert T.FREEZE == T.ROOT / "FREEZE.json" and not T.FREEZE.exists()
    assert not (T.LEDGER.exists() and "harness_commit test" in T.LEDGER.read_text())


def test_check_predictions_all_or_nothing(tmp_path):
    c = T.test_candidate_pairs()
    k0, k1 = c.kinase_idx.unique()[:2]
    full = c[c.kinase_idx == k0]
    part = c[c.kinase_idx == k1].iloc[:100]
    p = pd.concat([full, part]).assign(score=0.5)
    f = tmp_path / "p.parquet"
    p.to_parquet(f, index=False)
    r = T.check_predictions(f)
    assert r["kinases_full"] == 1 and r["kinases_partial"] == 1 and not r["ok"]
    full.assign(score=1.0).to_parquet(f, index=False)
    assert T.check_predictions(f)["ok"]
    non_test = int(np.setdiff1d(np.arange(1000), c.site_idx.unique())[0])                          # a non-TEST site
    bad = pd.concat([full, pd.DataFrame({"kinase_idx": [k0], "site_idx": [non_test]})]).assign(score=1.0)
    bad.to_parquet(f, index=False)
    r = T.check_predictions(f)
    assert r["rows_outside_candidates"] == 1 and not r["ok"]
