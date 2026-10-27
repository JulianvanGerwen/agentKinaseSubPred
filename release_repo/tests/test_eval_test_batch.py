"""G4 batch path (ksp.eval_test) on a toy 'TEST' = DEV CV fold 0. Every path is redirected to tmp: the real
FREEZE.json / results/test/LEDGER.md / sealed/ are never touched. Covers: known answers (perfect 1, reversed 0,
random ~0.5), the 9-test primary family incl. the unavailable-baseline contingency (p = 1) and Holm, ledger
STARTED/COMPLETED state machine with a simulated crash (attempt 2, .partial -> atomic rename), single-batch refusal,
schema / pin / pre-flight / real-policy / CLI guards, and the real prepare('test') branch."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from ksp import eval as E
from ksp import eval_test as T
from ksp import splits

pytestmark = pytest.mark.skipif(not (E.splits.SPLIT_DIR / "SEAL.json").exists(), reason="needs project data that is not part of this repository (frozen splits under KSP_DATA_DIR); see DATA_SOURCES.md")
RT = {"train_wall_s": 1.0, "inference_wall_s": 1.0}


def _write(df, path):
    df.to_parquet(path, index=False)
    return str(path)


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    w = tmp_path_factory.mktemp("g4")
    build = T.toy_data_builder(0)
    data = build()
    K = data.kinases.set_index("uniprot_acc")
    ks = [k for k in sorted(data.k_primary, key=lambda k: -len(data.positives[k])) if K.residue_class[k] == "ST"][:8]
    kidx = {a: i for i, a in enumerate(data.kinases.uniprot_acc)}
    uni = np.where(data.universe & data.sites.res.isin(["S", "T"]).values)[0]
    rng = np.random.default_rng(0)
    files = {}
    for name, f in [("perfect", lambda y, r: y * 2 + r), ("random", lambda y, r: r),
                    ("reversed", lambda y, r: -2 * y + r)]:
        rows = []
        for k in ks:
            y = np.isin(uni, data.positives[k]).astype(float)
            rows.append(pd.DataFrame({"kinase_idx": kidx[k], "site_idx": uni, "score": f(y, rng.random(len(uni)))}))
        p = pd.concat(rows, ignore_index=True)
        if name == "perfect":      # raw cross-kinase score for M2 (any finite) + A17 assignment probability in [0, 1]
            p["score_site"] = p.score * 10 - 3
            p["assignment_prob"] = (p.score - p.score.min()) / (p.score.max() - p.score.min())
        files[name] = _write(p, w / f"{name}.parquet")
    cfg = {"version": 1, "harness_commit": "test", "B": 2000,   # min p = 2/2001 -> Holm x9 < 0.05
           "configs": [{"name": "final", "role": "final", "predictions": files["perfect"], "runtime": RT},
                       {"name": "rnd", "role": "baseline", "baseline_id": "B1_arm_i", "predictions": files["random"],
                        "runtime": RT},
                       {"name": "rev", "role": "baseline", "baseline_id": "B1_arm_ii",
                        "predictions": files["reversed"], "runtime": RT},
                       {"name": "rnd3", "role": "baseline", "baseline_id": "B3_retrained",
                        "predictions": files["random"], "runtime": RT}],
           "comparisons": [{"a": "final", "b": "rnd", "endpoints": ["M1", "M3"], "primary": True},
                           {"a": "final", "b": "rev", "endpoints": ["M2"], "primary": True},
                           {"a": "final", "b": "rnd3", "endpoints": ["M1", "M2", "M3"], "primary": True},
                           {"a": "rnd", "b": "rev", "endpoints": ["M1"], "primary": False}],
           "unavailable_primary": [{"baseline_id": "B4_retrained", "reason": "toy"}]}
    T.pin_files(cfg)
    fz = w / "FREEZE.json"
    fz.write_text(json.dumps(cfg))
    T.write_toy_sealed(w / "sealed", w / "SEAL.json", 0)
    led = w / "LEDGER.md"
    kw = dict(out_root=w, mode="dry", data_builder=build, kinases=ks)
    # 1st attempt crashes inside the batch (after STARTED, with a .partial dir) ...
    real_run_batch = T.run_batch

    def crash(cfg_, data_, outd, *a, **k):
        outd.mkdir(parents=True, exist_ok=True)
        (outd / "junk").write_text("x")
        raise RuntimeError("simulated crash")
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(T, "run_batch", crash)
        with pytest.raises(RuntimeError, match="simulated crash"):
            T.main(fz, w / "SEAL.json", w / "sealed", led, **kw)
    sha12 = T._sha(fz)[:12]
    crashed = {"partial_exists": (w / f"{sha12}.partial").exists(), "final_exists": (w / sha12).exists(),
               "events": [e["event"] for e in T.ledger_events(led)]}
    assert T.run_batch is real_run_batch
    # ... the identical rerun is attempt 2 and completes
    summ = T.main(fz, w / "SEAL.json", w / "sealed", led, **kw)
    return {"w": w, "cfg": cfg, "summ": summ, "ks": ks, "fz": fz, "led": led, "sha12": sha12, "crashed": crashed,
            "build": build, "files": files}


def test_known_answers(toy):
    s = toy["summ"]
    assert s["final"]["M1"] == pytest.approx(1.0) and s["final"]["M3"] == pytest.approx(1.0)
    assert s["final"]["M2"] == pytest.approx(1.0) and s["final"]["M2_score_column"] == "score_site"
    assert s["rev"]["M1"] == pytest.approx(0.0) and s["rev"]["M3"] == pytest.approx(0.0)
    assert s["rev"]["M2"] < 0.1          # filtering can leave < 5 scored competitors at a few sites
    assert abs(s["rnd"]["M1"] - 0.5) < 0.1


def test_primary_family_holm_and_contingency(toy):
    C = pd.DataFrame(toy["summ"]["_comparisons"])
    prim = C[C.primary]
    assert len(prim) == T.N_PRIMARY
    nt = prim[prim.status.str.startswith("NOT TESTED")]
    assert len(nt) == 3 and (nt.baseline_id == "B4_retrained").all() and (nt.p_boot == 1.0).all()
    assert not nt.beats.any()
    assert np.allclose(prim.p_holm.values, E.holm(prim.p_boot.values))
    sec = C[~C.primary]
    assert np.allclose(sec.p_holm.values, E.holm(sec.p_boot.values))
    t = prim[prim.status == "tested"].set_index(["baseline_id", "endpoint"])
    assert bool(t.loc[("B1_arm_i", "M1"), "beats"]) and t.loc[("B1_arm_i", "M1"), "diff"] == pytest.approx(0.5, abs=0.1)
    for col in ["A_M1_n_kinases_defined", "matched_kinases", "frac_A_pairs_kept", "M2_score_column_A"]:
        assert col in t.columns
    assert t.loc[("B1_arm_i", "M1"), "matched_kinases"] == 8 and t.loc[("B1_arm_i", "M1"), "frac_A_pairs_kept"] == 1


def test_outputs_crash_safety_and_ledger(toy):
    w, sha12 = toy["w"], toy["sha12"]
    cr = toy["crashed"]
    assert cr["partial_exists"] and not cr["final_exists"] and cr["events"] == ["STARTED"]
    outd = w / sha12
    assert outd.exists() and not (w / f"{sha12}.partial").exists() and not (outd / "junk").exists()
    for n in ["final", "rnd", "rev", "rnd3"]:
        for suf in ["_per_kinase.tsv", "_site_level.tsv", "_ci.tsv", "_ci_kinases.tsv", "_strata.tsv",
                    "_ext_tiers.tsv", "_reliability.tsv"]:
            assert (outd / f"{n}{suf}").exists(), n + suf
    st = pd.read_csv(outd / "final_strata.tsv", sep="\t")
    assert {"manning_group", "residue_class", "train_bin", "priming", "pikk", "fam20c", "evidence", "temporal",
            "hub", "primed", "external_curated"} <= set(st.stratum)
    assert (st[st.stratum != "external_curated"].M1.dropna() == 1.0).all()
    ext = st[st.stratum == "external_curated"]
    assert ext.M4.notna().sum() == 1 and ext.M1.notna().sum() == 1      # K_any row carries M4 (minor fix)
    rel = pd.read_csv(outd / "final_reliability.tsv", sep="\t")
    assert {"labelled_TEST_sites"} <= set(rel.subset) and "ECE" in rel
    cov = toy["summ"]["final"]["A17_coverage"]
    assert cov and "coverage@0.5" in cov and "coverage@0.8_unlabelled" in cov
    assert "not head-to-head" in (outd / "summary.md").read_text().lower()
    ev = T.ledger_events(toy["led"])
    assert [e["event"] for e in ev] == ["STARTED", "STARTED", "COMPLETED"]
    assert ev[1]["attempt"] == 2 and ev[1]["batch_id"] == ev[0]["batch_id"] == ev[2]["batch_id"]
    assert "eval_test.py" in ev[0]["code_sha256"] and "strata.py" in ev[0]["code_sha256"]
    assert ev[0]["inputs"] == T.input_hashes(toy["cfg"])


def test_single_batch_refusal_any_completed(toy):
    with pytest.raises(T.FreezeError, match="already COMPLETED"):
        T.main(toy["fz"], toy["w"] / "SEAL.json", toy["w"] / "sealed", toy["led"], out_root=toy["w"], mode="dry",
               data_builder=toy["build"], kinases=toy["ks"])
    # a DIFFERENT FREEZE sha is refused too (single batch whatever the sha)
    with pytest.raises(T.FreezeError, match="already COMPLETED"):
        T.check_ledger(toy["led"], "f" * 64, {}, "x", tag="DRYRUN")
    assert T.check_ledger(toy["led"], "f" * 64, {}, "x", override="re-run approved", tag="DRYRUN") == 1


def test_ledger_state_machine(tmp_path):
    led = tmp_path / "L.md"
    inp, cid = {"a.parquet": "1" * 64}, "c" * 64
    assert T.check_ledger(led, "s" * 64, inp, cid) == 1
    T._ledger_append(led, "STARTED", {"event": "STARTED", "tag": "TEST", "batch_id": "b1", "attempt": 1,
                                      "freeze_sha": "s" * 64, "inputs": inp, "code_identity": cid})
    assert T.check_ledger(led, "s" * 64, inp, cid) == 2                       # identical rerun -> attempt 2
    for bad in [("t" * 64, inp, cid), ("s" * 64, {"a.parquet": "2" * 64}, cid), ("s" * 64, inp, "d" * 64)]:
        with pytest.raises(T.FreezeError, match="unfinished STARTED"):
            T.check_ledger(led, *bad)
    T._ledger_append(led, "COMPLETED", {"event": "COMPLETED", "tag": "TEST", "batch_id": "b1", "freeze_sha": "s" * 64})
    with pytest.raises(T.FreezeError, match="already COMPLETED"):
        T.check_ledger(led, "s" * 64, inp, cid)
    with pytest.raises(T.FreezeError):
        T.check_single_batch("z" * 64, led)
    # old-style text without markers does not count, marker-based blocks do
    led2 = tmp_path / "L2.md"
    led2.write_text("## TEST batch ... FREEZE.json sha256 " + "a" * 64 + "\n")
    T.check_single_batch("b" * 64, led2)


def test_validate_freeze_primary_structure(toy):
    base = toy["cfg"]

    def bad(mut, match):
        c = json.loads(json.dumps(base))
        mut(c)
        with pytest.raises(T.FreezeError, match=match):
            T.validate_freeze(c)
    T.validate_freeze(json.loads(json.dumps(base)))                         # the toy config itself is valid
    bad(lambda c: c["comparisons"][0].update(a="rnd"), "a = the final")
    bad(lambda c: c["comparisons"][0].update(endpoints=["M1", "M2"]), "not pre-declared")
    bad(lambda c: c["comparisons"][2].update(endpoints=["M1", "M2", "M4"]), "M4|not pre-declared")
    bad(lambda c: c["comparisons"].append({"a": "final", "b": "rnd", "endpoints": ["M1"], "primary": True}),
        "duplicate")
    bad(lambda c: c.update(unavailable_primary=[]), "exactly the 9")
    bad(lambda c: c["configs"][1].update(role="final"), "exactly one")
    bad(lambda c: c["configs"][1].pop("baseline_id"), "baseline_id")
    bad(lambda c: c["configs"][0].update(missing="Worst"), "missing must be 'worst'")
    bad(lambda c: c["configs"][0].update(predictions_sha256="0" * 64), "differs")
    bad(lambda c: c["configs"][0].pop("predictions_sha256"), "sha256")
    bad(lambda c: c["configs"][0].pop("runtime"), "runtime")
    # reviewer toy_checks #2: a wrong primary set with the right count is rejected
    c = json.loads(json.dumps(base))
    c["comparisons"] = [{"a": "rnd", "b": "rev", "primary": True}, {"a": "final", "b": "rnd3", "primary": True},
                        {"a": "final", "b": "rnd", "endpoints": ["M4", "M4", "M4"], "primary": True}]
    c.pop("unavailable_primary")
    with pytest.raises(T.FreezeError):
        T.validate_freeze(c)


@pytest.mark.skipif(not os.environ.get("KSP_ORIGINAL_REPO"),
                    reason="needs the original project repository: checks that the frozen harness files are unchanged "
                           "since git commit 974d520 (A28). This release changes their path handling and has no such "
                           "commit; set KSP_ORIGINAL_REPO=1 only in the original repository")
def test_real_policy(toy):
    c = json.loads(json.dumps(toy["cfg"]))
    with pytest.raises(T.FreezeError, match="B == 10000"):
        T.check_real_policy(c)
    c["B"] = 10000
    with pytest.raises(T.FreezeError, match="A28"):
        T.check_real_policy(c)
    c["harness_commit"] = T.A28_COMMIT
    T.check_real_policy(c)                                                 # frozen files unchanged since 974d520
    with pytest.raises(T.FreezeError, match="no test hooks"):
        T.main(toy["fz"], mode="real", data_builder=toy["build"])


def test_preflight_rejects(toy, tmp_path):
    data = toy["build"]()
    cfg0 = json.loads(json.dumps(toy["cfg"]))
    p = pd.read_parquet(toy["files"]["perfect"])

    def run(df, key="predictions"):
        c = json.loads(json.dumps(cfg0))
        f = tmp_path / f"x_{key}.parquet"
        df.to_parquet(f, index=False)
        c["configs"][0][key] = str(f)
        T.pin_files(c)
        return T.preflight(c, data.universe, data.k_any)
    with pytest.raises(T.FreezeError, match="non-finite score_site"):
        run(p.assign(score_site=np.where(np.arange(len(p)) == 0, np.nan, p.score_site)))
    with pytest.raises(T.FreezeError, match="non-finite assignment_prob"):
        run(p.assign(assignment_prob=np.where(np.arange(len(p)) == 0, np.nan, p.assignment_prob)))
    run(p.assign(score_site=p.score_site * 100))                              # any finite score_site is fine
    with pytest.raises(T.FreezeError, match="non-finite score"):
        run(p.assign(score=np.where(np.arange(len(p)) == 0, np.inf, p.score)))
    with pytest.raises(T.FreezeError, match="not exhaustive"):
        run(p.iloc[10:])                                                     # partial kinase
    with pytest.raises(T.FreezeError, match="out of range"):
        run(pd.concat([p, p.iloc[:1].assign(kinase_idx=503)]))              # canary index
    with pytest.raises(T.FreezeError, match=r"\[0, 1\]"):
        run(p.assign(assignment_prob=p.assignment_prob * 2))
    with pytest.raises(T.FreezeError, match="T3 predictions miss"):
        run(p, key="t3_predictions")                                         # covers 8 kinases, not all K_any


def test_cli_guards():
    assert T._cli([]) == 2
    with pytest.raises(SystemExit):
        T._cli(["--real", "--B", "100"])
    with pytest.raises(SystemExit):
        T._cli(["--dry-run", "--override-completed", "x"])


def test_real_prepare_branch(toy, tmp_path, monkeypatch):
    """Reviewer check_prepare_test: exercise ksp.eval.prepare('test', labels_df=...) (data_builder=None), with
    TEST proteins monkeypatched to fold-0 proteins and the toy sealed set. Never reads sealed/."""
    folds = splits.load_cv_folds()
    f0 = set(folds.index[folds == 0])
    monkeypatch.setattr(splits, "load_test_proteins", lambda: f0)
    c = {"version": 1, "harness_commit": "test", "B": 200,
         "configs": [{"name": "final", "role": "final", "predictions": toy["files"]["perfect"], "runtime": RT}]}
    T.pin_files(c)
    fz = tmp_path / "FREEZE.json"
    fz.write_text(json.dumps(c))
    s = T.main(fz, toy["w"] / "SEAL.json", toy["w"] / "sealed", tmp_path / "L.md", out_root=tmp_path, mode="dry",
               kinases=toy["ks"], allow_nonstandard_primary=True)
    assert s["final"]["M1"] == pytest.approx(1.0)


def test_check_predictions_all_or_nothing(tmp_path):
    c = T.test_candidate_pairs()
    k0, k1 = c.kinase_idx.unique()[:2]
    full = c[c.kinase_idx == k0]
    part = c[c.kinase_idx == k1].iloc[:100]
    f = tmp_path / "p.parquet"
    pd.concat([full, part]).assign(score=0.5).to_parquet(f, index=False)
    r = T.check_predictions(f)
    assert r["kinases_full"] == 1 and r["kinases_partial"] == 1 and not r["ok"]
    full.assign(score=1.0).to_parquet(f, index=False)
    assert T.check_predictions(f)["ok"]
    non_test = int(np.setdiff1d(np.arange(1000), c.site_idx.unique())[0])
    pd.concat([full, pd.DataFrame({"kinase_idx": [k0], "site_idx": [non_test]})]).assign(score=1.0) \
        .to_parquet(f, index=False)
    r = T.check_predictions(f)
    assert r["rows_outside_candidates"] == 1 and not r["ok"]


def test_real_paths_untouched():
    # Pre-G4 this asserted FREEZE.json did not exist. Since the G4 batch (FREEZE.json committed 5fc7870) the
    # invariant is: dry-run tests never write it, i.e. if present it is byte-identical to the committed version.
    import subprocess
    assert T.FREEZE == T.ROOT / "FREEZE.json"
    if T.FREEZE.exists():
        r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "FREEZE.json"], cwd=T.ROOT)
        assert r.returncode == 0, "FREEZE.json modified relative to HEAD"
    assert not T.LEDGER.exists() or "DRYRUN" not in T.LEDGER.read_text()


def test_side_files_merge_and_key_guards(toy, tmp_path):
    p = pd.read_parquet(toy["files"]["random"])
    main = _write(p[["kinase_idx", "site_idx", "score"]], tmp_path / "m.parquet")
    ap = _write(p[["kinase_idx", "site_idx"]].assign(assignment_prob=0.3), tmp_path / "ap.parquet")
    c = {"name": "x", "predictions": main, "assignment_prob_predictions": ap,
         "score_site_predictions": _write(p[["kinase_idx", "site_idx"]].assign(score_site=-5.0), tmp_path / "s.parquet")}
    q = T.load_predictions(c)
    assert (q.assignment_prob == 0.3).all() and (q.score_site == -5.0).all()
    from ksp import index
    S, K = index.load_sites(), index.load_kinases()
    mixed = pd.DataFrame({"kinase_acc": K.uniprot_acc.values[p.kinase_idx.values],
                          "site_id": S.site_id.values[p.site_idx.values], "assignment_prob": 0.3})
    c["assignment_prob_predictions"] = _write(mixed, tmp_path / "ap2.parquet")
    with pytest.raises(T.FreezeError, match="different key columns"):
        T.load_predictions(c)
    c["assignment_prob_predictions"] = ap
    c["predictions"] = _write(p.assign(assignment_prob=0.1), tmp_path / "m2.parquet")
    with pytest.raises(T.FreezeError, match="both inline and as a separate file"):
        T.load_predictions(c)
