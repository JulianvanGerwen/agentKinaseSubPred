"""ksp.deploy (G5): refuses without a COMPLETED TEST batch for the FREEZE sha; returns DEV + TEST labels (toy
sealed set on public TEST proteins - never the real sealed/), canary-free by default."""
import json

import numpy as np
import pandas as pd
import pytest

from ksp import deploy as D
from ksp import eval_test as ET
from ksp import index, labels, splits

pytestmark = pytest.mark.skipif(not labels.DEV_PAIRS.exists(), reason="needs project data that is not part of this repository (frozen splits under KSP_DATA_DIR); see DATA_SOURCES.md")


@pytest.fixture()
def toy(tmp_path):
    """Toy sealed dir: 30 fake label rows (DEV row schema) moved onto U_ms sites of TEST proteins + 1 canary row."""
    dev = labels.load_dev("all")
    S = index.load_sites()
    ts = S[S.in_U_ms & S.acc.isin(splits.load_test_proteins()) & S.res.isin(["S", "T"])].head(30)
    rows = dev[dev.provenance.eq("PSP") & dev.L2.astype(bool) & dev.res.isin(["S", "T"])].head(30).copy()
    rows["site_id"], rows["acc"], rows["pos"], rows["res"] = ts.site_id.values, ts.acc.values, ts.pos.values, \
        ts.res.values
    rows["in_U_ms"] = True
    can = rows.head(1).assign(kinase="CANARY_POS", kinase_acc="CANARY_POS", is_canary=True)
    sd = tmp_path / "sealed"
    sd.mkdir()
    pd.concat([rows.assign(is_canary=False), can], ignore_index=True).to_parquet(sd / "test_labels.parquet",
                                                                                 index=False)
    seal = tmp_path / "SEAL.json"
    seal.write_text(json.dumps({"files": {"test_labels.parquet": {"sha256": ET._sha(sd / "test_labels.parquet")}}}))
    fz = tmp_path / "FREEZE.json"
    fz.write_text('{"version": 1}')
    led = tmp_path / "LEDGER.md"
    kw = dict(freeze=fz, ledger=led, seal=seal, sealed_dir=sd, require_git=False)
    return {"kw": kw, "fz": fz, "led": led, "rows": rows, "sd": sd}


def _complete(led, sha, tag="TEST"):
    ET._ledger_append(led, "STARTED", {"event": "STARTED", "tag": tag, "batch_id": "b", "attempt": 1,
                                       "freeze_sha": sha, "inputs": {}, "code_identity": "x"})
    ET._ledger_append(led, "COMPLETED", {"event": "COMPLETED", "tag": tag, "batch_id": "b", "freeze_sha": sha})


def test_refuses_without_completed(toy):
    with pytest.raises(D.DeployError, match="no COMPLETED TEST batch"):
        D.load_deploy_labels(**toy["kw"])                                     # empty ledger
    ET._ledger_append(toy["led"], "STARTED", {"event": "STARTED", "tag": "TEST", "batch_id": "b", "attempt": 1,
                                              "freeze_sha": ET._sha(toy["fz"]), "inputs": {}, "code_identity": "x"})
    with pytest.raises(D.DeployError, match="no COMPLETED"):
        D.load_deploy_labels(**toy["kw"])                                     # STARTED only (batch running)
    _complete(toy["led"], "0" * 64)
    with pytest.raises(D.DeployError, match="no COMPLETED"):
        D.load_deploy_labels(**toy["kw"])                                     # COMPLETED for another FREEZE sha
    _complete(toy["led"], ET._sha(toy["fz"]), tag="DRYRUN")
    with pytest.raises(D.DeployError, match="no COMPLETED"):
        D.load_deploy_labels(**toy["kw"])                                     # a dry run does not count
    with pytest.raises(D.DeployError, match="not committed"):
        D.check_test_completed(toy["fz"], toy["led"], require_git=True)       # real call requires a tracked FREEZE


def test_superset_of_dev_when_completed(toy):
    _complete(toy["led"], ET._sha(toy["fz"]))
    for variant in ["L2", "all", "eval"]:
        dev = labels.load_dev(variant)
        dep = D.load_deploy_labels(variant, **toy["kw"])
        kd = set(zip(dev.kinase_acc, dev.site_id, dev.provenance))
        kp = set(zip(dep.kinase_acc, dep.site_id, dep.provenance))
        assert kd < kp, variant                                               # strict superset
        assert list(dep.columns) == list(dev.columns)
    dep = D.load_deploy_labels("L2", **toy["kw"])
    assert set(toy["rows"].site_id) <= set(dep.site_id) and not dep.is_canary.any()
    depc = D.load_deploy_labels("L2", with_canaries=True, **toy["kw"])
    assert depc.is_canary.any() and "CANARY_POS" in set(depc.kinase_acc)
    assert set(depc[depc.is_canary].kinase_acc) >= set(labels.load_dev("L2", with_canaries=True)
                                                       .query("is_canary").kinase_acc)


def test_refuses_tampered_seal(toy):
    _complete(toy["led"], ET._sha(toy["fz"]))
    f = toy["sd"] / "test_labels.parquet"
    f.write_bytes(f.read_bytes() + b"x")
    with pytest.raises(D.DeployError, match="SHA-256"):
        D.load_deploy_labels(**toy["kw"])


def test_refuses_protein_overlap(toy):
    _complete(toy["led"], ET._sha(toy["fz"]))
    dev = labels.load_dev("all")
    bad = dev.head(5).assign(is_canary=False)                                 # DEV proteins inside "TEST"
    f = toy["sd"] / "test_labels.parquet"
    bad.to_parquet(f, index=False)
    seal = toy["kw"]["seal"]
    seal.write_text(json.dumps({"files": {"test_labels.parquet": {"sha256": ET._sha(f)}}}))
    with pytest.raises(D.DeployError, match="share"):
        D.load_deploy_labels(**toy["kw"])
