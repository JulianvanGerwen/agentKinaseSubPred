"""ksp.eval on real frozen DEV data (CV fold 0): ID validation, missing predictions, oracle/random behaviour,
filtered site-centric ranks (A18), canary acceptance end-to-end, TEST refusal without FREEZE."""
import numpy as np
import pandas as pd
import pytest

from ksp import canary, eval as E, eval_test, manifest

pytestmark = pytest.mark.skipif(not (E.splits.SPLIT_DIR / "SEAL.json").exists(), reason="needs project data that is not part of this repository (frozen splits under KSP_DATA_DIR); see DATA_SOURCES.md")


@pytest.fixture(scope="module")
def data():
    return E.prepare("cv", fold=0)


def _preds(data, kinases, score_fn, seed=0):
    rng = np.random.default_rng(seed)
    S = data.sites
    uni = np.where(data.universe)[0]
    out = []
    K = data.kinases.set_index("uniprot_acc")
    for k in kinases:
        s = uni[E._res_ok(S.res.values[uni], K.residue_class[k])]
        sc = score_fn(k, s, rng)
        out.append(pd.DataFrame({"kinase_acc": k, "site_id": S.site_id.values[s], "score": sc}))
    return pd.concat(out, ignore_index=True)


def _top_kinases(data, n=8):
    return sorted((k for k in data.positives if not k.startswith("CANARY")),
                  key=lambda k: -len(data.positives[k]))[:n]


def test_unknown_ids_and_duplicates_raise(data):
    k = _top_kinases(data, 1)[0]
    p = _preds(data, [k], lambda k, s, r: r.random(len(s)))
    with pytest.raises(KeyError):
        E.evaluate(pd.concat([p, pd.DataFrame({"kinase_acc": ["NOTAKINASE"], "site_id": [p.site_id[0]],
                                                "score": [0.1]})]), data, kinases=[k])
    with pytest.raises(KeyError):
        E.evaluate(pd.concat([p, pd.DataFrame({"kinase_acc": [k], "site_id": ["P00000_S1"], "score": [0.1]})]),
                   data, kinases=[k])
    with pytest.raises(ValueError):
        E.evaluate(pd.concat([p, p.iloc[:1]]), data, kinases=[k])
    q = p.copy()
    q.loc[0, "score"] = np.nan
    with pytest.raises(ValueError):
        E.evaluate(q, data, kinases=[k])


def test_oracle_random_and_missing(data):
    ks = _top_kinases(data, 6)
    pos = {k: set(data.positives[k]) for k in ks}
    orc = _preds(data, ks, lambda k, s, r: r.random(len(s)) + np.isin(s, list(pos[k])))
    r = E.evaluate(orc, data, kinases=ks)
    assert (r.per_kinase.auroc_hard == 1).all() and (r.per_kinase.top1_recall == 1).all()
    rnd = _preds(data, ks, lambda k, s, r: r.random(len(s)))
    rr = E.evaluate(rnd, data, kinases=ks)
    assert abs(rr.per_kinase.auroc_all.mean() - 0.5) < 0.08
    # drop all positives of the first kinase -> counted; 'worst' ranks them last; 'error' raises
    k0 = ks[0]
    pid = set(data.sites.site_id.values[list(pos[k0])])
    part = orc[~((orc.kinase_acc == k0) & orc.site_id.isin(pid))]
    rm = E.evaluate(part, data, kinases=ks)
    assert rm.missing["per_kinase"][k0] == len(pid)
    assert rm.per_kinase.set_index("kinase_acc").auroc_hard[k0] == 0.0
    with pytest.raises(ValueError):
        E.evaluate(part, data, kinases=ks, missing="error")
    rd = E.evaluate(part, data, kinases=ks, missing="drop")
    assert rd.per_kinase.set_index("kinase_acc").n_pos[k0] == 0


def test_site_centric_filtering(data):
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY")]
    multi = pm.groupby("site_idx").kinase_acc.agg(list)
    multi = multi[multi.map(len) >= 2]
    s = multi.index[0]
    k_true, k_other = multi.iloc[0][0], multi.iloc[0][1]
    K = data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY")]
    res = data.sites.res.values[s]
    cand = K.uniprot_acc[E._res_ok(np.array([res] * len(K)), "ST") if res in "ST" else K.residue_class.ne("ST")]
    sid = data.sites.site_id.values[s]
    sc = {a: 0.0 for a in cand}
    sc[k_other] = 10.0     # other true kinase outranks: filtered away
    for j in data.par.get(k_true, ()):
        if j in sc:
            sc[j] = 9.0    # paralogs not labelled to s: excluded (A18)
    sc[k_true] = 5.0
    p = pd.DataFrame({"kinase_acc": list(sc), "site_id": sid, "score": list(sc.values())})
    sl = E._site_centric(E.map_predictions(p, data), data, "worst")
    row = sl[(sl.kinase_acc == k_true) & (sl.site_idx == s)].iloc[0]
    assert row.n_greater == 0 and row.hit1 == 1.0
    # separate score_site column is used for M2 when present (A20)
    p["score_site"] = p.score
    p.loc[p.kinase_acc == k_true, "score_site"] = -1.0
    sl2 = E._site_centric(E.map_predictions(p, data), data, "worst")
    assert sl2[(sl2.kinase_acc == k_true) & (sl2.site_idx == s)].iloc[0].hit1 < 1.0


def test_canaries_end_to_end(data):
    """Random scores pass the negative-canary threshold; a planted-motif score passes the positive one."""
    ks = [k for k in canary.CANARY_IDS if k in data.positives]
    W = data.sites.window15.values

    def fn(k, s, r):
        base = r.random(len(s))
        if k == canary.POS_ID:
            w = W[s]
            mot = np.array([x[4] == "R" and x[7] in "ST" and x[8] in "FILMV" for x in w])
            part = np.isin(data.sites.acc.values[s], canary.canary_ppi_partners())
            return base + mot + part
        return base
    r = E.evaluate(_preds(data, ks, fn), data, kinases=ks)
    assert r.canary["neg_pass"] and r.canary["pos_pass"], r.canary


def test_canary_v2_not_separable_by_label_free_or_labelledness_scores(data):
    """R3 audit of v1: -(#U_ms sites on the protein) separated v1 canaries (AUROC 0.76). Under v2 neither that
    label-free protein-size score nor a perfect 'is a labelled site' score may pass 0.57 on E-hard."""
    ks = [k for k in canary.NEG_IDS if k in data.positives]
    nprot = data.sites[data.sites.in_U_ms].groupby("acc").size()
    nsz = data.sites.acc.map(nprot).fillna(0).values
    lab = np.zeros(len(data.sites))
    lab[data.masks.hard_idx] = 1.0
    ev = E.labels.load_dev("eval")
    nlab = data.sites.acc.map(ev.groupby("acc").site_id.nunique()).fillna(0).values
    for fn in (lambda k, s, r: -nsz[s] + 1e-6 * r.random(len(s)), lambda k, s, r: lab[s] + 1e-6 * r.random(len(s)),
               lambda k, s, r: nlab[s] + 1e-6 * r.random(len(s))):
        r = E.evaluate(_preds(data, ks, fn), data, kinases=ks)
        assert r.canary["neg_mean_auroc_hard"] <= E.NEG_CANARY_MAX, r.canary
        assert r.canary["neg_pass"]


def test_bootstrap_and_paired(data):
    ks = _top_kinases(data, 5)
    pos = {k: set(data.positives[k]) for k in ks}
    a = _preds(data, ks, lambda k, s, r: r.random(len(s)) + 0.5 * np.isin(s, list(pos[k])))
    b = a.assign(score=np.random.default_rng(5).random(len(a)))
    ra, rb = E.evaluate(a, data, kinases=ks), E.evaluate(b, data, kinases=ks)
    ci = E.bootstrap(ra, B=100)
    assert set(ci.endpoint) == {"M1", "M2", "M3", "M4"}
    d = E.paired_bootstrap(ra, rb, B=200)
    m1 = d.set_index("endpoint").loc["M1"]
    assert m1["diff"] > 0 and m1.lo95 > 0
    # unmatched universe raises
    with pytest.raises(ValueError):
        E.paired_bootstrap(ra, E.evaluate(b.iloc[10:], data, kinases=ks, missing="drop"), B=10)
    x, y, rep = E.matched_universe(a, b.iloc[10:])
    assert rep["n_common"] == len(a) - 10


def test_eval_test_refuses_without_freeze(tmp_path):
    with pytest.raises(eval_test.FreezeError):
        eval_test.check_freeze(tmp_path / "FREEZE.json")


def test_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest, "RESULTS", tmp_path / "res")
    monkeypatch.setattr(manifest, "CACHE", tmp_path / "cache")
    f = tmp_path / "in.txt"
    f.write_text("x")
    with manifest.RunManifest("toy", {"a": 1}, seed=0, inputs=[f]) as run:
        run.count_fallback("nearest_matrix", 3)
        run.save_predictions(pd.DataFrame({"kinase_acc": ["A"], "site_id": ["B"], "score": [1.0]}))
        run.save_metrics({"M1": 0.5})
    import json
    m = json.loads((tmp_path / "res" / run.run_id / "manifest.json").read_text())
    assert m["status"] == "ok" and m["fallbacks"] == {"nearest_matrix": 3} and m["n_rows"]["predictions"] == 1
    assert m["inputs"][0]["sha256"]
    with pytest.raises(FileExistsError):
        manifest.RunManifest("toy", {}, 0, run_id=run.run_id)


def test_fold_offset_invariance_pooled_cv():
    """A24: adding a constant offset to one fold's scores leaves pooled M1/M3 unchanged under fold rank-normalisation
    (and changes them without it)."""
    d = E.prepare("cv", fold=None, with_canaries=False)
    ks = _top_kinases(d, 6)
    pos = {k: set(d.positives[k]) for k in ks}
    a = _preds(d, ks, lambda k, s, r: r.random(len(s)) + 0.7 * np.isin(s, list(pos[k])))
    b = a.copy()
    off = d.site_fold[pd.Series(np.arange(len(d.sites)), index=d.sites.site_id)[b.site_id].values] == 1
    b.loc[off, "score"] += 5.0
    ra, rb = E.evaluate(a, d, kinases=ks), E.evaluate(b, d, kinases=ks)
    assert ra.endpoints["fold_normalised"] and rb.endpoints["fold_normalised"]
    pa, pb = ra.per_kinase.set_index("kinase_acc"), rb.per_kinase.set_index("kinase_acc")
    assert np.allclose(pa.auroc_hard, pb.auroc_hard) and np.allclose(pa.top1_recall, pb.top1_recall)
    assert len(ra.endpoints["M1_per_fold"]) == 5
    rc = E.evaluate(b, d, kinases=ks, fold_normalise=False)
    assert not np.allclose(rc.per_kinase.set_index("kinase_acc").auroc_hard, pa.auroc_hard)
