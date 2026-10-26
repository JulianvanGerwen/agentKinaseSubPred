"""Toy-data tests for ksp.strata (canonical stratum definitions)."""
import numpy as np
import pandas as pd
import pytest

from ksp import strata as ST


def _kinases():
    return pd.DataFrame({
        "uniprot_acc": ["K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "K11"],
        "gene": ["GSK3B", "CSNK1D", "CSNK2A1", "PLK1", "ATM", "MTOR", "FAM20C", "SRC", "PRKACA", "TRRAP", "CSNK1G2"],
        "residue_class": ["ST", "ST", "ST", "ST", "ST", "ST", "ST", "Y", "ST", "ST", "ST"],
        "manning_group": ["CMGC", "CK1", "CMGC", "Other", "Atypical", "Atypical", "Atypical", "TK", "AGC", "Atypical",
                          "CK1"],
        "manning_family": ["GSK", "CK1", "CK2", "PLK", "PIKK", "PIKK", "FAM20", "Src", "PKA", "PIKK", "CK1"],
    })


def test_train_bin_edges():
    b = ST.train_bin([0, 1, 4, 5, 19, 20, 99, 100, 5000])
    assert list(b) == ["0", "1-4", "1-4", "5-19", "5-19", "20-99", "20-99", ">=100", ">=100"]


def test_priming_curated_offsets():
    p = ST.priming_curated(_kinases())
    assert p == {"K1": (4,), "K2": (-3,), "K3": (1, 3), "K4": (), "K11": (-3,)}


def _toy_matrices():
    cols = [f"{p}{a}" for p in ST._ST_POS for a in ST._AA23]
    M = pd.DataFrame(0.0, index=["GOOD", "NOTTOP", "WEAK", "OTHERIDX"], columns=cols)
    M.loc["GOOD", "4s"] = 3.0                   # pS at +4 is the top residue, >= 2 -> qualifies (+4)
    M.loc["GOOD", "-3t"] = 2.0                  # pT at -3 exactly 2.0 and top -> qualifies (-3)
    M.loc["NOTTOP", "-3s"] = 2.5
    M.loc["NOTTOP", "-3D"] = 3.0                # acidic beats pS -> not top -> no
    M.loc["WEAK", "1s"] = 1.5                   # below threshold
    M.loc["OTHERIDX", "4s"] = 5.0               # qualifies but accession not in the kinase index
    info = pd.DataFrame({"MATRIX_NAME": list(M.index), "UNIPROT_ID": ["K9", "K5", "K6", "ZZZ"]})
    return M, info


def test_cantley_phospho_preference_rule():
    M, info = _toy_matrices()
    d = ST.cantley_phospho_preference(matrices=M, info=info).set_index("matrix")
    assert d.loc["GOOD", "offsets"] == (-3, 4) and d.loc["GOOD", "qualifies"]
    assert not d.loc["NOTTOP", "qualifies"] and not d.loc["WEAK", "qualifies"]
    p = ST.priming_cantley(_kinases(), matrices=M, info=info)
    assert p == {"K9": (-3, 4)}                 # OTHERIDX dropped (not in index)


def test_kinase_table_columns():
    M, info = _toy_matrices()
    can = ST.priming_cantley(_kinases(), matrices=M, info=info)
    T = ST.kinase_table(_kinases(), train_counts=pd.Series({"K1": 300, "K5": 3}), cantley=can).set_index("gene")
    assert T.loc["GSK3B", "train_bin"] == ">=100" and T.loc["ATM", "train_bin"] == "1-4"
    assert T.loc["SRC", "train_bin"] == "0"
    assert set(T.index[T.pikk]) == {"ATM", "MTOR", "TRRAP"}
    assert set(T.index[T.pikk_sq]) == {"ATM"}
    assert set(T.index[T.secretory]) == {"FAM20C"} and set(T.index[T.fam20c]) == {"FAM20C"}
    assert set(T.index[T.priming]) == {"GSK3B", "CSNK1D", "CSNK2A1", "PLK1", "CSNK1G2"}
    assert set(T.index[T.priming_cantley]) == {"PRKACA"}
    assert T.loc["CSNK2A1", "priming_offsets"] == "+1,+3"
    # every kinase-level stratum is a partition of the kinases
    S = ST.kinase_sets(T.reset_index())
    for c in ST.KINASE_STRATA:
        parts = [v for (s, _), v in S.items() if s == c]
        assert sum(len(p) for p in parts) == len(T) and set().union(*parts) == set(T.kinase_acc)


def _sites():
    # protein P1 sequence positions 1..20: 'MSAAPSAAASAAAPEESTYG'
    seq = "MSAAPSAAASAAAPEESTYG"
    rows = []
    for i, r in enumerate(seq, start=1):
        if r in "STY":
            w = "".join(seq[j - 1] if 1 <= j <= len(seq) else "_" for j in range(i - 7, i + 8))
            rows.append({"site_id": f"P1_{r}{i}", "acc": "P1", "res": r, "pos": i, "window15": w})
    S = pd.DataFrame(rows)
    S["in_U_ms"] = S.site_id.isin(["P1_S2", "P1_S6", "P1_S10", "P1_S17"])
    return S


def test_offset_site_ids():
    S = _sites().set_index("site_id")
    o = ST.offset_site_ids(["P1_S6", "P1_S6", "P1_S2", "P1_S10"],
                           [S.loc["P1_S6", "window15"]] * 2 + [S.loc["P1_S2", "window15"], S.loc["P1_S10", "window15"]],
                           0)
    assert list(o) == ["P1_S6", "P1_S6", "P1_S2", "P1_S10"]
    o4 = ST.offset_site_ids(["P1_S6", "P1_S2"], [S.loc["P1_S6", "window15"], S.loc["P1_S2", "window15"]], 4)
    assert list(o4) == ["P1_S10", "P1_S6"]
    om = ST.offset_site_ids(["P1_S2"], [S.loc["P1_S2", "window15"]], -3)   # outside the protein ('_')
    assert om[0] is None


def test_primed_mask():
    S = _sites()
    offs = {"GSK": (4,), "CK1": (-3,), "PLK": ()}
    pos = pd.DataFrame({"kinase_acc": ["GSK", "GSK", "GSK", "CK1", "CK1", "PLK", "OTHER"],
                        "site_id": ["P1_S6", "P1_S2", "P1_S10", "P1_S17", "P1_S10", "P1_S6", "P1_S6"]})
    m = ST.primed_mask(pos, offs, sites=S)
    # GSK S6: +4 = S10 observed -> primed; GSK S2: +4 = S6 observed -> primed; GSK S10: +4 = A14 -> no
    # CK1 S17: -3 = P14 -> no; CK1 S10: -3 = A7 -> no; PLK has no offset; OTHER not a priming kinase
    assert list(m) == [True, True, False, False, False, False, False]
    # +4 residue exists and is S but not observed -> not primed
    S2 = S.assign(in_U_ms=S.site_id.isin(["P1_S6"]))
    assert list(ST.primed_mask(pos.iloc[:1], offs, sites=S2)) == [False]
    # CK1 at T18: -3 = E15 -> no; at S10 with -4 = S6 using a custom offset -> yes
    assert list(ST.primed_mask(pd.DataFrame({"kinase_acc": ["X"], "site_id": ["P1_S10"]}), {"X": (-4,)}, sites=S)) \
        == [True]


def test_hub_sites_counts_paralog_components():
    par = {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}, "D": set(), "E": set(), "F": set(), "G": set()}
    rows = [("s1", k) for k in "ABCDEF"]          # A-B-C one component (transitive) -> groups {ABC},D,E,F = 4
    rows += [("s2", k) for k in "ADEFG"]           # 5 groups -> hub
    rows += [("s3", k) for k in "DE"]
    rows += [("s2", "CANARY_NEG_1")]               # canaries never count
    pairs = pd.DataFrame(rows, columns=["site_id", "kinase_acc"])
    g = ST.site_kinase_groups(pairs, par)
    assert g["s1"] == 4 and g["s2"] == 5 and g["s3"] == 2
    assert ST.hub_sites(pairs, par) == {"s2"}
    assert ST.hub_sites(pairs, par, min_groups=4) == {"s1", "s2"}


def test_evidence_class():
    e = ST.evidence_class(pd.DataFrame({"evidence_native": ["both", "in_vivo_only", "none", "in_vitro_only"]}))
    assert list(e) == ["both", "in_vivo_only", "other", "in_vitro_only"]


def test_definitions_table():
    d = ST.definitions_table()
    assert {"priming", "priming_cantley", "pikk", "hub", "primed", "train_bin"} <= set(d.stratum)


def test_matches_eval_harness_on_real_index():
    """ksp.strata 'priming'/'pikk'/train_bin must equal the ksp.eval per-kinase attributes (same definitions)."""
    from ksp import eval as E
    from ksp import index
    try:
        K = index.load_kinases()
    except FileNotFoundError:
        pytest.skip("kinase index not built")
    T = ST.kinase_table(K, train_counts=pd.Series(dtype=float), cantley={})
    ev_prim = K.gene.isin(E.PRIMING) | K.manning_family.eq("CK1")
    assert np.array_equal(T.priming.values, ev_prim.values)
    assert np.array_equal(T.pikk.values, K.manning_family.eq("PIKK").values)
    assert ST.TRAIN_BINS == E.TRAIN_BINS and ST.TRAIN_BIN_LABELS == E.TRAIN_BIN_LABELS
