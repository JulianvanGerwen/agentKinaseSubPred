"""Toy-sequence tests for src/ksp/ids.py. Run: /workspace/.venv/bin/python -m pytest -q tests/test_ids.py
(or plain `python tests/test_ids.py`)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import ids  # noqa: E402

#            1234567890123456789012345678901234567890
SEQ_A = "MSTAKLRRASVGLPQESTWYLKKAAEGHHRRASVDDMPQT"  # 40 aa; 'RRASV' occurs twice (pos 10 S, pos 33 S)
SEQ_B = "MAAAAAAAAAYKKKKKKKKKKTEEEEE"               # 27 aa
SEQ_G = "MPPPPPPPPPPSQQQQQQQQQQ"                    # gene GENEG (Swiss-Prot), source uses TrEMBL acc


def _store():
    return ids.CanonicalStore(
        seqs={"P00001": SEQ_A, "P00002": SEQ_B, "P00003": SEQ_G},
        gene_to_accs={"GENEA": ["P00001"], "GENEB": ["P00002"], "GENEG": ["P00003"]},
        acc_to_gene={"P00001": "GENEA", "P00002": "GENEB", "P00003": "GENEG", "A0A000TR01": "GENEG"})


def test_site_id_roundtrip():
    s = ids.make_site_id("P05198", "s", 52)
    assert s == "P05198_S52"
    assert ids.parse_site_id(s) == ("P05198", "S", 52)
    for bad in ["P05198-2_S52", "P05198_K52", "P05198_S", "P05198S52"]:
        try:
            ids.parse_site_id(bad)
            raise AssertionError(bad)
        except ValueError:
            pass
    try:
        ids.make_site_id("P05198-2", "S", 52)
        raise AssertionError
    except ValueError:
        pass


def test_windows_padding():
    assert ids.window(SEQ_A, 2) == "______MSTAKLRRA"
    assert len(ids.window(SEQ_A, 2, 15)) == 31 and ids.window(SEQ_A, 2, 15).startswith("_" * 14 + "MS")
    w = ids.window(SEQ_A, 40)
    assert w == "SVDDMPQT" + "_" * 7
    assert ids.window(SEQ_A, 10) == "TAKLRRASVGLPQES"


def test_exact_at_stated_position():
    st = _store()
    m = st.map_window("P00001", 10, "TAKLRRASVGLPQES", "S")
    assert (m.method, m.site_id) == ("exact_pos", "P00001_S10")
    # padded window near N-terminus, lower-case input
    m = st.map_window("P00001", 2, "______mStaklrra", "S")
    assert m.site_id == "P00001_S2" and m.method == "exact_pos"
    # source padded where canonical has a residue (initiator Met removed in source numbering) -> wildcard
    m = st.map_window("P00001", 1, "_______STAKLRRA", "S")
    assert m.method == "shifted" and m.site_id == "P00001_S2"


def test_shifted_unique_and_ambiguous():
    st = _store()
    # correct window for S33 stated at wrong position 30 -> unique elsewhere
    m = st.map_window("P00001", 30, "EGHHRRASVDDMPQT", "S")
    assert (m.method, m.site_id) == ("shifted", "P00001_S33")
    # short-core window present twice ('RRASV' with padding wildcards) -> ambiguous
    m = st.map_window("P00001", 99, "____RRASV______", "S")
    assert m.method == "unmapped" and m.reason == "multiple_matches" and m.n_hits == 2
    # isoform accession: same window found in canonical
    m = st.map_window("P00001-2", 5, "EGHHRRASVDDMPQT", "S")
    assert m.site_id == "P00001_S33"


def test_unmapped_reasons():
    st = _store()
    assert st.map_window("P00001", 10, "TAKLRRASVGLPQE", "S").reason == "bad_window"   # 14-mer
    assert st.map_window("P00001", 10, "TAKLRRASVGLPQES", "T").reason == "bad_window"  # centre != residue
    assert st.map_window("Q99999", 10, "TAKLRRASVGLPQES", "S").reason == "no_candidate"
    assert st.map_window("P00001", 10, "TAKLRRWSVGLPQES", "S").reason == "no_match"


def test_gene_fallback_for_trembl_and_explicit_gene():
    st = _store()
    w = ids.window(SEQ_G, 12)
    m = st.map_window("A0A000TR01", 12, w, "S")         # TrEMBL acc -> gene GENEG -> P00003
    assert m.site_id == "P00003_S12" and m.method == "exact_pos"
    m = st.map_window("OBSOLETE1", 3, w, "S", gene="geneg")
    assert m.site_id == "P00003_S12" and m.method == "shifted"


def test_map_position():
    st = _store()
    assert st.map_position("P00002", 11, "Y").site_id == "P00002_Y11"
    assert st.map_position("P00002", 11, "S").reason == "residue_mismatch"
    assert st.map_position("P00002", 11, "Y", expected_length=30).reason == "length_changed"
    assert st.map_position("P00002", 99, "Y").reason == "pos_out_of_range"
    assert st.map_position("P00002-3", 11, "Y").reason == "isoform_accession"


def test_map_table_and_windows():
    st = _store()
    df = pd.DataFrame({"acc": ["P00001", "P00002", "P00001"], "pos": [10, 11, 99],
                       "window": ["TAKLRRASVGLPQES", ids.window(SEQ_B, 11), "XXXXXXXSXXXXXXX"]})
    o = ids.map_table(st, df)
    assert o.site_id.tolist()[:2] == ["P00001_S10", "P00002_Y11"]
    assert o.map_reason.iloc[2] == "multiple_matches"
    w = ids.add_windows(st, o.dropna(subset=["map_acc"]), "map_acc", "map_pos")
    assert w.window15.iloc[0] == "TAKLRRASVGLPQES" and len(w.window31.iloc[1]) == 31


def test_add_gene_synonyms():
    st = _store()
    assert st.add_gene_synonyms([("OLDA", "P00001"), ("GENEB", "P00001"), ("X", "NOTTHERE")]) == 1
    assert st.candidates(None, "olda") == ["P00001"]
    assert st.candidates(None, "GENEB") == ["P00002"]


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v()
            print("ok", k)
