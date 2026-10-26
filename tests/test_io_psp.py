"""Toy-data tests for src/ksp/io_psp.py (known answers). Run: PYTHONPATH=src python -m pytest -q tests/test_io_psp.py
(or plain `python tests/test_io_psp.py`)."""
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import io_psp  # noqa: E402

HDR = ("GENE\tKINASE\tKIN_ACC_ID\tKIN_ORGANISM\tSUBSTRATE\tSUB_GENE_ID\tSUB_ACC_ID\tSUB_GENE\tSUB_ORGANISM\t"
       "SUB_MOD_RSD\tSITE_GRP_ID\tSITE_+/-7_AA\tDOMAIN\tIN_VIVO_RXN\tIN_VITRO_RXN\tCST_CAT#")
ROWS = [
    # human kinase -> human substrate, both evidence, isoform substrate
    "AKT1\tAkt1\tP31749\thuman\tGSK3B\t2932\tP49841-2\tGSK3B\thuman\tS9\t100\tSGRPRTtsFAEsCKP\t\tX\tX\t",
    # mouse kinase (symbol Akt1) -> mouse substrate in same site group 100 -> maps to human P49841-2 S9
    "Akt1\tAkt1\tP31750\tmouse\tGsk3b\t56637\tQ9WV60\tGsk3b\tmouse\tS9\t100\tSGRPRTtsFAEsCKP\t\t \tX\t",
    # rat kinase, site group with no human member -> unmapped by SGID
    "Prkaca\tPKACA\tP27791\trat\tFoo\t1\tQ00001\tFoo\trat\tT4\t200\t____MAPtLLKKAAE\t\tX\t \t",
    # junk non-STY row
    "NME1\tNME1\tP15531\thuman\tNME1\t4830\tP15531\tNME1\thuman\tH118\t300\tQVGRNIIhGsDsVEs\tNDK\t \tX\t",
    # autophosphorylation
    "AKT1\tAkt1\tP31749\thuman\tAKT1\t207\tP31749\tAKT1\thuman\tT308\t400\tDGATMKtFCGTPEYL\t\tX\t \t",
]


def _toy_file():
    f = tempfile.NamedTemporaryFile("w", delete=False, encoding="latin1", suffix=".tsv")
    f.write("051925\nlicence text\n\n" + HDR + "\n" + "\n".join(ROWS) + "\n")
    f.close()
    return f.name


def test_parse():
    d = io_psp.load_psp(_toy_file())
    assert len(d) == 5
    r = d.iloc[0]
    assert (r.sub_acc, r.sub_isoform, r.residue, r.position) == ("P49841", "2", "S", 9)
    assert r.evidence == "both" and r.window == "SGRPRTTSFAESCKP"
    assert r.other_phos_offsets == (-1, 4)  # lower-case t at -1 and s at +4
    assert d.iloc[1].evidence == "in_vitro_only" and d.iloc[2].evidence == "in_vivo_only"
    assert d.is_sty.tolist() == [True, True, True, False, True]
    assert d.autophos.tolist() == [False, False, False, True, True]


def test_mapping():
    d = io_psp.filter_psp(io_psp.load_psp(_toy_file()))
    k = io_psp.map_kinases_to_human(d)
    assert k.kin_human_mapped.tolist() == [True, True, False, True]  # rat PRKACA has no human PSP kinase row
    m = io_psp.map_substrates_to_human(k)
    mm = m[m.sub_organism == "mouse"].iloc[0]
    assert mm.h_sub_acc_full == "P49841-2" and mm.h_position == 9
    assert not m[m.sub_organism == "rat"].iloc[0].sub_human_mapped


def test_window_tools():
    seq = "MAPTLLKKAAEEGGSRRSPQ"
    hit = io_psp._best_window_hit("____MAPTLLKKAAE", "T", seq)
    assert hit == (4, 11, 11)  # T at position 4, all 11 non-padding residues identical
    sites = pd.DataFrame({"sub_acc_full": ["X", "X", "Y"], "position": [4, 5, 1], "window": ["____MAPTLLKKAAE"] * 3})
    c = io_psp.check_windows_against_sequences(sites, {"X": seq})
    assert c.status.tolist() == ["exact", "centre_mismatch", "no_sequence"]
    assert c.shifted_pos.iloc[1] == 4


if __name__ == "__main__":
    test_parse(); test_mapping(); test_window_tools()
    print("all io_psp toy tests passed")
