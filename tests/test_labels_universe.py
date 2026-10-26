"""Real-data checks of the built kinase list, site universe and labels (hand-checked PSP rows -> expected site_ids).
Requires the outputs of scripts/build_kinase_list.py, build_site_universe.py and build_labels.py; tests are skipped
when they are missing. Run: /workspace/.venv/bin/python tests/test_labels_universe.py"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ksp import ids, kinases, labels  # noqa: E402

KL = ROOT / "results/labels/kinase_list.tsv"
UNI = ROOT / "cache/sites/universe.parquet"
KS = ROOT / "cache/labels/ks_pairs.parquet"


def _skip(p):
    if not p.exists():
        print(f"SKIP: {p} missing")
        return True
    return False


def test_kinase_list_and_synonyms():
    if _skip(KL):
        return
    r = kinases.load_resolver()
    kl = r.kl.set_index("gene")
    for old, new in [("ADRBK1", "GRK2"), ("ZAK", "MAP3K20"), ("MST4", "STK26"), ("KIAA0999", "SIK3"),
                     ("CABC1", "COQ8A"), ("GSG2", "HASPIN"), ("PKACa", "PRKACA"), ("ERK2", "MAPK1"), ("LKB1", "STK11")]:
        assert r.resolve(old)[0] == new, (old, r.resolve(old))
    assert r.resolve("PDK1")[0] == "PDK1"            # own symbol beats PDPK1 alias
    assert r.resolve(acc="O15530")[0] == "PDPK1"
    for fam in ["PKA", "CK2", "PKC", "AMPK"]:
        assert r.resolve_note(fam) == (None, "family_name")
    assert r.resolve_note("PKB/AKT1")[0] == "AKT1" and r.resolve_note("GSK3-beta")[0] == "GSK3B"
    for g in ["PKM", "PGK1", "WHR1", "PIK3CA", "CSNK2B", "BCR/ABL", "BRD4", "TAF1", "BAZ1B", "FAM20B"]:
        assert not kl.loc[g, "include"], g
    for g in ["FAM20C", "ATM", "MTOR", "PRKDC", "PDK1", "BCKDK", "TRPM7", "EEF2K", "RIOK1", "ALPK1", "COQ8A"]:
        assert kl.loc[g, "include"], g
    assert kl.loc["COQ8A", "activity_uncertain"] and kl.loc["FAM20A", "is_pseudokinase"]
    assert kl.loc["MAP2K1", "residue_class"] == "STY" and kl.loc["SRC", "residue_class"] == "Y"
    assert kl.loc["WEE1", "residue_class"] == "Y" and kl.loc["PRKACA", "residue_class"] == "ST"
    assert kl.loc["JAK2", "kinase_domain_name"] == "Protein kinase 2"  # JH1 (catalytic) chosen
    assert kl.loc["ATM", "kinase_domain_start"] == 2686
    assert kl.loc["SRC", "manning_group"] == "TK" and kl.loc["CDK1", "manning_family"] == "CDK"


def test_universe():
    if _skip(UNI):
        return
    u = pd.read_parquet(UNI)
    assert u.site_id.is_unique and u.res.isin(list("STY")).all()
    assert (u.window15.str[7] == u.res).all() and (u.window31.str[15] == u.res).all()
    assert (u.window31.str[8:23] == u.window15).all()
    r = u.set_index("site_id")
    assert r.loc["P05198_S52", "window15"] == "MILLSELSRRRIRSI" and r.loc["P05198_S52", "src_psp"]
    # PSP numbers histone H3 without the initiator Met: PSP S10 -> canonical S11
    assert r.loc["P68431_S11", "src_psp"] and r.loc["P68431_S11", "psp_map_method"] == "shifted"
    for sid in u.site_id.sample(200, random_state=0):
        a, res, pos = ids.parse_site_id(sid)
        assert ids.make_site_id(a, res, pos) == sid


def test_labels_known_pairs():
    if _skip(KS):
        return
    p = pd.read_parquet(KS)
    assert not p.duplicated(["kinase", "site_id", "provenance"]).any()
    psp = p[p.provenance.eq("PSP")].set_index(["kinase", "site_id"])
    # human-native rows (PSP: AKT1 -> GSK3B S9; PDPK1 -> AKT1 T308; EIF2AK3 -> EIF2S1 S52)
    for k, s in [("AKT1", "P49841_S9"), ("PDPK1", "P31749_T308"), ("EIF2AK3", "P05198_S52")]:
        r = psp.loc[(k, s)]
        assert r.L1 and r.L2 and r.native_human and r.in_vivo and r.in_vitro, (k, s)
    # human AURKB -> H3C1 'S10' (PSP numbering without Met) -> P68431_S11
    assert psp.loc[("AURKB", "P68431_S11"), "L1"]
    # CDK1 -> MAPT isoform P10636-8 T231 (2N4R tau numbering) -> canonical P10636 T548 (window relocated)
    r = psp.loc[("CDK1", "P10636_T548")]
    assert "T231" in r.psp_mod_rsd and r.L1
    # rat Camk2a -> rat Dlgap1 S389 (SGID 481511) -> human DLGAP1 O14490 S383: rodent-transferred (L2, not L1)
    r = psp.loc[("CAMK2A", "O14490_S383")]
    assert r.L2 and not r.L1 and r.species_origin == "rat>rat" and r.sub_map_methods == "sgid"


def test_label_variant_logic():
    if _skip(KS):
        return
    p = pd.read_parquet(KS)
    assert (p.L1 <= p.L2).all() and (p.L2 <= p.L3).all() and (p.L2_vitro <= p.L2).all()
    assert (p.L2_noauto <= p.L2).all() and not (p.L2_noauto & p.autophos).any()
    v = p.L1 | p.L2 | p.L3 | p.L2_vitro | p.L2_noauto | p.EXT
    assert (p.loc[v, "kinase_included"] & p.loc[v, "residue_in_class"]).all()
    assert not p.loc[p.provenance.ne("PSP"), ["L1", "L2", "L3"]].any().any()
    assert not p.loc[p.provenance.eq("PSP"), "EXT"].any()
    # residue class: no S/T labels for a Y-class kinase in any variant
    y = p[v & p.kinase_residue_class.eq("Y")]
    assert y.res.eq("Y").all()
    # OmniPath rows always carry >= 1 curated non-PSP source
    o = p[p.provenance.eq("OmniPath")]
    assert o.ext_sources.fillna("").ne("").all()
    assert all(set(s.split(";")) <= set(labels.OMNIPATH_CURATED + labels.OMNIPATH_OTHER_AGGREGATORS)
               for s in o.ext_sources)


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f()
            print("ok", k)
