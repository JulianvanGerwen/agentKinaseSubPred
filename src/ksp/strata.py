"""ksp.strata - canonical stratum definitions (DESIGN §6.4, §9, §17 A26). Owner: strata analyst.

Used on CV now (scripts/strata_cv_report.py) and later on TEST (ksp.eval_test). Every function takes DEV-safe inputs
only: the kinase index, the site index (window15, U_ms flag), Cantley matrices shipped with kinase-library 1.8.0,
DEV L2 training counts (ksp.splits.load_eval_kinase_sets().n_dev_L2) and, for the label-defined site strata (hub
sites, evidence class), a label table the CALLER passes in (CV: ksp.labels.load_dev('all'); TEST: whatever
ksp.eval_test already holds). Nothing here reads TEST labels itself.

Kinase-level strata (kinase_table(); one row per index kinase):
  manning_group     Manning group from cache/index/kinases.parquet ('(none)' if missing; same as ksp.eval).
  residue_class     ST / Y / STY (config/kinase_curation.tsv overrides; same as ksp.eval).
  train_bin         DEV L2 training-label count binned 0 / 1-4 / 5-19 / 20-99 / >=100 (= ksp.eval TRAIN_BINS).
  priming           PRIMARY priming stratum (curated, pre-registered in DESIGN §9; identical to ksp.eval 'priming'):
                    GSK3A/B, the Manning CK1 family (CSNK1A1/A1L/D/E/G1/G2/G3), CK2 alpha (CSNK2A1/A2/A3), PLK1.
  priming_cantley   SECONDARY, data-driven: S/T kinases whose kinase-library 1.8.0 log2 matrix has, at some scored
                    position p in -5..-1,+1..+4, a phospho-residue (pS or pT) that is the single most favoured of the
                    23 residues at p with log2 enrichment >= 2.0 (>= 4-fold over random). Rule fixed before any model
                    metric was looked at (docs/strata_cv.md §1). Tyrosine matrices are NOT used (many TK matrices
                    favour pY at -1/+1/+2, which reflects adjacent multisite phosphorylation, not priming).
  pikk              Manning family PIKK (ATM, ATR, PRKDC, MTOR, SMG1, TRRAP; same as ksp.eval 'pikk').
  pikk_sq           PIKKs with the S/T-Q motif (ATM, ATR, PRKDC, SMG1; i.e. pikk minus MTOR and pseudokinase TRRAP).
  secretory         Golgi/secretory-pathway kinases in the index: FAM20C, FAM20A (pseudokinase), PKDCC (VLK).
  fam20c            FAM20C alone.

Site/pair-level strata (masks over a positives table with columns kinase_acc, site_id):
  primed            positive (k, s) with k a priming kinase and the residue at the priming offset of k being an S/T
                    that is itself an observed phosphosite (U_ms). Offsets: GSK3 +4; CK1 -3; CK2 +1 or +3
                    (pS substitutes the acidic +1/+3 determinants; Meggio & Pinna 2003); PLK1: none (its priming is
                    PBD docking on a distal S-pT-P, not a fixed offset), so PLK1 positives are never 'primed'.
                    Cantley variant: offsets = the positions that qualified the kinase for priming_cantley.
  hub               site carrying evaluable-positive labels (human-native PSP, non-autophos) to >= 5 kinases that are
                    mutually non-paralogous, counted as connected components of the Par(k) graph among the site's
                    kinases (A26; e.g. GLO1 Y136).
  evidence          'in_vitro_only' / 'in_vivo_only' / 'both' (labels column evidence_native; same as ksp.eval).

Sources/justification for the curated priming list: GSK3 requires a phosphate at +4 (Fiol et al. 1987 JBC; Frame &
Cohen 2001 Biochem J); CK1 prefers pS/pT at -3 (Flotow et al. 1990 JBC; Knippschild 2005); CK2 is acidophilic and
pS/pT can replace the +1/+3 acidic determinants (Meggio & Pinna 2003 FASEB J); PLK1 substrate recruitment by the
polo-box domain requires a primed S-pT-P (Elia et al. 2003 Cell). DESIGN §9 lists exactly these four groups.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import index

ROOT = Path(__file__).resolve().parents[2]
KL_MATRICES = ROOT / "external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases"

TRAIN_BINS = [-1, 0, 4, 19, 99, 1e12]                      # identical to ksp.eval.TRAIN_BINS
TRAIN_BIN_LABELS = ["0", "1-4", "5-19", "20-99", ">=100"]

# curated priming kinases (primary): gene -> priming offsets relative to the acceptor (residue at pos+offset)
PRIMING_GENES = {"GSK3A": (4,), "GSK3B": (4,),
                 "CSNK2A1": (1, 3), "CSNK2A2": (1, 3), "CSNK2A3": (1, 3),
                 "PLK1": ()}
CK1_FAMILY_OFFSETS = (-3,)                                  # Manning family 'CK1' (CSNK1*), resolved from the index
PIKK_SQ = {"ATM", "ATR", "PRKDC", "SMG1"}
SECRETORY = {"FAM20C", "FAM20A", "PKDCC"}
CANTLEY_MIN_LOG2 = 2.0
HUB_MIN_GROUPS = 5
EVIDENCE_CLASSES = ("in_vitro_only", "in_vivo_only", "both")
_AA23 = "PGACSTVILMFYWHKRQNDEsty"
_ST_POS = [-5, -4, -3, -2, -1, 1, 2, 3, 4]


# ============================================================================================== kinase-level strata
def train_bin(counts) -> pd.Series:
    """Bin DEV L2 training-label counts (0 / 1-4 / 5-19 / 20-99 / >=100), as ksp.eval does."""
    c = pd.Series(counts).fillna(0)
    return pd.cut(c, TRAIN_BINS, labels=TRAIN_BIN_LABELS).astype(str)


def priming_curated(kinases: pd.DataFrame | None = None) -> dict:
    """Primary priming kinases: kinase_acc -> tuple of priming offsets (PLK1: empty tuple)."""
    K = kinases if kinases is not None else index.load_kinases()
    out = {}
    for a, g, fam in zip(K.uniprot_acc, K.gene, K.manning_family):
        if g in PRIMING_GENES:
            out[a] = PRIMING_GENES[g]
        elif fam == "CK1":
            out[a] = CK1_FAMILY_OFFSETS
    return out


def cantley_phospho_preference(min_log2: float = CANTLEY_MIN_LOG2, matrices: pd.DataFrame | None = None,
                               info: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per S/T Cantley matrix: max phospho-residue log2 enrichment, the qualifying offsets (pS/pT is the argmax of the
    23 residues at that position and >= min_log2) and the UniProt accession (kinome_information UNIPROT_ID; keyed by
    accession, never by name - LESSONS 2). `matrices`/`info` may be passed for tests (same layout as the shipped
    ser_thr_all_log2_matrices.tsv / kinome_information.tsv)."""
    if matrices is None:
        matrices = pd.read_csv(KL_MATRICES / "matrices/ser_thr/ser_thr_all_log2_matrices.tsv", sep="\t",
                               index_col=0)
    if info is None:
        info = pd.read_csv(KL_MATRICES / "kinase_data/kinome_information.tsv", sep="\t", dtype=str)
    info = info.set_index("MATRIX_NAME")
    cols = [f"{p}{a}" for p in _ST_POS for a in _AA23]
    L = matrices[cols].to_numpy(float).reshape(len(matrices), len(_ST_POS), len(_AA23))
    ph = [_AA23.index("s"), _AA23.index("t")]
    mph = L[:, :, ph].max(2)
    is_top = np.isin(L.argmax(2), ph)
    qual = is_top & (mph >= min_log2)
    rows = []
    for i, name in enumerate(matrices.index):
        offs = tuple(int(p) for p, q in zip(_ST_POS, qual[i]) if q)
        rows.append({"matrix": name, "uniprot_acc": info.UNIPROT_ID.get(name), "max_phospho_log2": float(mph[i].max()),
                     "best_offset": int(_ST_POS[int(mph[i].argmax())]), "offsets": offs, "qualifies": bool(offs)})
    return pd.DataFrame(rows)


def priming_cantley(kinases: pd.DataFrame | None = None, min_log2: float = CANTLEY_MIN_LOG2, **kw) -> dict:
    """Secondary (data-driven) priming kinases: kinase_acc -> offsets favoured as pS/pT by the Cantley matrix.
    Only ST/STY kinases of the index qualify; matrices whose accession is not in the index are ignored."""
    K = kinases if kinases is not None else index.load_kinases()
    ok = set(K.uniprot_acc[K.residue_class.isin(["ST", "STY"])])
    pref = cantley_phospho_preference(min_log2, **kw)
    out = {}
    for a, offs in zip(pref.uniprot_acc, pref.offsets):
        if offs and a in ok:
            out[a] = tuple(sorted(set(out.get(a, ())) | set(offs)))
    return out


def kinase_table(kinases: pd.DataFrame | None = None, train_counts: pd.Series | None = None,
                 cantley: dict | None = None) -> pd.DataFrame:
    """One row per kinase with every kinase-level stratum column. train_counts: kinase_acc -> DEV L2 count
    (default: ksp.splits.load_eval_kinase_sets().n_dev_L2 - the harness's train_counts)."""
    K = (kinases if kinases is not None else index.load_kinases()).copy()
    if train_counts is None:
        from . import splits
        ks = splits.load_eval_kinase_sets()
        train_counts = pd.Series(ks.n_dev_L2.values, index=ks.kinase_acc.values)
    cur = priming_curated(K)
    can = cantley if cantley is not None else priming_cantley(K)
    T = pd.DataFrame({"kinase_acc": K.uniprot_acc.values, "gene": K.gene.values})
    T["manning_group"] = K.manning_group.fillna("(none)").values
    T["manning_family"] = K.manning_family.values
    T["residue_class"] = K.residue_class.values
    T["n_train_dev_L2"] = T.kinase_acc.map(train_counts).fillna(0).astype(int)
    T["train_bin"] = train_bin(T.n_train_dev_L2).values
    T["priming"] = T.kinase_acc.isin(cur)
    T["priming_offsets"] = T.kinase_acc.map(lambda a: ",".join(f"{o:+d}" for o in cur.get(a, ())))
    T["priming_cantley"] = T.kinase_acc.isin(can)
    T["priming_cantley_offsets"] = T.kinase_acc.map(lambda a: ",".join(f"{o:+d}" for o in can.get(a, ())))
    T["pikk"] = K.manning_family.eq("PIKK").values
    T["pikk_sq"] = T.gene.isin(PIKK_SQ).values & T.pikk.values
    T["secretory"] = T.gene.isin(SECRETORY).values
    T["fam20c"] = T.gene.eq("FAM20C").values
    return T


KINASE_STRATA = ("manning_group", "residue_class", "train_bin", "priming", "priming_cantley", "pikk", "pikk_sq",
                 "secretory", "fam20c")


def kinase_sets(table: pd.DataFrame | None = None) -> dict:
    """{(stratum, value): set(kinase_acc)} for every kinase-level stratum value (booleans as 'True'/'False')."""
    T = table if table is not None else kinase_table()
    out = {}
    for c in KINASE_STRATA:
        for v, g in T.groupby(c):
            out[(c, str(v))] = set(g.kinase_acc)
    return out


# ============================================================================================== site/pair-level strata
def offset_site_ids(site_ids, windows15, offset: int) -> np.ndarray:
    """site_id of the residue at acceptor position + offset (None if outside the protein or not S/T)."""
    out = np.empty(len(site_ids), dtype=object)
    for i, (sid, w) in enumerate(zip(site_ids, windows15)):
        acc, rp = sid.rsplit("_", 1)
        pos = int(rp[1:])
        j = 7 + offset
        r = w[j] if 0 <= j < len(w) else "_"
        out[i] = f"{acc}_{r}{pos + offset}" if r in ("S", "T") else None
    return out


def primed_mask(pos: pd.DataFrame, offsets: dict, sites: pd.DataFrame | None = None,
                observed: set | None = None) -> np.ndarray:
    """Bool per row of pos[kinase_acc, site_id]: True iff the kinase has priming offsets and, for at least one
    offset, the residue there is S/T and an observed phosphosite (`observed`, default: U_ms site_ids of the index).
    Rows of kinases without offsets are False."""
    S = sites if sites is not None else index.load_sites()
    if observed is None:
        observed = set(S.site_id[S.in_U_ms.astype(bool)])
    win = pd.Series(S.window15.values, index=S.site_id.values)
    out = np.zeros(len(pos), bool)
    ka, sid = pos.kinase_acc.values, pos.site_id.values
    for off in sorted({o for v in offsets.values() for o in v}):
        sel = np.array([off in offsets.get(k, ()) for k in ka], bool) & ~out
        if not sel.any():
            continue
        oid = offset_site_ids(sid[sel], win.reindex(sid[sel]).fillna("").values, off)
        hit = np.array([o is not None and o in observed for o in oid], bool)
        idx = np.where(sel)[0]
        out[idx[hit]] = True
    return out


def _components(nodes: list, par: dict) -> int:
    nodes = list(nodes)
    seen, n = set(), 0
    ns = set(nodes)
    for a in nodes:
        if a in seen:
            continue
        n += 1
        stack = [a]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend((par.get(x, set()) & ns) - seen)
    return n


def site_kinase_groups(pairs: pd.DataFrame, paralogs: dict) -> pd.Series:
    """site_id -> number of mutually non-paralogous kinase groups (Par-graph components) labelled to the site.
    `pairs` must already be the label rows to count (e.g. evaluable positives; canaries are dropped here)."""
    p = pairs[~pairs.kinase_acc.astype(str).str.startswith("CANARY_")].drop_duplicates(["kinase_acc", "site_id"])
    return p.groupby("site_id").kinase_acc.agg(lambda s: _components(sorted(set(s)), paralogs))


def hub_sites(pairs: pd.DataFrame, paralogs: dict, min_groups: int = HUB_MIN_GROUPS) -> set:
    """A26 hub sites: labelled (rows of `pairs`) to >= min_groups mutually non-paralogous kinases."""
    g = site_kinase_groups(pairs, paralogs)
    return set(g.index[g >= min_groups])


def evidence_class(pos: pd.DataFrame) -> pd.Series:
    """Evidence class per positive pair (column evidence_native), restricted to EVIDENCE_CLASSES (else 'other')."""
    e = pos["evidence_native"].astype(str)
    return e.where(e.isin(EVIDENCE_CLASSES), "other")


# ============================================================================================== definitions table
DEFINITIONS = [
    ("manning_group", "kinase", "Manning group (cache/index/kinases.parquet)", "KinHub/Manning 2002"),
    ("residue_class", "kinase", "ST / Y / STY; config/kinase_curation.tsv overrides", "docs/core_infra.md §1"),
    ("train_bin", "kinase", "DEV L2 label count: 0 / 1-4 / 5-19 / 20-99 / >=100", "DESIGN §6.3"),
    ("priming", "kinase", "PRIMARY. GSK3A/B; Manning family CK1 (CSNK1*); CSNK2A1/A2/A3; PLK1 (= ksp.eval 'priming')",
     "DESIGN §9; Fiol 1987, Flotow 1990, Meggio & Pinna 2003, Elia 2003"),
    ("priming_cantley", "kinase", f"SECONDARY. S/T Cantley log2 matrix: pS/pT is the top residue at some position and "
     f"log2 >= {CANTLEY_MIN_LOG2}", "kinase-library 1.8.0 ser_thr_all_log2_matrices.tsv"),
    ("pikk", "kinase", "Manning family PIKK (= ksp.eval 'pikk')", "DESIGN §9"),
    ("pikk_sq", "kinase", "ATM, ATR, PRKDC, SMG1 (S/T-Q PIKKs; MTOR and pseudokinase TRRAP excluded)", "DESIGN §9"),
    ("secretory", "kinase", "FAM20C, FAM20A, PKDCC present in index", "DESIGN §9 / KB-11"),
    ("fam20c", "kinase", "FAM20C", "DESIGN §9"),
    ("primed", "pair", "positive of a priming kinase whose priming-offset residue (GSK3 +4, CK1 -3, CK2 +1/+3; PLK1 "
     "none) is an S/T U_ms site; Cantley variant uses the matrix offsets", "DESIGN §9"),
    ("hub", "site", f"site with evaluable-positive labels to >= {HUB_MIN_GROUPS} non-paralogous kinases "
     "(Par-graph components)", "DESIGN §17 A26"),
    ("evidence", "pair", "evidence_native: in_vitro_only / in_vivo_only / both", "DESIGN §6.4"),
]


def definitions_table() -> pd.DataFrame:
    return pd.DataFrame(DEFINITIONS, columns=["stratum", "level", "rule", "source"])
