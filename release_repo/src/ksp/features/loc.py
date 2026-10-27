"""F4 localisation: per-protein compartment vectors in a common coarse vocabulary + UniProt topology per site.

Coarse vocabulary (VOCAB, 15 terms; mapping tables below are the documentation of record, see docs/features_context.md):
  nucleoplasm, nucleolus, chromatin, nuclear_envelope, cytosol, plasma_membrane, cell_junction, mitochondrion, er,
  golgi, endolysosome_vesicle, peroxisome, centrosome_mt, actin_if, extracellular

Sources -> protein x VOCAB weight vectors (index = Swiss-Prot acc):
  hpa       HPA v25 subcellular IF (main location 1.0, additional 0.5) x per-location reliability
            (Enhanced 1.0, Supported 0.8, Approved 0.5, Uncertain 0.2); 'Predicted to be secreted' -> extracellular 0.5.
            Sperm-specific compartments (principal/mid/end/connecting piece, acrosome, ...) dropped.
  opencell  OpenCell live-cell annotations, grade 3/2/1 -> 1.0/0.6/0.3.
  hcm       Human Cell Map NMF rank profile (20 ranks, from Hein-repo copy of Go 2021 'BioID2021_NMF.xlsx', sheet
            'A) rank profiles'), projected to VOCAB through the rank -> MMF-localisation name mapping in preys-latest.txt.
  hein      Hein 2025 organelle IP-MS consensus cluster annotation (one-hot; 'unclassified' = missing).
  comp_exp  COMPARTMENTS experiments channel (confidence 0-5; NOTE: its only source is HPA, so it is largely
            redundant with 'hpa' and excluded from the combined vector).
  comp_know COMPARTMENTS knowledge channel (GO/UniProt curation) - ABLATION ONLY (imports curated annotations).
Native-space extras: hcm_nmf (20-dim NMF profile), hein_profile (z-scored enrichment over all organelle IPs).
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from .ppi_idmap import EXT, ROOT, MapLog, get_mapper

VOCAB = ["nucleoplasm", "nucleolus", "chromatin", "nuclear_envelope", "cytosol", "plasma_membrane", "cell_junction",
         "mitochondrion", "er", "golgi", "endolysosome_vesicle", "peroxisome", "centrosome_mt", "actin_if",
         "extracellular"]
VI = {v: i for i, v in enumerate(VOCAB)}
RAW = ROOT / "cache" / "features" / "context" / "raw"
HEIN_ZIP = EXT / "hein2025_zenodo_14000085" / "Organelle_IP_analyses_and_figures-2024-10-28.zip"
HEIN_P = "Organelle_IP_analyses_and_figures-main/"

HPA_MAP = {
    "Nucleoplasm": "nucleoplasm", "Nuclear bodies": "nucleoplasm", "Nuclear speckles": "nucleoplasm",
    "Nucleoli": "nucleolus", "Nucleoli fibrillar center": "nucleolus", "Nucleoli rim": "nucleolus",
    "Mitotic chromosome": "chromatin", "Kinetochore": "chromatin", "Nuclear membrane": "nuclear_envelope",
    "Cytosol": "cytosol", "Cytoplasmic bodies": "cytosol", "Aggresome": "cytosol", "Rods & Rings": "cytosol",
    "Plasma membrane": "plasma_membrane", "Cell Junctions": "cell_junction", "Focal adhesion sites": "cell_junction",
    "Mitochondria": "mitochondrion", "Endoplasmic reticulum": "er", "Golgi apparatus": "golgi",
    "Vesicles": "endolysosome_vesicle", "Endosomes": "endolysosome_vesicle", "Lysosomes": "endolysosome_vesicle",
    "Lipid droplets": "endolysosome_vesicle", "Peroxisomes": "peroxisome",
    "Centrosome": "centrosome_mt", "Centriolar satellite": "centrosome_mt", "Microtubules": "centrosome_mt",
    "Microtubule ends": "centrosome_mt", "Mitotic spindle": "centrosome_mt", "Cytokinetic bridge": "centrosome_mt",
    "Midbody": "centrosome_mt", "Midbody ring": "centrosome_mt", "Primary cilium": "centrosome_mt",
    "Primary cilium tip": "centrosome_mt", "Primary cilium transition zone": "centrosome_mt", "Basal body": "centrosome_mt",
    "Actin filaments": "actin_if", "Intermediate filaments": "actin_if", "Cleavage furrow": "actin_if",
}
HPA_DROP = {"Principal piece", "Mid piece", "End piece", "Connecting piece", "Acrosome", "Equatorial segment",
            "Perinuclear theca", "Calyx", "Flagellar centriole", "Annulus"}
HPA_REL = {"Enhanced": 1.0, "Supported": 0.8, "Approved": 0.5, "Uncertain": 0.2}

OPENCELL_MAP = {"cytoplasmic": {"cytosol": 1}, "nucleoplasm": {"nucleoplasm": 1}, "vesicles": {"endolysosome_vesicle": 1},
                "membrane": {"plasma_membrane": 1}, "nuclear_punctae": {"nucleoplasm": 1}, "er": {"er": 1},
                "chromatin": {"chromatin": 1}, "golgi": {"golgi": 1}, "nucleolus_gc": {"nucleolus": 1},
                "nucleolus_fc_dfc": {"nucleolus": 1}, "centrosome": {"centrosome_mt": 1}, "big_aggregates": {"cytosol": 1},
                "cell_contact": {"cell_junction": 1}, "cytoskeleton": {"actin_if": 0.5, "centrosome_mt": 0.5},
                "nuclear_membrane": {"nuclear_envelope": 1}, "focal_adhesions": {"cell_junction": 1},
                "mitochondria": {"mitochondrion": 1}}
OPENCELL_GRADE = {"annotations_grade_3": 1.0, "annotations_grade_2": 0.6, "annotations_grade_1": 0.3}

HCM_MAP = {"chromatin": {"chromatin": 1}, "plasma membrane": {"plasma_membrane": 1},
           "actin cytoskeleton, cytosol": {"actin_if": 0.5, "cytosol": 0.5}, "mitochondrial matrix": {"mitochondrion": 1},
           "nuclear body": {"nucleoplasm": 1}, "miscellaneous": {},
           "nuclear outer membrane-ER membrane network": {"er": 0.5, "nuclear_envelope": 0.5},
           "centrosome": {"centrosome_mt": 1}, "nucleoplasm": {"nucleoplasm": 1}, "cell junction": {"cell_junction": 1},
           "mitochondrial outer membrane, peroxisome": {"mitochondrion": 0.5, "peroxisome": 0.5},
           "endosome, lysosome": {"endolysosome_vesicle": 1}, "ER lumen": {"er": 1}, "nucleolus": {"nucleolus": 1},
           "Golgi apparatus": {"golgi": 1}, "early endosome, recycling endosome": {"endolysosome_vesicle": 1},
           "cytoplasmic ribonucleoprotein granule": {"cytosol": 1}, "ER membrane": {"er": 1},
           "mitochondrial inner membrane, mitochondrial intermembrane space": {"mitochondrion": 1},
           "microtubule cytoskeleton": {"centrosome_mt": 1}}

HEIN_MAP = {"cytosol": {"cytosol": 1}, "nucleus": {"nucleoplasm": 1}, "mitochondrion": {"mitochondrion": 1},
            "ER": {"er": 1}, "plasma_membrane": {"plasma_membrane": 1}, "actin_cytoskeleton": {"actin_if": 1},
            "nucleolus": {"nucleolus": 1}, "early_endosome": {"endolysosome_vesicle": 1}, "stress_granule": {"cytosol": 1},
            "Golgi": {"golgi": 1}, "lysosome": {"endolysosome_vesicle": 1}, "translation": {"cytosol": 1},
            "14-3-3_scaffold": {"cytosol": 1}, "trans-Golgi": {"golgi": 1}, "recycling_endosome": {"endolysosome_vesicle": 1},
            "ERGIC": {"er": 0.5, "golgi": 0.5}, "proteasome": {"cytosol": 1}, "peroxisome": {"peroxisome": 1},
            "centrosome": {"centrosome_mt": 1}, "p-body": {"cytosol": 1}}

# GO CC anchor terms for COMPARTMENTS (files are ancestor-propagated, so an anchor's score covers its descendants)
GO_ANCHORS = {
    "nucleoplasm": ["GO:0005654", "GO:0016604", "GO:0016607"], "nucleolus": ["GO:0005730"],
    "chromatin": ["GO:0000785", "GO:0005694", "GO:0000776"], "nuclear_envelope": ["GO:0005635"],
    "cytosol": ["GO:0005829"], "plasma_membrane": ["GO:0005886"], "cell_junction": ["GO:0030054"],
    "mitochondrion": ["GO:0005739"], "er": ["GO:0005783"], "golgi": ["GO:0005794"],
    "endolysosome_vesicle": ["GO:0031410", "GO:0005768", "GO:0005764", "GO:0005811"], "peroxisome": ["GO:0005777"],
    "centrosome_mt": ["GO:0005815", "GO:0015630", "GO:0005929", "GO:0030496"],
    "actin_if": ["GO:0015629", "GO:0045111", "GO:0005938"], "extracellular": ["GO:0005615"],
}
GO_NUCLEUS = "GO:0005634"  # used for 'nucleoplasm' only when no nuclear sub-compartment has >= that score


def _acc_vectors(rows: dict) -> pd.DataFrame:
    """rows: acc -> np.array(len(VOCAB)); returns DataFrame (acc x VOCAB), all-zero rows dropped."""
    if not rows:
        return pd.DataFrame(columns=VOCAB)
    df = pd.DataFrame.from_dict(rows, orient="index", columns=VOCAB)
    return df[df.sum(1) > 0]


def _add(rows, accs, vec):
    for a in accs:
        rows[a] = np.maximum(rows.get(a, np.zeros(len(VOCAB))), vec)


def load_hpa(log: MapLog) -> tuple[pd.DataFrame, pd.Series]:
    m = get_mapper()
    f = RAW / "hpa" / "subcellular_location.tsv"
    if not f.exists():
        f.parent.mkdir(parents=True, exist_ok=True)
        zipfile.ZipFile(EXT / "hpa" / "subcellular_location.tsv.zip").extractall(f.parent)
    h = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    rows, rel, n_ok, unknown = {}, {}, 0, set()
    for r in h.itertuples(index=False):
        accs = m.first((m.from_ensg, r.Gene), (m.from_symbol, r._1))
        if not accs:
            continue
        n_ok += 1
        loc_rel = {}
        for cat in HPA_REL:
            for loc in filter(None, getattr(r, cat).split(";")):
                loc_rel[loc] = max(loc_rel.get(loc, 0), HPA_REL[cat])
        v = np.zeros(len(VOCAB))
        for col, w in (("_3", 1.0), ("_4", 0.5)):
            for loc in filter(None, getattr(r, col).split(";")):
                if loc in HPA_DROP:
                    continue
                if loc not in HPA_MAP:
                    unknown.add(loc)
                    continue
                i = VI[HPA_MAP[loc]]
                v[i] = max(v[i], w * loc_rel.get(loc, HPA_REL.get(r.Reliability, 0.5)))
        if r._5 == "Predicted to be secreted":
            v[VI["extracellular"]] = max(v[VI["extracellular"]], 0.5)
        _add(rows, accs, v)
        for a in accs:
            rel[a] = max(rel.get(a, 0), HPA_REL.get(r.Reliability, 0))
    if unknown:
        raise ValueError(f"HPA locations without mapping: {sorted(unknown)}")
    log.add("hpa", "genes", len(h), n_ok)
    return _acc_vectors(rows), pd.Series(rel, name="hpa_reliability")


def load_opencell(log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    o = pd.read_csv(EXT / "opencell" / "opencell-localization-annotations.csv", dtype=str, keep_default_na=False)
    meta = pd.read_csv(EXT / "opencell" / "opencell-library-metadata.csv", dtype=str)
    uni = dict(zip(meta.ensg_id, meta.uniprot_id))
    rows, n_ok = {}, 0
    for r in o.itertuples(index=False):
        accs = m.first((m.from_uniprot, uni.get(r.ensg_id, "")), (m.from_ensg, r.ensg_id), (m.from_symbol, r.target_name))
        if not accs:
            continue
        n_ok += 1
        v = np.zeros(len(VOCAB))
        for col, g in OPENCELL_GRADE.items():
            for lab in filter(None, getattr(r, col).split(";")):
                for t, w in OPENCELL_MAP[lab].items():  # KeyError on unknown labels (never silent)
                    v[VI[t]] = max(v[VI[t]], g * w)
        _add(rows, accs, v)
    log.add("opencell_loc", "targets", len(o), n_ok)
    return _acc_vectors(rows)


def _hein_extract(name):
    out = RAW / "hein" / HEIN_P / name
    if not out.exists():
        zipfile.ZipFile(HEIN_ZIP).extract(HEIN_P + name, RAW / "hein")
    return out


def load_hcm(log: MapLog) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = get_mapper()
    nmf = pd.read_excel(_hein_extract("data/external/BioID2021_NMF.xlsx"), sheet_name="A) rank profiles")
    preys = pd.read_csv(EXT / "humancellmap" / "preys-latest.txt", sep="\t", dtype=str)
    rank_name = preys.groupby("NMF rank")["MMF localization"].agg(lambda s: s.value_counts().index[0])
    proj = np.zeros((20, len(VOCAB)))
    for rk, nm in rank_name.items():
        for t, w in HCM_MAP[nm].items():
            proj[int(rk) - 1, VI[t]] = w
    uni = dict(zip(preys.symbol, preys.UniProt))
    ent = dict(zip(preys.symbol, preys.Entrez))
    rows, nat, n_ok = {}, {}, 0
    X = nmf[[f"rank {i}" for i in range(1, 21)]].to_numpy(float)
    for g, x in zip(nmf.gene, X):
        accs = m.first((m.from_uniprot, uni.get(g, "")), (m.from_entrez, ent.get(g, "")), (m.from_symbol, g))
        if not accs or x.sum() <= 0:
            continue
        n_ok += 1
        _add(rows, accs, x @ proj)
        for a in accs:
            nat[a] = x
    log.add("hcm_loc", "preys", len(nmf), n_ok)
    return _acc_vectors(rows), pd.DataFrame.from_dict(nat, orient="index", columns=[f"rank{i}" for i in range(1, 21)])


def load_hein(log: MapLog) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = get_mapper()
    c = pd.read_csv(_hein_extract("data/labels/cluster_annotation_Dec6.csv"))
    rows, n_ok, n_cls = {}, 0, 0
    for ids, gname, ann in zip(c["Majority protein IDs"], c.Gene_name_canonical, c.cluster_annotation):
        if ann == "unclassified":
            continue
        n_cls += 1
        accs = set()
        for x in str(ids).split(";"):
            accs = m.from_uniprot(x)
            if accs:
                break
        accs = accs or m.from_symbol(gname)
        if not accs:
            continue
        n_ok += 1
        v = np.zeros(len(VOCAB))
        for t, w in HEIN_MAP[ann].items():
            v[VI[t]] = w
        _add(rows, accs, v)
    log.add("hein_loc", "classified_protein_groups", n_cls, n_ok)
    # native enrichment profile
    e = pd.read_csv(_hein_extract("notebooks/enrichment/output/enrichment_and_volcano_tables/"
                                  "2023-10-21-imp5-for-figures_enrichment_table_NOC_prop.csv"), header=[0, 1], index_col=0)
    samp = [col for col in e.columns if col[0] == "sample" and not col[1].startswith("NOC_")]
    Xs = e[samp].to_numpy(float)
    Xs = (Xs - np.nanmean(Xs, 0)) / np.nanstd(Xs, 0)
    nat, n2 = {}, 0
    for ids, gname, x in zip(e[("metadata", "Majority protein IDs")], e[("metadata", "Gene_name_canonical")], Xs):
        accs = set()
        for xx in str(ids).split(";"):
            accs = m.from_uniprot(xx)
            if accs:
                break
        accs = accs or m.from_symbol(str(gname))
        if not accs or not np.isfinite(x).all():
            continue
        n2 += 1
        for a in accs:
            nat.setdefault(a, x)  # first (majority) group wins
    log.add("hein_profile", "protein_groups", len(e), n2)
    return _acc_vectors(rows), pd.DataFrame.from_dict(nat, orient="index", columns=[c[1] for c in samp])


def load_compartments(channel: str, log: MapLog) -> pd.DataFrame:
    m = get_mapper()
    f = EXT / "compartments" / f"human_compartment_{channel}_full.tsv"
    d = pd.read_csv(f, sep="\t", header=None, dtype=str, usecols=[0, 1, 2, 6], names=["ensp", "sym", "go", "score"])
    anchors = {g: t for t, gs in GO_ANCHORS.items() for g in gs}
    d = d[d.go.isin(set(anchors) | {GO_NUCLEUS})]
    d["score"] = d.score.astype(float)
    acc_of, n_ok = {}, 0
    ids = d[["ensp", "sym"]].drop_duplicates()
    for p, s in zip(ids.ensp, ids.sym):
        acc_of[(p, s)] = m.first((m.from_ensp, p), (m.from_symbol, s))
    log.add(f"compartments_{channel}", "proteins", len(ids), sum(bool(v) for v in acc_of.values()))
    piv = d.pivot_table(index=["ensp", "sym"], columns="go", values="score", aggfunc="max").fillna(0)
    rows = {}
    for (p, s), r in piv.iterrows():
        accs = acc_of[(p, s)]
        if not accs:
            continue
        v = np.zeros(len(VOCAB))
        for g, t in anchors.items():
            if g in r.index:
                v[VI[t]] = max(v[VI[t]], r[g])
        nuc = r.get(GO_NUCLEUS, 0)
        if nuc > max(v[VI["nucleoplasm"]], v[VI["nucleolus"]], v[VI["chromatin"]], v[VI["nuclear_envelope"]]):
            v[VI["nucleoplasm"]] = nuc
        _add(rows, accs, v / 5.0)
    return _acc_vectors(rows)


# ---------------------------------------------------------------- similarity
def unit_rows(X: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(X, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return X / n


def cosine_kxp(V: pd.DataFrame, kin, prot) -> np.ndarray:
    """cosine similarity kinase x protein from per-acc vectors V; NaN if either is missing."""
    U = unit_rows(V.to_numpy(float))
    pos = {a: i for i, a in enumerate(V.index)}
    D = U.shape[1]
    Uk = np.full((len(kin), D), np.nan)
    Up = np.full((len(prot), D), np.nan)
    for i, a in enumerate(kin):
        if a in pos:
            Uk[i] = U[pos[a]]
    for j, a in enumerate(prot):
        if a in pos:
            Up[j] = U[pos[a]]
    return (Uk @ Up.T).astype(np.float32)


def corr_kxp(V: pd.DataFrame, kin, prot) -> np.ndarray:
    """Pearson correlation (row-centred cosine) kinase x protein."""
    X = V.to_numpy(float)
    X = X - X.mean(1, keepdims=True)
    return cosine_kxp(pd.DataFrame(X, index=V.index), kin, prot)


# ---------------------------------------------------------------- UniProt topology (site level)
TOPO_FILE = EXT / "uniprot_topology" / "UP000005640_reviewed_topology.tsv.gz"
LUMENAL_NOTES = {"Extracellular", "Lumenal", "Vesicular", "Lumenal, vesicle", "Perinuclear space", "Vacuolar",
                 "Lumenal, melanosome", "Intragranular", "Exoplasmic loop"}
CYTO_NOTES = {"Cytoplasmic", "Nuclear", "Mother cell cytoplasmic"}
MITO_NOTES = {"Mitochondrial intermembrane", "Mitochondrial matrix"}
PEROX_NOTES = {"Peroxisomal matrix", "Peroxisomal"}
_FT = re.compile(r"(TOPO_DOM|TRANSMEM|INTRAMEM|SIGNAL|TRANSIT) <?(\d+|\?)\.\.>?(\d+|\?);(?: /note=\"([^\"]*)\")?")

# site_topo codes
TOPO_CODES = {0: "no topology annotation, no signal/TM (default cytoplasmic/nuclear)", 1: "cytoplasmic/nuclear TOPO_DOM",
              2: "extracellular/lumenal TOPO_DOM", 3: "transmembrane/intramembrane segment",
              4: "signal or transit peptide", 5: "secreted/lumenal (signal peptide, no TM; after signal)",
              6: "mitochondrial (intermembrane/matrix TOPO_DOM, or after a mitochondrial transit peptide, no TM)",
              7: "peroxisomal TOPO_DOM", 8: "unknown (TM protein, position not covered by TOPO_DOM)"}


def protein_topology():
    """acc -> list of (kind, start, end, note)."""
    t = pd.read_csv(TOPO_FILE, sep="\t", dtype=str, keep_default_na=False)
    out = {}
    cols = ["Topological domain", "Transmembrane", "Intramembrane", "Signal peptide", "Transit peptide"]
    for entry, *vals in zip(t.Entry, *[t[c] for c in cols]):
        feats = []
        for col in vals:
            for mm in _FT.finditer(col):
                kind, s, e, note = mm.groups()
                if s == "?" or e == "?":
                    continue
                feats.append((kind, int(s), int(e), note or ""))
        out[entry] = feats
    return out


def site_topology(sites: pd.DataFrame) -> np.ndarray:
    """sites: DataFrame[acc, pos] -> int8 code per site (TOPO_CODES)."""
    topo = protein_topology()
    code = np.zeros(len(sites), np.int8)
    for n, (acc, pos) in enumerate(zip(sites.acc.to_numpy(), sites.pos.to_numpy())):
        feats = topo.get(acc, [])
        if not feats:
            continue
        has_tm = any(k in ("TRANSMEM",) for k, *_ in feats)
        c = None
        for k, s, e, note in feats:
            if s <= pos <= e:
                if k in ("TRANSMEM", "INTRAMEM"):
                    c = 3
                elif k in ("SIGNAL", "TRANSIT"):
                    c = 4
                elif k == "TOPO_DOM":
                    c = 2 if note in LUMENAL_NOTES else 1 if note in CYTO_NOTES else 6 if note in MITO_NOTES \
                        else 7 if note in PEROX_NOTES else 8
                if c in (3, 4):
                    break
        if c is None:
            sig = [e for k, s, e, nt in feats if k == "SIGNAL"]
            mito = [e for k, s, e, nt in feats if k == "TRANSIT" and nt.startswith("Mitochondri")]
            if has_tm:
                c = 8
            elif sig and pos > max(sig):
                c = 5
            elif mito and pos > max(mito):
                c = 6
            else:
                c = 0
        code[n] = c
    return code
