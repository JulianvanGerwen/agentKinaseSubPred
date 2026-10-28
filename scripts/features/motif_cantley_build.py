"""Build F1 Cantley motif feature arrays for all U_obs sites x all index kinases.

Outputs (cache/features/motif/), each float32 (n_kinases=503, n_sites=350,736), row/col order =
cache/index/{kinases,sites}.parquet, NaN where the kinase has no atlas matrix for the site's residue
type or the site residue is outside the kinase's residue class:
  cantley_{variant}_{metric}.npy  variant in {unprimed, primed_hard, primed_soft},
                                  metric in {score, pct_published, pct_universe}
  cantley_imputed_{variant}_{metric}.npy  same, only rows of kinases WITHOUT an own in-class matrix,
                                  filled from the nearest atlas kinase (kinase-domain identity).
  cantley_imputed_distance.npy    (n_kinases, 2) 1 - identity to the neighbour, columns [ST, Y].
  cantley_meta.json               description of every array.
Kinase-level table: results/features/motif/cantley_kinase_map.tsv.
Run: python -I scripts/features/motif_cantley_build.py
"""
import hashlib, json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.index import load_kinases, load_sites  # noqa: E402
from ksp.features import motif_cantley as mc  # noqa: E402
from ksp.uniprot import read_fasta  # noqa: E402
from ksp.store import save_matrix  # noqa: E402

OUT = ROOT / "cache/features/motif"
RES = ROOT / "results/features/motif"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)
VARIANTS = ("unprimed", "primed_hard", "primed_soft")
METRICS = ("score", "pct_published", "pct_universe")
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


kin = load_kinases()
sites = load_sites()
K, N = len(kin), len(sites)
kpos = {a: i for i, a in enumerate(kin.uniprot_acc)}
kclass = kin.residue_class.to_numpy()
log("index", K, N)

# ---- matrix -> UniProt mapping (kinome_information UNIPROT_ID), cross-checked with synonyms ----
info = mc.kinome_info()
syn = pd.read_csv(ROOT / "results/labels/kinase_synonyms.tsv", sep="\t", dtype=str)
syn = syn[syn.source.fillna("").str.contains("cantley_matrix_name")]
syn_map = dict(zip(syn.synonym_uc, syn.uniprot_acc))
disagree = [(r.MATRIX_NAME, r.UNIPROT_ID, syn_map.get(r.MATRIX_NAME.upper())) for r in info.itertuples()
            if r.MATRIX_NAME.upper() in syn_map and syn_map[r.MATRIX_NAME.upper()] != r.UNIPROT_ID]
assert not disagree, f"matrix-name mapping disagrees with kinase_synonyms.tsv: {disagree}"
assert info.UNIPROT_ID.isin(kpos).all(), "atlas kinase not in index"
assert info.set_index("MATRIX_NAME").loc["PDK1", "UNIPROT_ID"] == "O15530"  # PDK1 matrix = PDPK1
log("matrix mapping ok;", len(info), "matrices; synonyms cross-checked:", sum(n.upper() in syn_map for n in info.MATRIX_NAME))

# ---- priming weights ----
uni = pd.read_parquet(ROOT / "cache/sites/universe.parquet")
uni = uni[uni.in_U_obs]
w_soft_site = mc.ms_evidence_weight(uni)
w_hard_site = pd.Series(1.0, index=sites.site_id.to_numpy())          # every U_obs site
W = {"unprimed": None,
     "primed_hard": mc.flank_weights(sites, w_hard_site),
     "primed_soft": mc.flank_weights(sites, w_soft_site)}
log("flank weights: mean #primed flanks/site hard=%.3f soft=%.3f" % (W["primed_hard"].sum(1).mean(), W["primed_soft"].sum(1).mean()))

res = sites.res.to_numpy()
in_ms = sites.in_U_ms.to_numpy(bool)
arr = {}
for v in VARIANTS:
    for m in METRICS:
        arr[(v, m)] = np.full((K, N), np.nan, np.float32)

kmap_rows = []
own = {t: {} for t in mc.KIN_TYPES}   # acc -> atlas row index (in-class only)
for t in mc.KIN_TYPES:
    A = mc.load_atlas(t)
    ref = mc.load_reference(t)[A.names].to_numpy(float)
    cols = np.where(np.isin(res, mc.TYPE_RES[t]))[0]
    win = sites.window15.to_numpy()[cols]
    codes = mc.encode(win)
    centre = res[cols]
    rows = np.array([kpos[a] for a in A.uniprot])
    in_class = np.array([t in mc.CLASS_TYPES[kclass[r]] for r in rows])
    for i, (nm, a) in enumerate(zip(A.names, A.uniprot)):
        kmap_rows.append(dict(matrix_name=nm, kin_type=t, uniprot_acc=a, gene=kin.gene[kpos[a]],
                              residue_class=kclass[kpos[a]], used=bool(in_class[i]),
                              note="" if in_class[i] else "matrix type outside kinase residue class -> masked"))
        if in_class[i]:
            own[t][a] = i
    for v in VARIANTS:
        w = None if W[v] is None else W[v][cols]
        S = mc.score(A, codes, centre, phospho_w=w)                  # (n, K_t) float64, published rounding
        P_pub = mc.percentile(S, ref)
        P_uni = mc.percentile(S, S[in_ms[cols]], round_digits=None)
        for m, X in (("score", S), ("pct_published", P_pub), ("pct_universe", P_uni)):
            dst = arr[(v, m)]
            for i in np.where(in_class)[0]:
                dst[rows[i], cols] = X[:, i].astype(np.float32)
        log(t, v, "done", S.shape)
    del S, P_pub, P_uni

# ---- imputation for kinases lacking an own in-class matrix ----
fa = read_fasta()
seqs = dict(zip(fa[fa.isoform == ""].acc, fa[fa.isoform == ""].seq))
dom = mc.domain_sequences(kin, seqs)
imp = []
for t in mc.KIN_TYPES:
    need = [a for a, c in zip(kin.uniprot_acc, kclass) if t in mc.CLASS_TYPES[c] and a not in own[t]]
    nn = mc.nearest_atlas(need, list(own[t]), dom, n_jobs=12)
    nn["kin_type"] = t
    imp.append(nn)
    log("nearest atlas", t, len(need))
imp = pd.concat(imp, ignore_index=True)
imp["gene"] = imp.uniprot_acc.map(dict(zip(kin.uniprot_acc, kin.gene)))
imp["nn_gene"] = imp.nn_acc.map(dict(zip(kin.uniprot_acc, kin.gene)))
dist = np.full((K, 2), np.nan, np.float32)
type_cols = {t: np.where(np.isin(res, mc.TYPE_RES[t]))[0] for t in mc.KIN_TYPES}
for v in VARIANTS:
    for m in METRICS:
        src = arr[(v, m)]
        dst = np.full((K, N), np.nan, np.float32)
        for r in imp.itertuples():
            cols = type_cols[r.kin_type]
            dst[kpos[r.uniprot_acc], cols] = src[kpos[r.nn_acc], cols]
        save_matrix(f"motif/cantley_imputed_{v}_{m}", dst, "kinase", "site",
                    meta={"family": "F1", "variant": v, "metric": m, "imputed": True, "label_derived": False})
        del dst
for r in imp.itertuples():
    dist[kpos[r.uniprot_acc], 0 if r.kin_type == "ser_thr" else 1] = 1 - r.nn_identity
save_matrix("motif/cantley_imputed_distance", dist, "kinase", "none", col_ids=np.array(["ST", "Y"]),
            meta={"family": "F1", "what": "1 - kinase-domain identity to nearest atlas kinase", "label_derived": False})
log("imputed arrays written")

for (v, m), a in arr.items():
    save_matrix(f"motif/cantley_{v}_{m}", a, "kinase", "site",
                meta={"family": "F1", "variant": v, "metric": m, "imputed": False, "label_derived": False})
log("own-matrix arrays saved via ksp.store")
kmap = pd.DataFrame(kmap_rows)
imp_out = imp[["uniprot_acc", "gene", "kin_type", "nn_acc", "nn_gene", "nn_identity", "nn2_acc", "nn2_identity"]]
kmap.to_csv(RES / "cantley_kinase_map.tsv", sep="\t", index=False)
imp_out.to_csv(RES / "cantley_imputation_neighbours.tsv", sep="\t", index=False)

cov = {}
for v in VARIANTS[:1]:
    s = arr[(v, "score")]
    has = np.isfinite(s[:, :]).any(1)
    cov = {"kinases_with_own_matrix": int(has.sum()),
           "kinases_imputed_ST": int((imp.kin_type == "ser_thr").sum()),
           "kinases_imputed_Y": int((imp.kin_type == "tyrosine").sum())}
idx_manifest = json.loads((ROOT / "results/index_manifest.json").read_text())
meta = {
    "description": "F1 Cantley kinome-atlas motif features (kinase-library 1.8.0 matrices; reimplemented scoring). "
                   "Rows = cache/index/kinases.parquet order, columns = cache/index/sites.parquet order.",
    "shape": [K, N], "dtype": "float32",
    "index_manifest": idx_manifest,
    "source": "external/kinase_library/kinase_library-1.8.0 databases/matrices/*_all_log2_matrices.tsv, st_favorability.tsv, "
              "kinome_information.tsv (UNIPROT_ID), substrates/ochoa/scored_phosprots/*.parquet",
    "nan_rule": "NaN if kinase has no matrix for the site's residue type (S/T -> ser_thr, Y -> tyrosine) or the matrix type "
                "is outside the kinase's residue_class (index residue_class at build time; ST-class kinases with "
                "non-canonical *_TYR matrices are masked on Y sites).",
    "variants": {
        "unprimed": "Published definition: flanks upper-cased (kinase-library pp=False), st_fav=True, 3-dp rounding.",
        "primed_hard": "Flank S/T/Y that is itself a U_obs site scored with the phospho column (s/t/y), i.e. kinase-library pp=True "
                       "with those residues lower-cased.",
        "primed_soft": "Flank S/T/Y contributes (1-w)*log2M[upper] + w*log2M[lower], w = fraction of 4 MS evidence channels "
                       "(PSP MS_LIT+MS_CST>=1, compendium, CPTAC, Ochoa 2020) observing that flank site; 0 if not a U_obs site.",
    },
    "metrics": {
        "score": "log2 Cantley score.",
        "pct_published": "Percentile (0-100, 2 dp) vs the kinase's Ochoa-2020 reference distribution shipped with kinase-library "
                         "(82,735 S/T / 5,431 Y sites; unprimed scores). For primed variants this equals kinase-library pp=True behaviour.",
        "pct_universe": "Label-free universe-recalibrated percentile: 100*#{U_ms sites of same residue type (S/T or Y) with "
                        "score <= s}/n, computed within the same variant (primed scores vs primed U_ms distribution). Unrounded.",
    },
    "imputed": "cantley_imputed_* arrays hold ONLY rows of kinases lacking an own in-class matrix, copied from the nearest atlas "
               "kinase by kinase-domain identity (results/features/motif/cantley_imputation_neighbours.tsv); "
               "cantley_imputed_distance.npy (K x 2: ST, Y) = 1 - identity. Not part of the as-published arm.",
    "coverage": cov,
    "label_derived": False,
    "files": {},
}
meta["load"] = "ksp.store.load_matrix('motif/<file stem>') (verifies canonical index hash)"
for f in sorted(OUT.glob("cantley_*.npy")):
    meta["files"][f.name] = {"bytes": f.stat().st_size}
(OUT / "cantley_meta.json").write_text(json.dumps(meta, indent=1))
log("done", json.dumps(cov))
