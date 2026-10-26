"""SELPHI2.0 feature table, regenerated on our IDs (DESIGN §8.1 B4, PRIOR-4).

The released SELPHI2.0 feature table (external/selphi2/selphi_2.0/data/psp/pos_set_features_20240910.tsv) has
44 model columns: `residue` + 43 features. Their definitions were traced in the SELPHI code (file names relative
to external/selphi2/selphi_2.0/src/feature_table_generation/):

  site-level, Ochoa 2020 MOESM4/5 (sign-test-set.r, write-large-pos-set.r)          -> SITE_COLS (31 incl. residue)
  kinase x protein expression (coexpr.r, write-large-pos-set.r)                     -> gtex, rna.tiss, rna.cell
  kinase / substrate tissue skewness (coexpr.r: skewness of HPA tissue profile)     -> sel.1 (kinase), sel.2 (substrate)
  kinase x site Cantley PSSM, min-max per kinase (norm-pssm.r)                      -> pssm_ser_thr, pssm_tyr
  kinase x site PSP-derived PSSM with family fallback (pssm.py)        LABEL-DERIVED -> score
  KSEA kinase activity (PSP regulon) x site Spearman, -log10 p
     (calc-cuttilas-act.r, cutillas_act-cor.r; calc-cptac-act.r)        LABEL-DERIVED -> coreg293, coreg.ntera/hl/mcf7

Imputation follows their code: per-column replace_na defaults (sign-test-set.r), then every remaining NA -> 0
(merge_kin_sites.r). Non-numeric columns are integer-coded like R `as.integer(as.factor(x))` (levels sorted,
1-based; logicals FALSE=1, TRUE=2). `residue` uses learn.py's coding {T: 0, Y: 1, S: 2}.

Label-derived columns are ksp.crossfit families (`SelphiPSSM`, `SelphiKSEACoreg`); they never see labels other than
the `train_pairs` handed to `fit`.
"""
from __future__ import annotations

import json
import zipfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ksp.crossfit import LabelDerivedFeature, register
from ksp.index import ROOT, load_kinases, load_proteins, load_sites

CACHE = ROOT / "cache" / "selphi2"
EXT = ROOT / "external"
OCHOA4 = EXT / "ochoa2020" / "41587_2019_344_MOESM4_ESM.xlsx"
OCHOA5 = EXT / "ochoa2020" / "41587_2019_344_MOESM5_ESM.xlsx"
GTEX = EXT / "gtex_v8" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
HPA_CELL = EXT / "hpa" / "rna_celline.tsv.zip"
HPA_TISS = EXT / "hpa" / "rna_tissue_consensus.tsv.zip"

# --- the 44 model columns, in the order of SELPHI's released table -------------------------------------------------
SELPHI_COLS = ["residue", "score", "coreg293", "functional_score", "gtex", "rna.tiss", "rna.cell", "sel.1", "sel.2",
               "coreg.ntera", "coreg.hl", "coreg.mcf7", "is_disopred", "disopred_score", "exp3d_ala_ddG_effect",
               "exp3d_acid_ddG_effect", "log10_hotspot_pval_min", "isHotspot", "isInterface", "adj_ptms_w21",
               "netpho_max_all", "netpho_max_KIN", "paxdb_abundance_log10", "w0_mya", "w3_mya", "quant_top1",
               "quant_top5", "PWM_max_mss", "ACCpro", "SSpro", "SSpro8", "sift_min_score", "sift_mean_score",
               "sift_ala_score", "sift_acid_score", "isProteinDomain", "isProteinKinaseDomain", "isUniprotRegion",
               "isCytoplasmic", "isELMLinearMotif", "isEV_ala_prediction_epistatic5", "isKinaseCoreg",
               "pssm_ser_thr", "pssm_tyr"]
SITE_COLS = ["residue", "functional_score", "is_disopred", "disopred_score", "exp3d_ala_ddG_effect",
             "exp3d_acid_ddG_effect", "log10_hotspot_pval_min", "isHotspot", "isInterface", "adj_ptms_w21",
             "netpho_max_all", "netpho_max_KIN", "paxdb_abundance_log10", "w0_mya", "w3_mya", "quant_top1",
             "quant_top5", "PWM_max_mss", "ACCpro", "SSpro", "SSpro8", "sift_min_score", "sift_mean_score",
             "sift_ala_score", "sift_acid_score", "isProteinDomain", "isProteinKinaseDomain", "isUniprotRegion",
             "isCytoplasmic", "isELMLinearMotif", "isEV_ala_prediction_epistatic5", "isKinaseCoreg"]
EXPR_COLS = ["gtex", "rna.tiss", "rna.cell", "sel.1", "sel.2"]
PSSM_COLS = ["pssm_ser_thr", "pssm_tyr"]
LABEL_COLS = ["score", "coreg293", "coreg.ntera", "coreg.hl", "coreg.mcf7"]
# features derived from known KS annotations / another predictor (PRIOR-4 leakage ablation arm)
LEAKY_COLS = ["functional_score", "PWM_max_mss", "netpho_max_KIN", "isKinaseCoreg"]
assert set(SITE_COLS) | set(EXPR_COLS) | set(PSSM_COLS) | set(LABEL_COLS) == set(SELPHI_COLS)
assert len(SELPHI_COLS) == 44 and len(set(SELPHI_COLS)) == 44

RES_CODE = {"T": 0, "Y": 1, "S": 2}          # learn.py `aad`


def _rfactor(s: pd.Series) -> np.ndarray:
    """R as.integer(as.factor(x)): sorted levels, 1-based codes."""
    s = s.astype(str)
    lev = sorted(s.unique())
    return s.map({v: i + 1 for i, v in enumerate(lev)}).to_numpy(np.float32)


def _rbool(b: pd.Series) -> np.ndarray:
    return np.where(b.astype(bool).to_numpy(), 2.0, 1.0).astype(np.float32)


def _is_true(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "1.0"])


def _nonempty(s: pd.Series) -> pd.Series:
    return s.notna() & s.astype(str).ne("") & s.astype(str).ne("NA") & s.astype(str).ne("nan")


def read_ochoa() -> pd.DataFrame:
    """Ochoa 2020 annotated phosphoproteome (MOESM4) + functional score (MOESM5), keyed by our canonical site_id via
    the 15-mer mapping in cache/sites/mapping_log.parquet (source 'ochoa')."""
    p = CACHE / "ochoa_moesm4_5.parquet"
    if p.exists():
        return pd.read_parquet(p)
    a = pd.read_excel(OCHOA4, sheet_name="annotated_phosphoproteome", na_values=["NA"])
    f = pd.read_excel(OCHOA5, sheet_name="functional_score", na_values=["NA"])
    a = a.merge(f, on=["uniprot", "position"], how="left")
    a["src_id"] = a.uniprot.astype(str) + "_" + a.position.astype(int).astype(str)
    m = pd.read_parquet(ROOT / "cache" / "sites" / "mapping_log.parquet")
    m = m[m.source.eq("ochoa")][["src_id", "site_id"]]
    a = a.merge(m, on="src_id", how="left")
    for c in a.columns:
        if a[c].dtype == object:
            a[c] = a[c].astype(str).where(a[c].notna(), None)
    a.to_parquet(p)
    return a


def build_site_table(all_sites: bool = False) -> pd.DataFrame:
    """SITE_COLS per canonical site (index = site_id). Default: Ochoa-covered sites only (SELPHI's universe). With
    all_sites=True every index site gets a row; sites absent from Ochoa receive SELPHI's own NA defaults."""
    a = read_ochoa()
    n_raw = len(a)
    a = a[a.site_id.notna()].drop_duplicates("site_id")
    S = load_sites()[["site_id", "res"]]
    a = S.merge(a, on="site_id", how="left" if all_sites else "inner")
    d = pd.DataFrame(index=a.site_id.values)
    d["residue"] = a.res.map(RES_CODE).to_numpy(np.float32)
    med = lambda c: a[c].fillna(a[c].median())  # noqa: E731
    d["functional_score"] = med("functional_score").to_numpy(np.float32)
    d["is_disopred"] = _rbool(_is_true(a.is_disopred) | a.is_disopred.isna())          # replace_na(TRUE)
    d["disopred_score"] = med("disopred_score").to_numpy(np.float32)
    for c in ("exp3d_ala_ddG_effect", "exp3d_acid_ddG_effect"):
        d[c] = _rfactor(a[c].fillna("unknown"))
    hp = a.log10_hotspot_pval_min.replace([np.inf, -np.inf], np.nan)
    d["log10_hotspot_pval_min"] = a.log10_hotspot_pval_min.replace(np.inf, hp.max()).fillna(0).to_numpy(np.float32)
    d["isHotspot"] = _rbool(_is_true(a.isHotspot))
    d["isInterface"] = _rbool(_is_true(a.isInterface))
    d["adj_ptms_w21"] = a.adj_ptms_w21.fillna(0).to_numpy(np.float32)
    for c in ("netpho_max_all", "netpho_max_KIN", "paxdb_abundance_log10", "PWM_max_mss"):
        d[c] = med(c).to_numpy(np.float32)
    for c in ("w0_mya", "w3_mya"):
        d[c] = a[c].fillna(0).to_numpy(np.float32)
    for c in ("quant_top1", "quant_top5"):
        d[c] = a[c].fillna(0).to_numpy(np.float32)                                         # counts
    for c in ("ACCpro", "SSpro", "SSpro8"):
        d[c] = _rfactor(a[c].fillna("-"))
    for c in ("sift_min_score", "sift_mean_score", "sift_ala_score", "sift_acid_score"):
        d[c] = a[c].fillna(0).to_numpy(np.float32)                                        # NA kept, later NA -> 0
    d["isProteinDomain"] = _rbool(_nonempty(a.domain_uniprot))
    d["isProteinKinaseDomain"] = _rbool(a.domain_uniprot.astype(str).str.contains("inase", na=False))
    d["isUniprotRegion"] = _rbool(_nonempty(a.region_uniprot))
    d["isCytoplasmic"] = _rbool(a.topology_uniprot.astype(str).eq("Cytoplasmic"))
    elm = np.zeros(len(a), bool)
    for c in ("elm_CLV", "elm_DEG", "elm_DOC", "elm_LIG", "elm_MOD", "elm_TRG", "isELMkinaseMotif"):
        elm |= _is_true(a[c]).to_numpy()
    d["isELMLinearMotif"] = _rbool(pd.Series(elm))
    d["isEV_ala_prediction_epistatic5"] = _rbool(a.EV_ala_prediction_epistatic.lt(-5).fillna(False))
    d["isKinaseCoreg"] = _rbool(_nonempty(a.ptmdb_maxcoreg_kinase))
    d = d[SITE_COLS].astype(np.float32).replace([np.inf, -np.inf], np.nan).fillna(0)
    d.attrs["n_ochoa_rows"] = n_raw
    d.attrs["n_mapped_unique_sites"] = int(a.site_id.nunique())
    return d


# ----------------------------------------------------------------------------------------------- expression features
def _map_gene_profiles(genes: pd.DataFrame, X: np.ndarray, ens_col: str, sym_col: str) -> pd.DataFrame:
    """Gene rows -> Swiss-Prot accessions (ENSG first, then symbol); several genes per acc -> median (as SELPHI's
    aggregate(..., FUN=median) by gene name)."""
    from ksp.features.ppi_idmap import get_mapper
    m = get_mapper()
    acc_rows, idx = [], []
    for i, (g, s) in enumerate(zip(genes[ens_col], genes[sym_col])):
        for a in m.first((m.from_ensg, g), (m.from_symbol, s)):
            acc_rows.append(a)
            idx.append(i)
    P = pd.DataFrame(X[np.array(idx)], index=acc_rows)
    out = P.groupby(level=0).median()
    out.attrs["n_genes"] = len(genes)
    out.attrs["n_genes_mapped"] = len(set(idx))
    return out


def gtex_profiles() -> pd.DataFrame:
    g = pd.read_csv(GTEX, sep="\t", skiprows=2)
    return _map_gene_profiles(g[["Name", "Description"]], g.iloc[:, 2:].to_numpy(float), "Name", "Description")


def _hpa_long(path: Path, col: str, val: str) -> tuple[pd.DataFrame, np.ndarray, list]:
    with zipfile.ZipFile(path) as z:
        d = pd.read_csv(z.open(z.namelist()[0]), sep="\t", usecols=["Gene", "Gene name", col, val])
    X = d.pivot_table(index=["Gene", "Gene name"], columns=col, values=val, aggfunc="mean").fillna(0)
    genes = X.index.to_frame(index=False)
    return genes, X.to_numpy(float), list(X.columns)


def hpa_cell_profiles() -> pd.DataFrame:
    genes, X, _ = _hpa_long(HPA_CELL, "Cell line", "nTPM")
    return _map_gene_profiles(genes, X, "Gene", "Gene name")


def hpa_tissue_profiles_raw() -> pd.DataFrame:
    genes, X, _ = _hpa_long(HPA_TISS, "Tissue", "nTPM")
    return _map_gene_profiles(genes, X, "Gene", "Gene name")


def spearman_kxp(P: pd.DataFrame, kin_acc, prot_acc) -> np.ndarray:
    """Spearman (rank, then Pearson) kinase x protein; NaN where either is absent or has zero variance."""
    R = P.rank(axis=1).to_numpy(float)
    R = R - R.mean(1, keepdims=True)
    nrm = np.sqrt((R ** 2).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        R = R / nrm[:, None]
    pos = pd.Series(np.arange(len(P)), index=P.index)
    ki = pos.reindex(kin_acc).to_numpy()
    pi = pos.reindex(prot_acc).to_numpy()
    out = np.full((len(kin_acc), len(prot_acc)), np.nan, np.float32)
    km, pm = ~np.isnan(ki), ~np.isnan(pi)
    out[np.ix_(km, pm)] = (R[ki[km].astype(int)] @ R[pi[pm].astype(int)].T).astype(np.float32)
    return out


def build_expression(verbose: bool = True) -> dict:
    """gtex / rna.cell (kinase x protein Spearman), sel (protein tissue skewness) -> cache/selphi2/expr.npz.
    rna.tiss is read from our feature store (context/coexpr/coexpr_hpa_tissue_spearman; same definition)."""
    K = load_kinases().uniprot_acc.to_numpy()
    Pr = load_proteins().acc.to_numpy()
    info = {}
    G = gtex_profiles()
    info["gtex"] = dict(n_genes=G.attrs["n_genes"], n_genes_mapped=G.attrs["n_genes_mapped"], n_tissues=G.shape[1])
    gt = spearman_kxp(G, K, Pr)
    C = hpa_cell_profiles()
    info["rna_cell"] = dict(n_genes=C.attrs["n_genes"], n_genes_mapped=C.attrs["n_genes_mapped"], n_lines=C.shape[1])
    rc = spearman_kxp(C, K, Pr)
    T = hpa_tissue_profiles_raw()
    sk = pd.Series(stats.skew(T.to_numpy(float), axis=1, bias=True), index=T.index)
    info["sel"] = dict(n_tissues=T.shape[1], n_proteins=len(T))
    np.savez(CACHE / "expr.npz", kinase_acc=K.astype(str), protein_acc=Pr.astype(str), gtex=gt, rna_cell=rc,
             sel_protein=sk.reindex(Pr).to_numpy(np.float32), sel_kinase=sk.reindex(K).to_numpy(np.float32))
    (CACHE / "expr_info.json").write_text(json.dumps(info, indent=1))
    if verbose:
        print("expression:", info)
    return info


# ----------------------------------------------------------------------------------------------- loaders (by ID)
@lru_cache(maxsize=1)
def load_site_table() -> pd.DataFrame:
    return pd.read_parquet(CACHE / "site_features.parquet")


@lru_cache(maxsize=1)
def load_expression() -> dict:
    from ksp.store import load_matrix
    z = np.load(CACHE / "expr.npz")
    K = load_kinases().uniprot_acc.to_numpy().astype(str)
    Pr = load_proteins().acc.to_numpy().astype(str)
    if not (np.array_equal(z["kinase_acc"], K) and np.array_equal(z["protein_acc"], Pr)):
        raise RuntimeError("cache/selphi2/expr.npz is not aligned to the current kinase/protein index; rebuild it")
    rt, _ = load_matrix("context/coexpr/coexpr_hpa_tissue_spearman")
    return {"gtex": z["gtex"], "rna.cell": z["rna_cell"], "rna.tiss": np.asarray(rt, np.float32),
            "sel.1": z["sel_kinase"], "sel.2": z["sel_protein"]}


@lru_cache(maxsize=1)
def load_cantley_norm() -> tuple[np.ndarray, np.ndarray]:
    """Cantley (kinase-library 1.8.0) unprimed log2 score from our store, min-max normalised per kinase and residue
    type over SELPHI's universe (U_ms sites with Ochoa features), as norm-pssm.r does over their scored sites.
    Returns (st, y) arrays K x |site index| float32 with NaN where no matrix."""
    from ksp.store import load_matrix
    p = CACHE / "cantley_norm.npz"
    S = load_sites()
    if p.exists():
        z = np.load(p)
        if np.array_equal(z["site_id"], S.site_id.to_numpy().astype(str)):
            return z["st"], z["y"]
    M, _ = load_matrix("motif/cantley_unprimed_score")
    uni = set(load_site_table().index)
    in_u = S.site_id.isin(uni).to_numpy() & S.in_U_ms.to_numpy(bool)
    isY = S.res.eq("Y").to_numpy()
    out = {}
    for name, rmask in (("st", ~isY), ("y", isY)):
        A = np.full(M.shape, np.nan, np.float32)
        cols = np.where(rmask)[0]
        ref = np.where(rmask & in_u)[0]
        for i in range(M.shape[0]):
            row = np.asarray(M[i], np.float32)
            r = row[ref]
            r = r[np.isfinite(r)]
            if len(r) == 0:
                continue
            lo, hi = r.min(), r.max()
            v = row[cols]
            A[i, cols] = 1.0 if hi == lo else (v - lo) / (hi - lo)
        out[name] = A
    np.savez(p, site_id=S.site_id.to_numpy().astype(str), st=out["st"], y=out["y"])
    return out["st"], out["y"]


# ----------------------------------------------------------------------------------------------- label-derived
@register
class SelphiPSSM(LabelDerivedFeature):
    """`score`: SELPHI's PSP-derived kinase PSSM (pssm.py: per-kinase PSSM from known substrates, family PSSM for
    kinases with few substrates). Substituted by our F2a learned PSSM (ksp.features.motif_learned.fit_pssm, family
    prior, no phospho channel, alpha/beta tuned by grouped internal CV) - same information, cross-fitted."""
    name = "b4_selphi_pssm"

    def fit(self, train_pairs: pd.DataFrame, ctx=None) -> "SelphiPSSM":
        from ksp.features.motif_learned import fit_pssm
        groups = None
        if ctx is not None and getattr(ctx, "clusters", None) is not None and len(ctx.clusters):
            acc = train_pairs.site_id.str.rsplit("_", n=1).str[0]
            groups = pd.Series(acc.map(ctx.clusters).fillna(acc).astype(str).to_numpy(), index=train_pairs.site_id)
            groups = groups[~groups.index.duplicated()]
        self.sites = load_sites()
        self.model = fit_pssm(train_pairs[["kinase_acc", "site_id"]], self.sites, load_kinases(), "family", False,
                              groups=groups, seed=getattr(ctx, "seed", 0) if ctx is not None else 0)
        return self

    def transform(self, query: pd.DataFrame, ctx=None) -> pd.DataFrame:
        from ksp.features.motif_learned import score_pairs
        v = score_pairs(self.model, query.kinase_idx.to_numpy(int), query.site_idx.to_numpy(int), self.sites)
        return pd.DataFrame({"score": v})


def _coreg_blocks() -> list[tuple[str, np.ndarray, np.ndarray]]:
    """(column name, site-index rows, X conditions) for the four SELPHI co-regulation data sets:
    coreg293 <- CPTAC BRCA (substitute for Mertins 2016 BRCA; our BCM v1.2 block, protein-residualised);
    coreg.{ntera,hl,mcf7} <- Hijazi 2020 (PMID 31959955) 60-min inhibitor log2FC per cell line (compendium block,
    robust-SD scaling undone)."""
    from ksp.features.phospho import load_blocks, load_conditions
    bl = {b.block_id: b for b in load_blocks(kinds=("compendium", "cptac"))
          if b.block_id in ("pmid:31959955", "cptac:BRCA")}
    out = []
    b = bl["cptac:BRCA"]
    out.append(("coreg293", b.rows, np.asarray(b.X, np.float32)))
    h = bl["pmid:31959955"]
    cond = load_conditions().set_index("conditionID")
    sc = cond.scale.reindex(h.cols).to_numpy(np.float32)
    samp = cond.Sample.reindex(h.cols).to_numpy()
    X = np.asarray(h.X, np.float32) * sc[None, :]
    for name, cl in (("coreg.ntera", "NTERA2"), ("coreg.hl", "HL60"), ("coreg.mcf7", "MCF7")):
        keep = samp == cl
        Xc = X[:, keep]
        ok = ~np.isnan(Xc).all(1)
        out.append((name, h.rows[ok], Xc[ok]))
    return out


def _row_ranks(X: np.ndarray) -> np.ndarray:
    """Average ranks per row over observed entries (NaN kept)."""
    R = stats.rankdata(np.where(np.isnan(X), np.inf, X), axis=1, method="average").astype(np.float32)
    R[np.isnan(X)] = np.nan
    return R


def masked_spearman_logp(Rx: np.ndarray, Ra: np.ndarray, min_n: int = 5) -> np.ndarray:
    """-log10 two-sided p of the Pearson correlation of pre-ranked rows over pairwise-complete columns
    (approximate Spearman: ranks are taken over each row's own observed set). Rx: S x C, Ra: K x C -> S x K."""
    Mx, Ma = (~np.isnan(Rx)).astype(np.float32), (~np.isnan(Ra)).astype(np.float32)
    X0, A0 = np.nan_to_num(Rx), np.nan_to_num(Ra)
    n = Mx @ Ma.T
    sx, sa = X0 @ Ma.T, Mx @ A0.T
    sxx, saa = (X0 ** 2) @ Ma.T, Mx @ (A0 ** 2).T
    sxa = X0 @ A0.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = n * sxa - sx * sa
        var = (n * sxx - sx ** 2) * (n * saa - sa ** 2)
        r = cov / np.sqrt(var)
        r = np.clip(r, -0.999999, 0.999999)
        df = n - 2
        t = r * np.sqrt(df / (1 - r ** 2))
        lp = -(stats.t.logsf(np.abs(t), df) + np.log(2)) / np.log(10)
    lp[(n < min_n) | ~np.isfinite(lp)] = np.nan
    return lp.astype(np.float32)


@lru_cache(maxsize=1)
def _coreg_data():
    blocks = _coreg_blocks()
    out = []
    for name, rows, X in blocks:
        mu = np.nanmean(X, 0)
        sd = np.nanstd(X, 0)
        out.append((name, rows, X, mu, sd, _row_ranks(X)))
    return out


@register
class SelphiKSEACoreg(LabelDerivedFeature):
    """coreg293 / coreg.ntera / coreg.hl / coreg.mcf7: -log10 p of the Spearman correlation between a kinase's KSEA
    activity profile (regulon = the kinase's TRAINING substrates observed in the data set) and the site's
    phosphorylation profile across conditions of that data set (their cutillas_act-cor.r / calc-cptac-act.r).
    KSEA activity: z = (mean regulon log2FC - mean all) * sqrt(m) / sd all per condition (Casado 2013 z-score;
    substitutes their permutation-based ksea_batchKinases, trial=1000). NaN (kinase without regulon, < 5 shared
    conditions, site not quantified) -> 0 later, as SELPHI's NA -> 0."""
    name = "b4_selphi_coreg"
    MIN_N = 5

    def fit(self, train_pairs: pd.DataFrame, ctx=None) -> "SelphiKSEACoreg":
        S = load_sites()
        spos = pd.Series(np.arange(len(S)), index=S.site_id)
        K = load_kinases()
        kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc)
        tp = train_pairs[train_pairs.kinase_acc.isin(kpos.index) & train_pairs.site_id.isin(spos.index)]
        ki = kpos[tp.kinase_acc].to_numpy()
        si = spos[tp.site_id].to_numpy()
        self.nK = len(K)
        self.act = []
        for name, rows, X, mu, sd, R in _coreg_data():
            rpos = pd.Series(np.arange(len(rows)), index=rows)
            hit = rpos.reindex(si).to_numpy()
            ok = ~np.isnan(hit)
            A = np.full((len(K), X.shape[1]), np.nan, np.float32)
            if ok.any():
                kk, rr = ki[ok], hit[ok].astype(int)
                df = pd.DataFrame({"k": kk, "r": rr}).drop_duplicates()
                for k, g in df.groupby("k"):
                    sub = X[g.r.to_numpy()]
                    m = (~np.isnan(sub)).sum(0)
                    with np.errstate(invalid="ignore"):
                        z = (np.nanmean(sub, 0) - mu) * np.sqrt(m) / sd
                    z[m == 0] = np.nan
                    A[k] = z
            self.act.append((name, rows, R, A))
        return self

    def transform(self, query: pd.DataFrame, ctx=None) -> pd.DataFrame:
        ki = query.kinase_idx.to_numpy(int)
        si = query.site_idx.to_numpy(int)
        out = {}
        for name, rows, R, A in self.act:
            v = np.full(len(query), np.nan, np.float32)
            rpos = pd.Series(np.arange(len(rows)), index=rows)
            hit = rpos.reindex(si).to_numpy()
            has_k = ~np.isnan(A).all(1)
            ok = ~np.isnan(hit) & has_k[ki]
            if ok.any():
                us, inv_s = np.unique(hit[ok].astype(int), return_inverse=True)
                uk, inv_k = np.unique(ki[ok], return_inverse=True)
                L = masked_spearman_logp(R[us], _row_ranks(A[uk]), self.MIN_N)
                v[ok] = L[inv_s, inv_k]
            out[name] = v
        return pd.DataFrame(out)
