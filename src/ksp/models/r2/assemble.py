"""R2 pair-feature assembler (DESIGN §7, §10 CMP-1: gather per kinase chunk; never materialise all pairs).

Feature groups (names used by the ablation curve; see docs/model_r2.md):
  F1   Cantley (label-free): unprimed / primed-soft / primed-hard score + universe percentile, own matrix else the
       nearest-atlas imputed matrix, + imputation distance (0 = own matrix).
  F0   site-intrinsic: residue code, AlphaFold structure (18), ESM-2 650M PCA-64.
  F2a  learned PSSM (cross-fitted, xfit 'f2a').
  F3   PPI screens: per-source edge flags, any-edge, n_studies, RWR, RWR-dc (+pct), hop distance; + cross-fitted
       substrate neighbourhood (xfit 'f3nb').
  F4   localisation cosines (HPA, OpenCell, HCM, Hein, COMPARTMENTS-experiments, combined, NMF, Hein profile) +
       site topology (topo code, lumenal flag). COMPARTMENTS knowledge channel only with comp_know=True (ablation).
  F5   co-expression / co-dependency (HPA tissue, DepMap, DepMap-pc5).
  F6a  co-regulation (cross-fitted, xfit 'f6a' or 'f6a_clean').
  F6b  inhibitor response (ridge beta / z / contrast, kinase max-cos); clean=True -> *_f6_clean stores.
  F6c  CPTAC kinase protein / A-loop correlation.
  F7   kinase descriptors: residue class, Manning group, family (categorical), ESM-2 kinase-domain PCA-16, Cantley
       profile PCA (ST 8, Y 8). Canaries get their random ksp.canary descriptors in these slots.
  F8   Sugiyama 2019 (not available at build time -> group is empty; recorded).
Forbidden (A11 / A6 / DL2-7): no centre-site observation flags, no observation counts, no HPA reliability, no
loc_n_sources, no COMPARTMENTS knowledge (except the ablation), no PSP-derived anything outside xfit.
Canary kinases (index 503..515): kinase x site / kinase x protein stores have no rows for them -> NaN, except the
planted PPI flag of CANARY_POS (ppi_edge_any = 1 on ksp.canary.canary_ppi_partners(), as the canary design asks).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ksp import canary, index, store
from ksp.index import ROOT
from ksp.models.r2 import xfit

N_REAL = 503
GROUP_ORDER = ["F1", "F0", "F2a", "F3", "F4", "F6a", "F6b", "F6c", "F5", "F8", "F7", "Fknn"]
DESC_CACHE = ROOT / "cache" / "r2" / "kinase_desc.parquet"


def _ks(name):
    return ("ks", name)


class Assembler:
    def __init__(self, groups, fold: int | None, xfit_tag: str | None = None, clean: bool = False,
                 comp_know: bool = False, esm: str = "esm2_650m", skip_xfit: bool = False,
                 f2a_family: str = "f2a", desc_v: int = 1):
        self.groups = list(groups)
        unknown = set(self.groups) - set(GROUP_ORDER)
        if unknown:
            raise ValueError(f"unknown feature groups {unknown}")
        self.fold, self.clean, self.comp_know = fold, clean, comp_know
        self.xfit_tag = xfit_tag
        self.skip_xfit = skip_xfit
        self.f2a_family, self.desc_v = f2a_family, desc_v
        self.S = index.load_sites()
        self.K = canary.kinases_with_canaries()
        self.nK = len(self.K)
        P = index.load_proteins()
        self.site_prot = pd.Index(P.acc).get_indexer(self.S.acc)
        if (self.site_prot < 0).any():
            raise KeyError("sites on proteins missing from the protein index")
        self.cols: list[str] = []
        self.cat_cols: list[str] = []
        self.src: list[tuple] = []          # (kind, payload) per column
        self.fallbacks: dict = {}
        self._build()

    # ----------------------------------------------------------------------------------------- column registry
    def _add(self, name, kind, payload, cat=False):
        if kind == "k":
            payload = np.asarray(payload, np.float32)
            if len(payload) == N_REAL:                    # kinase vectors of real kinases only -> NaN for canaries
                payload = np.r_[payload, np.full(self.nK - N_REAL, np.nan, np.float32)]
            if len(payload) != self.nK:
                raise ValueError(f"{name}: kinase vector length {len(payload)}")
        self.cols.append(name)
        self.src.append((kind, payload))
        if cat:
            self.cat_cols.append(name)

    def _m(self, name):
        arr, desc = store.load_matrix(name)
        return arr

    def _build(self):
        g = set(self.groups)
        if "F1" in g:
            dist = self._m("motif/cantley_imputed_distance")        # 503 x 2 [ST, Y]
            for v in ["unprimed", "primed_soft", "primed_hard"]:
                for met in ["score", "pct_universe"]:
                    own = self._m(f"motif/cantley_{v}_{met}")
                    imp = self._m(f"motif/cantley_imputed_{v}_{met}")
                    self._add(f"f1_{v}_{met}", "ks_own_imp", (own, imp))
            own0 = self._m("motif/cantley_unprimed_score")
            imp0 = self._m("motif/cantley_imputed_unprimed_score")
            self._add("f1_impute_dist", "f1_dist", (own0, imp0, np.asarray(dist)))
        if "F0" in g:
            self._add("f0_res", "site_res", None, cat=True)
            st = self._m("structure/site_structure")
            scols = np.load(ROOT / "cache" / "features" / "structure" / "site_structure.cols.npy", allow_pickle=True)
            stv = np.asarray(st)
            for j, c in enumerate(scols):
                self._add(f"f0_{c}", "site", stv[:, j])
            esm = np.asarray(self._m("esm2_650m/site_pca64"))
            for j in range(esm.shape[1]):
                self._add(f"f0_esm{j:02d}", "site", esm[:, j])
        if "F2a" in g:
            self._xfit(self.f2a_family)
        if "F3" in g:
            for s in ["bioplex_293t", "bioplex_hct116", "huri", "opencell", "buljan", "hcm"]:
                self._add(f"f3_edge_{s}", "kp", self._m(f"context/ppi/ppi_edge_{s}"))
            self._add("f3_edge_any", "kp_edge_any", self._m("context/ppi/ppi_edge_any"))
            for c in ["ppi_n_studies", "ppi_rwr", "ppi_rwr_dc", "ppi_rwr_dc_pct", "ppi_dist"]:
                self._add(f"f3_{c[4:]}", "kp", self._m(f"context/ppi/{c}"))
            self._xfit("f3nb")
        if "F4" in g:
            for c in ["loc_cos_hpa", "loc_cos_opencell", "loc_cos_hcm", "loc_cos_hein", "loc_cos_comp_exp",
                      "loc_cos_combined", "loc_hcm_nmf_cos", "loc_hein_profile_corr"]:
                self._add(f"f4_{c}", "kp", self._m(f"context/loc/{c}"))
            if self.comp_know:
                for c in ["loc_cos_comp_know", "loc_cos_combined_know"]:
                    self._add(f"f4_{c}", "kp", self._m(f"context/loc/{c}"))
            self._add("f4_site_topo", "site", np.asarray(self._m("context/loc/site_topo_code")).astype(np.float32),
                      cat=True)
            self._add("f4_site_lumenal", "site", np.asarray(self._m("context/loc/site_lumenal")).astype(np.float32))
        if "F5" in g:
            for c in ["coexpr_hpa_tissue_spearman", "coexpr_depmap_pearson", "coexpr_depmap_pearson_pc5"]:
                self._add(f"f5_{c}", "kp", self._m(f"context/coexpr/{c}"))
        if "F6a" in g:
            self._xfit("f6a_clean" if self.clean else "f6a")
        if "F6b" in g:
            sfx = "_f6_clean" if self.clean else ""
            for c in ["f6b_ridge_beta", "f6b_ridge_z", "f6b_contrast"]:
                self._add(c, "ks", self._m(f"phospho/{c}{sfx}"))
            self._add("f6b_kinase_maxcos", "k", np.asarray(self._m(f"phospho/f6b_kinase_maxcos{sfx}")))
        if "F6c" in g:
            for c in ["f6c_prot", "f6c_aloop"]:
                self._add(c, "ks", self._m(f"phospho/{c}"))
        if "F8" in g:
            self.fallbacks["F8_unavailable"] = 1      # Sugiyama 2019 store not built; group contributes nothing
        if "Fknn" in g:
            self._xfit("knn")       # knn_sim (feat-knn), cross-fitted through ksp.crossfit by ksp.models.r2.xfit
        if "F7" in g:
            D = kinase_descriptors(version=self.desc_v)
            for c in D.columns:
                self._add(f"f7_{c}", "k", D[c].to_numpy(np.float32), cat=c in ("rclass", "group", "family"))

    def _xfit(self, fam):
        if self.skip_xfit:
            self.fallbacks[f'skip_xfit_{fam}'] = 1
            return
        arr, meta = xfit.load(fam, self.fold, tag=self.xfit_tag)
        if meta["n_kinases"] != self.nK or meta["n_sites"] != len(self.S):
            raise RuntimeError(f"xfit {fam} shape mismatch")
        for ci, c in enumerate(meta["columns"]):
            self._add(c if not fam.endswith("_clean") else c, "xfit", (arr, ci))

    # ----------------------------------------------------------------------------------------------- gathering
    def kinase_block(self, k: int, site_idx: np.ndarray, site_cache: dict | None = None) -> np.ndarray:
        """(len(site_idx), n_features) float32 for kinase index k (kinases_with_canaries order)."""
        site_idx = np.asarray(site_idx, np.int64)
        n = len(site_idx)
        X = np.empty((n, len(self.cols)), np.float32)
        prot = self.site_prot[site_idx]
        real = k < N_REAL
        for j, (kind, p) in enumerate(self.src):
            if kind == "ks":
                X[:, j] = p[k][site_idx] if real else np.nan
            elif kind == "ks_own_imp":
                if real:
                    a = np.asarray(p[0][k][site_idx], np.float32)
                    b = np.asarray(p[1][k][site_idx], np.float32)
                    X[:, j] = np.where(np.isnan(a), b, a)
                else:
                    X[:, j] = np.nan
            elif kind == "f1_dist":
                if real:
                    a = np.asarray(p[0][k][site_idx])
                    b = np.asarray(p[1][k][site_idx])
                    isY = self.S.res.values[site_idx] == "Y"
                    d = np.where(isY, p[2][k, 1], p[2][k, 0])
                    X[:, j] = np.where(~np.isnan(a), 0.0, np.where(~np.isnan(b), d, np.nan))
                else:
                    X[:, j] = np.nan
            elif kind == "kp":
                X[:, j] = p[k][prot] if real else np.nan
            elif kind == "kp_edge_any":
                if real:
                    X[:, j] = p[k][prot]
                else:
                    X[:, j] = 0.0
                    if self.K.uniprot_acc.values[k] == canary.POS_ID:
                        X[:, j] = np.isin(self.S.acc.values[site_idx], self._canary_partners()).astype(np.float32)
            elif kind == "site":
                X[:, j] = p[site_idx]
            elif kind == "site_res":
                r = self.S.res.values[site_idx]
                X[:, j] = np.where(r == "S", 0, np.where(r == "T", 1, 2))
            elif kind == "k":
                X[:, j] = p[k]
            elif kind == "xfit":
                arr, ci = p
                X[:, j] = arr[ci, k][site_idx]
            else:
                raise ValueError(kind)
        return X

    def _canary_partners(self):
        if not hasattr(self, "_cpp"):
            self._cpp = np.array(canary.canary_ppi_partners())
        return self._cpp

    def rows(self, kinase_idx, site_idx) -> np.ndarray:
        """Feature matrix for arbitrary (kinase_idx, site_idx) rows, gathered per kinase."""
        kinase_idx = np.asarray(kinase_idx, np.int64)
        site_idx = np.asarray(site_idx, np.int64)
        X = np.empty((len(kinase_idx), len(self.cols)), np.float32)
        order = np.argsort(kinase_idx, kind="stable")
        ks = kinase_idx[order]
        bounds = np.searchsorted(ks, np.arange(self.nK + 1))
        import os
        from concurrent.futures import ThreadPoolExecutor

        def one(k):
            a, b = bounds[k], bounds[k + 1]
            if b > a:
                X[order[a:b]] = self.kinase_block(k, site_idx[order[a:b]])
        with ThreadPoolExecutor(int(os.environ.get("R2_GATHER_THREADS", "6"))) as ex:
            list(ex.map(one, range(self.nK)))
        return X


# --------------------------------------------------------------------------------------------- F7 descriptors
def kinase_descriptors(n_esm: int = 16, n_cant: int = 8, rebuild: bool = False, version: int = 1) -> pd.DataFrame:
    """Label-free kinase descriptors for all kinases incl. canaries (index = kinases_with_canaries order).
    version 1: cache built with canary v1 metadata (all canaries in group CANARY); version 2: current ksp.canary
    (A25: one Manning group per canary)."""
    cache = DESC_CACHE if version == 1 else DESC_CACHE.with_name(f"kinase_desc_v{version}.parquet")
    if cache.exists() and not rebuild:
        return pd.read_parquet(cache)
    if version == 1:
        raise FileNotFoundError("descriptor v1 cache missing; it cannot be rebuilt after A25")
    from sklearn.decomposition import PCA
    K = canary.kinases_with_canaries()
    S = index.load_sites()
    out = pd.DataFrame(index=np.arange(len(K)))
    out["rclass"] = K.residue_class.map({"ST": 0, "Y": 1, "STY": 2}).astype(float).values
    grp = K.manning_group.fillna("(none)")
    fam = K.manning_family.fillna("(none)").where(~K.manning_family.fillna("").str.endswith("Unique"),
                                                   "(unique)")
    out["group"] = pd.factorize(grp, sort=True)[0].astype(float)
    out["family"] = pd.factorize(fam, sort=True)[0].astype(float)
    rng = np.random.default_rng(0)
    cdesc = canary.canary_descriptors()
    can_rows = np.arange(N_REAL, len(K))
    # ESM-2 kinase-domain PCA (fit on real kinases with a domain; others NaN)
    E = np.asarray(store.load_matrix("esm2_650m/kinase_domain")[0], np.float64)
    ok = np.isfinite(E).all(1)
    pca = PCA(n_esm, random_state=0).fit(E[ok])
    Z = np.full((len(K), n_esm), np.nan)
    Z[:N_REAL][ok] = pca.transform(E[ok])
    sd = Z[:N_REAL][ok].std(0)
    Z[can_rows] = cdesc.values[:, :n_esm] * sd
    for j in range(n_esm):
        out[f"esm_pc{j:02d}"] = Z[:, j]
    # Cantley score-profile PCA (own else imputed matrix), ST and Y separately, on random U_ms reference sites
    own = store.load_matrix("motif/cantley_unprimed_score")[0]
    imp = store.load_matrix("motif/cantley_imputed_unprimed_score")[0]
    off = n_esm
    for rt, res, nref in [("st", ("S", "T"), 20000), ("y", ("Y",), 8000)]:
        ref = np.where(S.in_U_ms.values & S.res.isin(res).values)[0]
        ref = np.sort(rng.choice(ref, min(nref, len(ref)), replace=False))
        M = np.asarray(own[:, ref], np.float64)
        Mi = np.asarray(imp[:, ref], np.float64)
        M = np.where(np.isnan(M), Mi, M)
        okr = np.isfinite(M).all(1)
        Mz = (M[okr] - M[okr].mean(1, keepdims=True)) / (M[okr].std(1, keepdims=True) + 1e-9)
        p = PCA(n_cant, random_state=0).fit(Mz)
        Zc = np.full((len(K), n_cant), np.nan)
        Zc[:N_REAL][okr] = p.transform(Mz)
        sdc = Zc[:N_REAL][okr].std(0)
        Zc[can_rows] = cdesc.values[:, off:off + n_cant] * sdc
        off += n_cant
        for j in range(n_cant):
            out[f"cant_{rt}_pc{j}"] = Zc[:, j]
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(cache)
    return out
