#!/usr/bin/env python
"""Export F0/F7 features into the ID-aligned feature store (src/ksp/store.py; canonical orders cache/index/).

Writes (cache/features/...):
  structure/site_structure        site x n_features float32 (col ids stored; col_axis='none')
  esm2_650m/site_pca64            site x 64  PCA of ESM-2 650M per-residue embeddings (PCA fit on U_ms sites only)
  esm2_3b/site_pca64              same for ESM-2 3B (if cache/esm2_3b exists)
  esm2_650m/protein_mean          protein x 1280 mean-pooled full-protein embedding (row_axis='protein')
  esm2_{650m,3b}/kinase_domain    kinase x D: primary kinase domain mean (domain-in-context; NaN if no domain)
  esm2_{650m,3b}/kinase_full      kinase x D: full-protein mean
No label-derived information is used. Usage: python -I scripts/export_structure_esm_features.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.features.esm_embed import load_sty  # noqa: E402
from ksp.index import load_kinases, load_proteins, load_sites  # noqa: E402
from ksp.store import reorder_to_index, save_matrix  # noqa: E402

NUM = ["plddt", "plddt_s5", "plddt_s10", "rsa", "sasa", "pse_12_70", "pse_12_180", "pse_24_70", "pse_24_180",
       "pse_24_180_s10", "idr", "high_acc", "disorder_plddt"]


def structure(sites):
    st = pd.read_parquet(ROOT / "cache/structure/site_structure.parquet")
    st = st[st.site_id.isin(set(sites.site_id))]
    X = st[NUM].astype("float64").to_numpy(np.float32)
    cols = list(NUM)
    for c in "HEC":
        X = np.c_[X, np.where(st.ss3.isna(), np.nan, (st.ss3 == c).astype(np.float32))]
        cols.append(f"ss3_{c}")
    X = np.c_[X, (st.af_status == "ok").astype(np.float32)]
    cols.append("af_ok")
    X = np.c_[X, (st.af_n_frags.fillna(0) > 1).astype(np.float32)]
    cols.append("af_fragmented")
    A = reorder_to_index(X, st.site_id.values, "site", 0)  # sites absent from the structure table -> NaN
    p = save_matrix("structure/site_structure", A, row_axis="site", col_axis="none", col_ids=cols,
                    meta=dict(source="cache/structure/site_structure.parquet (AFDB v6)",
                              doc="docs/features_structure_esm.md", label_derived=False))
    print(p, A.shape, "NaN frac", float(np.isnan(A).mean()))


def esm_sites(sites, model_dir, name, k=64):
    from sklearn.decomposition import PCA
    idx, M = load_sty(ROOT / "cache" / model_dir, sites.site_id.values)
    assert (idx.site_id.values == sites.site_id.values).all()
    M = M.astype(np.float32)
    ms = sites.in_U_ms.fillna(False).to_numpy(bool)
    pca = PCA(n_components=k, random_state=0).fit(M[ms])
    Z = pca.transform(M).astype(np.float32)
    p = save_matrix(f"{name}/site_pca64", Z, row_axis="site", col_axis="none", col_ids=[f"pc{i}" for i in range(k)],
                    meta=dict(model=model_dir, fit_on="U_ms sites (sites.in_U_ms)", n_fit=int(ms.sum()),
                              explained_variance=float(pca.explained_variance_ratio_.sum()), label_derived=False))
    np.save(ROOT / "cache/features" / name / "site_pca64.components.npy", pca.components_)
    np.save(ROOT / "cache/features" / name / "site_pca64.mean.npy", pca.mean_)
    print(p, Z.shape, "expl var", pca.explained_variance_ratio_.sum())


def protein_mean(model_dir, name):
    prots = load_proteins()
    pi = pd.read_parquet(ROOT / "cache" / model_dir / "protein_mean_index.parquet")
    V = np.load(ROOT / "cache" / model_dir / "protein_mean.npy")
    keep = pi.acc.isin(set(prots.acc)).values
    A = reorder_to_index(V[pi.row.values[keep]], pi.acc.values[keep], "protein", 0)
    p = save_matrix(f"{name}/protein_mean", A.astype(np.float32), row_axis="protein", col_axis="none",
                    col_ids=[f"d{i}" for i in range(A.shape[1])], meta=dict(model=model_dir, label_derived=False))
    print(p, A.shape, "NaN rows", int(np.isnan(A[:, 0]).sum()))


def kinase(model_dir, name):
    kin = load_kinases()
    kd = pd.read_parquet(ROOT / "cache" / model_dir / "kinases.parquet")
    kd = kd[kd.acc.isin(set(kin.uniprot_acc))]
    full = kd.drop_duplicates("acc")
    D = len(full.emb_full.iloc[0])
    F = reorder_to_index(np.stack(full.emb_full.values).astype(np.float32), full.acc.values, "kinase", 0)
    prim = kd[kd.is_primary].drop_duplicates("acc")
    K = reorder_to_index(np.stack(prim.emb_domain_ctx.values).astype(np.float32), prim.acc.values, "kinase", 0)
    cols = [f"d{i}" for i in range(D)]
    for nm, A, what in (("kinase_full", F, "full-protein mean"),
                        ("kinase_domain", K, "primary kinase domain mean, domain-in-context")):
        p = save_matrix(f"{name}/{nm}", A, row_axis="kinase", col_axis="none", col_ids=cols,
                        meta=dict(model=model_dir, what=what, domain_source="results/labels/kinase_list.tsv",
                                  label_derived=False))
        print(p, A.shape, "NaN rows", int(np.isnan(A[:, 0]).sum()))


def main():
    sites = load_sites()
    structure(sites)
    for model_dir, name in (("esm2_650m", "esm2_650m"), ("esm2_3b", "esm2_3b")):
        if (ROOT / "cache" / model_dir / "sty_index.parquet").exists():
            esm_sites(sites, model_dir, name)
            protein_mean(model_dir, name)
            kinase(model_dir, name)


if __name__ == "__main__":
    main()
