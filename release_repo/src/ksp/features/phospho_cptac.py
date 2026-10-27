"""CPTAC (BCM harmonised v1.2, LinkedOmics S3 mirror) tumour phospho blocks and F6c kinase-state correlation.

Residualisation (per cohort, per site, tumour samples only):
    site_log2 ~ 1 + host_protein_log2 (same ENSG, same sample) + confounders
Confounders available (PDC pan-cancer Clinical_meta_data_v1; the S3 *_meta.txt files have none):
    ESTIMATE_StromalScore, ESTIMATE_ImmuneScore (RNA-based), WES_purity (WES; missing for some samples -> cohort
    median imputed, plus a missingness indicator if >0 and <100 % missing). No TMT plex / batch and no ischemia time are
    in either metadata file (documented in docs/features_phospho.md).
Sites whose host protein is not quantified in the cohort are dropped (NaN), so every CPTAC value is protein-residualised.

F6c (not label-derived): per cohort Pearson r between a site's residualised profile and
    (i) the kinase's protein abundance (residualised on confounders only), and
    (ii) the kinase's activation-loop phospho activity profile = mean of the z-scored residualised profiles of the
         observed phosphosites in the A-loop (residues strictly between DFG and APE of the UniProt kinase domain).
Fisher z averaged over cohorts (weight n-3) with >= min_n paired tumours per cohort. Sites on the kinase's own protein
are set to NaN (trivial/self).
"""
from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from .phospho import Block, ROOT, masked_corr, fisher_z, to_dev

EXT = ROOT / "external"
S3 = EXT / "cptac_pancan_v1.2"
PDC_META = EXT / "cptac_pancancer_pdc" / "Clinical_meta_data_v1.zip"
COHORTS = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
META_NAME = {"OV": "HGSC", "CCRCC": "ccRCC"}
CONFOUNDERS = ["ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore", "WES_purity"]


def read_matrix(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    d = d.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    return d.mask(~np.isfinite(d))  # the BCM phospho files contain 'Inf' entries (e.g. 551 in BRCA) -> missing


def read_confounders(cohort: str) -> tuple[pd.DataFrame, dict]:
    """Numeric confounder table indexed by case id (tumour rows); report dict of availability."""
    z = zipfile.ZipFile(PDC_META)
    name = f"Clinical_meta_data_v1/{META_NAME.get(cohort, cohort)}_meta.txt"
    d = pd.read_csv(io.BytesIO(z.read(name)), sep="\t", low_memory=False, dtype=str)
    d = d[d.idx != "data_type"].set_index("idx")
    rep = {}
    out = {}
    for c in CONFOUNDERS:
        if c not in d.columns:
            rep[c] = "absent"
            continue
        v = pd.to_numeric(d[c], errors="coerce")
        rep[c] = f"{v.notna().mean():.2f} non-missing"
        out[c] = v
    return pd.DataFrame(out), rep


def design_confounders(conf: pd.DataFrame, samples: list[str]) -> tuple[np.ndarray, list[str]]:
    """Standardised confounder design (n_samples x q) for the given samples; median-imputed, + missingness flags."""
    C = conf.reindex(samples)
    cols, names = [], []
    for c in C.columns:
        v = C[c].to_numpy(dtype=float)
        miss = np.isnan(v)
        if miss.all() or np.nanstd(v) == 0:
            continue
        v = np.where(miss, np.nanmedian(v), v)
        cols.append((v - v.mean()) / (v.std() + 1e-12))
        names.append(c)
        if 0 < miss.mean() < 1 and miss.sum() >= 3:
            cols.append(miss.astype(float) - miss.mean())
            names.append(c + "_missing")
    if not cols:
        return np.zeros((len(samples), 0)), []
    return np.stack(cols, 1), names


def residualise(Y: np.ndarray, P: np.ndarray | None, C: np.ndarray, min_obs: int = 10, device="cuda",
                ridge: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
    """Per-row OLS residuals of Y (s x n, NaN missing) on [1, P_row (s x n, optional), C (n x q)].
    Rows with < max(min_obs, n_params + 5) usable samples -> all NaN. Returns (residuals float32, n_used)."""
    import torch
    dev = torch.device(device)
    Yt = torch.as_tensor(Y, dtype=torch.float32, device=dev)
    s, n = Yt.shape
    Ct = torch.as_tensor(C, dtype=torch.float32, device=dev)
    parts = [torch.ones((s, n, 1), device=dev)]
    M = ~torch.isnan(Yt)
    if P is not None:
        Pt = torch.as_tensor(P, dtype=torch.float32, device=dev)
        M &= ~torch.isnan(Pt)
        Pc = torch.nan_to_num(Pt - torch.nanmean(Pt, 1, keepdim=True), nan=0.0)
        parts.append(Pc[:, :, None])
    if Ct.shape[1]:
        parts.append(Ct[None].expand(s, n, Ct.shape[1]))
    D = torch.cat(parts, 2)  # s x n x p
    p = D.shape[2]
    Mf = M.float()
    Y0 = torch.nan_to_num(Yt, nan=0.0)
    Dm = D * Mf[:, :, None]
    XtX = Dm.transpose(1, 2) @ D + ridge * torch.eye(p, device=dev)[None]
    Xty = (Dm.transpose(1, 2) @ Y0[:, :, None])
    beta = torch.linalg.solve(XtX, Xty)
    R = Y0 - (D @ beta)[:, :, 0]
    nused = Mf.sum(1)
    need = max(min_obs, p + 5)
    R = torch.where(M & (nused[:, None] >= need), R, torch.full_like(R, float("nan")))
    return R.cpu().numpy().astype(np.float32), nused.cpu().numpy().astype(np.int32)


def cohort_site_map() -> pd.DataFrame:
    m = pd.read_parquet(ROOT / "cache" / "sites" / "mapping_log.parquet", columns=["source", "src_id", "site_id"])
    m = m[(m.source == "cptac") & m.site_id.notna()]
    return m[["src_id", "site_id"]]


def build_cohort_block(cohort: str, site_pos: dict, idmap: pd.DataFrame, min_obs: int = 10, device="cuda",
                       residualise_protein: bool = True) -> tuple[Block, dict, pd.DataFrame, np.ndarray, list]:
    """Residualised tumour phospho block for one cohort.
    Returns (block, report, protein matrix (ENSG x samples, raw log2), confounder design, sample list)."""
    ph = read_matrix(S3 / cohort / f"{cohort}_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt")
    pr = read_matrix(S3 / cohort / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
    samples = [s for s in ph.columns if s in pr.columns]
    ph, pr = ph[samples], pr[samples]
    conf, rep_conf = read_confounders(cohort)
    C, cnames = design_confounders(conf, samples)
    n_conf_samples = int(conf.reindex(samples).notna().any(axis=1).sum())
    # map rows -> canonical site; average duplicates (several BCM rows -> one canonical site)
    sid = idmap.set_index("src_id").site_id.reindex(ph.index)
    ok = sid.notna().to_numpy()
    ph = ph[ok]
    sid = sid[ok]
    ensg = pd.Series([i.split("|")[0] for i in ph.index], index=ph.index)
    ph = ph.assign(_site=sid.to_numpy(), _ensg=ensg.to_numpy())
    g = ph.groupby(["_site", "_ensg"], sort=False)[samples].mean()
    g = g.reset_index()
    # one ENSG per site (keep the first if a site maps from two genes)
    g = g.drop_duplicates("_site")
    Y = g[samples].to_numpy(np.float32)
    P = pr.reindex(g._ensg.to_numpy()).to_numpy(np.float32)
    has_prot = np.isfinite(P).sum(1) >= min_obs
    Y, P, g = Y[has_prot], P[has_prot], g[has_prot]
    R, nused = residualise(Y, P if residualise_protein else None, C, min_obs=min_obs, device=device)
    keep = np.isfinite(R).sum(1) >= min_obs
    rows = g._site.map(site_pos).to_numpy()
    keep &= ~pd.isna(rows)
    rows = rows[keep].astype(np.int32)
    R = R[keep]
    order = np.argsort(rows)
    blk = Block(f"cptac:{cohort}", "cptac", rows[order], R[order], samples, (), cohort)
    rep = dict(cohort=cohort, n_samples=len(samples), n_samples_with_confounders=n_conf_samples,
               confounders_used=";".join(cnames), confounder_availability=str(rep_conf),
               n_phospho_rows=int(len(ok)), n_rows_mapped=int(ok.sum()), n_sites=int(len(has_prot)),
               n_sites_with_host_protein=int(has_prot.sum()), n_sites_kept=int(keep.sum()),
               median_obs=float(np.median(np.isfinite(R).sum(1))) if len(R) else 0.0,
               plex_or_batch="not available", ischemia="not available")
    return blk, rep, pr, C, samples


# ----------------------------------------------------------------------------------------------- A-loop
ALOOP_START = [r"DFG", r"D[FLWY]G", r"D[A-Z]G"]
ALOOP_END = [r"APE", r"[ASPTGC]PE", r"[A-Z]PE", r"[A-Z]P[A-Z]"]


def activation_loop(seq: str, dom_start: int, dom_end: int) -> tuple[int, int, str] | None:
    """(start, end, how) 1-based residue range strictly between the DFG and APE motifs inside the kinase domain.
    DFG searched in the domain (first strict 'DFG', else first D[FLWY]G, else D.G) after the first 35 % of the domain
    (avoids N-lobe false hits); APE = first match of the most specific pattern within 12..45 residues after DFG."""
    if not seq or not dom_start or not dom_end:
        return None
    s0, e0 = int(dom_start) - 1, int(dom_end)
    dom = seq[s0:e0]
    off = int(len(dom) * 0.35)
    for i, pat in enumerate(ALOOP_START):
        m = re.search(pat, dom[off:])
        if not m:
            continue
        d = s0 + off + m.start()  # 0-based index of D
        for j, pe in enumerate(ALOOP_END):
            win = seq[d + 3 + 12: min(len(seq), d + 3 + 45)]
            m2 = re.search(pe, win)
            if m2:
                a = d + 3 + 12 + m2.start()  # 0-based start of APE
                return d + 4, a, f"{pat}..{pe}"  # 1-based: residue after G .. residue before A
        return None
    return None


def zrows(X: np.ndarray) -> np.ndarray:
    mu = np.nanmean(X, 1, keepdims=True)
    sd = np.nanstd(X, 1, keepdims=True)
    return (X - mu) / np.where(sd > 0, sd, np.nan)


def f6c_correlations(blocks: list[Block], kin_profiles: dict, n_kinases: int, n_sites: int, own_protein_rows: dict,
                     min_n: int = 20, device="cuda", chunk: int = 32768) -> tuple[np.ndarray, np.ndarray]:
    """Mean Fisher z (weight n-3) over cohorts of corr(site residual profile, kinase profile).
    kin_profiles: {cohort: (kinase_idx array, profile matrix n_k x n_samples)}. own_protein_rows: {kinase_idx: site
    idx array on that kinase's protein} -> NaN. Returns (mean_z (n_kinases x n_sites) float32, n_cohorts uint8)."""
    import torch
    dev = torch.device(device)
    sw = np.zeros((n_kinases, n_sites), np.float32)
    szw = np.zeros((n_kinases, n_sites), np.float32)
    nc = np.zeros((n_kinases, n_sites), np.uint8)
    for b in blocks:
        if b.cohort not in kin_profiles:
            continue
        kidx, K = kin_profiles[b.cohort]
        if len(kidx) == 0:
            continue
        Kt = to_dev(K, dev)
        MB = (~torch.isnan(Kt)).float()
        B0 = torch.nan_to_num(Kt, nan=0.0)
        prep = (MB, B0, B0 * B0)
        for c0 in range(0, len(b.rows), chunk):
            A = to_dev(b.X[c0:c0 + chunk], dev)
            r, n = masked_corr(A, None, min_n=min_n, prepared=prep)
            z = fisher_z(r)
            w = (n - 3).clamp(min=0)
            ok = torch.isfinite(z)
            zw = torch.where(ok, z * w, torch.zeros_like(z)).T.cpu().numpy()
            ww = torch.where(ok, w, torch.zeros_like(w)).T.cpu().numpy()
            okn = ok.T.cpu().numpy().astype(np.uint8)
            rows = b.rows[c0:c0 + chunk]
            szw[np.ix_(kidx, rows)] += zw
            sw[np.ix_(kidx, rows)] += ww
            nc[np.ix_(kidx, rows)] += okn
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(nc > 0, szw / np.maximum(sw, 1e-9), np.nan).astype(np.float32)
    for k, rows in own_protein_rows.items():
        out[k, rows] = np.nan
    return out, nc
