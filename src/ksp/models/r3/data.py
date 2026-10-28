"""R3 inputs (DESIGN §8.4, §5, §16 A6). Owner: R3-neural.

Site axis  : U_ms sites in cache/r3/win_sites.parquet order (explicit site_id array; matched to cache/index/sites.parquet
             by ID, never by position). Window embeddings: cache/r3/win_emb.npy rows via win_rows.npy (+-10; -1 = beyond
             a terminus). Built by scripts/models/r3_window_embed.py.
Kinase axis: ksp.canary.kinases_with_canaries() order (503 real kinases = cache/index/kinases.parquet, then 13 canaries)
             = the kinase_idx order of ksp.eval.map_predictions.
Targets    : int8 (n_sites, n_kinases): 1 = training positive (L2, U_ms, training proteins), 0 = candidate negative
             (residue class ok, not masked), -1 = ignore (masked by KP(k) / Hamming<=2 / Par(k) via
             ksp.labels.NegativeMasks built from TRAINING-protein labels only, out of class, or not a training row).
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ksp import canary, index, labels, splits, store

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
R3CACHE = ROOT / "cache" / "r3"
FLANK_MAX = 10
ESM_CACHES = {"650m": (ROOT / "cache" / "r3", 10), "3b": (ROOT / "cache" / "r3" / "esm3b", 7),
              # G5 deployment: identical 650M settings, windows for all U_obs sites (training rows stay U_ms)
              "650m_uobs": (ROOT / "cache" / "r3" / "uobs", 10)}
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}      # 0 = pad / unknown(X etc. -> 21)
N_AA_TOK = 22                                       # pad, 20 aa, other


@dataclass
class SiteInputs:
    site_ids: np.ndarray        # (N,) str
    site_idx: np.ndarray        # (N,) int64 position in cache/index/sites.parquet
    res_code: np.ndarray        # (N,) 0=S,1=T,2=Y
    acc: np.ndarray             # (N,) protein accession
    cv_fold: np.ndarray         # (N,) CV fold of the protein, -1 = TEST protein (never used here)
    rows: np.ndarray            # (N, 21) row into emb, -1 pad
    emb: np.ndarray             # (R, 1280) float16 (memmap)
    res_aa: np.ndarray          # (R,) aa token per embedding row
    res_phos: np.ndarray        # (R,) 1 if that residue is itself a U_ms phosphosite
    prot_row: np.ndarray        # (N,) row into prot_mean
    prot_mean: np.ndarray       # (P, 1280) float32
    info: dict
    in_ums: np.ndarray | None = None   # (N,) bool; all True except for the '650m_uobs' deployment cache


def load_sites(esm: str = "650m") -> SiteInputs:
    """esm: '650m' (cache/r3, +-10) or '3b' (cache/r3/esm3b, +-7; ESM-2 3B per-residue windows)."""
    R3CACHE, flank_max = ESM_CACHES[esm]
    ws = pd.read_parquet(R3CACHE / "win_sites.parquet").site_id.values
    ids_arr = np.load(R3CACHE / "win_sites.site_ids.npy")
    if not (np.array(list(ws), dtype="U") == ids_arr).all():
        raise RuntimeError("win_sites.parquet and win_sites.site_ids.npy disagree")
    rows = np.load(R3CACHE / "win_rows.npy")
    emb = np.load(R3CACHE / "win_emb.npy", mmap_mode="r")
    res = pd.read_parquet(R3CACHE / "win_res.parquet")
    if len(res) != emb.shape[0] or rows.shape[0] != len(ws) or rows.shape[1] != 2 * flank_max + 1:
        raise RuntimeError("R3 window cache shapes inconsistent; rebuild with scripts/models/r3_window_embed.py")
    S = index.load_sites()
    spos = pd.Series(np.arange(len(S)), index=S.site_id.values)
    sidx = spos.reindex(ws)
    if sidx.isna().any():
        raise KeyError(f"{sidx.isna().sum()} R3 window sites not in canonical site index")
    sidx = sidx.values.astype(np.int64)
    Sx = S.iloc[sidx]
    if esm.endswith("_uobs"):
        if len(ws) != len(S) or len(np.unique(sidx)) != len(S):
            raise RuntimeError(f"U_obs window cache has {len(ws)} sites but U_obs has {len(S)}; rebuild")
    else:
        if not Sx.in_U_ms.all():
            raise RuntimeError("R3 window sites must all be U_ms")
        n_ums = int(S.in_U_ms.sum())
        if len(ws) != n_ums:
            raise RuntimeError(f"R3 window cache has {len(ws)} sites but U_ms has {n_ums}; rebuild")
    res_code = Sx.res.map({"S": 0, "T": 1, "Y": 2}).values.astype(np.int64)
    acc = Sx.acc.values
    folds = splits.load_cv_folds()
    cvf = pd.Series(acc).map(folds).fillna(-1).astype(np.int64).values
    test_prot = splits.load_test_proteins()
    if ((cvf == -1) != pd.Series(acc).isin(test_prot).values).any():
        raise RuntimeError("protein without CV fold that is not a TEST protein")
    res_aa = np.array([AA_IDX.get(a, 21) for a in res.aa.values], dtype=np.int64)
    Su = S[S.in_U_ms]                           # phospho channel = flanking U_ms phosphosites (not U_obs)
    ums_key = set(zip(Su.acc.values, Su.pos.values))
    res_phos = np.fromiter(((a, p) in ums_key for a, p in zip(res.acc.values, res.pos.values)), bool,
                           len(res)).astype(np.int64)
    pm, desc = store.load_matrix("esm2_650m/protein_mean", mmap=False)
    pacc = index.load_proteins().acc.values
    prow = pd.Series(np.arange(len(pacc)), index=pacc).reindex(acc)
    if prow.isna().any():
        raise KeyError("site protein missing from protein index")
    meta = json.loads((R3CACHE / "win_embed_meta.json").read_text()) if (R3CACHE / "win_embed_meta.json").exists() \
        else {}
    return SiteInputs(ws.astype(str), sidx, res_code, acc, cvf, rows, emb, res_aa, res_phos,
                      prow.values.astype(np.int64), np.nan_to_num(np.asarray(pm, np.float32)),
                      {"n_sites": len(ws), "n_residues": len(res), "embed_meta": meta, "esm": esm,
                       "flank_max": flank_max,
                       "protein_mean_nan_rows": int(np.isnan(np.asarray(pm)).any(1).sum())},
                      Sx.in_U_ms.values.astype(bool))


# ------------------------------------------------------------------------------------------------ kinase descriptors
@dataclass
class KinaseInputs:
    kinases: pd.DataFrame       # kinases_with_canaries()
    esm: np.ndarray             # (K, 1280) standardised, 0 where missing
    esm_missing: np.ndarray     # (K,)
    cantley: np.ndarray         # (K, 2*n_pc) PCA scores (ST | Y), 0 where missing
    cantley_missing: np.ndarray  # (K, 2) missing flags (ST, Y)
    cantley_imputed: np.ndarray  # (K, 2) 1 = copied from nearest atlas kinase (distance in cantley_dist)
    cantley_dist: np.ndarray    # (K, 2)
    group: np.ndarray           # (K,) group id
    family: np.ndarray          # (K,) family id (-1 = none / unique)
    group_names: list
    family_names: list
    class_ok: np.ndarray        # (3, K) bool: residue code (S,T,Y) x kinase
    is_canary: np.ndarray
    info: dict


def _cantley_profiles(n_ref_st: int = 3000, n_ref_y: int = 1000, n_pc: int = 32, seed: int = 0):
    """Kinase Cantley descriptor = the kinase's (own, else nearest-atlas-imputed) unprimed log2 Cantley score profile
    over a fixed random reference set of U_ms sites (label-free; seed 0), z-scored per kinase, then PCA (fit over
    kinases having that block). Built from the F1 store only (motif/cantley_unprimed_score + imputed)."""
    S = index.load_sites()
    rng = np.random.default_rng(seed)
    ums = S.in_U_ms.values
    st = np.where(ums & S.res.isin(["S", "T"]).values)[0]
    yy = np.where(ums & (S.res.values == "Y"))[0]
    ref_st = np.sort(rng.choice(st, n_ref_st, replace=False))
    ref_y = np.sort(rng.choice(yy, n_ref_y, replace=False))
    own, _ = store.load_matrix("motif/cantley_unprimed_score")
    imp, _ = store.load_matrix("motif/cantley_imputed_unprimed_score")
    dist, ddesc = store.load_matrix("motif/cantley_imputed_distance")
    cols = list(np.load(store.STORE / "motif" / "cantley_imputed_distance.cols.npy"))
    out, miss, impf, dd = [], [], [], []
    for b, ref in enumerate([ref_st, ref_y]):
        A = np.asarray(own[:, ref], np.float32)
        B = np.asarray(imp[:, ref], np.float32)
        own_ok = np.isfinite(A).mean(1) > 0.5
        imp_ok = np.isfinite(B).mean(1) > 0.5
        X = np.where(own_ok[:, None], A, B)
        ok = own_ok | imp_ok
        X = np.nan_to_num(X, nan=0.0)
        mu = X.mean(1, keepdims=True)
        sd = X.std(1, keepdims=True) + 1e-6
        Z = (X - mu) / sd
        Zc = Z[ok] - Z[ok].mean(0)
        U, s, Vt = np.linalg.svd(Zc, full_matrices=False)
        pc = (Z - Z[ok].mean(0)) @ Vt[:n_pc].T
        pc /= (pc[ok].std(0) + 1e-6)
        pc[~ok] = 0
        out.append(pc.astype(np.float32))
        miss.append(~ok)
        impf.append(~own_ok & imp_ok)
        dd.append(np.nan_to_num(np.asarray(dist[:, cols.index(["ST", "Y"][b])], np.float32), nan=0.0))
    return (np.concatenate(out, 1), np.stack(miss, 1).astype(np.float32), np.stack(impf, 1).astype(np.float32),
            np.stack(dd, 1), {"n_ref_st": n_ref_st, "n_ref_y": n_ref_y, "n_pc": n_pc,
                              "n_with_st": int((~miss[0]).sum()), "n_with_y": int((~miss[1]).sum()),
                              "n_imputed_st": int(impf[0].sum()), "n_imputed_y": int(impf[1].sum())})


def load_kinases() -> KinaseInputs:
    K = canary.kinases_with_canaries().reset_index(drop=True)
    real = index.load_kinases()
    nR = len(real)
    if not (K.uniprot_acc.values[:nR] == real.uniprot_acc.values).all():
        raise RuntimeError("kinases_with_canaries() does not start with the canonical kinase index")
    isc = canary.is_canary(K.uniprot_acc.values)
    E, desc = store.load_matrix("esm2_650m/kinase_domain", mmap=False)
    E = np.asarray(E, np.float32)
    em = np.isnan(E).any(1)
    mu, sd = E[~em].mean(0), E[~em].std(0) + 1e-6
    Ez = np.where(em[:, None], 0.0, (E - mu) / sd).astype(np.float32)
    cd = canary.canary_descriptors(dim=E.shape[1])            # random-but-fixed canary descriptors
    cd = cd.loc[K.uniprot_acc.values[nR:]].values.astype(np.float32)
    esm = np.concatenate([Ez, cd], 0)
    esm_missing = np.r_[em, np.zeros(len(cd), bool)].astype(np.float32)
    C, Cm, Ci, Cd, cinfo = _cantley_profiles()
    nC = len(K) - nR
    C = np.concatenate([C, np.zeros((nC, C.shape[1]), np.float32)])
    Cm = np.concatenate([Cm, np.ones((nC, 2), np.float32)])     # canaries: no Cantley matrix
    Ci = np.concatenate([Ci, np.zeros((nC, 2), np.float32)])
    Cd = np.concatenate([Cd, np.zeros((nC, 2), np.float32)])
    grp = K.manning_group.fillna("(none)").astype(str).values
    fam = np.array([labels._real_family(f) or "" for f in K.manning_family.values], dtype=object)
    fam = np.where(isc, K.manning_family.astype(str).values, fam)
    gnames = sorted(set(grp))
    fnames = sorted(set(f for f in fam if f))
    gid = np.array([gnames.index(g) for g in grp])
    fid = np.array([fnames.index(f) if f else -1 for f in fam])
    cls = K.residue_class.values
    class_ok = np.stack([[c in ("ST", "STY") for c in cls], [c in ("ST", "STY") for c in cls],
                         [c in ("Y", "STY") for c in cls]])
    if not set(cls) <= {"ST", "Y", "STY"}:
        raise ValueError(f"unexpected residue classes {set(cls)}")
    return KinaseInputs(K, esm, esm_missing, C, Cm, Ci, Cd, gid, fid, gnames, fnames, class_ok, isc,
                        {"esm_missing": int(em.sum()), "cantley": cinfo, "n_kinases": len(K), "n_canaries": nC,
                         "n_groups": len(gnames), "n_families": len(fnames)})


# ------------------------------------------------------------------------------------------------ targets
def build_targets(si: SiteInputs, ki: KinaseInputs, train_sites_mask: np.ndarray, train_kinases: np.ndarray | None,
                  variant: str = "L2", par: dict | None = None, log: dict | None = None,
                  loader=None) -> np.ndarray:
    """int8 target matrix (N_sites, K) for one training run.
    train_sites_mask: bool (N,) rows usable for training (incl. inner validation rows).
    train_kinases: bool (K,) kinases whose labels may be used (LKO); None = all.
    Positives = `variant` DEV labels (+canaries) on training rows. Masks = NegativeMasks over the unfiltered DEV label
    table ('all' + canaries) restricted to training proteins (and training kinases), so no eval-fold label informs
    the training targets. Kinases with no training positive get -1 everywhere (not trained)."""
    log = log if log is not None else {}
    K = ki.kinases
    kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc.values)
    train_prot = set(si.acc[train_sites_mask])
    loader = loader or (lambda v: labels.load_dev(v, with_canaries=True))
    allp = loader("all")
    allp = allp[allp.acc.isin(train_prot)]
    if train_kinases is not None:
        keep_k = set(K.uniprot_acc.values[train_kinases])
        allp = allp[allp.kinase_acc.isin(keep_k)]
    masks = labels.NegativeMasks(allp, paralogs=par if par is not None else labels.paralog_sets())
    # map index positions -> R3 site rows
    nS_index = len(index.load_sites())
    to_r3 = np.full(nS_index, -1, np.int64)
    to_r3[si.site_idx] = np.arange(len(si.site_idx))
    T = np.where(ki.class_ok[si.res_code], 0, -1).astype(np.int8)          # (N, K)
    T[~train_sites_mask] = -1
    pos = loader(variant)
    pos = pos[pos.acc.isin(train_prot)]
    if train_kinases is not None:
        pos = pos[pos.kinase_acc.isin(keep_k)]
    pr = pd.Series(np.arange(len(si.site_ids)), index=si.site_ids).reindex(pos.site_id)
    pk = kpos.reindex(pos.kinase_acc)
    ok = pr.notna().values & pk.notna().values
    log["pos_pairs_loaded"] = int(len(pos))
    log["pos_pairs_not_in_r3_sites"] = int((~pr.notna()).sum())
    log["pos_pairs_unknown_kinase"] = int((~pk.notna()).sum())
    pr, pk = pr.values[ok].astype(np.int64), pk.values[ok].astype(np.int64)
    trained = np.zeros(len(K), bool)
    trained[np.unique(pk)] = True
    n_masked = 0
    for k in np.where(trained)[0]:
        ex = to_r3[masks.excluded(K.uniprot_acc.values[k])]
        ex = ex[ex >= 0]
        n_masked += int((T[ex, k] == 0).sum())
        T[ex, k] = -1
    out_cls = ~ki.class_ok[si.res_code[pr], pk]
    log["pos_pairs_out_of_class"] = int(out_cls.sum())
    T[pr[~out_cls], pk[~out_cls]] = 1
    T[:, ~trained] = -1
    log["n_trained_kinases"] = int(trained.sum())
    log["n_pos_entries"] = int((T == 1).sum())
    log["n_masked_negatives"] = n_masked
    return T
