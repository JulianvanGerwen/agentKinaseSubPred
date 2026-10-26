"""F1 / B1: Cantley kinome-atlas motif scores (Johnson 2023 Ser/Thr; Yaron-Barir 2024 Tyr).

Clean-room reimplementation of the kinase-library 1.8.0 score definition. Only the *data*
files shipped in external/kinase_library (matrices, S/T favourability, kinome_information,
Ochoa scored reference parquets) are read; no package code is imported or executed.

Score definition (kinase-library 1.8.0 `PhosphoProteomics.score`, `Substrate.score`; Johnson 2023
Suppl. Note 2 in log2 space):
  * window: 15-mer, phosphoacceptor at index 7 (0-based); out-of-protein positions are '_'.
  * scored positions: S/T matrices -5..-1,+1..+4 (9 rows); Tyr matrices -5..-1,+1..+5 (10 rows).
    Positions -7,-6,(+5,)+6,+7 and the centre are not scored.
  * log2 score = sum_p log2M[p, aa_p], where log2M = log2(norm-matrix / random-aa frequency) is the
    shipped `*_all_log2_matrices.tsv`. Columns: 20 aa + 's','t','y' (phospho-residues). Any symbol
    not among these (padding '_', 'X', 'U', ...) contributes 0 (= random). Unprimed scoring
    (pp=False, the published default) upper-cases every flank, so flank S/T/Y use the S/T/Y columns.
  * the dot product is rounded to 3 decimals; for S/T kinases log2(st_fav[kinase, centre S|T]) is
    added (st_fav=True default) and the result rounded to 3 decimals again.
  * percentile ("as published") = 100 * #{reference scores <= score} / n_ref, rounded to 2 decimals,
    reference = kinase's scores over the Ochoa 2020 phosphoproteome shipped with kinase-library
    (82,735 S/T and 5,431 Tyr sites, unprimed, st_fav=True).

Extensions (not part of the published method; see docs/features_motif.md):
  * `phospho_w` (N x 15 weights in [0,1]): flank S/T/Y at position p contributes
    (1-w) * log2M[p, upper] + w * log2M[p, lower]. w in {0,1} = hard priming (kinase-library pp=True
    with those flanks lower-cased); fractional w = soft (expected) priming.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
KL_DB = ROOT / "external/kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases"

AA23 = "PGACSTVILMFYWHKRQNDEsty"          # kinase-library column order (aa_phos)
NULL = 23                                  # code for padding / unknown -> contributes 0
POSITIONS = {"ser_thr": [-5, -4, -3, -2, -1, 1, 2, 3, 4],
             "tyrosine": [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]}
KIN_TYPES = ("ser_thr", "tyrosine")
TYPE_RES = {"ser_thr": ("S", "T"), "tyrosine": ("Y",)}
CLASS_TYPES = {"ST": ("ser_thr",), "Y": ("tyrosine",), "STY": ("ser_thr", "tyrosine")}

_CODE = np.full(256, NULL, dtype=np.int8)
for _i, _a in enumerate(AA23[:20]):
    _CODE[ord(_a)] = _i
_LOWER_OF = np.arange(24, dtype=np.int8)   # code -> phospho code (S->s, T->t, Y->y), identity else
_LOWER_OF[AA23.index("S")] = 20
_LOWER_OF[AA23.index("T")] = 21
_LOWER_OF[AA23.index("Y")] = 22


@dataclass
class Atlas:
    kin_type: str
    names: list                      # matrix names (kinase-library MATRIX_NAME)
    uniprot: list                    # UniProt accession per matrix (kinome_information UNIPROT_ID)
    L: np.ndarray                    # K x P x 24 float64 log2 values; [..., 23] = 0
    st_fav_log2: np.ndarray | None   # K x 2 (S, T) or None for tyrosine
    positions: list = field(default_factory=list)

    def subset(self, idx) -> "Atlas":
        idx = list(idx)
        return Atlas(self.kin_type, [self.names[i] for i in idx], [self.uniprot[i] for i in idx],
                     self.L[idx], None if self.st_fav_log2 is None else self.st_fav_log2[idx], self.positions)


def kinome_info() -> pd.DataFrame:
    return pd.read_csv(KL_DB / "kinase_data/kinome_information.tsv", sep="\t", dtype=str)


def load_atlas(kin_type: str) -> Atlas:
    """Read the shipped `<type>_all_log2_matrices.tsv` (+ st_favorability.tsv) into an Atlas."""
    assert kin_type in KIN_TYPES
    m = pd.read_csv(KL_DB / f"matrices/{kin_type}/{kin_type}_all_log2_matrices.tsv", sep="\t", index_col=0)
    pos = POSITIONS[kin_type]
    cols = [f"{p}{a}" for p in pos for a in AA23]
    assert set(cols) <= set(m.columns), "unexpected matrix layout"
    L = np.zeros((len(m), len(pos), 24))
    L[:, :, :23] = m[cols].to_numpy(float).reshape(len(m), len(pos), 23)
    assert np.isfinite(L).all()
    info = kinome_info().set_index("MATRIX_NAME")
    names = list(m.index)
    missing = [n for n in names if n not in info.index]
    assert not missing, f"matrices without kinome_information row: {missing}"
    assert (info.loc[names, "TYPE"] == kin_type).all()
    uni = list(info.loc[names, "UNIPROT_ID"])
    stf = None
    if kin_type == "ser_thr":
        f = pd.read_csv(KL_DB / "matrices/ser_thr/st_favorability.tsv", sep="\t", index_col=0)
        f.columns = [c.upper() for c in f.columns]
        stf = np.log2(f.loc[names, ["S", "T"]].to_numpy(float))
    return Atlas(kin_type, names, uni, L, stf, pos)


def load_reference(kin_type: str) -> pd.DataFrame:
    """Published percentile reference: Ochoa-phosphoproteome scores shipped with kinase-library
    (index = 15-mer with lower-case centre, columns = matrix names, log2 scores rounded to 3 dp)."""
    return pd.read_parquet(KL_DB / f"substrates/ochoa/scored_phosprots/{kin_type}_phosphoproteome_scored.parquet")


def encode(windows, keep_lower_flanks: bool = False) -> np.ndarray:
    """15-mer strings -> N x 15 int8 codes into AA23 (NULL for '_', X, U, ...).
    Upper-cases everything (kinase-library pp=False) unless keep_lower_flanks, in which case flank
    's','t','y' keep their phospho codes (kinase-library pp=True input convention)."""
    w = np.asarray(windows, dtype="U15")
    assert all(len(x) == 15 for x in w[:1000]), "windows must be 15-mers"
    b = np.frombuffer(np.char.encode(w, "ascii"), dtype=np.uint8).reshape(len(w), 15)
    codes = _CODE[np.where((b >= 97) & (b <= 122), b - 32, b)]
    if keep_lower_flanks:
        for ch, c in (("s", 20), ("t", 21), ("y", 22)):
            codes[b == ord(ch)] = c
        codes[:, 7] = _CODE[b[:, 7] - 32 * ((b[:, 7] >= 97) & (b[:, 7] <= 122))]
    return codes


def center_residue(windows) -> np.ndarray:
    return np.char.upper(np.asarray([x[7] for x in windows], dtype="U1"))


def score(atlas: Atlas, codes: np.ndarray, center: np.ndarray, phospho_w: np.ndarray | None = None,
          st_fav: bool = True, published_rounding: bool = True, chunk: int = 50_000) -> np.ndarray:
    """N x K log2 scores (float64). `center` = array of 'S'/'T'/'Y' (must match atlas type).
    phospho_w: optional N x 15 float weights for treating flank S/T/Y as phospho-residues."""
    N, K = len(codes), len(atlas.names)
    out = np.empty((N, K))
    lo = _LOWER_OF[codes.astype(np.int64)]
    for a in range(0, N, chunk):
        b = min(N, a + chunk)
        acc = np.zeros((b - a, K))
        for pi, p in enumerate(atlas.positions):
            col = 7 + p
            up = atlas.L[:, pi, :][:, codes[a:b, col]].T          # (n, K)
            if phospho_w is not None:
                w = phospho_w[a:b, col][:, None]
                if np.any(w > 0):
                    low = atlas.L[:, pi, :][:, lo[a:b, col]].T
                    up = up + w * (low - up)
            acc += up
        out[a:b] = acc
    if published_rounding:
        out = np.round(out, 3)
    if atlas.kin_type == "ser_thr" and st_fav:
        c = np.asarray(center)
        assert np.isin(c, ["S", "T"]).all(), "ser_thr atlas applied to non-S/T centre"
        out += np.where((c == "S")[:, None], atlas.st_fav_log2[:, 0][None, :], atlas.st_fav_log2[:, 1][None, :])
        if published_rounding:
            out = np.round(out, 3)
    elif atlas.kin_type == "tyrosine":
        assert np.isin(np.asarray(center), ["Y"]).all(), "tyrosine atlas applied to non-Y centre"
    return out


def percentile(scores: np.ndarray, ref: np.ndarray, round_digits: int | None = 2) -> np.ndarray:
    """Column-wise percentile of scores (N x K) against reference (R x K):
    100 * #{ref <= s} / R (kinase-library: searchsorted side='right')."""
    assert scores.shape[1] == ref.shape[1]
    out = np.empty(scores.shape)
    for j in range(ref.shape[1]):
        r = np.sort(ref[:, j][np.isfinite(ref[:, j])])
        out[:, j] = np.searchsorted(r, scores[:, j], side="right") / len(r) * 100
    out[~np.isfinite(scores)] = np.nan
    return np.round(out, round_digits) if round_digits is not None else out


# ----------------------------------------------------------------------------------------------
# priming weights
# ----------------------------------------------------------------------------------------------
def flank_weights(sites: pd.DataFrame, site_weight: pd.Series) -> np.ndarray:
    """N x 15 matrix: weight of window position (offset -7..+7) being a phosphosite, looked up as
    site_weight[f'{acc}_{res}{pos+off}'] for flank residues S/T/Y; 0 elsewhere and at the centre."""
    wmap = site_weight.to_dict()
    W = np.zeros((len(sites), 15), dtype=np.float32)
    acc = sites["acc"].to_numpy()
    pos = sites["pos"].to_numpy()
    win = sites["window15"].to_numpy()
    for j, off in enumerate(range(-7, 8)):
        if off == 0:
            continue
        ch = np.array([w[j] for w in win])
        m = np.isin(ch, ["S", "T", "Y"])
        idx = np.where(m)[0]
        keys = [f"{acc[i]}_{ch[i]}{pos[i] + off}" for i in idx]
        W[idx, j] = np.array([wmap.get(k, 0.0) for k in keys], dtype=np.float32)
    return W


def ms_evidence_weight(universe: pd.DataFrame) -> pd.Series:
    """Soft priming weight per U_obs site = fraction of the 4 independent MS evidence channels that
    observed it: PSP MS (MS_LIT+MS_CST>=1), compendium, CPTAC, Ochoa 2020. Label-free (no KS info,
    no low-throughput literature channel). Returns Series indexed by site_id, values in {0,.25,.5,.75,1}."""
    u = universe
    psp_ms = (u.psp_ms_lit.fillna(0) + u.psp_ms_cst.fillna(0)) > 0
    n = psp_ms.astype(int) + u.src_compendium.astype(int) + u.src_cptac.astype(int) + u.src_ochoa.astype(int)
    return pd.Series((n / 4.0).to_numpy(), index=u.site_id.to_numpy())


# ----------------------------------------------------------------------------------------------
# nearest atlas kinase for kinases without a matrix (imputation; NOT part of the published arm)
# ----------------------------------------------------------------------------------------------
def domain_sequences(kinases: pd.DataFrame, seqs: dict) -> dict:
    """acc -> kinase-domain sequence (UniProt domain coordinates in the index; full length if absent)."""
    out = {}
    for r in kinases.itertuples():
        s = seqs[r.uniprot_acc]
        if pd.notna(r.kinase_domain_start) and pd.notna(r.kinase_domain_end):
            s = s[int(r.kinase_domain_start) - 1:int(r.kinase_domain_end)]
        out[r.uniprot_acc] = s
    return out


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10.0, -0.5
    al.end_gap_score = 0.0           # free end gaps (domain boundaries are imprecise)
    return al


def pair_identity(a: str, b: str, al=None) -> float:
    """Identical aligned residues / length of the shorter sequence (global, free end gaps, BLOSUM62)."""
    al = al or _aligner()
    aln = al.align(a, b)[0]
    ident = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        ident += sum(x == y for x, y in zip(a[s1:e1], b[s2:e2]))
    return ident / min(len(a), len(b))


def _nearest_one(args):
    q, qseq, cands, top = args
    al = _aligner()
    sc = [(al.score(qseq, cs), c) for c, cs in cands if c != q]
    sc.sort(reverse=True)
    best = [(pair_identity(qseq, dict(cands)[c], al), c) for _, c in sc[:top]]
    best.sort(reverse=True)
    return q, best[0][1], best[0][0], best[1][1] if len(best) > 1 else None, best[1][0] if len(best) > 1 else np.nan


def nearest_atlas(query_accs, atlas_accs, dom: dict, n_jobs: int = 12, top: int = 8) -> pd.DataFrame:
    """For each query kinase, the atlas kinase with highest kinase-domain identity. Candidates are
    pre-ranked by BLOSUM62 global alignment score; identity is computed for the `top` best."""
    from concurrent.futures import ProcessPoolExecutor
    cands = [(c, dom[c]) for c in atlas_accs]
    jobs = [(q, dom[q], cands, top) for q in query_accs]
    with ProcessPoolExecutor(n_jobs) as ex:
        res = list(ex.map(_nearest_one, jobs, chunksize=2))
    return pd.DataFrame(res, columns=["uniprot_acc", "nn_acc", "nn_identity", "nn2_acc", "nn2_identity"])
