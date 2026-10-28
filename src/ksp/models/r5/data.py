"""R5 inputs (DESIGN §8.4, §5, §16 A6; R5 = fine-tuned ESM-2 late-interaction kinase-site model).

Reuses R3's data code wherever it already does the right thing:
  * site axis, CV folds, frozen ESM-2 650M window embeddings, protein means : ksp.models.r3.data.load_sites
  * kinase descriptor blocks (ESM-2 kinase-domain mean, Cantley PCA, Manning group/family), canaries, class_ok
                                                                              : ksp.models.r3.data.load_kinases
  * training-target matrix (positives / candidate negatives / masks)          : ksp.models.r3.data.build_targets

R5 adds two things R3 does not have:
  1. SUBSTRATE TOKENS. The +-15 window (`window31` of cache/index/sites.parquet) as ESM-2 token ids, plus a
     phospho-channel flag per position (1 = that flanking residue is itself a U_ms phosphosite). The centre site's
     own status is NEVER used (forced 0, A6/DL2-7). Terminus padding ('_' in window31) is masked out.
  2. KINASE-DOMAIN TOKENS. Cached per-residue ESM-2 650M embeddings of each kinase catalytic domain
     (cache/r5/kin_tok.npy, built by scripts/models/r5_kinase_embed.py). Canaries get random-but-fixed pseudo-token
     blocks (same role as ksp.canary.canary_descriptors for R3): no real sequence information.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ksp import canary, index
from ksp.models.r3 import data as r3data

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
R5CACHE = ROOT / "cache" / "r5"
FLANK = 15                      # +-15 window -> 31 residues
N_CANARY_TOK = 16               # pseudo-domain tokens per canary
CANARY_TOK_SEED = 20261026 + 5

# ESM-2 vocabulary (identical for every esm2_t*_UR50D checkpoint); ids taken from the tokenizer at build time and
# asserted against it in build_window_tokens().
PAD_ID, CLS_ID, EOS_ID, UNK_ID = 1, 0, 2, 3


@dataclass
class R5Sites:
    si: r3data.SiteInputs       # R3 site inputs (site_idx, res_code, acc, cv_fold, rows, emb, prot_mean, ...)
    tok: np.ndarray             # (N, 33) int16 ESM token ids: CLS + 31 residues + EOS
    attn: np.ndarray            # (N, 33) bool
    phos: np.ndarray            # (N, 33) int8, 1 = flanking residue is a U_ms site (centre forced 0)
    info: dict


def _esm_vocab():
    from transformers import AutoTokenizer
    tk = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D")
    v = tk.get_vocab()
    if (v["<pad>"], v["<cls>"], v["<eos>"], v["<unk>"]) != (PAD_ID, CLS_ID, EOS_ID, UNK_ID):
        raise RuntimeError(f"unexpected ESM-2 special token ids: {v['<pad>'], v['<cls>'], v['<eos>'], v['<unk>']}")
    return v


def build_window_tokens(site_ids: np.ndarray) -> dict:
    """ESM token ids + phospho flags for the +-15 window of `site_ids` (cached, keyed by the site-id array)."""
    f = R5CACHE / "win31_tokens.npz"
    if f.exists():
        z = np.load(f, allow_pickle=False)
        want = np.asarray(site_ids, dtype="U")
        if len(z["site_ids"]) == len(want) and (z["site_ids"] == want).all():
            return {"tok": z["tok"], "attn": z["attn"], "phos": z["phos"],
                    "meta": json.loads((R5CACHE / "win31_tokens.json").read_text())}
    v = _esm_vocab()
    S = index.load_sites()
    pos_of = pd.Series(np.arange(len(S)), index=S.site_id.values)
    row = pos_of.reindex(site_ids)
    if row.isna().any():
        raise KeyError("R5: site_id missing from cache/index/sites.parquet")
    Sx = S.iloc[row.values.astype(np.int64)]
    w31 = Sx.window31.values.astype(str)
    L = 2 * FLANK + 1
    # token ids
    lut = np.full(256, UNK_ID, np.int16)
    for aa, i in v.items():
        if len(aa) == 1 and aa.isalpha():
            lut[ord(aa)] = i
    lut[ord("_")] = PAD_ID
    W = np.frombuffer("".join(w31).encode("ascii"), dtype=np.uint8).reshape(len(w31), L)
    body = lut[W]
    attn_body = W != ord("_")
    tok = np.c_[np.full(len(w31), CLS_ID, np.int16), body, np.full(len(w31), EOS_ID, np.int16)]
    attn = np.c_[np.ones(len(w31), bool), attn_body, np.ones(len(w31), bool)]
    # phospho channel: flanking residue is itself a U_ms site
    ums = S[S.in_U_ms]
    pacc = index.load_proteins().acc.values
    pidx = pd.Series(np.arange(len(pacc)), index=pacc)
    key_ums = np.sort(pidx.reindex(ums.acc.values).values.astype(np.int64) * 100_000 + ums.pos.values.astype(np.int64))
    centre = Sx.pos.values.astype(np.int64)
    base = pidx.reindex(Sx.acc.values).values.astype(np.int64) * 100_000
    offs = np.arange(-FLANK, FLANK + 1)
    keys = base[:, None] + (centre[:, None] + offs[None, :])
    j = np.searchsorted(key_ums, keys)
    hit = (j < len(key_ums)) & (key_ums[np.minimum(j, len(key_ums) - 1)] == keys)
    hit &= attn_body
    hit[:, FLANK] = False                                     # never the centre site's own status (A6 / DL2-7)
    phos = np.c_[np.zeros(len(w31), np.int8), hit.astype(np.int8), np.zeros(len(w31), np.int8)]
    meta = {"n_sites": int(len(w31)), "flank": FLANK, "L_tokens": int(tok.shape[1]),
            "n_padded_positions": int((~attn_body).sum()), "n_phospho_flags": int(hit.sum()),
            "frac_sites_with_flank_phos": float((hit.any(1)).mean())}
    R5CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(f, site_ids=np.asarray(site_ids, dtype="U"), tok=tok, attn=attn, phos=phos)
    (R5CACHE / "win31_tokens.json").write_text(json.dumps(meta, indent=1))
    return {"tok": tok, "attn": attn, "phos": phos, "meta": meta}


def load_sites(esm: str = "650m") -> R5Sites:
    si = r3data.load_sites(esm)
    w = build_window_tokens(si.site_ids)
    return R5Sites(si, w["tok"], w["attn"], w["phos"], {"window": w["meta"], **si.info})


@dataclass
class R5Kinases:
    ki: r3data.KinaseInputs     # R3 descriptor blocks (+ canaries, class_ok, group/family)
    tok: np.ndarray             # (K, L_KIN, 1280) float16 per-residue kinase-domain embeddings
    mask: np.ndarray            # (K, L_KIN) bool
    info: dict


def load_kinases() -> R5Kinases:
    ki = r3data.load_kinases()
    E = np.load(R5CACHE / "kin_tok.npy")
    M = np.load(R5CACHE / "kin_mask.npy")
    meta = json.loads((R5CACHE / "kin_meta.json").read_text())
    nR = int((~ki.is_canary).sum())
    if E.shape[0] != nR:
        raise RuntimeError(f"cache/r5/kin_tok.npy has {E.shape[0]} rows, expected {nR} real kinases; "
                           "rebuild with scripts/models/r5_kinase_embed.py")
    acc = [m["acc"] for m in meta["kinases"]]
    if acc != list(ki.kinases.uniprot_acc.values[:nR]):
        raise RuntimeError("cache/r5/kin_tok.npy kinase order does not match cache/index/kinases.parquet")
    nC = len(ki.kinases) - nR
    rng = np.random.default_rng(CANARY_TOK_SEED)
    cE = np.zeros((nC, E.shape[1], E.shape[2]), np.float16)
    cM = np.zeros((nC, E.shape[1]), bool)
    scale = float(np.asarray(E[M], np.float32).std())
    cE[:, :N_CANARY_TOK] = (rng.standard_normal((nC, N_CANARY_TOK, E.shape[2])) * scale).astype(np.float16)
    cM[:, :N_CANARY_TOK] = True
    return R5Kinases(ki, np.concatenate([E, cE], 0), np.concatenate([M, cM], 0),
                     {"L_kin": int(E.shape[1]), "d_kin_tok": int(E.shape[2]), "n_real": nR, "n_canary": nC,
                      "n_truncated": meta["n_truncated"], "n_no_domain": meta["n_no_domain"],
                      "canary_tokens": N_CANARY_TOK, "canary_token_scale": round(scale, 4),
                      "kin_embed_model": meta["model"]})


build_targets = r3data.build_targets
