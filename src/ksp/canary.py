"""Canary pseudo-kinases (DESIGN §11.3, AIFM-3; PI r2 update). Owner: core-infra.

* 12 NEGATIVE canaries CANARY_NEG_00..11: labels are random U_ms sites, CLUSTERED BY PROTEIN like real kinases
  (random proteins, then 1-4 random U_ms sites of the canary's residue class on each), sizes drawn from the label
  counts of real kinases with 10-100 L2 labels; residue classes mixed (6 ST, 4 Y, 2 STY); descriptors are random but
  fixed (canary_descriptors). Any model that scores them above chance on held-out data is leaking.
* 1 POSITIVE canary CANARY_POS: planted motif R-x-x-[ST]-[FILMV] (R at -3, hydrophobic at +1) on a fixed random set
  of "PPI partner" proteins (canary_ppi_partners(): the planted PPI flag; feature stores should expose it as an
  ordinary PPI edge so the full pipeline sees it). Labels = 1-2 motif U_ms sites per partner protein.

Canary labels are generated over ALL proteins (deterministic, seed CANARY_SEED) BEFORE the split is frozen, so the
split optimiser balances them like real kinases; scripts/build_splits.py seals their TEST part with the TEST labels.
`canary_labels()` returns the DEV part only (used by ksp.labels.load_dev(with_canaries=True)).
Acceptance (checked in ksp.eval.canary_check): mean negative-canary AUROC <= 0.57 and positive >= 0.90, on both
E-hard-type and E-all-type negatives. Canaries are excluded from every reported metric.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

CANARY_SEED = 20261026
N_NEG = 12
NEG_CLASSES = ["ST"] * 6 + ["Y"] * 4 + ["STY"] * 2
POS_ID = "CANARY_POS"
NEG_IDS = [f"CANARY_NEG_{i:02d}" for i in range(N_NEG)]
CANARY_IDS = NEG_IDS + [POS_ID]
POS_N_PARTNERS = 120          # planted PPI partner proteins
DESC_DIM = 64


def is_canary(kinase_acc) -> np.ndarray:
    return pd.Series(list(kinase_acc) if not isinstance(kinase_acc, str) else [kinase_acc],
                     dtype=object).astype(str).str.startswith("CANARY_").values


def canary_kinases() -> pd.DataFrame:
    """Rows in the format of cache/index/kinases.parquet; appended AFTER the 503 real kinases (index 503..515)."""
    cls = NEG_CLASSES + ["ST"]
    return pd.DataFrame({"uniprot_acc": CANARY_IDS, "gene": CANARY_IDS, "residue_class": cls,
                         "manning_group": "CANARY", "manning_family": [f"CANARY{i}" for i in range(len(cls))],
                         "manning_subfamily": None, "is_pseudokinase": False,
                         "kinase_domain_start": np.nan, "kinase_domain_end": np.nan})


def kinases_with_canaries() -> pd.DataFrame:
    from . import index
    return pd.concat([index.load_kinases(), canary_kinases()], ignore_index=True)


def canary_descriptors(dim: int = DESC_DIM) -> pd.DataFrame:
    """Random-but-fixed kinase descriptor vectors (standard normal) for all canaries, rows = CANARY_IDS."""
    rng = np.random.default_rng(CANARY_SEED + 1)
    return pd.DataFrame(rng.standard_normal((len(CANARY_IDS), dim)).astype(np.float32), index=CANARY_IDS,
                        columns=[f"d{i}" for i in range(dim)])


def _motif_mask(w15: pd.Series) -> np.ndarray:
    w = w15.values.astype(str)
    return np.array([x[4] == "R" and x[7] in "ST" and x[8] in "FILMV" for x in w])


@lru_cache(maxsize=1)
def _generate(size_source: tuple) -> tuple[pd.DataFrame, tuple]:
    """size_source: tuple of real-kinase L2 label counts (10-100) to draw negative canary sizes from."""
    from . import index
    rng = np.random.default_rng(CANARY_SEED)
    s = index.load_sites()
    s = s[s.in_U_ms].reset_index(drop=True)
    sizes = rng.choice(np.array(size_source), size=N_NEG, replace=True)
    rows = []
    for cid, cls, n in zip(NEG_IDS, NEG_CLASSES, sizes):
        pool = s[s.res.isin(list("ST") if cls == "ST" else list(cls))]
        by_p = pool.groupby("acc").indices
        prots = np.array(sorted(by_p))
        rng.shuffle(prots)
        chosen = []
        for p in prots:
            if len(chosen) >= n:
                break
            idx = by_p[p]
            k = min(int(rng.integers(1, 5)), len(idx), n - len(chosen))
            chosen.extend(rng.choice(idx, size=k, replace=False).tolist())
        for j in chosen:
            rows.append((cid, cls, pool.site_id.values[j]))
    # positive canary
    st = s[s.res.isin(["S", "T"])]
    mot = st[_motif_mask(st.window15)]
    by_p = mot.groupby("acc").indices
    prots = np.array(sorted(by_p))
    partners = np.sort(rng.choice(prots, size=min(POS_N_PARTNERS, len(prots)), replace=False))
    for p in partners:
        idx = by_p[p]
        k = min(int(rng.integers(1, 3)), len(idx))
        for j in rng.choice(idx, size=k, replace=False):
            rows.append((POS_ID, "ST", mot.site_id.values[j]))
    d = pd.DataFrame(rows, columns=["kinase", "kinase_residue_class", "site_id"])
    return d, tuple(partners)


def _size_source() -> tuple:
    """Label counts of real kinases with 10-100 L2 labels (from the pre-split table during build, else DEV)."""
    from . import labels
    if labels.FULL_PAIRS.exists():
        p = pd.read_parquet(labels.FULL_PAIRS, columns=["kinase_acc", "site_id", "L2"])
        n = p[p.L2].groupby("kinase_acc").site_id.nunique()
    else:
        raise FileNotFoundError("canary sizes are frozen in cache/labels/canary_all.parquet; it is missing")
    return tuple(sorted(int(x) for x in n[(n >= 10) & (n <= 100)]))


def _as_pairs(d: pd.DataFrame) -> pd.DataFrame:
    out = d.copy()
    out["kinase_acc"] = out.kinase
    out["acc"] = out.site_id.str.rsplit("_", n=1).str[0]
    out["res"] = out.site_id.str.rsplit("_", n=1).str[1].str[0]
    out["pos"] = out.site_id.str.rsplit("_", n=1).str[1].str[1:].astype(int)
    out["provenance"] = "PSP"          # behaves like a human-native PSP pair everywhere
    for c in ["L1", "L2", "L3", "L2_noauto", "native_human", "hmr", "sub_native_human", "in_vivo", "in_vitro",
              "in_vivo_native", "in_vitro_native", "in_vivo_hmr", "in_vitro_hmr", "kinase_included",
              "residue_in_class", "in_universe"]:
        out[c] = True
    for c in ["L2_vitro", "EXT", "EXT_not_in_psp", "autophos", "in_psp_ks"]:
        out[c] = False
    out["L2_vitro"] = True
    out["evidence_native"] = "both"
    out["evidence_hmr"] = "both"
    out["species_origin"] = "canary"
    out["is_canary"] = True
    return out


def generate_all(size_source: tuple | None = None) -> tuple[pd.DataFrame, list]:
    """All canary label pairs (DEV + TEST proteins) and the planted PPI partner list. Used by build_splits only."""
    d, partners = _generate(size_source if size_source is not None else _size_source())
    return _as_pairs(d), list(partners)


CANARY_FILE_DEV = "canary_dev.parquet"


def canary_labels() -> pd.DataFrame:
    """DEV-protein canary label pairs (frozen at split time in cache/labels/canary_dev.parquet)."""
    from . import labels
    f = labels.DEV_PAIRS.parent / CANARY_FILE_DEV
    if not f.exists():
        raise FileNotFoundError(f"{f} missing: run scripts/build_splits.py")
    return pd.read_parquet(f)


def canary_ppi_partners() -> list:
    """Planted PPI partner proteins of CANARY_POS (all splits; frozen in results/splits/canary_ppi_partners.tsv)."""
    from . import splits
    return pd.read_csv(splits.SPLIT_DIR / "canary_ppi_partners.tsv", sep="\t").acc.tolist()
