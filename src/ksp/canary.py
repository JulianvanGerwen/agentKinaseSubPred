"""Canary pseudo-kinases (DESIGN §11.3, AIFM-3; PI r2 + canary v2 decision). Owner: core-infra.

CANARY v2 (current; PI decision after the R3 audit of v1):
* 12 NEGATIVE canaries CANARY_NEG_00..11: labels are REAL DEV-labelled sites (evaluation-positive sites of real
  kinases: human-native PSP, non-autophos, in U_ms, pooled across kinases), protein-clustered and site-uniform:
  proteins are drawn without replacement with weight n_p / E[min(K, n_p)], then min(K, n_p) of their labelled sites,
  K ~ U{1..4}, until the canary's size is reached (every labelled site has ~equal inclusion probability; ~70% of
  canary sites have a same-protein sibling). Sizes are the v1 sizes (drawn from real kinases with 10-100 L2 labels);
  residue classes 6 ST / 4 Y / 2 STY with the natural residue mix of the labelled pool. A model that learned generic 'labelled-site-ness' cannot separate them from E-hard negatives
  (other labelled sites), so the acceptance gate is mean negative-canary AUROC on E-HARD <= 0.57 (E-all reported,
  not gated). KP masks of canaries = their own sites only. Because the pool is DEV-only, v2 negative canaries have no
  TEST labels (TEST share 0; canaries are evaluated on CV only).
* 1 POSITIVE canary CANARY_POS (unchanged from v1): planted motif R-x-x-[ST]-[FILMV] on a fixed random set of
  "PPI partner" proteins (canary_ppi_partners(): the planted PPI flag). Gate: AUROC >= 0.90 on E-hard and E-all.

CANARY v1 (provenance; used for the split freeze at 14:34): negative labels were random U_ms sites on uniformly
chosen proteins (1-4 per protein) - biased toward proteins with few U_ms sites (a label-free -#U_ms-sites score gave
AUROC 0.76). Kept in cache/labels/canary_dev_v1.parquet. Canary sites inherit DEV/TEST from their protein; the split
never changes. `canary_labels()` returns the DEV part (used by ksp.labels.load_dev(with_canaries=True)).
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
CANARY_VERSION = 2


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


def _generate_v2_negatives(sizes: dict, dev_pool: pd.DataFrame) -> pd.DataFrame:
    """v2 negative canaries from real DEV-labelled sites (see module docstring). dev_pool: unique sites with columns
    site_id, acc, res (evaluation-positive sites of real kinases on DEV proteins, in U_ms, non-autophos).
    Site-uniform protein clustering: proteins are drawn without replacement with weight n_p / E[min(K, n_p)]
    (n_p = labelled pool sites on the protein, K ~ U{1..4}), then min(K, n_p) of their sites uniformly; every
    labelled site thus has (approximately) equal inclusion probability, so protein-level label counts or protein
    size do not separate canary sites from other labelled sites."""
    rng = np.random.default_rng(CANARY_SEED + 2)
    pool_all = dev_pool.drop_duplicates("site_id").sort_values("site_id").reset_index(drop=True)
    rows = []
    for cid, cls in zip(NEG_IDS, NEG_CLASSES):
        n = int(sizes[cid])
        pool = pool_all[pool_all.res.isin(list("ST") if cls == "ST" else list(cls))].reset_index(drop=True)
        by_p = pool.groupby("acc").indices
        prots = np.array(sorted(by_p))
        npp = np.array([len(by_p[p]) for p in prots], dtype=float)
        e_take = np.array([np.mean([min(k, m) for k in (1, 2, 3, 4)]) for m in npp])
        w = npp / e_take
        avail = np.ones(len(prots), bool)
        chosen = []
        while len(chosen) < n and avail.any():
            ww = np.where(avail, w, 0.0)
            i = rng.choice(len(prots), p=ww / ww.sum())
            avail[i] = False
            idx = by_p[prots[i]]
            k = min(int(rng.integers(1, 5)), len(idx), n - len(chosen))
            chosen += rng.choice(idx, size=k, replace=False).tolist()
        rows += [(cid, cls, pool.site_id.values[x]) for x in chosen]
    return pd.DataFrame(rows, columns=["kinase", "kinase_residue_class", "site_id"])


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
    out["in_index"] = True
    out["in_U_ms"] = True
    out["new_since_2021"] = False
    out["ext_test"] = False
    out["ext_tier"] = ""
    out["mask_only"] = False
    out["mask_sensitivity_only"] = False
    return out


def generate_v1(size_source: tuple | None = None) -> tuple[pd.DataFrame, list]:
    """v1 canaries (provenance; the split freeze used these)."""
    d, partners = _generate(size_source if size_source is not None else _size_source())
    return _as_pairs(d), list(partners)


def generate_all(size_source: tuple | None = None, dev_pool: pd.DataFrame | None = None
                 ) -> tuple[pd.DataFrame, list]:
    """Current (v2) canary label pairs and the planted PPI partner list: v1 positive canary (identical) + v2 negative
    canaries drawn from `dev_pool` (required; built by scripts/build_splits.py --relabel from DEV labels)."""
    if dev_pool is None:
        raise ValueError("canary v2 needs dev_pool (DEV evaluation-positive sites); see scripts/build_splits.py")
    d1, partners = _generate(size_source if size_source is not None else _size_source())
    sizes = d1[d1.kinase.ne(POS_ID)].groupby("kinase").size().to_dict()
    neg = _generate_v2_negatives(sizes, dev_pool)
    d = pd.concat([neg, d1[d1.kinase.eq(POS_ID)]], ignore_index=True)
    out = _as_pairs(d)
    out["canary_version"] = CANARY_VERSION
    return out, list(partners)


CANARY_FILE_DEV = "canary_dev.parquet"
CANARY_FILE_DEV_V1 = "canary_dev_v1.parquet"


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
