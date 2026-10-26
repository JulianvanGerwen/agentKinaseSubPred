"""Training-label variants for the regime comparison (DESIGN §3 ablations). Wrapper only: ksp.labels is frozen (A28).

Variants (training labels only; evaluation positives are ALWAYS ksp.labels.evaluable_positive via ksp.eval):
  L2        default (human-native + mouse/rat PSP transferred to human)            -> ksp.labels 'L2'
  L1        human-native PSP rows only                                             -> ksp.labels 'L1'
  L2_vitro  L2 pairs with in vitro evidence (HMR rows; includes 'both')            -> ksp.labels 'L2_vitro'
  L2_vivo   L2 pairs with in vivo evidence ONLY (= L2 minus L2_vitro; evidence_hmr == 'in_vivo_only').
            ASSUMPTION (stated): "in-vivo-only" as in the brief, i.e. the complement of L2_vitro inside L2, so that
            L2 = L2_vitro (disjoint union) L2_vivo.
  L2ext     L2 + the filtered external curated pairs on DEV proteins: provenance OmniPath/UniProt rows with
            ext_test == True (DL-5 near-copy filter + A9 PhosphoSite-source drop + 20/PMID cap, tiers A and B;
            computed by ksp.labels.external_test_filter/ext_test_tiers at split freeze), U_ms only.
Canary labels: identical in every variant (all 428 canary rows are L1/L2/L2_vitro-flagged), so the canary gate is
comparable across variants.

`patch()` replaces the attribute ksp.labels.load_dev in THIS process by `load_dev_variant` (pass-through for every
built-in variant incl. 'all'/'eval'), so ksp.models.r2.sample.FoldData and ksp.crossfit.outer_split see the new
variants without editing frozen / R2-owned code.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ksp import labels
from ksp.index import ROOT

EXTRA = ("L2_vivo", "L2ext")
ALL_VARIANTS = ("L2", "L1", "L2_vitro", "L2_vivo", "L2ext")
XFIT_ROOT = ROOT / "cache" / "regime" / "xfit"

_ORIG = getattr(labels.load_dev, "_regime_orig", labels.load_dev)


def load_dev_variant(variant: str = "L2", with_canaries: bool = False, path=None, ums_only=None) -> pd.DataFrame:
    if variant not in EXTRA:
        return _ORIG(variant, with_canaries=with_canaries, path=path, ums_only=ums_only)
    d = _ORIG("L2", with_canaries=with_canaries, path=path, ums_only=ums_only)
    canm = d.is_canary.astype(bool)
    if variant == "L2_vivo":
        keep = canm | ~d.L2_vitro.astype(bool)
        out = d[keep].reset_index(drop=True)
        real = out[~out.is_canary.astype(bool)]
        if not (real.evidence_hmr == "in_vivo_only").all():
            raise AssertionError("L2_vivo contains pairs that are not in-vivo-only")
        return out
    # L2ext
    a = _ORIG("all", with_canaries=False, path=path, ums_only=False)
    e = a[a.provenance.isin(["OmniPath", "UniProt"]) & a.ext_test.fillna(False).astype(bool)]
    if ums_only is None or ums_only:
        e = e[e.in_U_ms.astype(bool)]
    e = e.sort_values(["kinase_acc", "site_id", "provenance"]).drop_duplicates(["kinase_acc", "site_id"])
    have = set(zip(d.kinase_acc, d.site_id))
    e = e[[(k, s) not in have for k, s in zip(e.kinase_acc, e.site_id)]].copy()
    e["is_canary"] = False
    return pd.concat([d, e], ignore_index=True)


load_dev_variant._regime_orig = _ORIG


def patch() -> None:
    labels.load_dev = load_dev_variant


def xfit_tag(variant: str) -> str:
    """Absolute path used as ksp.models.r2.xfit 'tag' (XFIT_DIR / <abs path> == <abs path>): our own cache dir."""
    return str(XFIT_ROOT / variant)


def variant_summary() -> pd.DataFrame:
    rows = []
    for v in ALL_VARIANTS:
        d = load_dev_variant(v)
        rows.append({"variant": v, "pairs": len(d), "kinases": d.kinase_acc.nunique(), "sites": d.site_id.nunique(),
                     "kinases_ge20": int((d.groupby("kinase_acc").size() >= 20).sum())})
    return pd.DataFrame(rows)
