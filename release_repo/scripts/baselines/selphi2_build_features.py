"""Build the label-free SELPHI2.0 feature inputs (B4): Ochoa site table, expression matrices, normalised Cantley.

Outputs (cache/selphi2/): ochoa_moesm4_5.parquet, site_features.parquet (Ochoa-covered sites),
site_features_all.parquet (every index site, SELPHI NA defaults), expr.npz + expr_info.json, cantley_norm.npz;
results/baselines/selphi2/feature_build.json (counts).
Usage: python -I scripts/baselines/selphi2_build_features.py
"""
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

from ksp.baselines.selphi2 import features as F  # noqa: E402
from ksp.index import load_sites  # noqa: E402


def main():
    F.CACHE.mkdir(parents=True, exist_ok=True)
    out = ROOT / "results" / "baselines" / "selphi2"
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rep = {}
    st = F.build_site_table(all_sites=False)
    st.to_parquet(F.CACHE / "site_features.parquet")
    S = load_sites().set_index("site_id")
    ums = S.in_U_ms.reindex(st.index).fillna(False).astype(bool)
    rep["site_table"] = {"ochoa_rows": st.attrs["n_ochoa_rows"], "sites_mapped": len(st),
                         "sites_mapped_in_U_ms": int(ums.sum()),
                         "by_residue_in_U_ms": S.res.reindex(st.index[ums.values]).value_counts().to_dict(),
                         "U_ms_sites_total": int(S.in_U_ms.sum()), "seconds": round(time.time() - t0, 1)}
    print(rep["site_table"], flush=True)
    sa = F.build_site_table(all_sites=True)
    sa.to_parquet(F.CACHE / "site_features_all.parquet")
    t1 = time.time()
    rep["expression"] = F.build_expression()
    rep["expression"]["seconds"] = round(time.time() - t1, 1)
    t2 = time.time()
    F.load_site_table.cache_clear()
    stn, yn = F.load_cantley_norm()
    rep["cantley_norm"] = {"seconds": round(time.time() - t2, 1)}
    rep["total_seconds"] = round(time.time() - t0, 1)
    (out / "feature_build.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()
