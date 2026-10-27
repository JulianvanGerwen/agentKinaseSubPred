"""Minimal example: read the released predictions (no project data, GPU or ksp install needed).

Usage (from the repository root):
    python examples/read_predictions.py                    # summary + a few look-ups
    python examples/read_predictions.py --site P04637_S15  # top-5 kinases for one site
    python examples/read_predictions.py --kinase ATM       # top-20 sites for one kinase

Requires only pandas. The files are derived from PhosphoSitePlus labels (CC BY-NC-SA 3.0): non-commercial use,
share-alike, attribution required (predictions/LICENSE-PREDICTIONS.md).
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

PRED = Path(__file__).resolve().parents[1] / "predictions"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_checksums() -> None:
    manifest = json.loads((PRED / "RELEASE_MANIFEST.json").read_text())
    for key in ("kinase_top200", "site_top5"):
        name = Path(manifest[key]["path"]).name
        ok = sha256(PRED / name) == manifest[key]["sha256"]
        print(f"{name}: sha256 {'OK' if ok else 'MISMATCH'}")
        if not ok:
            raise SystemExit(f"{name} does not match RELEASE_MANIFEST.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="P04637_S15", help="site_id = <UniProt accession>_<residue><position>")
    ap.add_argument("--kinase", default="ATM", help="kinase gene symbol or UniProt accession")
    ap.add_argument("--no-verify", action="store_true", help="skip the SHA-256 check")
    a = ap.parse_args()
    if not a.no_verify:
        verify_checksums()

    site5 = pd.read_csv(PRED / "site_top5.tsv.gz", sep="\t")
    kin200 = pd.read_csv(PRED / "kinase_top200.tsv.gz", sep="\t")
    print(f"site_top5: {len(site5):,} rows, {site5.site_id.nunique():,} sites, {site5.kinase_acc.nunique()} kinases")
    print(f"kinase_top200: {len(kin200):,} rows, {kin200.kinase_acc.nunique()} kinases")

    top1 = site5[site5["rank"] == 1]
    print(f"sites whose rank-1 assignment_prob >= 0.5: {(top1.assignment_prob >= 0.5).sum():,}; "
          f">= 0.8: {(top1.assignment_prob >= 0.8).sum():,}")

    print(f"\nTop-5 kinases for {a.site} (ranked by assignment_prob; treat as a shortlist, not a call):")
    rows = site5[site5.site_id == a.site]
    print(rows.to_string(index=False) if len(rows) else "  site not in the released universe")

    sel = kin200[(kin200.gene == a.kinase) | (kin200.kinase_acc == a.kinase)]
    print(f"\nTop-20 sites for {a.kinase} (ranked by fused score):")
    print(sel.head(20).to_string(index=False) if len(sel) else "  kinase not in the released set")


if __name__ == "__main__":
    main()
