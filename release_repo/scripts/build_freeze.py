"""Assemble FREEZE.json from the owners' pinned FREEZE entries and validate it (PI, G4; DESIGN §4.3, A19, A23, A29).

Usage: python -I scripts/build_freeze.py [--out FREEZE.json] [--unavailable B4_retrained:"reason"]
Reads only entry JSONs listed in ENTRY_FILES (no labels). Re-verifies every pin and runs ksp.eval_test.validate_freeze.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval_test as ET  # noqa: E402

# Owner-delivered entry files (paths relative to ROOT). Each holds one config dict or {"configs": [...]}.
ENTRY_FILES = [
    "results/freeze_entries/final.json",                                 # R2+R3 fusion (R2 owner)
    "results/freeze_entries/r2_alone.json",                              # R2 component (secondary)
    "results/freeze_entries/r3.json",                                     # R3 component (secondary)
    "results/baselines/simple/freeze_entries/ALL.json",                  # B1 i/ii, B2, B7, B0 floors
    "results/freeze_entries/b3.json",                                    # B3 Phosformer retrained
    "results/freeze_entries/b4.json",                                    # B4 SELPHI2.0 retrained
]

PRIMARY = [
    {"a": "final", "b": "B1_i", "endpoints": ["M1", "M3"], "primary": True},
    {"a": "final", "b": "B1_ii", "endpoints": ["M2"], "primary": True},
    {"a": "final", "b": "B3", "endpoints": ["M1", "M2", "M3"], "primary": True},
    {"a": "final", "b": "B4", "endpoints": ["M1", "M2", "M3"], "primary": True},
]
SECONDARY = [
    {"a": "final", "b": "R2", "endpoints": ["M1", "M2", "M3", "M4"], "primary": False},
    {"a": "final", "b": "R3", "endpoints": ["M1", "M2", "M3", "M4"], "primary": False},
    {"a": "final", "b": "B7", "endpoints": ["M1", "M2", "M3", "M4"], "primary": False},
    {"a": "final", "b": "B2", "endpoints": ["M1", "M2", "M3", "M4"], "primary": False},
    {"a": "final", "b": "B1_i", "endpoints": ["M4"], "primary": False},
]


def load_entries() -> list[dict]:
    out = []
    for rel in ENTRY_FILES:
        p = ROOT / rel
        if not p.exists():
            print(f"[missing] {rel}")
            continue
        d = json.loads(p.read_text())
        out.extend(d["configs"] if "configs" in d else [d])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "FREEZE.json"))
    ap.add_argument("--unavailable", action="append", default=[],
                    help="baseline_id:reason for a primary baseline not delivered (A29 contingency)")
    a = ap.parse_args()
    configs = load_entries()
    names = {c["name"] for c in configs}
    unavailable = [dict(zip(("baseline_id", "reason"), u.split(":", 1))) for u in a.unavailable]
    gone = {u["baseline_id"] for u in unavailable}
    comps = [c for c in PRIMARY + SECONDARY
             if c["b"] in names and not any(x.get("baseline_id") in gone for x in configs if x["name"] == c["b"])]
    cfg = {"version": 1, "harness_commit": "974d520", "B": 10000, "B_secondary": 1000,
           "notes": "Single TEST batch (G4). Final = R2+R3 late fusion (DESIGN A13, docs/model_r2.md §6). "
                    "Primary family per §6.8/A19/A23; contingency per A29; B1 per A30; assignment_prob per A31(b).",
           "configs": configs, "comparisons": comps}
    if unavailable:
        cfg["unavailable_primary"] = unavailable
    ET.validate_freeze(cfg)  # raises on any schema / pin / primary-family violation
    Path(a.out).write_text(json.dumps(cfg, indent=1) + "\n")
    print(f"wrote {a.out}: {len(configs)} configs, {len(comps)} comparisons, unavailable={sorted(gone)}")


if __name__ == "__main__":
    main()
