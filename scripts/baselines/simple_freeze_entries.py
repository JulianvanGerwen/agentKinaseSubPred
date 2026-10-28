"""Collect the simple-baseline FREEZE entries (results/baselines/simple/freeze_entries/<method>.json) into
freeze_entries/ALL.json, re-verifying every SHA-256 pin with ksp.eval_test.pin_files (raises on any mismatch)."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

from ksp import eval_test as T  # noqa: E402

FE = ROOT / "results" / "baselines" / "simple" / "freeze_entries"
ORDER = ["b1_arm_i", "b1_arm_ii", "b2_pssm", "b7_lr", "b0_prior", "b0_pop", "b0_miss"]


def main():
    cfgs, pre = [], {}
    for m in ORDER:
        f = FE / f"{m}.json"
        if not f.exists():
            print("MISSING", m)
            continue
        d = json.loads(f.read_text())
        e = d["freeze_entry"]
        chk = T.pin_files({"configs": [copy.deepcopy(e)]})["configs"][0]
        for k in T.FILE_KEYS:
            if e.get(k) and chk[f"{k}_sha256"] != e[f"{k}_sha256"]:
                raise RuntimeError(f"{m}: {k} changed since pinning")
        cfgs.append(e)
        pre[e["name"]] = {"main": d["preflight"], "t3": d.get("t3_preflight")}
    comps = [{"a": "final", "b": "B1_i", "endpoints": ["M1", "M3"], "primary": True},
             {"a": "final", "b": "B1_ii", "endpoints": ["M2"], "primary": True}]
    comps += [{"a": "final", "b": c["name"], "endpoints": ["M1", "M2", "M3", "M4"], "primary": False}
              for c in cfgs if c["name"] in ("B2", "B7", "B0_prior", "B0_pop", "B0_miss")]
    out = {"configs": cfgs, "suggested_comparisons": comps, "preflight": pre}
    (FE / "ALL.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
