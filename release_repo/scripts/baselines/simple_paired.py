"""Paired cluster bootstrap (ksp.eval.paired_bootstrap, B=1000) between simple baselines on CV.
    .venv/bin/python -I scripts/baselines/simple_paired.py b7_lr b1_published [b2_pssm b1_published ...]
Two views: 'full' (all candidates, missing='worst') and 'matched' (PRIOR-1: pairs both methods score, missing='drop').
Writes results/baselines/simple/paired.json (appends/updates entries)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import pandas as pd  # noqa: E402

from ksp import eval as kev  # noqa: E402

OUT = ROOT / "results" / "baselines" / "simple"


def preds(m):
    rid = json.loads((OUT / f"{m}.json").read_text())["run_id"]
    src = json.loads((ROOT / "results" / "runs" / rid / "manifest.json").read_text()).get("reevaluated_from") or rid
    return pd.read_parquet(ROOT / "cache" / "runs" / src / "predictions.parquet"), rid


def main():
    a = sys.argv[1:]
    pairs = list(zip(a[::2], a[1::2]))
    data = kev.prepare("cv")
    todo = [k for k in sorted(data.positives) if not k.startswith("CANARY_")]
    S = data.sites.site_id.values
    K = data.kinases.uniprot_acc.values
    f = OUT / "paired.json"
    allres = json.loads(f.read_text()) if f.exists() else {}
    for x, y in pairs:
        px, rx = preds(x)
        py, ry = preds(y)
        rA, rB = kev.evaluate(px, data, kinases=todo), kev.evaluate(py, data, kinases=todo)
        full = kev.paired_bootstrap(rA, rB, B=1000, seed=0)
        ax = px.assign(kinase_acc=K[px.kinase_idx.values], site_id=S[px.site_idx.values]).drop(columns=["kinase_idx", "site_idx"])
        ay = py.assign(kinase_acc=K[py.kinase_idx.values], site_id=S[py.site_idx.values]).drop(columns=["kinase_idx", "site_idx"])
        mx, my, rep = kev.matched_universe(ax, ay)
        kk = sorted(set(mx.kinase_acc) & set(todo))
        mA = kev.evaluate(mx, data, kinases=kk, missing="drop")
        mB = kev.evaluate(my, data, kinases=kk, missing="drop")
        matched = kev.paired_bootstrap(mA, mB, B=1000, seed=0)
        key = f"{x} - {y}"
        allres[key] = {"runs": [rx, ry], "full_missing_worst": full.to_dict(orient="records"),
                       "matched_universe": matched.to_dict(orient="records"), "matched_report": rep,
                       "matched_endpoints": {x: mA.endpoints, y: mB.endpoints}}
        print(key, "\nfull\n", full.round(4).to_string(), "\nmatched\n", matched.round(4).to_string(), "\n", rep)
    f.write_text(json.dumps(allres, indent=1, default=kev._json_default))


if __name__ == "__main__":
    main()
