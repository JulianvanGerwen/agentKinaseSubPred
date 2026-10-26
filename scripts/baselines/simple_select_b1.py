"""Choose the official B1 'Cantley best arm' on CV M1 (DESIGN §16 A19) and write results/baselines/simple/b1_best_arm.json.

Pre-declared rule: candidates = the three DESIGN arms (i) b1_published, (ii) b1_universe, (iii) b1_primed_hard /
b1_primed_soft (universe-recalibrated percentiles). Highest pooled-CV M1 (K_cv, E-hard, missing='worst') wins;
tie-break M2. Supplementary variants (primed + published percentile, nearest-atlas imputed) are NOT candidates.
Also: paired cluster bootstrap (ksp.eval.paired_bootstrap, B=1000) of best vs each other candidate.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from ksp import eval as kev  # noqa: E402

OUT = ROOT / "results" / "baselines" / "simple"
CANDS = ["b1_published", "b1_universe", "b1_primed_hard", "b1_primed_soft"]


def main():
    rows = []
    for m in CANDS:
        d = json.loads((OUT / f"{m}.json").read_text())
        rows.append({"method": m, "store": d["config"]["store"], "run_id": d["run_id"],
                     "M1": d["endpoints"]["M1_auroc_hard"], "M2": d["endpoints"]["M2_hit5"],
                     "M3": d["endpoints"]["M3_top1_recall"], "M4": d["endpoints"]["M4_auroc_hard_Kany"]})
    t = pd.DataFrame(rows).sort_values(["M1", "M2"], ascending=False).reset_index(drop=True)
    best = t.iloc[0]
    data = kev.prepare("cv")
    todo = [k for k in sorted(data.positives) if not k.startswith("CANARY_")]
    res = {}
    for r in t.itertuples():
        p = pd.read_parquet(ROOT / "cache" / "runs" / r.run_id / "predictions.parquet")
        res[r.method] = kev.evaluate(p, data, kinases=todo, missing="worst")
    paired = []
    for m in t.method[1:]:
        pb = kev.paired_bootstrap(res[best.method], res[m], B=1000, seed=0)
        pb.insert(0, "comparison", f"{best.method} - {m}")
        paired.append(pb)
    paired = pd.concat(paired, ignore_index=True)
    out = {"rule": "max pooled-CV M1 (K_cv, E-hard, missing=worst) among DESIGN arms (i),(ii),(iii); tie-break M2",
           "best_method": best.method, "store": best.store, "run_id": best.run_id,
           "candidates": t.to_dict(orient="records"), "paired_bootstrap_best_vs_other": paired.to_dict(orient="records")}
    (OUT / "b1_best_arm.json").write_text(json.dumps(out, indent=1, default=kev._json_default))
    print(t.to_string())
    print(paired.to_string())
    print("BEST:", best.method, best.store)


if __name__ == "__main__":
    main()
