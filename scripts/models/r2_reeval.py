"""Re-evaluate pooled R2 runs with the CURRENT ksp.eval (A24 within-(fold,kinase) rank normalisation for pooled CV)
and collect per-fold M1 from the fold runs. Appends to results/models/r2/cv_table.tsv.
Usage: r2_reeval.py --runs <pooled run_id> ..."""
import argparse, json, time
import numpy as np, pandas as pd
from ksp import eval as kev
from ksp.index import ROOT

ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+"); ap.add_argument("--label", default="")
a = ap.parse_args()
data = kev.prepare("cv", fold=None)
out = ROOT / "results" / "models" / "r2" / "cv_table.tsv"
for rid in a.runs:
    m = json.load(open(ROOT / "results" / "runs" / rid / "manifest.json"))
    P = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
    res = kev.evaluate(P, data, missing="error")
    E, C = res.endpoints, res.canary
    bs = kev.bootstrap(res, B=1000, seed=0).set_index("endpoint")
    fm = []
    for fr in m.get("fold_runs", []):
        mm = json.load(open(ROOT / "results" / "runs" / fr / "metrics.json"))
        fm.append(mm["endpoints"]["M1_auroc_hard"])
    d = ROOT / "results" / "models" / "r2" / "reeval_A24"; d.mkdir(parents=True, exist_ok=True)
    res.save_json(d / f"{rid}.json")
    res.per_kinase.to_csv(d / f"{rid}.per_kinase.tsv", sep="\t", index=False)
    row = {"time": time.strftime("%H:%M"), "run_id": rid, "name": m["name"], "groups": "+".join(m["config"]["groups"]),
           "M1": E["M1_auroc_hard"], "M1_lo": bs.loc["M1", "lo95"], "M1_hi": bs.loc["M1", "hi95"],
           "M2": E["M2_hit5"], "M3": E["M3_top1_recall"], "M3_lo": bs.loc["M3", "lo95"], "M3_hi": bs.loc["M3", "hi95"],
           "M4": E["M4_auroc_hard_Kany"], "perfold_M1_mean": float(np.mean(fm)) if fm else np.nan,
           "perfold_M1": ";".join(f"{x:.4f}" for x in fm), "neg_canary_hard": C.get("neg_mean_auroc_hard"),
           "neg_canary_all": C.get("neg_mean_auroc_all"), "pos_canary_hard": C.get("pos_auroc_hard"),
           "canary_neg_pass": C.get("neg_pass"), "canary_pos_pass": C.get("pos_pass"), "eval": "A24 rank-normalised"}
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False, mode="a", header=not out.exists())
    print(json.dumps(row, default=str), flush=True)
