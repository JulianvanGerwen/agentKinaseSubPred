"""Backfill results/cv/REGISTRY.tsv with R2 per-fold component runs (run with --no-register by r2_cv5.sh) and a
stage-(a) selection-context note (audit 17:30). Idempotent: skips run_ids already in the registry."""
import json, glob, time
import pandas as pd
from ksp.index import ROOT
REG = ROOT / "results" / "cv" / "REGISTRY.tsv"
reg = pd.read_csv(REG, sep="\t")
have = set(reg.run_id.dropna())
cols = list(reg.columns)
rows = []
pooled = {}
for p in glob.glob(str(ROOT / "results" / "runs" / "r2_*_cv5__*" / "manifest.json")):
    m = json.load(open(p))
    for fr in m.get("fold_runs", []):
        pooled[fr] = m["run_id"]
for p in sorted(glob.glob(str(ROOT / "results" / "runs" / "r2_*" / "manifest.json"))):
    m = json.load(open(p))
    rid = m["run_id"]
    if rid in have or m.get("status") != "ok" or "fold_runs" in m or not (ROOT / "results" / "runs" / rid / "metrics.json").exists():
        continue
    if m["name"].startswith(("r2_fusion", "r2_partial", "r2_combo")):
        continue
    E = json.load(open(ROOT / "results" / "runs" / rid / "metrics.json"))["endpoints"]
    fold = m["notes"].split("folds=[")[1].split("]")[0] if "folds=[" in m["notes"] else ""
    note = (f"per-fold component of pooled {pooled[rid]}" if rid in pooled else
            "per-fold run NOT pooled (superseded: stale canary stores / killed)") + \
           f"; M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f}; fold-level eval (no A24 needed within one fold)"
    rows.append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R2", "run_id": rid, "split": "cv",
                 "fold": fold, "config_name": m["name"], "config_json": json.dumps({k: v for k, v in m["config"].items() if k != "groups"} | {"groups": "+".join(m["config"]["groups"])}),
                 "selection_metric": "M1_auroc_hard(K_cv fold)", "value": round(E["M1_auroc_hard"], 5),
                 "n_configs_tried_model": "", "notes": "backfill(audit): " + note})
rows.append({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "model": "R2", "run_id": "", "split": "cv", "fold": "0",
             "config_name": "stage_a_selection_context", "config_json": "", "selection_metric": "", "value": "",
             "n_configs_tried_model": "", "notes": "Stage (a)-(d) were selected on fold-0 M1 with canary-v1/intermediate "
             "xfit stores before A24; 19 fold-0 configs + 1 label-free control; chosen leaves127 lr0.03 a0.5 1:50 h0.5 binary; "
             "final v3 adds f2a_sep (canaries fitted separately), descriptor v2 (A25), ES on real kinases only"})
pd.DataFrame(rows)[cols].to_csv(REG, sep="\t", index=False, header=False, mode="a")
print(f"appended {len(rows)} rows")
