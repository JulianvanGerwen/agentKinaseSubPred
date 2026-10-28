"""R5: collate run summaries into the tables used in docs/model_r5.md.

  python -I scripts/models/r5_report.py folds            # every R5 CV training run, newest first
  python -I scripts/models/r5_report.py pooled <run_id>  # one pooled run: endpoints, bins, canaries
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
RUNS = ROOT / "results" / "runs"


def _summary(d: Path):
    f = d / "summary.json"
    return json.loads(f.read_text()) if f.exists() else None


def folds():
    rows = []
    for d in sorted(RUNS.glob("r5_*"), key=lambda p: p.name):
        s = _summary(d)
        if s is None:
            continue
        m = json.loads((d / "manifest.json").read_text())
        ti = s.get("train_info", {}) or {}
        c = s.get("canary") or {}
        rows.append({
            "run": d.name, "config": m["config"].get("trunk", ""), "fold": m["config"].get("_fold", ""),
            "M1": s.get("M1_auroc_hard"), "M2": s.get("M2_hit5"), "M3": s.get("M3_top1_recall"),
            "M4": s.get("M4_auroc_hard_Kany"), "innerval": ti.get("ensemble_val_auroc_hard")
            or ti.get("best_val_auroc_hard"),
            "neg_can": c.get("neg_mean_auroc_hard"), "pos_can": c.get("pos_auroc_hard"),
            "t_train_s": ti.get("t_train_s"), "t_score_s": ti.get("t_score_s"), "t_eval_s": s.get("t_eval_s"),
            "t_total_s": ti.get("t_total_s"), "peak_gb": ti.get("peak_gpu_mem_gb"), "n_pred": ti.get("n_pred_rows"),
        })
    hdr = list(rows[0]) if rows else []
    print("\t".join(hdr))
    for r in rows:
        print("\t".join("" if r[h] is None else (f"{r[h]:.4f}" if isinstance(r[h], float) else str(r[h]))
                        for h in hdr))


def pooled(run_id: str):
    s = _summary(RUNS / run_id)
    print(json.dumps({k: s[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany",
                                        "M5_auroc_hard_Kcv", "sec_auroc_all", "sec_mrr", "sec_hit1", "sec_hit10",
                                        "n_K_primary_evaluated", "n_K_any_evaluated", "n_missing_pairs",
                                        "n_pred_rows", "t_eval_s"] if k in s}, indent=1))
    print("canary:", json.dumps(s.get("canary"), indent=1))
    print("M5 by train bin:")
    for b in s.get("M5_by_train_bin", []):
        print(" ", b.get("value"), b.get("n_kinases"), round(b.get("auroc_hard", float("nan")), 4))
    for b in s.get("bootstrap", []) or []:
        print("boot:", json.dumps(b))


if __name__ == "__main__":
    if sys.argv[1] == "folds":
        folds()
    else:
        pooled(sys.argv[2])
