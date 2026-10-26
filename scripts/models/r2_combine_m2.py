"""Evaluate a binary R2 run (score -> M1/M3/M4) together with a site-grouped R2 run (score_site -> M2), A20.
Both runs must cover the same folds and identical (kinase, site) rows.
Usage: r2_combine_m2.py --score-run <run_id> --site-run <run_id> [--fold 0] [--name r2_combo]"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, manifest
from ksp.index import ROOT
from ksp.models.r2 import train

ap = argparse.ArgumentParser()
ap.add_argument("--score-run", required=True); ap.add_argument("--site-run", required=True)
ap.add_argument("--fold", type=int, default=None); ap.add_argument("--name", default="r2_combo_m2")
ap.add_argument("--stage", default="d")
a = ap.parse_args()


def load(r):
    m = json.load(open(ROOT / "results" / "runs" / r / "manifest.json"))
    return pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"]), m


pa, ma = load(a.score_run)
pb, mb = load(a.site_run)
pa = pa.sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
pb = pb.sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
if not (np.array_equal(pa.kinase_idx.values, pb.kinase_idx.values) and np.array_equal(pa.site_idx.values, pb.site_idx.values)):
    raise ValueError("runs do not cover identical rows")
P = pa.assign(score_site=pb.score.values)
cfg = {"score_run": a.score_run, "site_run": a.site_run, "score_cfg": ma["config"], "site_cfg": mb["config"]}
with manifest.RunManifest(a.name, config=cfg, seed=ma["config"]["seed"], split="cv",
                          notes=f"folds={'all' if a.fold is None else a.fold}; score from {a.score_run}, score_site from {a.site_run}") as run:
    data = kev.prepare("cv", fold=a.fold)
    res = kev.evaluate(P, data, missing="error")
    run.save_metrics(res)
    E = res.endpoints
    print(f"[{a.name}] M1={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} "
          f"hit1={E['sec_hit1']:.4f} mrr={E['sec_mrr']:.4f}")
    train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)" if a.fold is None else a.fold,
                    "config_name": a.name, "config_json": json.dumps({"score_run": a.score_run, "site_run": a.site_run}),
                    "selection_metric": "M2_hit5(K_cv)", "value": round(E["M2_hit5"], 5),
                    "notes": f"{a.stage}; combination (no new model); M1={E['M1_auroc_hard']:.4f}"})
