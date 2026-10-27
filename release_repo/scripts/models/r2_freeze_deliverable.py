"""Assemble the final R2+R3 fusion deliverable for G4 (docs/G4_RUNBOOK.md): one run directory holding
test_predictions.parquet (score, score_site, assignment_prob) and test_t3_predictions.parquet, with train/inference
wall-clock = sum over the component runs (R2 5-seed DEV models, R3 DEV model, fusion and assignment steps), and the
component run ids + sha256 pins. Prints the FREEZE.json entry.
Usage: r2_freeze_deliverable.py --test-run <assign run> --t3-run <fused t3 run> --r2-test <run> --r2-t3 <run> --r3-run <run>"""
import argparse, json, shutil
import pandas as pd
from ksp import manifest
from ksp.index import ROOT

ap = argparse.ArgumentParser()
for k in ["test-run", "t3-run", "r2-test", "r2-t3", "r3-run", "fuse-test"]:
    ap.add_argument(f"--{k}", required=True)
a = ap.parse_args()
M = lambda r: json.load(open(ROOT / "results" / "runs" / r / "manifest.json"))
comp = {"r2_test": M(a.r2_test), "r2_t3": M(a.r2_t3), "r3": M(a.r3_run), "fuse_test": M(a.fuse_test),
        "assign_test": M(a.test_run), "fuse_t3": M(a.t3_run)}
train_s = comp["r2_test"]["train_wall_s"] + comp["r3"]["train_wall_s"]
inf_s = (comp["r2_test"]["inference_wall_s"] + comp["r3"]["inference_wall_s"] + comp["fuse_test"]["runtime_s"]
         + comp["assign_test"]["runtime_s"])
src_test = ROOT / comp["assign_test"]["outputs"]["test_predictions"]["path"]
src_t3 = ROOT / comp["fuse_t3"]["outputs"]["test_t3_predictions"]["path"]
cfg = {"model": "R2 (v3, 5-seed ensemble) + R3 late fusion", "w_score": 0.6, "w_score_site": 0.5,
       "components": {k: v["run_id"] for k, v in comp.items()},
       "component_sha256": {k: {n: o["sha256"] for n, o in v["outputs"].items()} for k, v in comp.items()}}
with manifest.RunManifest("r2r3_final", config=cfg, seed=0, split="final_dev", inputs=[src_test, src_t3],
                          notes="G4 deliverable: final fused model; TEST + T3 predictions") as run:
    run.save_predictions(pd.read_parquet(src_test), name="test_predictions")
    run.save_predictions(pd.read_parquet(src_t3), name="test_t3_predictions")
    run.set("train_wall_s", round(train_s, 1)); run.set("inference_wall_s", round(inf_s, 1))
    run.set("t3_note", f"T3: R2 part {comp['r2_t3']['config'].get('seeds')} seeds per LKO fold, fused with R3 T3 file")
    entry = {"name": "final", "role": "final", "run_id": run.run_id,
             "predictions": f"cache/runs/{run.run_id}/test_predictions.parquet",
             "score_site_predictions": None, "assignment_prob_predictions": None,
             "t3_predictions": f"cache/runs/{run.run_id}/test_t3_predictions.parquet",
             "runtime": {"train_wall_s": round(train_s, 1), "inference_wall_s": round(inf_s, 1)}}
    run.set("freeze_entry", entry)
print(json.dumps(entry, indent=1))
