"""Run one R2 configuration on CV folds (train, exhaustive OOF predict, ksp.eval, REGISTRY row).
Usage: r2_run.py --name r2_base --folds 0 [--set num_leaves=127 alpha=1.0 groups='["F1","F0"]'] [--stage a]"""
import argparse, json, sys
from ksp.models.r2 import train

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True)
ap.add_argument("--folds", nargs="+", type=int, default=[0])
ap.add_argument("--set", nargs="*", default=[])
ap.add_argument("--stage", default="")
ap.add_argument("--notes", default="")
ap.add_argument("--no-register", action="store_true")
a = ap.parse_args()
cfg = train.default_config()
for kv in a.set:
    k, v = kv.split("=", 1)
    if k not in cfg and k not in ("site_max_cand",):
        sys.exit(f"unknown config key {k}")
    try:
        cfg[k] = json.loads(v)
    except json.JSONDecodeError:
        cfg[k] = v
print("config:", json.dumps(cfg), flush=True)
out = train.run(cfg, folds=a.folds, name=a.name, config_name=a.name, stage=a.stage, notes=a.notes,
                register_row=not a.no_register)
print("RESULT", json.dumps({"run_id": out["run_id"], "endpoints": out.get("endpoints"), "canary": out.get("canary")},
                           default=str), flush=True)
