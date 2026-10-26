"""Pool per-fold R2 runs (each trained/predicted on one CV fold) into one 5-fold OOF run: concatenate predictions,
evaluate with ksp.eval.prepare('cv') (pooled), save via RunManifest, register in REGISTRY.tsv.
Usage: r2_pool.py --name r2_x_cv5 --runs <run_id fold0> ... <run_id fold4> [--stage s] [--notes ..]"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, manifest
from ksp.index import ROOT
from ksp.models.r2 import train

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--runs", nargs="+", required=True)
ap.add_argument("--stage", default=""); ap.add_argument("--notes", default="")
ap.add_argument("--no-register", action="store_true")
a = ap.parse_args()
ms = [json.load(open(ROOT / "results" / "runs" / r / "manifest.json")) for r in a.runs]
_dflt = train.default_config()
cfgs = [{**_dflt, **m["config"]} for m in ms]          # keys added later to default_config keep their defaults
cfg = cfgs[0]
for c in cfgs[1:]:
    if c != cfg:
        raise ValueError(f"pooled runs have different configs: {set(k for k in c if c[k] != cfg.get(k))}")
folds = sorted(int(m["notes"].split("folds=[")[1].split("]")[0]) for m in ms)
if folds != [0, 1, 2, 3, 4]:
    raise ValueError(f"need folds 0..4 exactly, got {folds}")
inputs = [ROOT / m["outputs"]["predictions"]["path"] for m in ms]
with manifest.RunManifest(a.name, config=cfg, seed=cfg["seed"], split="cv", inputs=inputs,
                          notes=f"pooled 5-fold OOF of {a.runs}; {a.notes}") as run:
    P = pd.concat([pd.read_parquet(p) for p in inputs], ignore_index=True)
    run.set("fold_runs", a.runs)
    run.set("fold_logs", {m["run_id"]: m.get("fold_logs") for m in ms})
    run.save_predictions(P)
    data = kev.prepare("cv", fold=None)
    res = kev.evaluate(P, data, missing="error")
    run.save_metrics(res)
    res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
    bs = kev.bootstrap(res, B=1000, seed=0)
    bs.to_csv(run.dir / "bootstrap.tsv", sep="\t", index=False)
    st = kev.strata(res); st.to_csv(run.dir / "strata.tsv", sep="\t", index=False)
    E = res.endpoints
    print(f"[{a.name}] M1/M5={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} "
          f"M4={E['M4_auroc_hard_Kany']:.4f} canary={res.canary}", flush=True)
    print(bs.to_string())
    if not a.no_register:
        train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)", "config_name": a.name,
                        "config_json": json.dumps({k: v for k, v in cfg.items() if k != "groups"} | {"groups": "+".join(cfg["groups"])}),
                        "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(E["M1_auroc_hard"], 5),
                        "notes": f"{a.stage}; M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} M4={E['M4_auroc_hard_Kany']:.4f}; "
                                 f"canary_neg={res.canary.get('neg_mean_auroc_hard'):.3f}/{res.canary.get('neg_mean_auroc_all'):.3f} "
                                 f"pos={res.canary.get('pos_auroc_hard'):.3f}/{res.canary.get('pos_auroc_all'):.3f}; {a.notes}"})
    print("RESULT", json.dumps({"run_id": run.run_id, "endpoints": E, "canary": res.canary}, default=str))
