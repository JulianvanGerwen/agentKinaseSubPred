"""Label-variant run through the frozen final R2 config (v3, seed 0): per-fold train + exhaustive OOF predict, then
pooled ksp.eval (A24) + REGISTRY row (model=REGIME).
Usage: regime_run.py --variant L1 --feat refit|L2fixed [--folds 0 1 2 3 4] [--seed 0] [--fold0-eval]
In refit mode each fold waits until its variant xfit stores exist (built by regime_xfit.py).
Resumable: per-fold run ids in results/models/regime/state/<name>.json."""
import argparse, json, os, sys, time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--variant", required=True)
ap.add_argument("--feat", required=True, choices=["refit", "L2fixed"])
ap.add_argument("--folds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--fold0-eval", action="store_true")
ap.add_argument("--name", default=None)
a = ap.parse_args()
from ksp.index import ROOT
from ksp.models.regime import run as R, variants as V

name = a.name or f"regime_{a.variant}_{a.feat}_s{a.seed}"
st_p = ROOT / "results" / "models" / "regime" / "state" / f"{name}.json"
st_p.parent.mkdir(parents=True, exist_ok=True)
state = json.loads(st_p.read_text()) if st_p.exists() else {"folds": {}}
t0 = time.time()
for f in a.folds:
    if str(f) in state["folds"]:
        print(f"fold {f} done earlier: {state['folds'][str(f)]}", flush=True)
        continue
    if a.feat == "refit" and a.variant != "L2":
        need = [Path(V.xfit_tag(a.variant)) / f"fold{f}" / f"{fam}.json" for fam in ("f2a_sep", "f3nb", "f6a")]
        while not all(p.exists() for p in need):
            time.sleep(30)
    out = R.run_fold(a.variant, a.feat, f, f"{name}_f{f}", seed=a.seed, evaluate_fold=(a.fold0_eval and f == 0))
    state["folds"][str(f)] = out["run_id"]
    state.setdefault("logs", {})[str(f)] = out["log"]
    if "summary" in out:
        state["fold0_summary"] = out["summary"]
        R.register({"run_id": out["run_id"], "split": "cv", "fold": "0", "config_name": f"{name}_fold0",
                    "config_json": json.dumps({"variant": a.variant, "feat": a.feat, "seed": a.seed,
                                               "base": "r2 final v3"}),
                    "selection_metric": "M1_auroc_hard(K_cv)", "value": round(out["summary"]["M1"], 5),
                    "notes": f"fold-0 pilot (timing); M2={out['summary']['M2']:.4f}"})
    st_p.write_text(json.dumps(state, indent=1, default=str))
if all(str(f) in state["folds"] for f in range(5)) and "pooled" not in state:
    cfg = {**R.final_config(a.seed), "variant": a.variant, "feat": a.feat}
    po = R.pool(f"{name}_cv5", [state["folds"][str(f)] for f in range(5)], cfg,
                notes=f"label variant {a.variant}, features {a.feat}; train wall {time.time() - t0:.0f}s")
    state["pooled"] = po["run_id"]
    state["summary"] = po["summary"]
    st_p.write_text(json.dumps(state, indent=1, default=str))
print("DONE", name, time.time() - t0, flush=True)
