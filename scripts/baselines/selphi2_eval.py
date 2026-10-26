"""Evaluate B4 SELPHI2.0 OOF predictions with ksp.eval (DESIGN §6; matched-universe rule §6.8 / PRIOR-1).

Reported per arm:
  coverage            : ksp.eval.coverage (required pairs = positives + E-all candidates of every kinase with a positive)
  full   (missing='worst'): unscored pairs ranked last (deployment view, secondary)
  matched (missing='drop'): only pairs B4 scores = U_B4 (U_ms sites with Ochoa features) - the claim-bearing view
  reference on U_B4   : B1 Cantley (universe percentile, unprimed; label-free) on exactly the same pairs, so the reader
                        can see what U_B4 restriction does to a comparator.
Writes results/baselines/selphi2/<arm>/metrics_{full,matched}.json, per_kinase_matched.tsv, summary.json and a
ksp.manifest run (cache/runs/<run_id>/predictions.parquet + results/runs/<run_id>/).
Usage: python -I scripts/baselines/selphi2_eval.py --arm faithful [--B 1000]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ksp import eval as E  # noqa: E402
from ksp.baselines.selphi2 import features as F  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402
from ksp.manifest import RunManifest  # noqa: E402


def cantley_on(preds: pd.DataFrame) -> pd.DataFrame:
    from ksp.store import load_matrix
    M, _ = load_matrix("motif/cantley_unprimed_pct_universe")
    kp = pd.Series(np.arange(len(load_kinases())), index=load_kinases().uniprot_acc)
    sp = pd.Series(np.arange(len(load_sites())), index=load_sites().site_id)
    ki, si = kp[preds.kinase_acc].to_numpy(), sp[preds.site_id].to_numpy()
    v = np.empty(len(ki), np.float32)
    o = np.argsort(ki, kind="stable")
    for k in np.unique(ki):
        sel = o[ki[o] == k]
        v[sel] = np.asarray(M[k])[si[sel]]
    out = preds[["kinase_acc", "site_id"]].assign(score=v)
    return out[np.isfinite(out.score)].reset_index(drop=True)


def summarise(res: E.EvalResult, boot: pd.DataFrame | None) -> dict:
    d = {k: v for k, v in res.endpoints.items() if not isinstance(v, (list, dict))}
    if boot is not None:
        d["bootstrap"] = boot.to_dict(orient="records")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="faithful")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--no-cantley", action="store_true")
    a = ap.parse_args()
    out = ROOT / "results" / "baselines" / "selphi2" / a.arm
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    files = sorted((F.CACHE / a.arm).glob("oof_fold*.parquet"))
    if len(files) != 5:
        raise FileNotFoundError(f"expected 5 OOF fold files for arm {a.arm}, found {len(files)}")
    preds = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    fold_json = [json.loads((out / f"fold{i}.json").read_text()) for i in range(5)]
    cfg = fold_json[0]["config"]
    with RunManifest(f"b4_selphi2_{a.arm}", config={"arm": a.arm, **cfg,
                                                     "best_params_per_fold": [f["best_params"] for f in fold_json]},
                     seed=cfg.get("seed", 0), inputs=files, split="cv",
                     notes="B4 SELPHI2.0 retrained, OOF over DEV 5-fold CV; universe U_B4 = U_ms with Ochoa features")\
            as run:
        run.save_predictions(preds)
        data = E.prepare("cv")
        cov = E.coverage(preds, data)
        run.set("coverage", cov)
        print("coverage", cov, flush=True)
        res_m = E.evaluate(preds, data, missing="drop")
        res_f = E.evaluate(preds, data, missing="worst")
        boot_m = E.bootstrap(res_m, B=a.B)
        boot_f = E.bootstrap(res_f, B=a.B)
        res_m.save_json(out / "metrics_matched.json")
        res_f.save_json(out / "metrics_full.json")
        run.save_metrics(res_m, "metrics_matched")
        run.save_metrics(res_f, "metrics_full")
        res_m.per_kinase.to_csv(out / "per_kinase_matched.tsv", sep="\t", index=False)
        summ = {"arm": a.arm, "coverage": cov, "matched_U_B4": summarise(res_m, boot_m),
                "full_worst": summarise(res_f, boot_f),
                "strata_matched": E.strata(res_m).to_dict(orient="records"),
                "canary_matched": res_m.canary}
        if not a.no_cantley:
            cp = cantley_on(preds)
            run.count_fallback("cantley_ref_pairs_without_matrix", len(preds) - len(cp))
            res_c = E.evaluate(cp, data, missing="drop")
            summ["cantley_B1_recal_on_U_B4"] = summarise(res_c, E.bootstrap(res_c, B=a.B))
            summ["cantley_B1_recal_on_U_B4"]["n_pairs"] = len(cp)
            try:
                pa, pb, rep = E.matched_universe(preds, cp)
                ra, rb = E.evaluate(pa, data, missing="drop"), E.evaluate(pb, data, missing="drop")
                summ["paired_selphi_minus_cantley_on_common"] = {
                    "universe": rep, "selphi": summarise(ra, None), "cantley": summarise(rb, None),
                    "paired": E.paired_bootstrap(ra, rb, B=a.B).to_dict(orient="records")}
            except Exception as e:   # recorded, not hidden
                summ["paired_selphi_minus_cantley_on_common"] = {"error": repr(e)}
        summ["eval_seconds"] = round(time.time() - t0, 1)
        summ["train_timing_s_per_fold"] = [f["timing_s"] for f in fold_json]
        summ["run_id"] = run.run_id
        run.save_metrics(summ, "summary")
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=E._json_default))
    print(json.dumps({k: summ[k] for k in ("coverage", "matched_U_B4", "full_worst")}, indent=1,
                     default=E._json_default))


if __name__ == "__main__":
    main()
