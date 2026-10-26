"""Run one simple baseline on the 5-fold DEV CV, evaluate ONLY with ksp.eval, save via ksp.manifest.

    OMP_NUM_THREADS=8 .venv/bin/python -I scripts/baselines/simple_run.py --method b1_universe [--canary]

Outputs
  cache/runs/<run_id>/predictions.parquet      OOF predictions (kinase_idx, site_idx, score) over all eval candidates
  results/runs/<run_id>/{manifest,metrics,bootstrap,extra}.json
  results/baselines/simple/<method>[_canary].json   compact copy (endpoints, CIs, coverage, runtime, run_id)
  results/baselines/simple/<method>_per_kinase.tsv  per-kinase metrics
Methods: see ksp.baselines.simple. --canary trains with load_dev(with_canaries=True) and scores canary kinases only.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp import canary as kcan  # noqa: E402
from ksp import eval as kev  # noqa: E402
from ksp import manifest  # noqa: E402
from ksp.baselines import simple as B  # noqa: E402

OUT = ROOT / "results" / "baselines" / "simple"
REGISTRY = ROOT / "results" / "cv" / "REGISTRY.tsv"
SEED = 0


def registry_append(model, run_id, config_name, config, value, notes):
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    new = not REGISTRY.exists()
    n_prev = 0
    if not new:
        r = pd.read_csv(REGISTRY, sep="\t")
        n_prev = int((r.model == model).sum())
    with open(REGISTRY, "a") as fh:
        if new:
            fh.write("timestamp\tmodel\trun_id\tsplit\tfold\tconfig_name\tconfig_json\tselection_metric\tvalue\t"
                     "n_configs_tried_model\tnotes\n")
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{model}\t{run_id}\tcv\tpooled\t{config_name}\t"
                 f"{json.dumps(config)}\tM1_auroc_hard(K_cv)\t{value:.5f}\t{n_prev + 1}\t{notes}\n")


def best_cantley_arm() -> str:
    f = OUT / "b1_best_arm.json"
    if not f.exists():
        raise FileNotFoundError("run scripts/baselines/simple_select_b1.py first (B1 best arm chosen on CV M1)")
    return json.loads(f.read_text())["store"]


def compute(method: str, cand: pd.DataFrame, data, canary: bool, log: dict) -> tuple[np.ndarray, dict]:
    cfg = {"method": method, "seed": SEED}
    if method == "b0_prior":
        return B.b0_prior(cand, data), cfg
    if method == "b0_pop":
        return B.b0_pop(cand, data), cfg
    if method == "b0_miss":
        cfg.update(neg_ratio=50, hard_frac=0.5, lgb="lr0.05 leaves31 minleaf100 es50 on inner-cluster 20% split")
        return B.b0_miss(cand, data, seed=SEED, log=log), cfg
    if method in B.B1_ARMS:
        st = B.B1_ARMS[method][0]
        cfg["store"] = st
        return B.cantley(cand, st), cfg
    if method.startswith("b1_imputed_"):
        st = B.B1_ARMS["b1_" + method[len("b1_imputed_"):]][0]
        cfg.update(store=st, imputed=True)
        return B.cantley(cand, st, imputed=True), cfg
    if method in ("b2_pssm", "b2_pssm_ph"):
        cfg.update(prior="family", phospho=method.endswith("_ph"), score="universe percentile")
        return B.b2_pssm(cand, data, phospho=method.endswith("_ph"), with_canaries=canary, seed=SEED, log=log), cfg
    if method == "b7_lr":
        st = best_cantley_arm()
        cfg.update(cantley_store=st, C=1.0, neg_ratio=50, hard_frac=0.5,
                   features=["logit Cantley best arm", "logit B2 pct (cross-fitted)", "ppi_edge_any",
                             "loc_cos_combined", "+ missing indicators"], kinase_intercept="L2 one-hot, IPW weights")
        return B.b7_lr(cand, data, st, with_canaries=canary, seed=SEED, log=log), cfg
    raise ValueError(f"unknown method {method}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True)
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--B", type=int, default=1000)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tag = a.method + ("_canary" if a.canary else "")
    t_all = time.time()
    data = kev.prepare("cv")
    if a.canary:
        todo = [k for k in kcan.CANARY_IDS if k in data.positives]
        cand = B.candidate_pairs(data, kinases=todo, include_m2=False)
    else:
        todo = [k for k in sorted(data.positives) if not k.startswith("CANARY_")]
        cand = B.candidate_pairs(data, kinases=todo)
    log: dict = {}
    inputs = [ROOT / "cache" / "labels" / "dev_pairs.parquet", ROOT / "results" / "splits" / "cv_folds.tsv"]
    with manifest.RunManifest(f"simple_{tag}", config={"method": a.method, "canary": a.canary}, seed=SEED,
                              inputs=inputs, split="cv") as run:
        t0 = time.time()
        score, cfg = compute(a.method, cand, data, a.canary, log)
        t_fit = time.time() - t0
        site_score = None
        if isinstance(score, dict):       # separate within-kinase score and site-centric (M2) score
            score, site_score = score["score"], score["score_site"]
        ok = np.isfinite(score)
        run.count_fallback("candidate_pairs_not_scored", int((~ok).sum()))
        preds = pd.DataFrame({"kinase_idx": cand.kinase_idx.values[ok].astype(np.int16),
                              "site_idx": cand.site_idx.values[ok].astype(np.int32),
                              "score": score[ok].astype(np.float32)})
        if site_score is not None:
            if not np.isfinite(site_score[ok]).all():
                raise ValueError("score_site non-finite where score is finite")
            preds["score_site"] = site_score[ok].astype(np.float32)
        run.m["config"].update(cfg)
        run.save_predictions(preds)
        t1 = time.time()
        res = kev.evaluate(preds, data, kinases=todo, missing="worst")
        extra = {"method": a.method, "canary": a.canary, "run_id": run.run_id, "config": cfg,
                 "runtime_fit_predict_s": round(t_fit, 1), "n_candidates": int(len(cand)),
                 "n_scored": int(ok.sum()), "missing": {k: v for k, v in res.missing.items() if k != "per_kinase"},
                 "log": log}
        if a.canary:
            extra["canary"] = res.canary
            extra["canary_per_kinase"] = res.per_kinase[["kinase_acc", "n_pos", "n_neg_hard", "n_neg_all",
                                                         "auroc_hard", "auroc_all", "top1_recall"]] \
                .to_dict(orient="records")
        else:
            boot = kev.bootstrap(res, B=a.B, seed=SEED)
            pk = res.per_kinase
            scored_k = set(data.kinases.uniprot_acc.values[np.unique(preds.kinase_idx.values)])
            kcv = set(data.k_primary)
            extra["coverage"] = {
                "kinases_with_positives": len(todo), "kinases_scored": len(set(todo) & scored_k),
                "K_cv": len(kcv), "K_cv_scored": len(kcv & scored_k),
                "real_kinases_scored_any": len(scored_k - set(kcan.CANARY_IDS)),
                "frac_required_pairs_scored": 1 - res.missing["pairs_missing"] / max(res.missing["pairs_required"], 1)}
            extra["endpoints"] = res.endpoints
            extra["bootstrap"] = boot.to_dict(orient="records")
            # covered-only view (missing='drop'): only for methods with incomplete coverage
            if res.missing["pairs_missing"] > 0:
                rd = kev.evaluate(preds, data, kinases=sorted(set(todo) & scored_k), missing="drop",
                                  keep_items=False)
                extra["endpoints_covered_only"] = rd.endpoints
            pk.to_csv(OUT / f"{a.method}_per_kinase.tsv", sep="\t", index=False)
            extra["strata"] = kev.strata(res).to_dict(orient="records")
            extra["strata_K_cv"] = kev.strata(res, kinase_set="K_primary").to_dict(orient="records")
            extra["pair_strata"] = res.pair_strata.to_dict(orient="records")
            run.save_metrics(boot.to_dict(orient="records"), "bootstrap")
        extra["runtime_eval_s"] = round(time.time() - t1, 1)
        extra["runtime_total_s"] = round(time.time() - t_all, 1)
        run.save_metrics(res)
        run.save_metrics(extra, "extra")
        (OUT / f"{tag}.json").write_text(json.dumps(extra, indent=1, default=kev._json_default))
        if not a.canary:
            registry_append(f"simple:{a.method}", run.run_id, a.method, cfg, res.endpoints["M1_auroc_hard"],
                            "simple baseline (no tuning)")
        E = res.endpoints if not a.canary else res.canary
        print(json.dumps({"tag": tag, "run_id": run.run_id, "fit_s": round(t_fit, 1),
                          **({k: E[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall",
                                                "M4_auroc_hard_Kany"]} if not a.canary else E)},
                         default=kev._json_default))


if __name__ == "__main__":
    main()
