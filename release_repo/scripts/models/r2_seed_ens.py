"""Average the pooled OOF logits of several seed runs of one configuration (seed ensemble), evaluate (ksp.eval A24),
save via RunManifest, register; also report per-seed endpoints mean +- sd.
Usage: r2_seed_ens.py --name r2_final_ens --runs <pooled run ids>"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, manifest
from ksp.index import ROOT
from ksp.models.r2 import train

ap = argparse.ArgumentParser(); ap.add_argument("--name", required=True); ap.add_argument("--runs", nargs="+", required=True)
a = ap.parse_args()
ms = [json.load(open(ROOT / "results" / "runs" / r / "manifest.json")) for r in a.runs]
inputs = [ROOT / m["outputs"]["predictions"]["path"] for m in ms]
with manifest.RunManifest(a.name, config={"seed_runs": a.runs, "base_config": {k: v for k, v in ms[0]["config"].items() if k != "seed"}},
                          seed=-1, split="cv", inputs=inputs, notes="seed ensemble (mean logit)") as run:
    P = None
    for p in inputs:
        x = pd.read_parquet(p).sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
        if P is None:
            P = x[["kinase_idx", "site_idx"]].copy(); acc = x.score.values.astype(np.float64)
        else:
            if not (np.array_equal(P.kinase_idx.values, x.kinase_idx.values) and np.array_equal(P.site_idx.values, x.site_idx.values)):
                raise ValueError("seed runs cover different rows")
            acc += x.score.values
        del x
    P["score"] = (acc / len(inputs)).astype(np.float32)
    run.save_predictions(P)
    data = kev.prepare("cv", fold=None)
    res = kev.evaluate(P, data, missing="error")
    run.save_metrics(res)
    res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
    kev.bootstrap(res, B=1000, seed=0).to_csv(run.dir / "bootstrap.tsv", sep="\t", index=False)
    kev.strata(res).to_csv(run.dir / "strata.tsv", sep="\t", index=False)
    per_seed = []
    for r in a.runs:
        E = json.load(open(ROOT / "results" / "runs" / r / "metrics.json"))["endpoints"]
        per_seed.append({k: E[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany"]} |
                        {"M1_fold_mean": E.get("M1_fold_mean")})
    S = pd.DataFrame(per_seed, index=a.runs)
    S.to_csv(run.dir / "per_seed.tsv", sep="\t")
    E = res.endpoints
    print(S.round(4).to_string()); print("mean", S.mean().round(4).to_dict(), "sd", S.std().round(4).to_dict())
    print(f"[{a.name}] M1={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} M4={E['M4_auroc_hard_Kany']:.4f} "
          f"fold-mean M1={E.get('M1_fold_mean')} per-fold={E.get('M1_per_fold')} canary={res.canary}")
    train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)", "config_name": a.name,
                    "config_json": json.dumps({"seed_runs": a.runs}), "selection_metric": "M1_auroc_hard(K_cv) [=M5]",
                    "value": round(E["M1_auroc_hard"], 5),
                    "notes": f"seed ensemble of {len(a.runs)} (no new config); M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f}"})
