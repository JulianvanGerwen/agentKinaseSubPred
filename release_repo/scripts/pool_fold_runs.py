#!/usr/bin/env python
"""Pool per-fold prediction runs into ONE run (new run dir via ksp.manifest) and evaluate it with ksp.eval (A24).
Owner: core-infra.

  CV  : 5 CV-fold runs (each = OOF predictions for one DEV fold) -> pooled OOF; fold_normalise=True (A24).
  LKO : 5 leave-kinase-out runs (each = held-out kinases x DEV sites) -> pooled T2; fold_normalise=False
        (one model per kinase).
Checks: every source run exists with cache/runs/<id>/predictions.parquet; no duplicate (kinase, site) across sources;
CV sources must not overlap in substrate fold (each source's sites lie in one CV fold, distinct across sources);
LKO sources must not overlap in kinases.
Writes cache/runs/<new>/predictions.parquet, results/runs/<new>/{manifest,metrics,metrics_a24}.json and
bootstrap.tsv (cluster bootstrap, --B resamples) so scripts/reeval_cv_runs.py lists it without recomputation.
Usage: pool_fold_runs.py --name r3_final_distill_cv_pooled --mode cv --runs id0 id1 id2 id3 id4 [--B 1000]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import index, splits  # noqa: E402
from ksp.manifest import RunManifest  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--mode", choices=["cv", "lko"], required=True)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--notes", default="")
    a = ap.parse_args()
    t0 = time.time()
    srcs = [ROOT / "cache" / "runs" / r / "predictions.parquet" for r in a.runs]
    for s in srcs:
        if not s.exists():
            sys.exit(f"missing {s}")
    S = index.load_sites()
    pos = pd.Series(np.arange(len(S)), index=S.site_id)
    sfold = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values
    parts, info = [], {}
    for r, s in zip(a.runs, srcs):
        cols = pd.read_parquet(s, columns=None).columns if False else None  # noqa: F841 (schema read below)
        p = pd.read_parquet(s)
        if {"kinase_idx", "site_idx"} <= set(p.columns):
            p = p[["kinase_idx", "site_idx", "score"]]
        else:
            kpos = {k: i for i, k in enumerate(E._canary.kinases_with_canaries().uniprot_acc)}
            p = pd.DataFrame({"kinase_idx": p.kinase_acc.map(kpos).values, "site_idx": pos[p.site_id].values,
                              "score": p.score.values})
        f = np.unique(sfold[p.site_idx.values])
        info[r] = {"n_rows": int(len(p)), "cv_folds_of_sites": f.tolist(),
                   "n_kinases": int(p.kinase_idx.nunique())}
        parts.append(p.assign(_src=len(parts)))
    if a.mode == "cv":
        fs = [tuple(v["cv_folds_of_sites"]) for v in info.values()]
        if any(len(x) != 1 for x in fs) or len(set(fs)) != len(fs):
            sys.exit(f"CV sources must each cover exactly one distinct CV fold: {fs}")
    else:
        ks = [set(p.kinase_idx.unique()) for p in parts]
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                if ks[i] & ks[j]:
                    sys.exit(f"LKO sources {a.runs[i]} and {a.runs[j]} share {len(ks[i] & ks[j])} kinases")
    P = pd.concat(parts, ignore_index=True)
    key = P.kinase_idx.values.astype(np.int64) * len(S) + P.site_idx.values
    if len(np.unique(key)) != len(key):
        sys.exit("duplicate (kinase, site) across source runs")
    P = P.drop(columns="_src")
    src_man = {}
    for r in a.runs:
        m = ROOT / "results" / "runs" / r / "manifest.json"
        if m.exists():
            mm = json.loads(m.read_text())
            src_man[r] = {k: mm.get(k) for k in ("config", "started", "runtime_s", "git_hash")}
    cfg = {"pooled_from": a.runs, "mode": a.mode, "source_info": info, "source_manifests": src_man,
           "fold_normalise": a.mode == "cv", "bootstrap_B": a.B}
    with RunManifest(a.name, cfg, seed=0, split="cv" if a.mode == "cv" else "lko_T2",
                     notes=a.notes or f"pooled {a.mode} predictions of {len(a.runs)} source runs") as run:
        run.set("source_runs", a.runs)
        run.save_predictions(P)
        data = E.prepare("cv", fold=None)
        res = E.evaluate(P, data, missing="worst", fold_normalise=a.mode == "cv")
        run.save_metrics(res)
        j = res.to_json()
        st = E.strata(res, kinase_set="K_any")
        pk = res.per_kinase
        keep = [c for c in ["kinase_acc", "gene", "n_pos", "n_pos_required", "n_neg_hard", "n_neg_all", "n_missing",
                            "auroc_hard", "auroc_all", "top1_recall", "hit5", "in_K_primary", "in_K_any", "train_bin",
                            "residue_class"] if c in pk]
        j["per_kinase"] = pk[keep].to_dict(orient="records")
        ci = E.bootstrap(res, B=a.B, seed=0) if a.B > 0 else pd.DataFrame()
        ci.to_csv(run.dir / "bootstrap.tsv", sep="\t", index=False)
        res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
        j["a24"] = {"source_run": run.run_id, "canary_scored": bool((P.kinase_idx >= len(index.load_kinases())).any()),
                    "fold_normalise": a.mode == "cv", "missing_policy": "worst",
                    "evaluated": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "n_pred_rows": int(len(P)),
                    "canary_labels_version": 2, "strata_K_any": st.to_dict(orient="records"),
                    "bootstrap": ci.to_dict(orient="records")}
        (run.dir / "metrics_a24.json").write_text(json.dumps(j, indent=1, default=E._json_default))
        e = res.endpoints
        print(json.dumps({"run_id": run.run_id, "M1": e["M1_auroc_hard"], "M1_fold_mean": e.get("M1_fold_mean"),
                          "M2": e["M2_hit5"], "M3": e["M3_top1_recall"], "M4": e["M4_auroc_hard_Kany"],
                          "canary": res.canary, "ci": ci.to_dict(orient="records"),
                          "seconds": round(time.time() - t0)}, indent=1, default=E._json_default))


if __name__ == "__main__":
    main()
