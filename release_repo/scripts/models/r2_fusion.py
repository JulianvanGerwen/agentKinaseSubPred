"""Late fusion R2 + R3 (DESIGN §16 A13): blended logit = w * z(R2) + (1 - w) * z(R3), where z() standardises each
model's raw logits by their SD within CV fold (one scale per fold, keeps cross-kinase comparability for M2).
w chosen on pooled CV M1 (ksp.eval, A24 rank-normalised); M2 uses the blended logit (cross-kinase).
Usage: r2_fusion.py --r2 <pooled R2 run> [--r3 <pooled R3 run>] [--ws 0 0.2 0.4 0.5 0.6 0.8 1]"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, manifest, splits
from ksp.index import ROOT
from ksp.models.r2 import train

ap = argparse.ArgumentParser()
ap.add_argument("--r2", required=True); ap.add_argument("--r3", default=None)
ap.add_argument("--ws", nargs="+", type=float, default=[0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
a = ap.parse_args()
R3F = json.load(open(ROOT / "results" / "models" / "r3" / "FINAL_OOF.json"))
r3 = a.r3 or R3F["pooled_cv_run"]


def load(r):
    m = json.load(open(ROOT / "results" / "runs" / r / "manifest.json"))
    if "predictions" not in m["outputs"]:          # R3 pooled run stores no predictions: use its 5 fold runs
        return pd.concat([load(x) for x in R3F["cv_fold_runs"]], ignore_index=True)
    p = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"], columns=["kinase_idx", "site_idx", "score"])
    return p.astype({"kinase_idx": np.int64, "site_idx": np.int64})


A, B = load(a.r2), load(r3)
M = A.merge(B, on=["kinase_idx", "site_idx"], how="inner", suffixes=("_r2", "_r3"))
rep = {"n_r2": len(A), "n_r3": len(B), "n_common": len(M)}
if len(M) != len(A):
    raise ValueError(f"R2/R3 rows differ: {rep}")
data = kev.prepare("cv", fold=None)
fold = data.sites.acc.map(splits.load_cv_folds()).values[M.site_idx.values]
for c in ["score_r2", "score_r3"]:
    z = M[c].values.astype(np.float64)
    for f in range(5):
        mk = fold == f
        z[mk] = (z[mk] - z[mk].mean()) / z[mk].std()
    M[c + "_z"] = z
rows = []
with manifest.RunManifest("r2_fusion_r3", config={"r2": a.r2, "r3": r3, "ws": a.ws}, seed=0, split="cv",
                          notes="late fusion sweep (A13); score = blend with w chosen on M1, score_site = blend with w "
                                "chosen on M2 (A20: each component selected on its own CV endpoint)") as run:
    def blend(w):
        return (w * M.score_r2_z.values + (1 - w) * M.score_r3_z.values).astype(np.float32)
    for w in a.ws:
        P = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values, "score": blend(w)})
        res = kev.evaluate(P, data, missing="error", keep_items=False)
        E = res.endpoints
        r = {"w_r2": w, "M1": E["M1_auroc_hard"], "M2": E["M2_hit5"], "M3": E["M3_top1_recall"],
             "M4": E["M4_auroc_hard_Kany"], "mrr": E["sec_mrr"], "M1_fold_mean": E.get("M1_fold_mean"),
             "neg_canary_hard": res.canary.get("neg_mean_auroc_hard"), "pos_canary_hard": res.canary.get("pos_auroc_hard")}
        rows.append(r); print(r, flush=True)
    R = pd.DataFrame(rows)
    tag = a.r2.split("__")[0]
    R.to_csv(ROOT / "results" / "models" / "r2" / f"fusion_r3_sweep__{tag}.tsv", sep="\t", index=False)
    w1 = float(R.loc[R.M1.idxmax(), "w_r2"]); w2 = float(R.loc[R.M2.idxmax(), "w_r2"])
    P = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values, "score": blend(w1),
                      "score_site": blend(w2)})
    run.save_predictions(P)
    res = kev.evaluate(P, data, missing="error")
    run.save_metrics(res)
    res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
    kev.bootstrap(res, B=1000, seed=0).to_csv(run.dir / "bootstrap.tsv", sep="\t", index=False)
    kev.strata(res).to_csv(run.dir / "strata.tsv", sep="\t", index=False)
    run.set("weights", {"w_r2_score(M1)": w1, "w_r2_score_site(M2)": w2, "rows": rep})
    E = res.endpoints
    train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)", "config_name": "r2_fusion_r3",
                    "config_json": json.dumps({"r2": a.r2, "r3": r3, "ws": a.ws, "w_M1": w1, "w_M2": w2}),
                    "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(E["M1_auroc_hard"], 5),
                    "notes": f"late fusion sweep ({len(a.ws)} weights = 1 registry row); w_M1={w1} w_M2={w2}; "
                             f"M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} fold-mean M1={E.get('M1_fold_mean')}; "
                             f"canary={res.canary.get('neg_mean_auroc_hard'):.3f}/{res.canary.get('pos_auroc_hard'):.3f}"})
    print("FINAL", {"w_M1": w1, "w_M2": w2, **{k: E[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany"]},
                    "M1_fold_mean": E.get("M1_fold_mean"), "M1_per_fold": E.get("M1_per_fold"), "canary": res.canary}, flush=True)
    print("RUN", run.run_id)
