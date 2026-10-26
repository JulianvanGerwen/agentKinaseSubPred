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
r3 = a.r3 or json.load(open(ROOT / "results" / "models" / "r3" / "FINAL_OOF.json"))["pooled_cv_run"]


def load(r):
    m = json.load(open(ROOT / "results" / "runs" / r / "manifest.json"))
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
rows, best = [], None
with manifest.RunManifest("r2_fusion_r3", config={"r2": a.r2, "r3": r3, "ws": a.ws}, seed=0, split="cv",
                          notes="late fusion sweep (A13)") as run:
    for w in a.ws:
        P = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                          "score": (w * M.score_r2_z.values + (1 - w) * M.score_r3_z.values).astype(np.float32)})
        res = kev.evaluate(P, data, missing="error", keep_items=False)
        E = res.endpoints
        r = {"w_r2": w, "M1": E["M1_auroc_hard"], "M2": E["M2_hit5"], "M3": E["M3_top1_recall"],
             "M4": E["M4_auroc_hard_Kany"], "mrr": E["sec_mrr"], "neg_canary_hard": res.canary.get("neg_mean_auroc_hard"),
             "pos_canary_hard": res.canary.get("pos_auroc_hard")}
        rows.append(r); print(r, flush=True)
        if best is None or r["M1"] > best[0]["M1"]:
            best = (r, P)
    R = pd.DataFrame(rows)
    R.to_csv(ROOT / "results" / "models" / "r2" / "fusion_r3_sweep.tsv", sep="\t", index=False)
    run.save_predictions(best[1])
    run.save_metrics({"sweep": rows, "best": best[0], "rows": rep})
    train.register({"run_id": run.run_id, "split": "cv", "fold": "all(pooled)", "config_name": "r2_fusion_r3",
                    "config_json": json.dumps({"r2": a.r2, "r3": r3, "ws": a.ws}),
                    "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(best[0]["M1"], 5),
                    "notes": f"late fusion sweep ({len(a.ws)} weights = 1 registry row); best w_r2={best[0]['w_r2']} "
                             f"M2={best[0]['M2']:.4f} M3={best[0]['M3']:.4f}"})
    print("BEST", best[0])
