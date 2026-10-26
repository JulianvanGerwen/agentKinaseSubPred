"""LKO-T2 evaluation of the late-fusion final model: R2 LKO (r2_final_lko) + R3 LKO (FINAL_OOF lko_T2_runs), same
standardisation (per CV fold of the site) and the CV-chosen weights (w 0.6 score, 0.5 score_site). Evaluated like the
LKO pools: prepare('cv'), kinases = held-out kinases with positives, missing='drop'.
Usage: r2_fusion_lko.py --r2 <r2_final_lko run>"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, manifest, splits, canary
from ksp.index import ROOT
from ksp.models.r2 import train
ap = argparse.ArgumentParser(); ap.add_argument("--r2", required=True); ap.add_argument("--w", type=float, default=0.6)
ap.add_argument("--w-site", type=float, default=0.5); a = ap.parse_args()
R3F = json.load(open(ROOT / "results" / "models" / "r3" / "FINAL_OOF.json"))
def load(r):
    m = json.load(open(ROOT / "results" / "runs" / r / "manifest.json"))
    return pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"], columns=["kinase_idx", "site_idx", "score"]).astype(
        {"kinase_idx": np.int64, "site_idx": np.int64})
A = load(a.r2); B = pd.concat([load(r) for r in R3F["lko_T2_runs"]], ignore_index=True)
M = A.merge(B, on=["kinase_idx", "site_idx"], suffixes=("_r2", "_r3"))
rep = {"n_r2": len(A), "n_r3": len(B), "n_common": len(M),
       "r2_duplicates": int(A.duplicated(["kinase_idx", "site_idx"]).sum()),
       "r3_duplicates": int(B.duplicated(["kinase_idx", "site_idx"]).sum())}
print(rep, flush=True)
if not (len(M) == len(A) == len(B)) or rep["r2_duplicates"] or rep["r3_duplicates"]:      # audit r2 (a)
    raise ValueError(f"R2/R3 LKO rows are not aligned 1:1: {rep}")
data = kev.prepare("cv", fold=None)
fold = data.sites.acc.map(splits.load_cv_folds()).values[M.site_idx.values]
for c in ["score_r2", "score_r3"]:
    z = M[c].values.astype(np.float64)
    for f in range(5):
        mk = fold == f; z[mk] = (z[mk] - z[mk].mean()) / z[mk].std()
    M[c + "_z"] = z
K = canary.kinases_with_canaries()
with manifest.RunManifest("r2r3_fusion_lko", config=vars(a) | {"r3_lko_runs": R3F["lko_T2_runs"]}, seed=0, split="lko_T2",
                          notes="LKO T2 of the late-fusion final model") as run:
    P = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                      "score": (a.w * M.score_r2_z + (1 - a.w) * M.score_r3_z).astype(np.float32).values,
                      "score_site": (a.w_site * M.score_r2_z + (1 - a.w_site) * M.score_r3_z).astype(np.float32).values})
    run.set("row_alignment", rep | {"aligned_1to1": True})
    run.save_predictions(P)
    held = sorted(set(K.uniprot_acc.values[np.unique(P.kinase_idx.values)]) & set(data.positives))
    res = kev.evaluate(P, data, kinases=held, missing="drop")
    run.save_metrics(res); res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
    E = res.endpoints
    print(f"[fusion LKO] M1={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} M4={E['M4_auroc_hard_Kany']:.4f} rows={rep}")
    train.register({"run_id": run.run_id, "split": "lko_T2", "fold": "all(pooled)", "config_name": "r2r3_fusion_lko",
                    "config_json": json.dumps(vars(a)), "selection_metric": "M1_auroc_hard(LKO T2)", "value": round(E["M1_auroc_hard"], 5),
                    "notes": f"LKO of final fusion (no selection); M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f}"})
