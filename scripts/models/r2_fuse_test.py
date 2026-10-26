"""Late fusion R2 + R3 for the final TEST batch (A13), with the CV-chosen weights (r2_fusion_r3__20261026-191901__s0:
w_R2 = 0.6 for score [M1/M3/M4], 0.5 for score_site [M2]). Each model's logits are standardised over its own TEST
candidate rows (one model = one 'fold'), as on CV where they were standardised within each CV fold.
Usage: r2_fuse_test.py --name r2r3_final_test --r2 <R2 parquet> --r3 <R3 parquet> [--w 0.6 --w-site 0.5]
       [--kind test|t3]   (output file name: test_predictions / test_t3_predictions)"""
import argparse, json
import numpy as np, pandas as pd
from ksp import manifest
from ksp.index import ROOT

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--r2", required=True); ap.add_argument("--r3", required=True)
ap.add_argument("--w", type=float, default=0.6); ap.add_argument("--w-site", type=float, default=0.5)
ap.add_argument("--kind", choices=["test", "t3"], default="test")
a = ap.parse_args()


def load(p):
    d = pd.read_parquet(p)
    if "kinase_idx" not in d:
        raise ValueError(f"{p}: need kinase_idx/site_idx")
    return d[["kinase_idx", "site_idx", "score"]].astype({"kinase_idx": np.int64, "site_idx": np.int64})


with manifest.RunManifest(a.name, config=vars(a), seed=0, split="final_dev", inputs=[a.r2, a.r3],
                          notes="late fusion of final R2 and R3 TEST predictions with CV-chosen weights") as run:
    A, B = load(a.r2), load(a.r3)
    M = A.merge(B, on=["kinase_idx", "site_idx"], how="inner", suffixes=("_r2", "_r3"))
    rep = {"n_r2": len(A), "n_r3": len(B), "n_common": len(M)}
    if not (len(M) == len(A) == len(B)):
        raise ValueError(f"R2/R3 TEST rows differ: {rep}")
    z2 = (M.score_r2 - M.score_r2.mean()) / M.score_r2.std()
    z3 = (M.score_r3 - M.score_r3.mean()) / M.score_r3.std()
    out = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                        "score": (a.w * z2 + (1 - a.w) * z3).astype(np.float32).values,
                        "score_site": (a.w_site * z2 + (1 - a.w_site) * z3).astype(np.float32).values})
    if not np.isfinite(out[["score", "score_site"]].values).all():
        raise RuntimeError("non-finite fused scores")
    run.save_predictions(out, name="test_predictions" if a.kind == "test" else "test_t3_predictions")
    run.set("rows", rep)
    print("RUN", run.run_id, rep)
