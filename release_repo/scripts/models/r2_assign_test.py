"""A17 assignment probability for the final TEST predictions (PI 22:12; audit round 2 fix R2-M1).

Pipeline (every map is fitted on DEV OOF only, never on TEST):
  1. calibrated logit  = a_k + b_k * score            (score = fused M1 score; tiered calibrators fitted on ALL DEV
                                                         OOF of the CV fusion run: results/models/r2/calib/<run>__calibrated/calibrators.tsv)
  2. p(k|s)            = softmax over real kinases of matching class at site s, temperature T = mean of the 5
                         cross-fitted CV temperatures (0.921-0.932)
  3. isotonic map      : p_top -> P(top kinase is a curated kinase), fitted SEPARATELY for ST and Y sites on all
                         labelled DEV sites' OOF p_top (cache/r2/assign_<run>__calibrated.parquet; auditor's fix,
                         docs/audit/audit_cv_r2_scripts/calib_check3.py, here without fold cross-fitting = "map fitted on
                         all DEV OOF" for TEST)
  4. assignment_prob   : top kinase -> iso(p_top); other kinases -> p(k|s) * (1 - iso) / (1 - p_top)
                         (probabilities stay in [0, 1] and sum to 1 per site; ordering within a site unchanged)
Writes assignment_prob into a copy of the fused predictions -> cache/runs/<new run>/test_predictions.parquet.
Usage: r2_assign_test.py --name r2r3_final_test_assign --pred <fused test_predictions.parquet>"""
import argparse, json
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
from ksp import index, manifest
from ksp.index import ROOT
from ksp.models.r2 import calib

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--pred", required=True)
ap.add_argument("--cv-run", default="r2_fusion_r3__20261026-191901__s0")
a = ap.parse_args()
cdir = ROOT / "results" / "models" / "r2" / "calib" / f"{a.cv_run}__calibrated"
C = pd.read_csv(cdir / "calibrators.tsv", sep="\t", index_col=0)
T = float(np.mean(list(json.load(open(cdir / "assignment.json"))["temperature_by_fold"].values())))
A = pd.read_parquet(ROOT / "cache" / "r2" / f"assign_{a.cv_run}__calibrated.parquet")
S = index.load_sites()
A["res"] = np.where(S.res.values[A.site_idx.values] == "Y", "Y", "ST")
iso = {}
for r in ("ST", "Y"):
    L = A[A.labelled & (A.res == r)]
    iso[r] = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(L.p_top.values, L.top_is_true.values.astype(float))
with manifest.RunManifest(a.name, config={"pred": a.pred, "cv_run": a.cv_run, "temperature": T}, seed=0, split="final_dev",
                          inputs=[a.pred, cdir / "calibrators.tsv"], notes="A17 assignment_prob for TEST (isotonic ST/Y)") as run:
    P = pd.read_parquet(a.pred)
    if (P.kinase_idx >= len(C)).any():
        raise ValueError("kinase outside calibrators")
    lc = C.a.values[P.kinase_idx.values] + C.b.values[P.kinase_idx.values] * P.score.values.astype(np.float64)
    o = np.lexsort((P.kinase_idx.values, P.site_idx.values))
    s = P.site_idx.values[o]
    p = calib.site_softmax(s, lc[o], T)
    b = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
    cnt = np.diff(np.r_[b, len(s)])
    ptop = np.maximum.reduceat(p, b)
    res = np.where(S.res.values[s[b]] == "Y", "Y", "ST")
    q = np.empty(len(b))
    for r in ("ST", "Y"):
        m = res == r
        q[m] = iso[r].predict(ptop[m])
    ptop_r, q_r = np.repeat(ptop, cnt), np.repeat(q, cnt)
    is_top = p >= ptop_r                                   # ties: all tied tops get iso value (rare)
    with np.errstate(divide="ignore", invalid="ignore"):
        other = np.where(ptop_r < 1, p * (1 - q_r) / (1 - ptop_r), 0.0)
    ap_ = np.where(is_top, q_r, other)
    ap_ = np.clip(ap_, 0, 1)
    out = P.copy()
    out["assignment_prob"] = np.empty(len(P), np.float32)
    out.loc[o, "assignment_prob"] = ap_.astype(np.float32)   # o indexes rows of P (RangeIndex)
    if not np.isfinite(out.assignment_prob).all():
        raise RuntimeError("non-finite assignment_prob")
    run.save_predictions(out, name="test_predictions")
    summ = {"temperature": T, "n_sites": int(len(b)), "coverage@0.5": int((q >= 0.5).sum()), "coverage@0.8": int((q >= 0.8).sum()),
            "coverage@0.5_by_res": {r: int(((q >= 0.5) & (res == r)).sum()) for r in ("ST", "Y")},
            "iso_max": {r: float(iso[r].predict([1.0])[0]) for r in ("ST", "Y")}}
    run.set("summary", summ)
    print("RUN", run.run_id, json.dumps(summ))
