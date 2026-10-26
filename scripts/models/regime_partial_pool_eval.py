"""Evaluate the partial-pooling lambda grid: --stage fold0 -> fold-0 M1 per lambda (each lambda = one REGISTRY row),
selected lambda = argmax fold-0 M1 (pre-declared rule; ties -> larger lambda). --stage cv5 -> pooled 5-fold for the
selected lambda (registered) and, for transparency only, every lambda (not used for selection, not registered).
Also writes per training-count bin M1 per lambda ('where per-kinase models start to pay off')."""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev
from ksp.index import ROOT
from ksp.models.regime import run as R

ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["fold0", "cv5"], required=True)
a = ap.parse_args()
LAMS = [0.0, 3.0, 10.0, 30.0, 100.0, 300.0, float("inf")]
OUTD = ROOT / "results" / "models" / "regime"
sel_p = OUTD / "ppool_selected_lambda.json"
if a.stage == "fold0":
    P = pd.read_parquet(ROOT / "cache" / "regime" / "ppool" / "fold0.parquet")
    data = kev.prepare("cv", fold=0)
    rows = []
    for lam in LAMS:
        res = kev.evaluate(R.to_eval(P.assign(score=P[f"score_lam{lam:g}"])), data, missing="error")
        s = R.summarise(res)
        pk = res.per_kinase[res.per_kinase.in_K_primary]
        r = {"lambda": lam, **{k: s[k] for k in ("M1", "M2", "M3", "M4", "canary_neg_hard", "canary_pos")}}
        for b, g in pk.groupby("train_bin"):
            r[f"M1_bin_{b}"] = g.auroc_hard.mean()
        rows.append(r)
        R.register({"run_id": f"regime_ppool_lam{lam:g}_fold0", "split": "cv", "fold": "0",
                    "config_name": f"regime_ppool_lam{lam:g}",
                    "config_json": json.dumps({"regime": "partial pooling ridge residual", "lambda": lam, "mu": 10,
                                               "base": R.BASE_RUN}),
                    "selection_metric": "M1_auroc_hard(K_cv)", "value": round(s["M1"], 5),
                    "notes": f"lambda grid fold-0 selection; M2={s['M2']:.4f} M3={s['M3']:.4f}"})
        print(r, flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(OUTD / "ppool_lambda_fold0.tsv", sep="\t", index=False)
    best = T.M1.max()
    lam_sel = float(T[T.M1 >= best - 1e-12]["lambda"].max())
    sel_p.write_text(json.dumps({"lambda": lam_sel, "fold0_M1": best}))
    print("selected lambda", lam_sel)
else:
    lam_sel = json.loads(sel_p.read_text())["lambda"]
    P = pd.concat([pd.read_parquet(ROOT / "cache" / "regime" / "ppool" / f"fold{f}.parquet") for f in range(5)],
                  ignore_index=True)
    data = kev.prepare("cv", fold=None)
    rows = []
    for lam in LAMS:
        Pl = P.assign(score=P[f"score_lam{lam:g}"])
        if lam == lam_sel:
            po = R.pool(f"regime_ppool_lam{lam:g}_cv5", [], {"regime": "partial pooling", "lambda": lam, "mu": 10,
                                                             "base": R.BASE_RUN, "seed": 0},
                        notes="lambda selected on fold-0 M1", preds=R.to_eval(Pl))
            res, s = po["res"], po["summary"]
        else:
            res = kev.evaluate(R.to_eval(Pl), data, missing="error")
            s = R.summarise(res)
        pk = res.per_kinase[res.per_kinase.in_K_primary]
        r = {"lambda": lam, "selected": lam == lam_sel, **{k: s[k] for k in s if not k.startswith("pair_")}}
        for b, g in pk.groupby("train_bin"):
            r[f"M1Kcv_bin_{b}"] = g.auroc_hard.mean()
        rows.append(r)
        print(r, flush=True)
    pd.DataFrame(rows).to_csv(OUTD / "ppool_lambda_cv5.tsv", sep="\t", index=False)
