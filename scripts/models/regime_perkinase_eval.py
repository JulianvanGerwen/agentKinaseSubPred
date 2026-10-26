"""Combine per-kinase predictions with the global R2 OOF fallback and evaluate (fold 0 or pooled 5-fold).
Variants written: pk_lgb, pk_lr (score = per-kinase model where available, else global logit), each with
score_site = qmap onto the global scores of the same (fold, kinase) (label-free; for M2), plus raw-score M2.
Usage: regime_perkinase_eval.py --folds 0 | --folds 0 1 2 3 4"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, splits
from ksp.index import ROOT
from ksp.models.regime import perkinase as PK, run as R

ap = argparse.ArgumentParser()
ap.add_argument("--folds", nargs="+", type=int, default=[0])
ap.add_argument("--n-min", type=int, default=20)
a = ap.parse_args()
base = json.loads((ROOT / "results" / "runs" / R.BASE_RUN / "manifest.json").read_text())
G = pd.read_parquet(ROOT / base["outputs"]["predictions"]["path"])
from ksp import index
S = index.load_sites()
fold_of = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values
G["fold"] = fold_of[G.site_idx.values]
G = G[G.fold.isin(a.folds)].reset_index(drop=True)
PKs = pd.concat([pd.read_parquet(ROOT / "cache" / "regime" / "perkinase" / f"fold{f}.parquet") for f in a.folds])
M = G.merge(PKs, on=["kinase_idx", "site_idx"], how="left")
has = M.score_lgb.notna().values
fold_label = "0" if a.folds == [0] else "all(pooled)"
data = kev.prepare("cv", fold=None if len(a.folds) == 5 else a.folds[0])
out = {}
for mdl in ("lgb", "lr"):
    sc = np.where(has, M[f"score_{mdl}"].values, M.score.values)
    ss = sc.copy()
    hidx = np.flatnonzero(has)
    for (k, f), ix in M[has].groupby(["kinase_idx", "fold"]).indices.items():
        rows = hidx[ix]
        ss[rows] = PK.qmap(sc[rows], M.score.values[rows])
    for m2 in ("qmap", "raw"):
        P = M[["kinase_idx", "site_idx"]].assign(score=sc.astype(np.float64), score_site=(ss if m2 == "qmap" else sc))
        name = f"regime_perkinase_{mdl}_N{a.n_min}_{m2}" + ("_fold0" if fold_label == "0" else "_cv5")
        cfg = {"regime": "per-kinase", "model": mdl, "n_min": a.n_min, "fallback": "global R2 " + R.BASE_RUN,
               "m2_score": m2, "seed": 0}
        if len(a.folds) == 5:
            done = sorted((ROOT / "results" / "runs").glob(f"{name}__*/summary.json"))
            if done:          # already pooled earlier (process was OOM-killed after it): re-evaluate, don't re-register
                res = kev.evaluate(R.to_eval(P), data, missing="error")
                s = R.summarise(res)
                pk = res.per_kinase.set_index("kinase_acc")
                out[(mdl, m2)] = (s, pk)
                print(name, "(re-evaluated; pooled run exists)", flush=True)
                continue
            po = R.pool(name, [], cfg, notes=f"per-kinase {mdl} N>={a.n_min}, else global; M2 score {m2}",
                        preds=R.to_eval(P), register_row=(m2 == "qmap"))
            res, s = po["res"], po["summary"]
        else:
            res = kev.evaluate(R.to_eval(P), data, missing="error")
            s = R.summarise(res)
            if m2 == "qmap":
                R.register({"run_id": name, "split": "cv", "fold": "0", "config_name": name,
                            "config_json": json.dumps(cfg), "selection_metric": "M1_auroc_hard(K_cv)",
                            "value": round(s["M1"], 5), "notes": f"fold-0 pilot; M2={s['M2']:.4f}"})
        pk = res.per_kinase.set_index("kinase_acc")
        out[(mdl, m2)] = (s, pk)
        print(name, json.dumps({k: v for k, v in s.items() if not k.startswith("pair_")}, default=str), flush=True)
# subset view: kinases with an own model in every evaluated fold vs the global run
from ksp import canary
K = canary.kinases_with_canaries()
own = PKs.groupby("kinase_idx").site_idx.size().index
nf = PKs.drop_duplicates("kinase_idx").shape[0]
cnt = pd.Series(0, index=range(len(K)))
for f in a.folds:
    t = pd.read_csv(ROOT / "cache" / "regime" / "perkinase" / f"fold{f}_models.tsv", sep="\t")
    cnt[t.kinase_idx.values] += 1
sub = set(K.uniprot_acc.values[cnt.values == len(a.folds)])
Gres = kev.evaluate(R.to_eval(G[["kinase_idx", "site_idx", "score"]]), data, missing="error").per_kinase.set_index("kinase_acc")
kp = [k for k in data.k_primary if k in sub and k in Gres.index]
rows = [{"model": "global_R2", "subset": f">=N in all folds (n={len(kp)})", "M1": Gres.loc[kp, "auroc_hard"].mean(),
         "M3": Gres.loc[kp, "top1_recall"].mean(), "M2": Gres.loc[kp, "hit5"].mean()}]
rest = [k for k in data.k_primary if k not in sub and k in Gres.index]
rows.append({"model": "global_R2", "subset": f"<N in some fold (n={len(rest)})", "M1": Gres.loc[rest, "auroc_hard"].mean(),
             "M3": Gres.loc[rest, "top1_recall"].mean(), "M2": Gres.loc[rest, "hit5"].mean()})
for (mdl, m2), (s, pk) in out.items():
    rows.append({"model": f"pk_{mdl}_{m2}", "subset": f">=N in all folds (n={len(kp)})", "M1": pk.loc[kp, "auroc_hard"].mean(),
                 "M3": pk.loc[kp, "top1_recall"].mean(), "M2": pk.loc[kp, "hit5"].mean()})
    rows.append({"model": f"pk_{mdl}_{m2}", "subset": "all K_cv", "M1": s["M1"], "M3": s["M3"], "M2": s["M2"]})
T = pd.DataFrame(rows)
tag = "fold0" if fold_label == "0" else "cv5"
T.to_csv(ROOT / "results" / "models" / "regime" / f"perkinase_subset_{tag}.tsv", sep="\t", index=False)
print(T.to_string())
