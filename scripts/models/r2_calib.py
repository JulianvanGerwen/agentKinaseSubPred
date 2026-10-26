"""R2 tiered calibration (A16) + site-level assignment probability and coverage (A17) on a pooled 5-fold OOF run.
Everything is cross-fitted over CV folds (calibrator / temperature fitted on folds != f, applied to fold f).
Usage: r2_calib.py --run <pooled run_id> [--cantley <B1 run_id>]
Outputs: results/models/r2/calib/<run_id>/{calibrators.tsv, reliability_by_bin.tsv, assignment.json, ...}"""
import argparse, json, time
import numpy as np, pandas as pd
from ksp import eval as kev, splits, canary, index
from ksp.index import ROOT
from ksp.models.r2 import calib

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--cantley", default=None, help="B1 run id for the Cantley-only coverage baseline")
ap.add_argument("--min-pos", type=int, default=30)
a = ap.parse_args()
t0 = time.time()
out = ROOT / "results" / "models" / "r2" / "calib" / a.run
out.mkdir(parents=True, exist_ok=True)
m = json.load(open(ROOT / "results" / "runs" / a.run / "manifest.json"))
P = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
data = kev.prepare("cv", fold=None)
S = data.sites
fold_of_site = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values
K = data.kinases
real = ~K.uniprot_acc.str.startswith("CANARY_").values
# ---------------------------------------------------------------- calibration table + cross-fitted calibrators
T = calib.table_from_eval(P, data)
T["fold"] = fold_of_site[T.site_idx.values]
T["p_cf"] = np.nan
Cs = {}
for f in range(5):
    C = calib.fit_calibrator(T[T.fold != f], K, a.min_pos)
    Cs[f] = C
    mk = T.fold.values == f
    T.loc[mk, "p_cf"] = 1 / (1 + np.exp(-calib.apply(C, T.kinase_idx.values[mk], T.score.values[mk])))
C_all = calib.fit_calibrator(T, K, a.min_pos)
C_all.assign(gene=K.gene.values[C_all.index]).to_csv(out / "calibrators.tsv", sep="\t")
tc = data.train_counts
T["n_train"] = pd.Series(K.uniprot_acc.values[T.kinase_idx.values]).map(tc).fillna(0).values
T["bin"] = pd.cut(T.n_train, kev.TRAIN_BINS, labels=kev.TRAIN_BIN_LABELS).astype(str)
rel_rows, ece = [], {}
for b, g in list(T.groupby("bin")) + [("all", T)]:
    R, e = calib.reliability(g.p_cf.values, g.y.values, g.w.values)
    R.insert(0, "train_bin", b)
    rel_rows.append(R)
    ece[b] = e
pd.concat(rel_rows).to_csv(out / "reliability_by_bin.tsv", sep="\t", index=False)
tiers = C_all.tier.value_counts().to_dict()
# expected vs observed positives (calibration-in-the-large) per bin, out-of-sample
itl = T.groupby("bin").apply(lambda g: pd.Series({"expected": float((g.p_cf * g.w).sum()), "observed": int(g.y.sum())}))
itl.to_csv(out / "calibration_in_the_large.tsv", sep="\t")
print("ECE by train bin (cross-fitted):", ece, "\ntiers:", tiers, "\n", itl, flush=True)

# ---------------------------------------------------------------- assignment probability (A17)
# calibrated logits for every OOF row of real kinases (cross-fitted calibrator by the site's fold)
P = P[real[P.kinase_idx.values]].reset_index(drop=True)
pf = fold_of_site[P.site_idx.values]
lc = np.empty(len(P), np.float64)
for f in range(5):
    mk = pf == f
    lc[mk] = calib.apply(Cs[f], P.kinase_idx.values[mk], P.score.values[mk])
o = np.lexsort((P.kinase_idx.values, P.site_idx.values))
s_sorted = P.site_idx.values[o]
k_sorted = P.kinase_idx.values[o]
l_sorted = lc[o]
# true kinases per site (DEV eval positives) -> is_true flags
pos_pairs = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
kpos = {acc: i for i, acc in enumerate(K.uniprot_acc)}
key_true = set(zip(pos_pairs.kinase_acc.map(kpos).values, pos_pairs.site_idx.values))
is_true = np.fromiter(((k, s) in key_true for k, s in zip(k_sorted, s_sorted)), bool, len(k_sorted))
sfold = fold_of_site[s_sorted]
Tf, p_assign = {}, np.empty(len(s_sorted))
for f in range(5):
    tr = sfold != f
    Tf[f], nll = calib.fit_temperature(s_sorted[tr], l_sorted[tr], is_true[tr])
    te = ~tr
    p_assign[te] = calib.site_softmax(s_sorted[te], l_sorted[te], Tf[f])
b = np.flatnonzero(np.r_[True, s_sorted[1:] != s_sorted[:-1]])
cnt = np.diff(np.r_[b, len(s_sorted)])
top_p = np.maximum.reduceat(p_assign, b)
arg = b + np.array([np.argmax(p_assign[x:x + c]) for x, c in zip(b, cnt)])
site_ids = s_sorted[b]
labelled = np.add.reduceat(is_true.astype(int), b) > 0
top_true = is_true[arg]
rel = []
for lo, hi in [(0, .2), (.2, .5), (.5, .8), (.8, 1.01)]:
    mk = labelled & (top_p >= lo) & (top_p < hi)
    rel.append({"p_top_bin": f"[{lo},{hi})", "n_labelled_sites": int(mk.sum()),
                "mean_p_top": float(top_p[mk].mean()) if mk.any() else None,
                "frac_top_is_true": float(top_true[mk].mean()) if mk.any() else None})
taus = (0.5, 0.8)
cov = {f"coverage@{t}": int((top_p >= t).sum()) for t in taus}
cov.update({f"coverage@{t}_unlabelled": int(((top_p >= t) & ~labelled).sum()) for t in taus})
cov["n_sites"] = int(len(site_ids)); cov["n_labelled_sites"] = int(labelled.sum())
# random-score baseline: same calibrators/temperature machinery on random scores is degenerate -> use uniform
# random logits with the fitted temperature distribution (conservative): p_top under random N(0,1) logits.
rng = np.random.default_rng(0)
rl = rng.standard_normal(len(s_sorted))
p_rand = calib.site_softmax(s_sorted, rl, float(np.median(list(Tf.values()))))
top_rand = np.maximum.reduceat(p_rand, b)
cov_rand = {f"coverage@{t}": int((top_rand >= t).sum()) for t in taus}
res = {"temperature_by_fold": Tf, "assignment_reliability_labelled_sites": rel, "coverage_DEV_Ums_OOF": cov,
       "coverage_random_scores": cov_rand, "top1_accuracy_labelled": float(top_true[labelled].mean()),
       "ece_by_bin": ece, "tiers": tiers, "runtime_s": round(time.time() - t0, 1),
       "note": "PU: calibrated probabilities estimate P(curated human-native PSP positive | score), a lower bound; "
               "coverage is conditional on the site being a target of a modelled kinase (A17)."}
json.dump(res, open(out / "assignment.json", "w"), indent=1, default=float)
pd.DataFrame({"site_idx": site_ids, "top_kinase_idx": k_sorted[arg], "p_top": top_p, "labelled": labelled,
              "top_is_true": top_true}).to_parquet(ROOT / "cache" / "r2" / f"assign_{a.run}.parquet")
# M2 with the assignment score (A20): score_site = cross-fitted calibrated logit (monotone in p(k|s) within site)
P2 = pd.DataFrame({"kinase_idx": k_sorted, "site_idx": s_sorted, "score": P.score.values[o], "score_site": l_sorted})
Pc = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
Pc = Pc[~real[Pc.kinase_idx.values]]
Pc["score_site"] = Pc.score
P2 = pd.concat([P2, Pc], ignore_index=True)
r2 = kev.evaluate(P2, data, missing="error")
res["endpoints_with_assignment_score_site"] = r2.endpoints
json.dump(res, open(out / "assignment.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
