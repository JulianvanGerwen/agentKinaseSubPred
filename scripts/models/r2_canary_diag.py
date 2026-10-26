"""Negative-canary diagnostic: is above-0.5 canary AUROC explained by the canary sampler's protein site-density bias?
AUROC (E-all) overall and stratified by deciles of #U_ms sites on the site's protein (within-stratum AUROC,
weighted by positives). Usage: r2_canary_diag.py --runs <run_id> ... [--fold 0]"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, canary
from ksp.index import ROOT

ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+"); ap.add_argument("--fold", type=int, default=None)
a = ap.parse_args()
data = kev.prepare("cv", fold=a.fold)
S = data.sites
nps = S[S.in_U_ms].groupby("acc").size()
dens = S.acc.map(nps).fillna(0).values
K = data.kinases
def auc(p, n):
    if len(p) == 0 or len(n) == 0: return np.nan
    r = pd.Series(np.r_[p, n]).rank().values
    return (r[:len(p)].sum() - len(p) * (len(p) + 1) / 2) / (len(p) * len(n))
rows = []
for rid in a.runs:
    m = json.load(open(ROOT / "results" / "runs" / rid / "manifest.json"))
    P = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
    for kacc in canary.NEG_IDS + [canary.POS_ID]:
        i = int(np.where(K.uniprot_acc.values == kacc)[0][0])
        pos = data.positives.get(kacc, np.zeros(0, int))
        x = P[P.kinase_idx.values == i]
        sc = np.full(len(S), np.nan); sc[x.site_idx.values] = x.score.values
        cand = data.universe & kev._res_ok(S.res.values, K.residue_class.values[i])
        excl = np.zeros(len(S), bool); excl[data.masks.excluded(kacc)] = True
        isp = np.zeros(len(S), bool); isp[pos] = True
        neg = np.flatnonzero(cand & ~excl & ~isp)
        pos = pos[np.isfinite(sc[pos])]; neg = neg[np.isfinite(sc[neg])]
        q = np.quantile(dens[neg], np.linspace(0, 1, 11)[1:-1])
        bp, bn = np.digitize(dens[pos], q), np.digitize(dens[neg], q)
        strat = [(auc(sc[pos[bp == b]], sc[neg[bn == b]]), (bp == b).sum()) for b in range(10)]
        strat = [(v, w) for v, w in strat if np.isfinite(v)]
        rows.append({"run": rid.split("__")[0], "canary": kacc, "n_pos": len(pos), "auroc_all": auc(sc[pos], sc[neg]),
                     "auroc_all_density_stratified": np.average([v for v, w in strat], weights=[w for v, w in strat]) if strat else np.nan,
                     "proxy_auroc_all": auc(-dens[pos], -dens[neg])})
R = pd.DataFrame(rows)
neg = R[R.canary.str.startswith("CANARY_NEG")]
print(neg.groupby("run")[["auroc_all", "auroc_all_density_stratified", "proxy_auroc_all"]].mean().round(3))
print(R[R.canary == canary.POS_ID].round(3).to_string())
out = ROOT / "results" / "models" / "r2" / "canary_density_diagnostic.tsv"
R.to_csv(out, sep="\t", index=False, mode="a", header=not out.exists())
