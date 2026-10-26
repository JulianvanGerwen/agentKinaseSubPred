"""T2 diagnostic: AUROC split by whether a site is a training positive of some trained (non-held) kinase ('seen')."""
import json, sys
sys.path.insert(0, "/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/audit2")
from common import *
from ksp import eval as kev, labels, canary
data = kev.prepare("cv", fold=None)
K = canary.kinases_with_canaries(); kpos = {a: i for i, a in enumerate(K.uniprot_acc)}
kf = splits.load_kinase_folds(); fold_of_k = dict(zip(kf.kinase_acc, kf.lko_fold))
S = data.sites; spos = pd.Series(np.arange(len(S)), index=S.site_id.values)
L2 = labels.load_dev("L2"); L2 = L2[L2.site_id.isin(spos.index)]
L2["si"] = spos[L2.site_id].values; L2["kf"] = L2.kinase_acc.map(fold_of_k)
nS = len(S)
seen = {}
for j in range(5):
    m = np.zeros(nS, bool); m[L2.si.values[(L2.kf != j).values]] = True; seen[j] = m
def per_kinase(P):
    P = P.sort_values("kinase_idx"); kb = np.searchsorted(P.kinase_idx.values, np.arange(len(K) + 1))
    rows = []
    for k in data.k_primary:
        if k not in fold_of_k: continue
        i = kpos[k]; j = fold_of_k[k]
        sc = np.full(nS, np.nan); sl = slice(kb[i], kb[i + 1]); sc[P.site_idx.values[sl]] = P.score.values[sl]
        pos = data.positives.get(k, np.zeros(0, int))
        cand = data.universe & kev._res_ok(S.res.values, K.residue_class.values[i])
        ex = np.zeros(nS, bool); ex[data.masks.excluded(k)] = True
        hard = np.zeros(nS, bool); hard[data.masks.hard_candidates(k)] = True
        pm = np.zeros(nS, bool); pm[pos] = True
        neg = np.where(cand & ~ex & ~pm & hard & ~np.isnan(sc))[0]
        pos = pos[~np.isnan(sc[pos])]
        sj = seen[j]
        r = {"k": k, "all": kev.auroc(sc[pos], sc[neg]),
             "pos_seen_vs_neg_seen": kev.auroc(sc[pos[sj[pos]]], sc[neg[sj[neg]]]),
             "pos_unseen_vs_neg_seen": kev.auroc(sc[pos[~sj[pos]]], sc[neg[sj[neg]]]),
             "pos_seen_vs_all": kev.auroc(sc[pos[sj[pos]]], sc[neg]), "pos_unseen_vs_all": kev.auroc(sc[pos[~sj[pos]]], sc[neg]),
             "frac_pos_seen": sj[pos].mean() if len(pos) else np.nan, "frac_neg_seen": sj[neg].mean() if len(neg) else np.nan}
        rows.append(r)
    return pd.DataFrame(rows)
out = {}
L = load(C + LKO + "/predictions.parquet")
d = per_kinase(L); d.to_csv("out/lko_strat_T2.tsv", sep="\t", index=False); out["T2"] = d.drop(columns="k").mean().to_dict(); del L
sf = site_fold()
E = load(C + ENS + "/predictions.parquet")
# CV control: rank-normalise within (fold, kinase) as A24 so pooled per-kinase AUROC is not offset-distorted
E["score"] = E.groupby([E.kinase_idx.values, sf[E.site_idx.values]]).score.rank(pct=True).values
d = per_kinase(E); d.to_csv("out/lko_strat_CV.tsv", sep="\t", index=False); out["CV_control"] = d.drop(columns="k").mean().to_dict()
json.dump(out, open("out/lko_strat.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
