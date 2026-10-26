"""Leave-one-fold-out selection of fusion weights from per-fold sweeps, then pooled A24 evaluation of the nested blend."""
import json, sys
sys.path.insert(0, "/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/audit2")
from common import *
from ksp import eval as kev
SW = pd.concat([pd.read_csv(f"out/sweep_fold{f}.tsv", sep="\t") for f in range(5)])
sel = {}
for f in range(5):
    o = SW[SW.fold != f].groupby("w")[["M1", "M2", "M3"]].mean()
    sel[f] = {"w_M1": float(o.M1.idxmax()), "w_M2": float(o.M2.idxmax())}
glob = SW.groupby("w")[["M1", "M2", "M3"]].mean()
res = {"selected_lofo": sel, "global_argmax_foldmean": {"w_M1": float(glob.M1.idxmax()), "w_M2": float(glob.M2.idxmax())},
       "foldmean_by_w": glob.round(5).to_dict("index")}
# per-fold nested estimate (fold-mean metrics)
nm1 = np.mean([SW[(SW.fold == f) & (SW.w == sel[f]["w_M1"])].M1.iloc[0] for f in range(5)])
nm2 = np.mean([SW[(SW.fold == f) & (SW.w == sel[f]["w_M2"])].M2.iloc[0] for f in range(5)])
nm3 = np.mean([SW[(SW.fold == f) & (SW.w == sel[f]["w_M1"])].M3.iloc[0] for f in range(5)])
res["nested_foldmean"] = {"M1": nm1, "M2": nm2, "M3": nm3}
res["inselect_foldmean_at_0.6_0.5"] = {"M1": float(glob.loc[0.6, "M1"]), "M2": float(glob.loc[0.5, "M2"]), "M3": float(glob.loc[0.6, "M3"])}
res["oracle_foldmean"] = {"M1": float(glob.M1.max()), "M2": float(glob.M2.max())}
print(json.dumps(res, indent=1)); json.dump(res, open("out/nested_sel.json", "w"), indent=1)
# pooled nested blend
sf = site_fold(); ens = load(C + ENS + "/predictions.parquet")
parts = []
for f in range(5):
    M, _ = blend_fold(f, sf, ens)
    a, b = sel[f]["w_M1"], sel[f]["w_M2"]
    parts.append(pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                               "score": a * M.score_r2_z.values + (1 - a) * M.score_r3_z.values,
                               "score_site": b * M.score_r2_z.values + (1 - b) * M.score_r3_z.values}))
    del M
del ens
P = pd.concat(parts, ignore_index=True); del parts
data = kev.prepare("cv", fold=None)
r = kev.evaluate(P, data, missing="error", keep_items=False)
E = r.endpoints
res["nested_pooled_A24"] = {k: E.get(k) for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M1_fold_mean"]}
res["nested_pooled_canary"] = r.canary
json.dump(res, open("out/nested_sel.json", "w"), indent=1, default=kev._json_default)
print(json.dumps(res["nested_pooled_A24"], default=kev._json_default), r.canary)
