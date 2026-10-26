"""Per-fold weight sweep for leave-one-fold-out selection of the fusion weight. Usage: fold_sweep.py <fold>..."""
import json, sys, time
sys.path.insert(0, "/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/audit2")
from common import *
from ksp import eval as kev
WS = [round(0.1 * i, 1) for i in range(11)]
sf = site_fold()
ens = load(C + ENS + "/predictions.parquet")
for f in map(int, sys.argv[1:]):
    t0 = time.time()
    M, info = blend_fold(f, sf, ens)
    data = kev.prepare("cv", fold=f)
    rows = []
    for w in WS:
        P = pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                          "score": (w * M.score_r2_z.values + (1 - w) * M.score_r3_z.values)})
        res = kev.evaluate(P, data, missing="error", keep_items=False)
        E = res.endpoints
        rows.append({"fold": f, "w": w, "M1": E["M1_auroc_hard"], "M2": E["M2_hit5"], "M3": E["M3_top1_recall"],
                     "M4": E["M4_auroc_hard_Kany"], "neg_canary": res.canary.get("neg_mean_auroc_hard")})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"out/sweep_fold{f}.tsv", sep="\t", index=False)
    json.dump(info | {"t_s": time.time() - t0}, open(f"out/sweep_fold{f}_info.json", "w"))
