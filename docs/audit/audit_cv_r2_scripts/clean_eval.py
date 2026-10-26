import json, sys
sys.path.insert(0, "/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/audit2")
from common import *
from ksp import eval as kev
sf = site_fold(); S = index.load_sites()
m = json.load(open(R + "r2_final_cleanbench_s0_cv5__20261026-203921__s0/manifest.json"))
parts, integ = [], {}
for fr in m["fold_runs"]:
    f = int(fr.split("_f")[-1][0]); p = load(C + fr + "/predictions.parquet")
    integ[fr] = bool((sf[p.site_idx.values] == f).all()); parts.append(p)
    fm = json.load(open(R + fr + "/manifest.json")); integ[fr + "_fallbacks"] = fm.get("fallbacks"); integ[fr + "_clean"] = fm["config"].get("clean")
P = pd.concat(parts, ignore_index=True)
data = kev.prepare("cv", fold=None)
r = kev.evaluate(P, data, missing="error", keep_items=False)
out = {"integrity": integ, "endpoints": {k: r.endpoints.get(k) for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M1_fold_mean"]}, "canary": r.canary}
json.dump(out, open("out/clean_eval.json", "w"), indent=1, default=kev._json_default); print(json.dumps(out, indent=1, default=kev._json_default))
