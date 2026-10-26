import sys, json
sys.path.insert(0, "/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/audit2")
from common import *
sf = site_fold(); S = index.load_sites(); ums = S.in_U_ms.values
m = json.load(open(R + ENS + "/manifest.json"))
bad = []; n = 0
for sr in m["config"]["seed_runs"]:
    for fr in json.load(open(R + sr + "/manifest.json"))["fold_runs"]:
        f = int(fr.split("_f")[-1][0])
        p = pd.read_parquet(C + fr + "/predictions.parquet", columns=["kinase_idx", "site_idx"])
        ok = (sf[p.site_idx.values] == f).all()
        cov = np.isin(np.where((sf == f) & ums)[0], p.site_idx.unique()).all()
        n += 1
        if not (ok and cov and p.kinase_idx.nunique() == 516):
            bad.append((fr, bool(ok), bool(cov), int(p.kinase_idx.nunique())))
print({"checked": n, "bad": bad})
