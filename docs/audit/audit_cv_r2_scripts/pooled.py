"""Pooled reproduction: fusion (saved), R2 ens5, rebuilt blend equality + alignment, R2 LKO raw vs A24."""
import json, sys, time
sys.path.insert(0, "<scratchpad>/audit2")
from common import *
from ksp import eval as kev
KEYS = ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M1_fold_mean", "M1_per_fold",
        "n_missing_pairs", "n_required_pairs", "fold_normalised", "n_K_primary_evaluated"]
out = {}
def ev(name, P, data, **kw):
    t = time.time(); res = kev.evaluate(P, data, keep_items=False, **kw)
    out[name] = {"endpoints": {k: res.endpoints.get(k) for k in KEYS}, "canary": res.canary, "t": round(time.time() - t)}
    res.per_kinase.to_csv(f"out/{name}_per_kinase.tsv", sep="\t", index=False)
    print(name, out[name], flush=True)
    json.dump(out, open("out/pooled.json", "w"), indent=1, default=kev._json_default)
sf = site_fold()
data = kev.prepare("cv", fold=None)
# 1. saved fusion predictions
F = pd.read_parquet(C + FUS + "/predictions.parquet")
F = F.astype({"kinase_idx": np.int64, "site_idx": np.int64})
ev("fusion_saved", F, data, missing="error")
# 3. rebuilt blend: alignment + equality
ens = load(C + ENS + "/predictions.parquet")
chk = {}
Fs = F.sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
del F
parts = []
for f in range(5):
    M, info = blend_fold(f, sf, ens)
    r3 = load(C + R3F[f] + "/predictions.parquet")
    info["r3_all_sites_in_fold"] = bool((sf[r3.site_idx.values] == f).all())
    info["r2_keys_eq_r3_keys"] = info["n_r2"] == info["n_r3"] == info["n_common"]
    chk[f] = info
    parts.append(pd.DataFrame({"kinase_idx": M.kinase_idx.values, "site_idx": M.site_idx.values,
                               "s06": 0.6 * M.score_r2_z.values + 0.4 * M.score_r3_z.values,
                               "s05": 0.5 * M.score_r2_z.values + 0.5 * M.score_r3_z.values}))
    del M, r3
Bl = pd.concat(parts, ignore_index=True).sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
same_keys = len(Bl) == len(Fs) and (Bl.kinase_idx.values == Fs.kinase_idx.values).all() and (Bl.site_idx.values == Fs.site_idx.values).all()
chk["same_keys_as_saved"] = bool(same_keys)
if same_keys:
    chk["max_abs_diff_score_w0.6"] = float(np.abs(Bl.s06.values - Fs.score.values).max())
    chk["max_abs_diff_score_site_w0.5"] = float(np.abs(Bl.s05.values - Fs.score_site.values).max())
chk["ens_sites_fold_ok"] = True
out["alignment"] = chk; print(chk, flush=True)
json.dump(out, open("out/pooled.json", "w"), indent=1, default=kev._json_default)
del Bl, Fs, parts
# 2. R2 ensemble alone
ev("ens5", ens, data, missing="error")
del ens
# 4. R2 LKO
L = load(C + LKO + "/predictions.parquet")
from ksp import canary
K = canary.kinases_with_canaries()
held = sorted(set(K.uniprot_acc.values[np.unique(L.kinase_idx.values)]) & set(data.positives))
out["lko_info"] = {"rows": len(L), "n_kinases": int(L.kinase_idx.nunique()), "n_held_with_pos": len(held),
                   "sites_all_dev": bool((sf[L.site_idx.values] >= 0).all())}
ev("lko_A24", L, data, kinases=held, missing="drop")
ev("lko_raw", L, data, kinases=held, missing="drop", fold_normalise=False)
print("DONE", flush=True)
