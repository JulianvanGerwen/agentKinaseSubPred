"""Audit pass 2: evaluate each spec with fold_normalise=False (harness as used for reported numbers) and True (A24)."""
import json, sys, time
sys.path.insert(0, "/workspace/src")
sys.path.insert(0, "/workspace/docs/audit")
import numpy as np, pandas as pd
from ksp import eval as kev, canary
import audit_cv_r1_recompute as rc

KEYS = ["M1_auroc_hard", "M2_hit5", "M3_top1_recall", "M4_auroc_hard_Kany", "M1_fold_mean", "n_missing_pairs",
        "n_required_pairs", "n_K_primary_evaluated", "n_K_any_evaluated", "sec_mrr", "sec_auroc_all"]

def main(name, out, modes):
    sp = rc.SPECS[name]; t0 = time.time()
    if "fold_runs" in sp:
        P = pd.concat([pd.read_parquet(rc.C + r + "/predictions.parquet", columns=["kinase_idx", "site_idx", "score"])
                       for r in sp["fold_runs"]], ignore_index=True)
    else:
        P = pd.read_parquet(sp["pooled"])
    P["kinase_idx"] = P.kinase_idx.astype(np.int64); P["site_idx"] = P.site_idx.astype(np.int64)
    data = kev.prepare("cv", fold=None)
    K = canary.kinases_with_canaries()
    log = {"spec": name}
    for fn in modes:
        if sp["mode"] == "lko":
            held = sorted(set(K.uniprot_acc.values[np.unique(P.kinase_idx.values)]) & set(data.positives))
            res = kev.evaluate(P, data, kinases=held, missing="drop", fold_normalise=fn, keep_items=False)
        elif sp["mode"] == "simple":
            real_k = sorted(k for k in data.positives if not k.startswith("CANARY_"))
            res = kev.evaluate(P, data, kinases=real_k, missing="worst", fold_normalise=fn, keep_items=False)
        else:
            res = kev.evaluate(P, data, missing="worst", fold_normalise=fn, keep_items=False)
        E = res.endpoints
        log[f"norm={fn}"] = {"endpoints": {k: E.get(k) for k in KEYS}, "canary": res.canary,
                             "missing_kinases": len(res.missing.get("per_kinase", {}))}
        res.per_kinase.to_csv(f"{out}/{name}_norm{int(fn)}_per_kinase.tsv", sep="\t", index=False)
    try:
        rep = json.load(open(sp["reported"])); log["reported"] = {k: rep.get("endpoints", {}).get(k) for k in KEYS}
        log["reported_canary"] = rep.get("canary")
    except Exception as e:
        log["reported_error"] = repr(e)
    log["t_s"] = round(time.time() - t0, 1)
    json.dump(log, open(f"{out}/{name}_v2.json", "w"), indent=1, default=kev._json_default)
    print("ok", name, log["t_s"])

if __name__ == "__main__":
    modes = [m == "1" for m in sys.argv[3].split(",")]
    main(sys.argv[1], sys.argv[2], modes)
