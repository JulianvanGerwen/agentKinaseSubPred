"""Paired comparison of two CV prediction runs on their MATCHED universe (PRIOR-1) with ksp.eval.paired_bootstrap.
Usage: r2_compare.py --a <run_id> --b <run_id> [--B 1000] [--fold N]   (default: pooled 5-fold CV)"""
import argparse, json
import numpy as np, pandas as pd
from ksp import eval as kev, canary, index
from ksp.index import ROOT

ap = argparse.ArgumentParser()
ap.add_argument("--a", required=True); ap.add_argument("--b", required=True)
ap.add_argument("--B", type=int, default=1000); ap.add_argument("--fold", type=int, default=None)
ap.add_argument("--out", default=None)
a = ap.parse_args()
K = canary.kinases_with_canaries(); S = index.load_sites()
kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc.values); spos = pd.Series(np.arange(len(S)), index=S.site_id.values)

def load(run_id):
    m = json.load(open(ROOT / "results" / "runs" / run_id / "manifest.json"))
    p = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"])
    if "kinase_idx" not in p:
        ki = p.kinase_acc.map(kpos); si = p.site_id.map(spos)
        if ki.isna().any() or si.isna().any():
            raise KeyError(f"{run_id}: unknown ids")
        p = pd.DataFrame({"kinase_idx": ki.values.astype(np.int64), "site_idx": si.values.astype(np.int64),
                          "score": p.score.values} | ({"score_site": p.score_site.values} if "score_site" in p else {}))
    p["key"] = p.kinase_idx.values.astype(np.int64) * len(S) + p.site_idx.values
    return p

pa, pb = load(a.a), load(a.b)
common = np.intersect1d(pa.key.values, pb.key.values)
rep = {"n_A": len(pa), "n_B": len(pb), "n_common": len(common)}
pa = pa[np.isin(pa.key.values, common)].sort_values("key").drop(columns="key")
pb = pb[np.isin(pb.key.values, common)].sort_values("key").drop(columns="key")
data = kev.prepare("cv", fold=a.fold)
ra = kev.evaluate(pa, data, missing="drop"); rb = kev.evaluate(pb, data, missing="drop")
pk_a = set(ra.items); pk_b = set(rb.items)
# kinases must have identical items; restrict to the kinase set both evaluate
pb_ = kev.paired_bootstrap(ra, rb, B=a.B, seed=0)
pb_["A"] = a.a; pb_["B"] = a.b
pb_["A_point"] = [ra.endpoints[{"M1": "M1_auroc_hard", "M2": "M2_hit5", "M3": "M3_top1_recall", "M4": "M4_auroc_hard_Kany"}[e]] for e in pb_.endpoint]
pb_["B_point"] = [rb.endpoints[{"M1": "M1_auroc_hard", "M2": "M2_hit5", "M3": "M3_top1_recall", "M4": "M4_auroc_hard_Kany"}[e]] for e in pb_.endpoint]
for k, v in rep.items():
    pb_[k] = v
pb_["n_kinases_A"] = ra.endpoints["n_K_primary_evaluated"]
print(pb_.to_string())
out = a.out or ROOT / "results" / "models" / "r2" / f"compare__{a.a.split('__')[0]}__vs__{a.b.split('__')[0]}.tsv"
pb_.to_csv(out, sep="\t", index=False)
print("written", out)
