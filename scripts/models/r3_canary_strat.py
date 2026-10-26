"""R3 canary diagnostic 2: negative-canary AUROC (E-all and E-hard) STRATIFIED by the protein's number of U_ms sites
(the cue the canary sampler leaks: canary sites are drawn uniformly over proteins, then 1-4 sites per protein, so they
are enriched on proteins with few U_ms sites - see scripts/models/r3_canary_diag.py). Within a density stratum that cue
is removed; a model that only exploits it gives ~0.5, while label leakage would remain visible.
Uses ksp.eval.prepare/masks and ksp.eval.auroc (no own metric code). Positives-weighted mean over strata with
>= 1 positive and >= 20 negatives.
  python -I scripts/models/r3_canary_strat.py --runs <cv fold run ids...> --tag <name>
Output: results/models/r3/canary_stratified_<tag>.tsv
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp import canary, index, manifest  # noqa: E402
from ksp import eval as kev  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--runs", nargs="+")
ap.add_argument("--tag", required=True)
a = ap.parse_args()
S = index.load_sites()
ums = S[S.in_U_ms]
n_prot = S.acc.map(ums.groupby("acc").size()).fillna(0).values.astype(int)
edges = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 23, 32, 45, 64, 100, 1e9])
stratum = np.searchsorted(edges, n_prot, side="right")
K = canary.kinases_with_canaries()
cidx = np.where(K.uniprot_acc.str.startswith("CANARY_").values)[0]
parts = []
for r in a.runs:
    p = pd.read_parquet(manifest.CACHE / r / "predictions.parquet", columns=["kinase_idx", "site_idx", "score"])
    parts.append(p[p.kinase_idx.isin(cidx)])
P = pd.concat(parts, ignore_index=True)
data = kev.prepare("cv", fold=None)
rows = []
for k in cidx:
    acc = K.uniprot_acc.values[k]
    sc = np.full(len(S), np.nan)
    x = P[P.kinase_idx == k]
    sc[x.site_idx.values] = x.score.values
    cand = data.universe & kev._res_ok(S.res.values, K.residue_class.values[k])
    excl = np.zeros(len(S), bool)
    excl[data.masks.excluded(acc)] = True
    pos = data.positives.get(acc, np.zeros(0, int))
    posm = np.zeros(len(S), bool)
    posm[pos] = True
    neg_all = cand & ~excl & ~posm & ~np.isnan(sc)
    hard = np.zeros(len(S), bool)
    hard[data.masks.hard_candidates(acc)] = True
    out = {"kinase_acc": acc, "n_pos": len(pos)}
    for nm, neg in [("all", neg_all), ("hard", neg_all & hard)]:
        out[f"auroc_{nm}_raw"] = kev.auroc(sc[pos], sc[neg])
        num = den = 0.0
        for b in np.unique(stratum[pos]):
            pb = pos[stratum[pos] == b]
            nb = neg & (stratum == b)
            if len(pb) and nb.sum() >= 20:
                num += len(pb) * kev.auroc(sc[pb], sc[nb])
                den += len(pb)
        out[f"auroc_{nm}_density_stratified"] = num / den if den else np.nan
        out[f"n_pos_used_{nm}"] = int(den)
    rows.append(out)
d = pd.DataFrame(rows)
neg = d[d.kinase_acc.str.startswith("CANARY_NEG")]
summ = {c: float(neg[c].mean()) for c in d.columns if c.startswith("auroc")}
d.to_csv(ROOT / f"results/models/r3/canary_stratified_{a.tag}.tsv", sep="\t", index=False)
print(d.to_string())
print("NEG mean:", summ)
