"""Shrink an exhaustive U_obs x kinase prediction parquet to the per-kinase top (2*|PSP_all(k)| + 200) rows, which
is all benchmarKIN's arms need (top-N with N <= |PSP|, 50, 100; '-PSP' arm ranks after removing PSP pairs). Ranking
within kinase is unchanged. Usage: r2_topk_filter.py --pred <in.parquet> --out <out.parquet>"""
import argparse
import numpy as np, pandas as pd
from ksp import index
from ksp.downstream import libraries as L
ap = argparse.ArgumentParser(); ap.add_argument("--pred"); ap.add_argument("--out"); a = ap.parse_args()
P = pd.read_parquet(a.pred, columns=["kinase_idx", "site_idx", "score"])
K = index.load_kinases(); S = index.load_sites()
sz = L.sizes(L.psp_all())
keep_n = (2 * K.uniprot_acc.map(sz).fillna(0).astype(int) + 200).values
o = np.lexsort((-P.score.values, P.kinase_idx.values))
k = P.kinase_idx.values[o]
start = np.searchsorted(k, np.arange(len(K) + 1))
rank = np.arange(len(k)) - start[k]
sel = o[rank < keep_n[k]]
Q = P.iloc[sel].reset_index(drop=True)
Q["kinase_acc"] = K.uniprot_acc.values[Q.kinase_idx.values]
Q["site_id"] = S.site_id.values[Q.site_idx.values]
Q.to_parquet(a.out, index=False)
print(len(P), "->", len(Q), "rows;", Q.kinase_idx.nunique(), "kinases")
