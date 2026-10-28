"""Stream an exhaustive U_obs x kinase predictions parquet (pyarrow iter_batches; never loaded whole) and
  (1) check that a truncated per-kinase top-K file contains every row each library arm needs, including ties at the
      truncation boundary (superset = all rows with score >= the kinase's minimum score in the truncated file);
  (2) write the per-kinase top-K rows WITHIN the benchmark-measured universe (measured-universe arms), K = 2*|PSP|+200.
Arms need: top-N (N <= 100 or |PSP|) and top-N after removing PSP_all pairs (<= |PSP| + N rows) -> K = 2|PSP| + 200.
Usage: pred_topk.py --pred full.parquet --trunc topk.parquet --tag r2_clean
Outputs: cache/downstream/<tag>/{uobs_superset,measured_topk}.parquet, results/downstream/<tag>_topk_check.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from ksp.downstream import bench_data  # noqa: E402
from ksp.downstream import libraries as L  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402


def topk_per_kinase(d: pd.DataFrame, k_of: np.ndarray) -> pd.DataFrame:
    d = d.sort_values(["kinase_idx", "score", "site_idx"], ascending=[True, False, True], kind="mergesort")
    r = d.groupby("kinase_idx").cumcount().to_numpy()
    return d[r < k_of[d.kinase_idx.to_numpy()]].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--trunc", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--batch", type=int, default=4_000_000)
    a = ap.parse_args()
    out = ROOT / "cache" / "downstream" / a.tag
    out.mkdir(parents=True, exist_ok=True)
    K, S = load_kinases(), load_sites()
    psp = L.psp_all()
    sz = K.uniprot_acc.map(L.sizes(psp)).fillna(0).astype(int).to_numpy()
    k_of = 2 * sz + 200
    T = pd.read_parquet(a.trunc, columns=["kinase_idx", "site_idx", "score"])
    cut = np.full(len(K), np.inf)
    g = T.groupby("kinase_idx").score.min()
    cut[g.index.to_numpy()] = g.to_numpy()
    meas = bench_data.measured_sites(bench_data.load_perturbation(), bench_data.load_tumour())
    is_meas = pd.Index(S.site_id).isin(meas)
    sup, mea = [], []
    n_rows, n_k = 0, np.zeros(len(K), int)
    f = pq.ParquetFile(a.pred)
    for b in f.iter_batches(batch_size=a.batch, columns=["kinase_idx", "site_idx", "score"]):
        ki = b.column(0).to_numpy()
        si = b.column(1).to_numpy()
        sc = b.column(2).to_numpy()
        n_rows += len(ki)
        n_k += np.bincount(ki, minlength=len(K))
        m1 = sc >= cut[ki]
        sup.append(pd.DataFrame({"kinase_idx": ki[m1], "site_idx": si[m1], "score": sc[m1]}))
        m2 = is_meas[si]
        d = pd.DataFrame({"kinase_idx": ki[m2], "site_idx": si[m2], "score": sc[m2]})
        mea.append(topk_per_kinase(d, k_of))           # per-batch top-K is a superset of the global top-K
    sup = pd.concat(sup, ignore_index=True)
    mea = topk_per_kinase(pd.concat(mea, ignore_index=True), k_of)
    for d in (sup, mea):
        d["kinase_acc"] = K.uniprot_acc.to_numpy()[d.kinase_idx.to_numpy()]
        d["site_id"] = S.site_id.to_numpy()[d.site_idx.to_numpy()]
    sup.to_parquet(out / "uobs_superset.parquet", index=False)
    mea.to_parquet(out / "measured_topk.parquet", index=False)
    # --- checks: (a) superset vs truncated row counts (ties at the cut); (b) arms identical from both files
    ns = sup.groupby("kinase_idx").size().reindex(range(len(K)), fill_value=0).to_numpy()
    nt = T.groupby("kinase_idx").size().reindex(range(len(K)), fill_value=0).to_numpy()
    chk = {"full_rows": int(n_rows), "kinases_in_full": int((n_k > 0).sum()), "trunc_rows": int(len(T)),
           "superset_rows": int(len(sup)), "kinases_with_ties_at_cut": int((ns > nt).sum()),
           "kinases_trunc_short_of_K": int(((nt < k_of) & (n_k > nt)).sum())}
    Tl = T.assign(kinase_acc=K.uniprot_acc.to_numpy()[T.kinase_idx], site_id=S.site_id.to_numpy()[T.site_idx])
    arms = {}
    for n in ("PSP", 50, 100):
        nn = L.sizes(psp) if n == "PSP" else n
        for tag, ex in (("", None), ("-PSP", psp)):
            x = L.topn_from_long(Tl[["kinase_acc", "site_id", "score"]], nn, exclude=ex)
            y = L.topn_from_long(sup[["kinase_acc", "site_id", "score"]], nn, exclude=ex)
            same = set(zip(x.kinase_acc, x.site_id)) == set(zip(y.kinase_acc, y.site_id))
            # rows needed by the arm: max per-kinase rank (in the full ranking) of a selected row
            arms[f"top{n}{tag}"] = {"identical_trunc_vs_superset": bool(same), "pairs": int(len(y))}
    # rank needed by the '-PSP' arms: position of the N-th non-PSP row
    key = set(zip(psp.kinase_acc, psp.site_id))
    s2 = sup.sort_values(["kinase_idx", "score", "site_id"], ascending=[True, False, True])
    s2["rank"] = s2.groupby("kinase_idx").cumcount() + 1
    s2["is_psp"] = [(k, s) in key for k, s in zip(s2.kinase_acc, s2.site_id)]
    s2["nonpsp_rank"] = (~s2.is_psp).groupby(s2.kinase_idx).cumsum()
    need = {}
    for n in ("PSP", 50, 100):
        nv = sz if n == "PSP" else np.full(len(K), n)
        r = s2[s2.nonpsp_rank.to_numpy() == nv[s2.kinase_idx.to_numpy()]].groupby("kinase_idx")["rank"].min()
        need[f"top{n}-PSP_max_rank_minus_K"] = int((r.to_numpy() - k_of[r.index.to_numpy()]).max())
    chk["arms"] = arms
    chk["minus_psp_rank_needed_vs_K"] = need
    chk["measured_topk_rows"] = int(len(mea))
    (ROOT / "results" / "downstream" / f"{a.tag}_topk_check.json").write_text(json.dumps(chk, indent=1))
    print(json.dumps(chk, indent=1))


if __name__ == "__main__":
    main()
