"""Sanity check of context features: enrichment of PSP/HMR KS protein pairs vs random pairs (NOT an official metric).

Positives: ksp.labels.load_dev('L2') (DEV proteins only) minus autophos, kinase in the 503-kinase index, substrate in the
protein index; collapsed to unique (kinase, protein) pairs.
Backgrounds (same kinase multiset as positives, seeds 0-4):
  RAND: kinase x uniformly random DEV protein, excluding any protein that is a known substrate (any row of
        load_dev('all'), any evidence level) of that kinase.
  SHUF: kinase x random L2 substrate protein of ANOTHER kinase (controls substrate study bias / degree).
For each kinase x protein feature in the descriptors (cache/features/context/*/descriptor.json):
  coverage = fraction of positives with a non-missing value; frac_pos>0 / frac_bg>0 (for count/flag features);
  AUROC_all (NaN imputed to -inf, i.e. 'absent' ranks lowest; ties averaged) and AUROC_obs (pairs non-missing only).
Flags any AUROC > 0.8 as potential leakage.
Run: .venv/bin/python -I scripts/features/context_sanity.py [--only ppi,loc,...]
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from ksp.index import load_kinases, load_proteins  # noqa: E402
from ksp.store import load_matrix  # noqa: E402

CTX = ROOT / "cache" / "features" / "context"
RES = ROOT / "results" / "features" / "context"
SEEDS = [0, 1, 2, 3, 4]


def positives():
    """Non-official sanity labels: ksp.labels.load_dev('L2') minus autophos (DEV proteins only; TEST is sealed).
    Returns (kinase idx, protein idx, known-pair set from load_dev('all'), K, P, DEV protein indices)."""
    from ksp import labels, splits
    kin = load_kinases().uniprot_acc.tolist()
    prot = load_proteins().acc.tolist()
    kpos = {a: i for i, a in enumerate(kin)}
    ppos = {a: i for i, a in enumerate(prot)}
    allp = labels.load_dev("all")[["kinase_acc", "acc"]].drop_duplicates()
    allp = allp[allp.kinase_acc.isin(kpos) & allp.acc.isin(ppos)]
    known = set(zip(allp.kinase_acc.map(kpos), allp.acc.map(ppos)))
    ks = labels.load_dev("L2")
    pos = ks[~ks.autophos.astype(bool)][["kinase_acc", "acc"]].drop_duplicates()
    pos = pos[pos.kinase_acc.isin(kpos) & pos.acc.isin(ppos)]
    pi = pos.kinase_acc.map(kpos).to_numpy()
    pj = pos.acc.map(ppos).to_numpy()
    dev = np.array(sorted(ppos[a] for a in splits.load_cv_folds().index if a in ppos))
    return pi, pj, known, len(kin), len(prot), dev


def backgrounds(pi, pj, known, dev, seed):
    rng = np.random.default_rng(seed)
    rand = np.empty_like(pj)
    shuf = np.empty_like(pj)
    pool = np.unique(pj)
    by_k = {}
    for k, j in zip(pi, pj):
        by_k.setdefault(k, set()).add(j)
    for n, k in enumerate(pi):
        while True:
            j = dev[rng.integers(len(dev))]
            if (k, j) not in known:
                rand[n] = j
                break
        while True:
            j = pool[rng.integers(len(pool))]
            if j not in by_k[k] and (k, j) not in known:
                shuf[n] = j
                break
    return rand, shuf


def auroc(pos, neg, impute=True):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    x = np.r_[pos, neg].astype(np.float64)
    if impute:
        lo = np.nanmin(x) if np.isfinite(x).any() else 0.0
        x = np.where(np.isnan(x), lo - 1.0, x)
    else:
        ok = ~np.isnan(x)
        y, x = y[ok], x[ok]
    if len(np.unique(y)) < 2 or len(x) < 20:
        return np.nan
    return float(roc_auc_score(y, x))


def evaluate_matrix(M, pi, pj, bgs):
    vp = M[pi, pj].astype(np.float64)
    row = {"coverage_pos": float(np.mean(~np.isnan(vp))), "frac_pos_gt0": float(np.nanmean(vp > 0) if (~np.isnan(vp)).any() else np.nan)}
    acc = {k: [] for k in ["frac_rand_gt0", "frac_shuf_gt0", "coverage_rand", "auc_rand_all", "auc_rand_obs",
                           "auc_shuf_all", "auc_shuf_obs"]}
    for rand, shuf in bgs:
        vr = M[pi, rand].astype(np.float64)
        vs = M[pi, shuf].astype(np.float64)
        acc["frac_rand_gt0"].append(np.mean(np.nan_to_num(vr, nan=0) > 0))
        acc["frac_shuf_gt0"].append(np.mean(np.nan_to_num(vs, nan=0) > 0))
        acc["coverage_rand"].append(np.mean(~np.isnan(vr)))
        acc["auc_rand_all"].append(auroc(vp, vr))
        acc["auc_rand_obs"].append(auroc(vp, vr, False))
        acc["auc_shuf_all"].append(auroc(vp, vs))
        acc["auc_shuf_obs"].append(auroc(vp, vs, False))
    for k, v in acc.items():
        row[k] = float(np.nanmean(v))
        if k.startswith("auc"):
            row[k + "_sd"] = float(np.nanstd(v))
    fp = np.mean(np.nan_to_num(vp, nan=0) > 0)
    row["ratio_gt0_pos_shuf"] = float(fp / row["frac_shuf_gt0"]) if row["frac_shuf_gt0"] > 0 else np.inf
    row["ratio_gt0_pos_rand"] = float(fp / row["frac_rand_gt0"]) if row["frac_rand_gt0"] > 0 else np.inf
    row["flag_gt0.8"] = bool(max(row["auc_rand_all"], row["auc_shuf_all"]) > 0.8)
    return row


def site_checks():
    """Site-level: fraction of L2 KS sites flagged extracellular/lumenal, FAM20C vs all other kinases vs all U_obs."""
    try:
        lum = np.asarray(load_matrix("context/loc/site_lumenal")[0], dtype=np.float32)
    except FileNotFoundError:
        return None
    from ksp.index import site_pos
    sp_ = site_pos()
    from ksp import labels
    ks = labels.load_dev("L2")
    ks = ks[ks.site_id.isin(sp_)][["kinase", "site_id"]].drop_duplicates()
    ks["lum"] = lum[ks.site_id.map(sp_).to_numpy()]
    g = ks.groupby("kinase").lum.agg(["mean", "count"])
    out = {"all_U_obs_sites_frac_lumenal": float(np.nanmean(lum)),
           "L2_sites_all_kinases_frac_lumenal": float(ks.lum.mean()),
           "top_kinases_by_frac_lumenal(min 10 sites)": g[g["count"] >= 10].sort_values("mean", ascending=False)
           .head(8).round(3).reset_index().to_dict("records")}
    if "FAM20C" in g.index:
        out["FAM20C"] = g.loc["FAM20C"].to_dict()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    pi, pj, known, K, P, dev = positives()
    bgs = [backgrounds(pi, pj, known, dev, s) for s in SEEDS]
    print(f"positives: {len(pi)} pairs, {len(np.unique(pi))} kinases, {len(np.unique(pj))} proteins", flush=True)
    mods = [d for d in sorted(CTX.glob("*/descriptor.json"))]
    if args.only:
        mods = [d for d in mods if d.parent.name in args.only.split(",")]
    rows = []
    for d in mods:
        desc = json.loads(d.read_text())
        for name, f in desc["features"].items():
            if f["shape"] != [K, P]:
                continue
            if "store_name" in f:
                M = load_matrix(f["store_name"])[0]
            else:
                M = np.load(ROOT / f["file"], mmap_mode="r")
            M = np.asarray(M, dtype=np.float32) * f.get("sign", 1)
            r = evaluate_matrix(M, pi, pj, bgs)
            r.update(module=d.parent.name, feature=name)
            rows.append(r)
            print(f"{name:40s} cov={r['coverage_pos']:.3f} AUCrand={r['auc_rand_all']:.3f}/{r['auc_rand_obs']:.3f} "
                  f"AUCshuf={r['auc_shuf_all']:.3f}/{r['auc_shuf_obs']:.3f} gt0 ratio shuf={r['ratio_gt0_pos_shuf']:.1f}", flush=True)
    df = pd.DataFrame(rows)
    cols = ["module", "feature"] + [c for c in df.columns if c not in ("module", "feature")]
    df = df[cols]
    tag = args.only.replace(",", "_") if args.only else "all"
    df.to_csv(RES / f"sanity_{tag}.tsv", sep="\t", index=False, float_format="%.4f")
    meta = {"n_pos_pairs": int(len(pi)), "n_kinases": int(len(np.unique(pi))), "n_proteins": int(len(np.unique(pj))),
            "seeds": SEEDS, "label_set": "load_dev(L2) & ~autophos (DEV only), protein-level, kinase in index, protein in index"}
    if not args.only or "loc" in args.only:
        meta["site_lumenal_check"] = site_checks()
        print(json.dumps(meta["site_lumenal_check"], indent=1))
    (RES / f"sanity_{tag}_meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
