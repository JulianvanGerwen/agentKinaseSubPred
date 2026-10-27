#!/usr/bin/env python
"""NON-OFFICIAL sanity check of F6 features (deliverable 5): single-feature AUROC of DEV L2 positives vs random
sites per kinase (macro), plus missingness (ascertainment) of labelled vs unlabelled sites.

* Labels: ksp.labels.load_dev('L2') (DEV proteins, U_ms). Never the sealed TEST labels.
* Negatives: per kinase, up to --n-neg random U_ms sites on DEV proteins of the kinase's residue class that carry no
  label of that kinase (seed 0). Not the official E-hard / E-all sets (no paralog masks) -> sanity only.
* F6a is label-derived: computed with a TEMPORARY protein-grouped 5-fold split of DEV proteins (seed 0); the
  features of a site in fold f use only labels on proteins outside fold f (plus within-call same-protein exclusion).
* Scores: F6b uses -beta, -z, -contrast (decrease on inhibition = substrate); F6a/F6c use +.
Outputs: results/features/phospho/sanity_auroc.tsv (per feature x kinase), sanity_summary.tsv, sanity_missingness.tsv
Usage: CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_sanity.py [--variant v]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import store  # noqa: E402
from ksp.labels import load_dev  # noqa: E402
from ksp.features import phospho as ph  # noqa: E402

RES_CLASS = {"ST": set("ST"), "Y": set("Y"), "STY": set("STY")}


def f6a_crossfit(lab, sites, kin, dev_mask, k=5, seed=0):
    accs = sites.acc[dev_mask].unique()
    rng = np.random.default_rng(seed)
    fold_of = pd.Series(rng.integers(0, k, len(accs)), index=accs)
    sfold = sites.acc.map(fold_of).fillna(-1).astype(int).to_numpy()
    lfold = lab.acc.map(fold_of).fillna(-1).astype(int).to_numpy()
    out = {f"f6a_{t}_{s}": np.full((len(kin), len(sites)), np.nan, np.float32)
           for t in ("comp", "cptac") for s in ("top", "med")}
    for f in range(k):
        q = np.nonzero(sfold == f)[0]
        tr = lab[lfold != f]
        for kind, tag in (("compendium", "comp"), ("cptac", "cptac")):
            r = ph.coreg_features(tr, q, ph.cached_blocks(kind), n_kinases=len(kin))
            out[f"f6a_{tag}_top"][:, q] = r["top"]
            out[f"f6a_{tag}_med"][:, q] = r["med"]
        print(f"  F6a fold {f}: {len(q)} query sites, {len(tr)} training pairs", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="")
    ap.add_argument("--n-neg", type=int, default=3000)
    ap.add_argument("--min-pos", type=int, default=5)
    ap.add_argument("--only", default="", help="comma list of non-F6a features to evaluate (default all)")
    ap.add_argument("--skip-f6a", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    sfx = f"_{a.variant}" if a.variant else ""
    sites, kin = ph.site_index(), ph.kinase_index()
    lab = load_dev("L2")
    lab = lab[lab.kinase_acc.isin(set(kin.uniprot_acc))]
    spos = {s: j for j, s in enumerate(sites.site_id)}
    kpos = {x: i for i, x in enumerate(kin.uniprot_acc)}
    lab = lab.assign(site_idx=lab.site_id.map(spos), kinase_idx=lab.kinase_acc.map(kpos)).dropna(
        subset=["site_idx"])
    lab["site_idx"] = lab.site_idx.astype(int)
    dev_prot = set(lab.acc)
    try:
        from ksp.splits import load_cv_folds
        dev_prot = set(load_cv_folds().index)
    except Exception as e:  # noqa: BLE001
        print("cv folds unavailable, DEV proteins = labelled proteins only:", e)
    dev_mask = (sites.acc.isin(dev_prot) & sites.in_U_ms.fillna(False).astype(bool)).to_numpy()
    feats = {}
    for name, sign in (("f6b_ridge_z", -1), ("f6b_ridge_beta", -1), ("f6b_contrast", -1), ("f6c_prot", 1),
                       ("f6c_aloop", 1)):
        if a.only and name not in a.only.split(","):
            continue
        M, _ = store.load_matrix(f"phospho/{name}{sfx}")
        feats[name] = (M, sign)
    cf = {} if a.skip_f6a else f6a_crossfit(lab[["kinase_acc", "site_id", "acc", "kinase_idx", "site_idx"]], sites,
                                           kin, dev_mask)
    for n, M in cf.items():
        store.save_matrix(f"phospho/sanity_{n}{sfx}", M.astype(np.float16), "kinase", "site",
                          meta=dict(non_official=True, desc="temporary protein-grouped 5-fold cross-fit, sanity only"))
        feats[n] = (M, 1)
    print(f"features ready {time.time()-t0:.0f}s", flush=True)
    rng = np.random.default_rng(0)
    res = sites.res.to_numpy()
    recs = []
    pos_by_k = lab.groupby("kinase_idx").site_idx.apply(lambda x: np.unique(x.to_numpy()))
    for k, pos in pos_by_k.items():
        if len(pos) < a.min_pos:
            continue
        cls = RES_CLASS[kin.residue_class[k]]
        cand = np.nonzero(dev_mask & np.isin(res, list(cls)))[0]
        cand = np.setdiff1d(cand, pos)
        neg = rng.choice(cand, size=min(a.n_neg, len(cand)), replace=False)
        for name, (M, sign) in feats.items():
            vp = np.asarray(M[k, pos], np.float32) * sign
            vn = np.asarray(M[k, neg], np.float32) * sign
            fp, fn = np.isfinite(vp), np.isfinite(vn)
            rec = dict(feature=name, kinase=kin.gene[k], n_pos=len(pos), cov_pos=fp.mean(), cov_neg=fn.mean())
            if fp.sum() >= a.min_pos and fn.sum() >= 20:
                y = np.r_[np.ones(fp.sum()), np.zeros(fn.sum())]
                rec["auroc_nonmissing"] = roc_auc_score(y, np.r_[vp[fp], vn[fn]])
                rec["n_pos_nonmissing"] = int(fp.sum())
            recs.append(rec)
    d = pd.DataFrame(recs)
    d.to_csv(ph.RESULTS / f"sanity_auroc{sfx}.tsv", sep="\t", index=False)
    summ = d.groupby("feature").agg(n_kinases_evaluable=("auroc_nonmissing", "count"),
                                    macro_auroc=("auroc_nonmissing", "mean"),
                                    median_auroc=("auroc_nonmissing", "median"),
                                    mean_cov_pos=("cov_pos", "mean"), mean_cov_neg=("cov_neg", "mean")).reset_index()
    summ.to_csv(ph.RESULTS / f"sanity_summary{sfx}.tsv", sep="\t", index=False)
    print(summ.to_string())
    # site-level missingness: labelled (any DEV L2 label) vs unlabelled DEV U_ms sites
    labelled = np.zeros(len(sites), bool)
    labelled[lab.site_idx.unique()] = True
    mrec = []
    for name, (M, _) in feats.items():
        anyk = np.zeros(len(sites), bool)
        for c0 in range(0, len(sites), 50000):
            anyk[c0:c0 + 50000] = np.isfinite(np.asarray(M[:, c0:c0 + 50000], np.float32)).any(0)
        mrec.append(dict(feature=name, frac_nonmissing_labelled=anyk[labelled & dev_mask].mean(),
                         frac_nonmissing_unlabelled=anyk[~labelled & dev_mask].mean(),
                         n_labelled=int((labelled & dev_mask).sum()), n_unlabelled=int((~labelled & dev_mask).sum())))
    m = pd.DataFrame(mrec)
    m.to_csv(ph.RESULTS / f"sanity_missingness{sfx}.tsv", sep="\t", index=False)
    print(m.to_string())
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
