"""B4 SELPHI2.0 retrained on our 5-fold DEV CV -> OOF predictions (DESIGN §8.1 B4).

Per outer fold f (results/splits/cv_folds.tsv):
  training labels = ksp.crossfit.outer_split(f, 'L2') (load_dev: DEV, U_ms);
  universe U_B4   = U_ms sites with Ochoa 2020 site features (SELPHI's feature-available set);
  positives       = training labels on U_B4 training sites;
  negatives       = per bag, neg_ratio x |positives| uniform random (kinase, site) pairs, masked (model.sample_negatives);
  label features  = ksp.crossfit (SelphiPSSM, SelphiKSEACoreg);
  RF              = grid search (GroupKFold(3) by substrate cluster) in cfg.grid_bags, then one RF per bag; mean prob;
  eval rows       = every real kinase x U_B4 fold-f site of its residue class.
Outputs: cache/selphi2/<arm>/oof_fold{f}.parquet [kinase_acc, site_id, score], results/baselines/selphi2/<arm>/
fold{f}.json (timings, grid table, counts).
Usage: python -I scripts/baselines/selphi2_cv.py --arm faithful --folds 0 1 2 3 4 [--pilot]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.ensemble import RandomForestClassifier  # noqa: E402

from ksp import crossfit, labels, splits  # noqa: E402
from ksp.baselines.selphi2 import features as F  # noqa: E402
from ksp.baselines.selphi2 import model as M  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402

ARMS = {
    # faithful: their features (incl. leaky ones), their RF grid, 10 bags (of 100), 10x negatives (paper)
    "faithful": dict(),
    # PRIOR-4 leakage ablation: drop functional_score, PWM_max_mss, netpho_max_KIN, isKinaseCoreg
    "noleak": dict(drop_cols=tuple(F.LEAKY_COLS)),
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def build_fold(fold: int, cfg: M.Config, clusters: pd.Series, pilot: bool = False):
    T = {}
    t0 = time.time()
    tr_lab, tr_prot, ev_prot = crossfit.outer_split(fold, cfg.label_variant)
    S = load_sites()
    K = load_kinases()
    kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc)
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    U = M.selphi_universe()
    acc = S.acc.to_numpy()
    train_mask = U & np.isin(acc, list(tr_prot))
    eval_mask = U & np.isin(acc, list(ev_prot))
    pos = tr_lab[tr_lab.kinase_acc.isin(kpos.index)][["kinase_acc", "site_id"]].drop_duplicates()
    pos = pos.assign(kinase_idx=kpos[pos.kinase_acc].to_numpy(), site_idx=spos.reindex(pos.site_id).to_numpy())
    pos = pos[pos.site_idx.notna()]
    pos["site_idx"] = pos.site_idx.astype(np.int64)
    pos = pos[train_mask[pos.site_idx.to_numpy()]].reset_index(drop=True)
    # masks from TRAINING-fold label rows only (no evaluation-fold label information enters training)
    allrows = labels.load_dev("all")
    allrows = allrows[allrows.acc.isin(tr_prot)]
    masks = labels.NegativeMasks(allrows)
    T["load_and_masks_s"] = round(time.time() - t0, 1)
    kins = sorted(pos.kinase_acc.unique())
    negs = []
    t1 = time.time()
    for b in range(cfg.n_bags):
        rng = np.random.default_rng(cfg.seed * 100000 + fold * 1000 + b)
        negs.append(M.sample_negatives(pos, train_mask, masks, kins, cfg.neg_ratio * len(pos), rng))
    T["neg_sampling_s"] = round(time.time() - t1, 1)
    trows = pd.concat([pos] + negs, ignore_index=True).drop_duplicates(["kinase_idx", "site_idx"]).reset_index(drop=True)
    erows = M.candidate_pairs(eval_mask, K)
    if pilot:
        erows = erows.sample(n=min(300_000, len(erows)), random_state=0).reset_index(drop=True)
    info = {"fold": fold, "n_train_labels": len(tr_lab), "n_pos": len(pos), "n_pos_kinases": len(kins),
            "n_train_sites_U_B4": int(train_mask.sum()), "n_eval_sites_U_B4": int(eval_mask.sum()),
            "n_train_rows_unique": len(trows), "n_eval_rows": len(erows),
            "neg_per_bag": [len(n) for n in negs]}
    log(f"fold {fold}: {info}")
    t2 = time.time()
    lab_tr, lab_ev = M.label_features(tr_lab, trows, erows, fold, clusters, seed=cfg.seed, log=log)
    T["crossfit_s"] = round(time.time() - t2, 1)
    info["label_feature_nonzero_frac_train"] = {c: float(np.nanmean(np.nan_to_num(lab_tr[c].to_numpy(float)) != 0))
                                                for c in F.LABEL_COLS}
    info["label_feature_nonzero_frac_eval"] = {c: float(np.nanmean(np.nan_to_num(lab_ev[c].to_numpy(float)) != 0))
                                               for c in F.LABEL_COLS}
    return pos, negs, trows, erows, lab_tr, lab_ev, info, T


def run_arm(fold, cfg, data, arm, clusters, out_cache, out_res, grid_override=None, pilot=False):
    pos, negs, trows, erows, lab_tr, lab_ev, info, T = data
    T = dict(T)
    cols = [c for c in F.SELPHI_COLS if c not in cfg.drop_cols]
    asm = M.Assembler()
    t0 = time.time()
    Xall = asm.matrix(trows.kinase_idx.to_numpy(), trows.site_idx.to_numpy(), lab_tr, cols)
    key = pd.Series(np.arange(len(trows)), index=pd.MultiIndex.from_arrays([trows.kinase_idx, trows.site_idx]))
    def idx(d):
        return key.reindex(pd.MultiIndex.from_arrays([d.kinase_idx, d.site_idx])).to_numpy().astype(np.int64)
    ip = idx(pos)
    groups_all = splits.site_cluster(trows.site_id, clusters)
    T["assemble_train_s"] = round(time.time() - t0, 1)
    best, grid_tab = grid_override, None
    models = []
    t1 = time.time()
    for b, nb in enumerate(negs):
        ii = np.r_[ip, idx(nb)]
        y = np.r_[np.ones(len(ip)), np.zeros(len(nb))]
        if best is None and b in cfg.grid_bags:
            tg = time.time()
            best, grid_tab = M.grid_search(Xall[ii], y, groups_all[ii], cfg.grid, cfg.n_jobs, cfg.seed, log=log)
            T["grid_search_s"] = round(time.time() - tg, 1)
            log(f"  grid best: {best}")
        ne = best["n_estimators"] if cfg.n_estimators_cap is None else min(best["n_estimators"], cfg.n_estimators_cap)
        tf = time.time()
        rf = RandomForestClassifier(n_estimators=ne, max_depth=best["max_depth"],
                                    min_samples_split=best["min_samples_split"], bootstrap=True, n_jobs=cfg.n_jobs,
                                    random_state=cfg.seed * 1000 + fold * 100 + b)
        rf.fit(Xall[ii], y)
        models.append(rf)
        log(f"  bag {b}: fit {time.time() - tf:.0f}s ({len(ii)} rows, {ne} trees)")
    T["fit_bags_s"] = round(time.time() - t1, 1)
    t2 = time.time()
    score = M.predict_chunked(models, asm, erows, lab_ev, cols)
    T["predict_s"] = round(time.time() - t2, 1)
    T["predict_rows_per_s"] = round(len(erows) / max(T["predict_s"], 1e-9))
    out = pd.DataFrame({"kinase_acc": erows.kinase_acc.to_numpy(), "site_id": erows.site_id.to_numpy(), "score": score})
    out_cache.mkdir(parents=True, exist_ok=True)
    out_res.mkdir(parents=True, exist_ok=True)
    tag = "pilot_" if pilot else ""
    out.to_parquet(out_cache / f"{tag}oof_fold{fold}.parquet")
    imp = pd.Series(np.mean([m.feature_importances_ for m in models], 0), index=cols).sort_values(ascending=False)
    rep = {"arm": arm, "config": {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.__dict__.items()},
           "best_params": best, "info": info, "timing_s": T, "feature_importance": imp.round(5).to_dict()}
    if grid_tab is not None:
        rep["grid"] = grid_tab.to_dict(orient="records")
    (out_res / f"{tag}fold{fold}.json").write_text(json.dumps(rep, indent=1, default=str))
    log(f"fold {fold} arm {arm} done: {T}")
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", default=["faithful"])
    ap.add_argument("--folds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    ap.add_argument("--n-bags", type=int, default=10)
    ap.add_argument("--neg-ratio", type=int, default=10)
    ap.add_argument("--grid", choices=["full", "pilot"], default="full")
    ap.add_argument("--grid-every-fold", action="store_true",
                    help="grid search in each fold (default: fold of first run only; winner reused)")
    ap.add_argument("--fixed-params", type=str, default=None, help="json dict: skip grid search")
    ap.add_argument("--n-estimators-cap", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=12)
    ap.add_argument("--pilot", action="store_true")
    a = ap.parse_args()
    clusters = splits.load_clusters()
    grid = dict(M.GRID) if a.grid == "full" else {"max_depth": [100], "min_samples_split": [8],
                                                   "n_estimators": [150, 500]}
    best = json.loads(a.fixed_params) if a.fixed_params else None
    for fold in a.folds:
        cfg0 = M.Config(neg_ratio=a.neg_ratio, n_bags=a.n_bags, grid=grid, n_jobs=a.n_jobs,
                        n_estimators_cap=a.n_estimators_cap)
        data = build_fold(fold, cfg0, clusters, pilot=a.pilot)
        fold_best = best if not a.grid_every_fold else None
        for arm in a.arms:       # the first arm grid-searches; later arms (ablations) reuse its per-fold winner
            cfg = M.Config(**{**cfg0.__dict__, **ARMS[arm]})
            b = run_arm(fold, cfg, data, arm, clusters, F.CACHE / arm, ROOT / "results" / "baselines" / "selphi2" / arm,
                        grid_override=fold_best, pilot=a.pilot)
            fold_best = b
            if best is None and not a.grid_every_fold:
                best = b


if __name__ == "__main__":
    main()
