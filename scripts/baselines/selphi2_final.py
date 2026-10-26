"""B4 SELPHI2.0-retrained final model (DESIGN §8.6; docs/G4_RUNBOOK.md): train on ALL DEV with the CV configuration,
score the TEST candidate set exhaustively.

Training = exactly the CV recipe of selphi2_cv.py (faithful arm): positives = DEV L2 labels (U_ms) on U_B4 sites
(U_ms with Ochoa features), 10 bags x 10x uniform random masked negatives, label-derived columns cross-fitted
(inner 5-fold over DEV substrate clusters), full SELPHI RF grid (GroupKFold(3) by cluster) on bag 0, 10 RFs, mean.
Inference = every TEST candidate pair (ksp.eval_test.test_candidate_pairs(): 503 kinases x 69,430 TEST U_ms sites,
residue-class matched). The runbook requires all-or-nothing per kinase, so sites WITHOUT Ochoa features receive
SELPHI's own NA defaults (cache/selphi2/site_features_all.parquet; how SELPHI itself scored PSP sites absent from the
Ochoa table) - counted as a manifest fallback. Label-derived TEST features: fit on all DEV labels (ksp.crossfit
eval path; LeakageError on any overlap). TEST labels are never read.

--mode cv0: sanity check of the imputation on DEV CV fold 0 (train on folds 1-4, score ALL fold-0 U_ms candidates,
evaluate with ksp.eval missing='error') -> results/baselines/selphi2/faithful_imputed_fold0.json.
Usage: python -I scripts/baselines/selphi2_final.py --mode test [--n-jobs 16]
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


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def train(tr_lab, tr_prot, cfg: M.Config, clusters, tag: int):
    """CV recipe on (tr_lab, tr_prot). tag = seed offset (CV fold number; 5 = final DEV model)."""
    S, K = load_sites(), load_kinases()
    kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc)
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    train_mask = M.selphi_universe() & np.isin(S.acc.to_numpy(), list(tr_prot))
    pos = tr_lab[tr_lab.kinase_acc.isin(kpos.index)][["kinase_acc", "site_id"]].drop_duplicates()
    pos = pos.assign(kinase_idx=kpos[pos.kinase_acc].to_numpy(), site_idx=spos.reindex(pos.site_id).to_numpy())
    pos = pos[pos.site_idx.notna()].astype({"site_idx": np.int64})
    pos = pos[train_mask[pos.site_idx.to_numpy()]].reset_index(drop=True)
    allrows = labels.load_dev("all")
    allrows = allrows[allrows.acc.isin(tr_prot)]
    masks = labels.NegativeMasks(allrows)
    kins = sorted(pos.kinase_acc.unique())
    negs = [M.sample_negatives(pos, train_mask, masks, kins, cfg.neg_ratio * len(pos),
                               np.random.default_rng(cfg.seed * 100000 + tag * 1000 + b)) for b in range(cfg.n_bags)]
    trows = pd.concat([pos] + negs, ignore_index=True).drop_duplicates(["kinase_idx", "site_idx"]).reset_index(drop=True)
    info = {"n_pos": len(pos), "n_pos_kinases": len(kins), "n_train_rows_unique": len(trows),
            "n_train_sites_U_B4": int(train_mask.sum()), "neg_per_bag": [len(n) for n in negs]}
    log("train set", info)
    lab_tr, _ = M.label_features(tr_lab, trows, None, None if tag == 5 else tag, clusters, seed=cfg.seed, log=log)
    asm = M.Assembler()
    X = asm.matrix(trows.kinase_idx.to_numpy(), trows.site_idx.to_numpy(), lab_tr)
    key = pd.Series(np.arange(len(trows)), index=pd.MultiIndex.from_arrays([trows.kinase_idx, trows.site_idx]))
    idx = lambda d: key.reindex(pd.MultiIndex.from_arrays([d.kinase_idx, d.site_idx])).to_numpy().astype(np.int64)  # noqa
    ip = idx(pos)
    groups = splits.site_cluster(trows.site_id, clusters)
    best, grid_tab, models = None, None, []
    for b, nb in enumerate(negs):
        ii = np.r_[ip, idx(nb)]
        y = np.r_[np.ones(len(ip)), np.zeros(len(nb))]
        if best is None:
            best, grid_tab = M.grid_search(X[ii], y, groups[ii], cfg.grid, cfg.n_jobs, cfg.seed, log=log)
            log("grid best", best)
        rf = RandomForestClassifier(n_estimators=best["n_estimators"], max_depth=best["max_depth"],
                                    min_samples_split=best["min_samples_split"], bootstrap=True, n_jobs=cfg.n_jobs,
                                    random_state=cfg.seed * 1000 + tag * 100 + b)
        rf.fit(X[ii], y)
        models.append(rf)
        log(f"bag {b} fitted")
    imp = pd.Series(np.mean([m.feature_importances_ for m in models], 0), index=F.SELPHI_COLS)
    info.update(best_params=best, grid=grid_tab.to_dict(orient="records"),
                feature_importance=imp.sort_values(ascending=False).round(5).to_dict())
    return models, info


def infer(models, tr_lab, cand: pd.DataFrame, clusters, cfg: M.Config, outer_fold, chunk_rows=6_500_000):
    """Score candidate pairs (kinase_idx, site_idx) in kinase-contiguous chunks; label features fit on all tr_lab."""
    S, K = load_sites(), load_kinases()
    sid, kacc = S.site_id.to_numpy(), K.uniprot_acc.to_numpy()
    asm = M.Assembler(impute_non_ochoa=True)
    cand = cand.sort_values(["kinase_idx", "site_idx"], kind="stable").reset_index(drop=True)
    kb = np.r_[0, np.where(np.diff(cand.kinase_idx.to_numpy()) != 0)[0] + 1, len(cand)]
    out = np.empty(len(cand), np.float64)
    a = 0
    while a < len(cand):
        b = kb[np.searchsorted(kb, a + chunk_rows, side="right") - 1]
        b = b if b > a else kb[np.searchsorted(kb, a, side="right")]
        q = cand.iloc[a:b].reset_index(drop=True)
        q = q.assign(kinase_acc=kacc[q.kinase_idx.to_numpy()], site_id=sid[q.site_idx.to_numpy()])
        q = q[["kinase_acc", "site_id", "kinase_idx", "site_idx"]]
        empty = q.iloc[:0]
        labs = []
        for fam in (F.SelphiPSSM, F.SelphiKSEACoreg):
            _, ev = crossfit.crossfit(fam, tr_lab, empty, q, seed=cfg.seed, outer_fold=outer_fold, clusters=clusters)
            labs.append(ev)
        lab = pd.concat(labs, axis=1)
        out[a:b] = M.predict_chunked(models, asm, q, lab, F.SELPHI_COLS)
        log(f"scored rows {a}..{b} of {len(cand)}")
        a = b
    n_imp = int((~asm.has_ochoa[cand.site_idx.to_numpy()]).sum())
    return cand.assign(score=out), n_imp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["test", "cv0"], default="test")
    ap.add_argument("--n-jobs", type=int, default=16)
    a = ap.parse_args()
    cfg = M.Config(neg_ratio=10, n_bags=10, grid=dict(M.GRID), n_jobs=a.n_jobs)
    clusters = splits.load_clusters()
    if a.mode == "cv0":
        from ksp import eval as E
        tr_lab, tr_prot, ev_prot = crossfit.outer_split(0, "L2")
        models, info = train(tr_lab, tr_prot, cfg, clusters, tag=0)
        S = load_sites()
        uni = S.in_U_ms.to_numpy() & np.isin(S.acc.to_numpy(), list(ev_prot))
        from ksp.eval_test import candidate_pairs
        pred, n_imp = infer(models, tr_lab, candidate_pairs(uni), clusters, cfg, 0)
        pred.to_parquet(F.CACHE / "faithful_imputed_fold0_predictions.parquet")
        data = E.prepare("cv", fold=0, with_canaries=False)
        p = pred[["kinase_idx", "site_idx", "score"]]
        p = p.assign(kinase_acc=load_kinases().uniprot_acc.to_numpy()[p.kinase_idx],
                     site_id=S.site_id.to_numpy()[p.site_idx])[["kinase_acc", "site_id", "score"]]
        res = E.evaluate(p, data, missing="error")
        ochoa = M.selphi_universe()
        r_m = E.evaluate(p[ochoa[pred.site_idx.to_numpy()]], data, missing="drop")
        oof = pd.read_parquet(F.CACHE / "faithful" / "oof_fold0.parquet")
        r_cv = E.evaluate(oof, data, missing="worst")
        keep = lambda e: {k: v for k, v in e.items() if not isinstance(v, (dict, list))}  # noqa
        rep = {"what": "fold-0 check of SELPHI-default imputation for non-Ochoa sites (train folds 1-4)",
               "n_pairs": len(p), "n_pairs_imputed_site": n_imp, "best_params": info["best_params"],
               "imputed_full_universe": keep(res.endpoints), "retrained_on_U_B4_only_drop": keep(r_m.endpoints),
               "cv_oof_U_B4_unscored_worst": keep(r_cv.endpoints)}
        (ROOT / "results" / "baselines" / "selphi2" / "faithful_imputed_fold0.json").write_text(
            json.dumps(rep, indent=1, default=E._json_default))
        print(json.dumps(rep, indent=1, default=E._json_default))
        return
    from ksp import eval_test as T
    from ksp.manifest import RunManifest
    inputs = [F.CACHE / "site_features.parquet", F.CACHE / "site_features_all.parquet", F.CACHE / "expr.npz",
              ROOT / "results" / "splits" / "cv_folds.tsv", ROOT / "results" / "splits" / "substrate_clusters.tsv",
              ROOT / "cache" / "labels" / "dev_pairs.parquet"]
    conf = {"baseline_id": "B4_retrained", "recipe": "SELPHI2.0 faithful arm of the CV (docs/baseline_selphi2.md)",
            **{k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.__dict__.items()},
            "train_universe": "U_B4 = U_ms sites with Ochoa 2020 features on DEV proteins",
            "test_scoring": "exhaustive over ksp.eval_test.test_candidate_pairs(); non-Ochoa sites get SELPHI NA defaults"}
    with RunManifest("b4_selphi2_final", config=conf, seed=cfg.seed, inputs=inputs, split="test",
                     notes="B4_retrained TEST batch deliverable; trained on all DEV; TEST labels never read") as run:
        t0 = time.time()
        tr_lab, tr_prot, _ = crossfit.outer_split(None, "L2")
        models, info = train(tr_lab, tr_prot, cfg, clusters, tag=5)
        train_wall = time.time() - t0
        log(f"train_wall_s {train_wall:.0f}")
        t1 = time.time()
        cand = T.test_candidate_pairs()
        pred, n_imp = infer(models, tr_lab, cand, clusters, cfg, None)
        out = pd.DataFrame({"kinase_idx": pred.kinase_idx.to_numpy(np.int64), "site_idx": pred.site_idx.to_numpy(np.int64),
                            "score": pred.score.to_numpy(np.float64)})
        f = run.save_predictions(out, name="test_predictions")
        inference_wall = time.time() - t1
        run.count_fallback("test_pairs_site_without_ochoa_features_SELPHI_NA_defaults", n_imp)
        run.set("train_wall_s", round(train_wall, 1))
        run.set("inference_wall_s", round(inference_wall, 1))
        run.set("train_info", info)
        run.set("n_test_pairs", int(len(out)))
        run.set("n_kinases_scored", int(out.kinase_idx.nunique()))
        chk = T.check_predictions(f)
        run.set("check_predictions", chk)
        log("check", chk)
        entry = {"name": "B4", "role": "baseline", "baseline_id": "B4_retrained",
                 "run_id": run.run_id, "predictions": str(f.relative_to(ROOT)), "score_site_predictions": None,
                 "t3_predictions": None, "runtime": {"train_wall_s": round(train_wall, 1),
                                                     "inference_wall_s": round(inference_wall, 1)},
                 "missing": "worst"}
        pinned = T.pin_files({"configs": [entry]})["configs"][0]
        run.set("freeze_entry", pinned)
        (ROOT / "results" / "baselines" / "selphi2" / "freeze_entry_B4.json").write_text(json.dumps(pinned, indent=1))
        print(json.dumps(pinned, indent=1))


if __name__ == "__main__":
    main()
