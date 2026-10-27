"""R2 LightGBM global kinase-site pair model: train per CV fold, exhaustive OOF prediction, ksp.eval evaluation.

    cfg = default_config(); cfg.update(...)
    res = run(cfg, folds=[0], name="r2_xxx")         # -> RunManifest dir, metrics, REGISTRY row

Training rows : ksp.models.r2.sample.FoldData(f).sample(ratio, hard, alpha, seed).
Early stopping: an inner validation split of TRAINING substrate clusters (ksp.splits.train_valid_split_inner, 20%),
                metric = macro per-kinase AUROC on the validation rows (binary / lambdarank) or NDCG@5 (site-grouped).
                The outer eval fold is never seen. With refit=True the model is refitted on all training rows with
                the early-stopped number of rounds.
Prediction    : EXHAUSTIVE over eval-fold U_ms sites x every kinase (503 real + 13 canaries) of matching residue
                class, gathered per kinase (CMP-1). Score = raw logit (margin).
Objectives    : 'binary' | 'lambdarank_kinase' (queries = random <=1000-row chunks within a kinase) |
                'site_grouped' (queries = training positive sites, candidates = kinases of matching class; positives =
                the site's training kinases; kinases labelled to the site in other sources / Par of a true kinase are
                removed, mirroring the M2 filter). Ranking objectives use no per-row weights (LightGBM ranking).
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd

from ksp import eval as kev
from ksp import labels, manifest, splits
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample

REGISTRY = ROOT / "results" / "cv" / "REGISTRY.tsv"
RESULTS = ROOT / "results" / "models" / "r2"
N_THREADS = int(os.environ.get("R2_THREADS", "16"))

MAIN_GROUPS = ["F1", "F0", "F2a", "F3", "F4", "F6a", "F6b", "F6c"]
FULL_GROUPS = MAIN_GROUPS + ["F5", "F7"]


def default_config() -> dict:
    return dict(groups=list(FULL_GROUPS), ratio=50, hard=0.5, alpha=0.5, objective="binary", num_leaves=63,
                learning_rate=0.05, min_data_in_leaf=100, feature_fraction=0.8, bagging_fraction=0.8,
                lambda_l2=1.0, max_bin=255, max_rounds=3000, es_rounds=100, refit=True, clean=False,
                comp_know=False, seed=0, variant="L2", miss_only=False, min_rounds=50, es_metric="mixed", es_real_only=False,
                f2a_family="f2a", desc_v=1)


# ---------------------------------------------------------------------------------------------------- metrics
def macro_auc_by_group(score: np.ndarray, y: np.ndarray, grp: np.ndarray) -> float:
    """Macro AUROC over groups having both classes (ties averaged)."""
    o = np.lexsort((score, grp))
    g, s, yy = grp[o], score[o], y[o]
    bounds = np.flatnonzero(np.r_[True, g[1:] != g[:-1], True])
    aucs = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        yb = yy[a:b]
        npos = yb.sum()
        nneg = len(yb) - npos
        if npos == 0 or nneg == 0:
            continue
        r = pd.Series(s[a:b]).rank().values
        aucs.append((r[yb == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))
    return float(np.mean(aucs)) if aucs else float("nan")


# ---------------------------------------------------------------------------------------------------- data
def _miss_view(A: assemble.Assembler, X: np.ndarray) -> tuple[np.ndarray, list]:
    """B0-miss: missingness indicators of the phospho (F6) columns only."""
    idx = [j for j, c in enumerate(A.cols) if c.startswith(("f6a_", "f6b_", "f6c_")) and c != "f6b_kinase_maxcos"]
    return np.isnan(X[:, idx]).astype(np.float32), [f"miss_{A.cols[j]}" for j in idx]


def site_grouped_rows(fd: sample.FoldData, seed: int, max_cand: int | None = None) -> pd.DataFrame:
    """Rows for the site-grouped ranking objective: every training positive site x candidate kinases of matching class
    (real kinases + canaries); label 1 for the site's training kinases; kinases with any other label row at the site
    (mask pairs) and paralogs of a true kinase that are not labelled to the site are removed."""
    rng = np.random.default_rng(seed)
    K = fd.K
    res = fd.S.res.values
    cls = K.residue_class.values
    pos = fd.pos
    by_site = pos.groupby("site_idx").kinase_idx.agg(set)
    mp = fd.mask_pairs
    spos = pd.Series(np.arange(len(fd.S)), index=fd.S.site_id)
    mp = mp[mp.site_id.isin(spos.index)]
    other = mp.assign(si=spos[mp.site_id].values).groupby("si").kinase_acc.agg(set)
    par = fd.masks.par
    out = []
    for s, ks in by_site.items():
        r = res[s]
        ok = np.array([sample.residue_ok(np.array([r]), c)[0] for c in cls])
        true_acc = {K.uniprot_acc.values[k] for k in ks}
        drop = set(other.get(s, set())) - true_acc
        for a in true_acc:
            drop |= set(par.get(a, set())) - true_acc
        for a in drop:
            if a in fd.kpos:
                ok[fd.kpos[a]] = False
        cand = np.flatnonzero(ok)
        negs = np.setdiff1d(cand, np.fromiter(ks, np.int64))
        if max_cand is not None and len(negs) > max_cand:
            negs = rng.choice(negs, max_cand, replace=False)
        kk = np.r_[np.fromiter(ks, np.int64), negs]
        out.append(pd.DataFrame({"kinase_idx": kk, "site_idx": s, "y": np.r_[np.ones(len(ks)), np.zeros(len(negs))]}))
    d = pd.concat(out, ignore_index=True)
    d["y"] = d.y.astype(int)
    d["w"] = 1.0
    return d


def build_train(fd: sample.FoldData, A: assemble.Assembler, cfg: dict):
    if cfg["objective"] == "site_grouped":
        d = site_grouped_rows(fd, cfg["seed"], cfg.get("site_max_cand"))
        stats = None
    else:
        d = fd.sample(cfg["ratio"], cfg["hard"], cfg["alpha"], cfg["seed"])
        stats = d.attrs["stats"]
    X = A.rows(d.kinase_idx.values, d.site_idx.values)
    cols = list(A.cols)
    if cfg["miss_only"]:
        X, cols = _miss_view(A, X)
    return d, X, cols, stats


def _groups_for(d: pd.DataFrame, obj: str, seed: int):
    """Sort order + query sizes for ranking objectives."""
    rng = np.random.default_rng(seed + 1)
    if obj == "lambdarank_kinase":
        key = np.zeros(len(d), np.int64)
        for k, ix in d.groupby("kinase_idx").indices.items():
            nq = max(1, int(np.ceil(len(ix) / 1000)))
            key[ix] = k * 100000 + rng.integers(0, nq, len(ix))
    elif obj == "site_grouped":
        key = d.site_idx.values.astype(np.int64)
    else:
        raise ValueError(obj)
    o = np.argsort(key, kind="stable")
    _, cnt = np.unique(key[o], return_counts=True)
    return o, cnt


def lgb_params(cfg: dict) -> dict:
    p = dict(learning_rate=cfg["learning_rate"], num_leaves=cfg["num_leaves"], min_data_in_leaf=cfg["min_data_in_leaf"],
             feature_fraction=cfg["feature_fraction"], bagging_fraction=cfg["bagging_fraction"], bagging_freq=1,
             lambda_l2=cfg["lambda_l2"], max_bin=cfg["max_bin"], num_threads=N_THREADS, seed=cfg["seed"],
             deterministic=True, force_col_wise=True, verbose=-1, max_cat_to_onehot=8, cat_smooth=10)
    obj = cfg["objective"]
    if obj == "binary":
        p.update(objective="binary", metric="None")
    elif obj == "lambdarank_kinase":
        p.update(objective="lambdarank", metric="None", lambdarank_truncation_level=1000)
    elif obj == "site_grouped":
        p.update(objective="lambdarank", metric="ndcg", eval_at=[5], lambdarank_truncation_level=50)
    else:
        raise ValueError(obj)
    return p


def train_fold(fd, A, cfg, log):
    t0 = time.time()
    d, X, cols, stats = build_train(fd, A, cfg)
    t_gather = time.time() - t0
    cat = [c for c in A.cat_cols if c in cols]
    va = splits.train_valid_split_inner(fd.S.site_id.values[d.site_idx.values], 0.2, cfg["seed"])
    obj = cfg["objective"]
    params = lgb_params(cfg)

    def mk(mask, ref=None):
        dd = d[mask].reset_index(drop=True)
        XX = X[mask]
        if obj == "binary":
            return lgb.Dataset(XX, dd.y.values, weight=dd.w.values, feature_name=cols, categorical_feature=cat,
                               reference=ref, free_raw_data=False), dd, XX
        o, cnt = _groups_for(dd, obj, cfg["seed"])
        return lgb.Dataset(XX[o], dd.y.values[o], group=cnt, feature_name=cols, categorical_feature=cat,
                           reference=ref, free_raw_data=False), dd.iloc[o].reset_index(drop=True), XX[o]
    dtr, dtr_rows, _ = mk(~va)
    dva, dva_rows, Xva = mk(va, dtr)
    hist = []

    es_metric = cfg.get("es_metric", "mixed")
    if obj != "site_grouped" and es_metric == "hard":
        # macro AUROC of positives vs HARD-type negatives only (aligned with M1 = E-hard); real kinases only
        keep_va = (dva_rows.neg_type.values != "all") & (dva_rows.kinase_idx.values < assemble.N_REAL)
    else:
        keep_va = np.ones(len(dva_rows), bool)
    if cfg.get("es_real_only"):        # PI 16:58: tune on real kinases only
        keep_va &= dva_rows.kinase_idx.values < assemble.N_REAL

    def feval(preds, ds):
        if obj == "site_grouped":
            return []
        auc = macro_auc_by_group(preds[keep_va], dva_rows.y.values[keep_va], dva_rows.kinase_idx.values[keep_va])
        return [("macro_auc_kinase", auc, True)]
    t1 = time.time()
    cb = [lgb.early_stopping(cfg["es_rounds"], first_metric_only=True, verbose=False), lgb.record_evaluation(hist_d := {})]
    bst = lgb.train(params, dtr, num_boost_round=cfg["max_rounds"], valid_sets=[dva], valid_names=["inner_val"],
                    feval=None if obj == "site_grouped" else feval, callbacks=cb)
    best_it = max(bst.best_iteration or cfg["max_rounds"], cfg["min_rounds"])
    best_val = None
    try:
        key = list(hist_d["inner_val"].keys())[0]
        best_val = float(hist_d["inner_val"][key][best_it - 1])
    except Exception as e:  # recorded, not hidden
        log.setdefault("warnings", []).append(f"best_val unavailable: {e!r}")
    t_es = time.time() - t1
    if cfg["refit"]:
        dall, _, _ = mk(np.ones(len(d), bool))
        t2 = time.time()
        bst = lgb.train(params, dall, num_boost_round=best_it)
        t_refit = time.time() - t2
    else:
        t_refit = 0.0
    log.update({"n_train_rows": int(len(d)), "n_pos": int(d.y.sum()), "n_val_rows": int(va.sum()),
                "best_iteration": int(best_it), "inner_val_metric": best_val, "t_gather_train_s": round(t_gather, 1),
                "t_train_es_s": round(t_es, 1), "t_refit_s": round(t_refit, 1), "n_features": len(cols)})
    return bst, cols, stats


def predict_fold(fd, A, bst, cfg, kinases=None):
    """Exhaustive OOF scores: eval-fold U_ms sites x kinases of matching class."""
    from concurrent.futures import ThreadPoolExecutor
    res = fd.S.res.values[fd.eval_sites]
    out_k, out_s, out_v = [], [], []
    ks = list(range(A.nK) if kinases is None else kinases)
    site_sets = {c: fd.eval_sites[sample.residue_ok(res, c)] for c in ("ST", "Y", "STY")}

    def gather(k):
        s = site_sets[A.K.residue_class.values[k]]
        X = A.kinase_block(k, s)
        if cfg["miss_only"]:
            X, _ = _miss_view(A, X)
        return k, s, X
    def bounded(ex, ks, width=24):
        for i in range(0, len(ks), width):
            yield from ex.map(gather, ks[i:i + width])
    with ThreadPoolExecutor(int(os.environ.get("R2_GATHER_THREADS", "6"))) as ex:
        for k, s, X in bounded(ex, ks):
            if not len(s):
                continue
            v = bst.predict(X, raw_score=True, num_threads=N_THREADS)
            out_k.append(np.full(len(s), k, np.int32))
            out_s.append(s.astype(np.int32))
            out_v.append(v.astype(np.float32))
    return pd.DataFrame({"kinase_idx": np.concatenate(out_k), "site_idx": np.concatenate(out_s),
                         "score": np.concatenate(out_v)})


# ---------------------------------------------------------------------------------------------------- driver
_EVAL_CACHE: dict = {}


def eval_data(fold):
    if fold not in _EVAL_CACHE:
        _EVAL_CACHE[fold] = kev.prepare("cv", fold=fold)
    return _EVAL_CACHE[fold]


def register(row: dict):
    cols = ["timestamp", "model", "run_id", "split", "fold", "config_name", "config_json", "selection_metric", "value",
            "n_configs_tried_model", "notes"]
    exists = REGISTRY.exists()
    n_prev = 0
    if exists:
        r = pd.read_csv(REGISTRY, sep="\t")
        n_prev = int((r.model == "R2").sum())
    row = {**row, "model": "R2", "n_configs_tried_model": n_prev + 1,
           "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    line = "\t".join(str(row.get(c, "")) for c in cols)
    with open(REGISTRY, "a") as fh:
        if not exists:
            fh.write("\t".join(cols) + "\n")
        fh.write(line + "\n")


def run(cfg: dict, folds=(0,), name: str = "r2", config_name: str = "", notes: str = "", save_preds: bool = True,
        evaluate: bool = True, register_row: bool = True, stage: str = ""):
    folds = list(folds)
    cfgc = dict(cfg)
    with manifest.RunManifest(name, config=cfgc, seed=cfg["seed"], split="cv",
                              notes=f"folds={folds}; {notes}") as run_:
        logs = {}
        preds = []
        fold_stats = {}
        for f in folds:
            t = time.time()
            fd = sample.FoldData(f, cfg["variant"])
            A = assemble.Assembler(cfg["groups"], f, clean=cfg["clean"], comp_know=cfg["comp_know"],
                                   f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
            for k_, v_ in A.fallbacks.items():
                run_.count_fallback(k_, v_)
            lg = {}
            bst, cols, stats = train_fold(fd, A, cfg, lg)
            t3 = time.time()
            p = predict_fold(fd, A, bst, cfg)
            lg["t_predict_s"] = round(time.time() - t3, 1)
            lg["n_pred_rows"] = int(len(p))
            lg["t_fold_s"] = round(time.time() - t, 1)
            logs[f] = lg
            preds.append(p)
            imp = pd.DataFrame({"feature": cols, "gain": bst.feature_importance("gain"),
                                "split": bst.feature_importance("split")})
            run_.dir.mkdir(parents=True, exist_ok=True)
            imp.to_csv(run_.dir / f"importance_fold{f}.tsv", sep="\t", index=False)
            if stats is not None:
                stats.to_csv(run_.dir / f"train_stats_fold{f}.tsv", sep="\t")
            print(f"[{name}] fold {f}: {lg}", flush=True)
        P = pd.concat(preds, ignore_index=True)
        run_.set("fold_logs", logs)
        if save_preds:
            run_.save_predictions(P)
        out = {"run_id": run_.run_id, "logs": logs}
        if evaluate:
            fe = folds[0] if len(folds) == 1 else None
            data = eval_data(fe)
            res = kev.evaluate(P, data, missing="error" if len(folds) == 5 or len(folds) == 1 else "drop")
            run_.save_metrics(res)
            res.per_kinase.to_csv(run_.dir / "per_kinase.tsv", sep="\t", index=False)
            out.update({"endpoints": res.endpoints, "canary": res.canary, "res": res})
            E = res.endpoints
            print(f"[{name}] M1/M5={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} "
                  f"M4={E['M4_auroc_hard_Kany']:.4f} canary={res.canary}", flush=True)
            if register_row:
                register({"run_id": run_.run_id, "split": "cv", "fold": "all" if len(folds) == 5 else ",".join(map(str, folds)),
                          "config_name": config_name or name, "config_json": json.dumps(
                              {k: v for k, v in cfg.items() if k != "groups"} | {"groups": "+".join(cfg["groups"])}),
                          "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(E["M1_auroc_hard"], 5),
                          "notes": f"{stage}; M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f}; "
                                   f"canary_neg={res.canary.get('neg_mean_auroc_hard', float('nan')):.3f}/"
                                   f"{res.canary.get('neg_mean_auroc_all', float('nan')):.3f} "
                                   f"pos={res.canary.get('pos_auroc_hard', float('nan')):.3f}; {notes}"})
        return out
