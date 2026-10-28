"""NON-OFFICIAL sanity check of knn_sim on CV fold 0 (not a project metric; ksp.eval is the only metric code).

* Labels: ksp.crossfit.outer_split(0) (DEV labels of outer-training folds; real kinases, no canaries).
* Eval rows: fold-0 proteins. Positives = DEV L2 labels (U_ms, in residue class) of kinases with >= 5 fold-0
  positives; negatives = 500 random fold-0 U_ms sites per kinase of in-class residue types carrying no DEV label
  (load_dev('all')) of that kinase. No paralog / Hamming masks.
* knn_sim and the learned PSSM (f2a_pssm, ksp.features.motif_learned) both go through ksp.crossfit.crossfit with
  the SAME rows, seed and inner folds (eval rows: randomly assigned inner-fold fit, A1 default).
* Real-data leak demonstration on training rows (outer-training positives + random negatives):
  leaky (fit on all labels, no own-cluster exclusion) vs transform-robust (fit on all, exclusion on) vs driver.
* Timing: one fit on all fold-0 training labels + dense transform of all 503 kinases x all U_ms sites.
Writes results/features/knn/sanity_fold0_{per_kinase.tsv,summary.tsv,summary.json}.
Usage: CUDA_VISIBLE_DEVICES=2 python -I scripts/features/knn_sanity.py
"""
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "6")
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from ksp import crossfit, labels, splits  # noqa: E402
from ksp.eval import auroc  # noqa: E402
from ksp.features import knn  # noqa: E402
from ksp.features import motif_learned as ml  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402
from ksp.store import load_matrix  # noqa: E402

torch.set_num_threads(6)
OUT = ROOT / "results" / "features" / "knn"
N_NEG, N_NEG_TR, MIN_POS, SEED = 500, 200, 5, 0
SPACES = tuple(a for a in sys.argv[1:] if not a.startswith("-")) or knn.SPACES


def class_ok(res, cls):
    isY = res == "Y"
    return np.where(isY, np.isin(cls, ["Y", "STY"]), np.isin(cls, ["ST", "STY"]))


def per_kinase_auc(df, cols, min_pos=MIN_POS):
    rec = []
    for k, idx in df.groupby("kinase_acc").indices.items():
        y = df.y.to_numpy()[idx]
        r = {"kinase_acc": k, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())}
        for c in cols:
            v = df[c].to_numpy()[idx].astype(float)
            m = ~np.isnan(v)
            r[c] = auroc(v[m][y[m] == 1], v[m][y[m] == 0]) if ((y[m] == 1).sum() >= min_pos and
                                                                (y[m] == 0).sum() >= 20) else np.nan
            r[c + "__nan_frac"] = float(1 - m.mean())
        rec.append(r)
    return pd.DataFrame(rec)


def sample_rows(pos_pairs, cand_sites, any_lab, kcls, res_of, n_neg, rng, min_pos):
    rows = []
    cres = cand_sites.res.to_numpy()
    for k, g in pos_pairs.groupby("kinase_acc"):
        res = np.array([res_of[s] for s in g.site_id])
        pos = np.unique(g.site_id.to_numpy()[class_ok(res, np.full(len(res), kcls[k]))])
        if len(pos) < min_pos:
            continue
        cand = cand_sites.site_id.to_numpy()[class_ok(cres, np.full(len(cres), kcls[k]))]
        cand = cand[~np.isin(cand, list(any_lab.get(k, ())))]
        neg = rng.choice(cand, min(n_neg, len(cand)), replace=False)
        rows.append(pd.DataFrame({"kinase_acc": k, "site_id": np.r_[pos, neg],
                                  "y": np.r_[np.ones(len(pos)), np.zeros(len(neg))].astype(int)}))
    return pd.concat(rows, ignore_index=True)


def main():
    t00 = time.time()
    rng = np.random.default_rng(SEED)
    S, K = load_sites(), load_kinases()
    reps = knn.get_reps(SPACES)
    cl = splits.load_clusters()
    kset = set(K.uniprot_acc)
    kcls = dict(zip(K.uniprot_acc, K.residue_class))
    res_of = dict(zip(S.site_id, S.res))
    tr, tr_prot, ev_prot = crossfit.outer_split(0)
    tr = tr[["kinase_acc", "site_id"]].drop_duplicates()
    tr = tr[tr.kinase_acc.isin(kset)].reset_index(drop=True)
    dev = labels.load_dev("L2")[["kinase_acc", "site_id"]].drop_duplicates()
    acc = splits.site_protein(dev.site_id)
    ev_pos = dev[pd.Series(acc).isin(ev_prot).values & dev.kinase_acc.isin(kset).values & dev.site_id.isin(reps.pos)]
    allk = labels.load_dev("all")[["kinase_acc", "site_id"]]
    any_lab = allk.groupby("kinase_acc").site_id.apply(set).to_dict()
    ums = S[S.in_U_ms.values & S.site_id.isin(reps.pos).values]
    E = sample_rows(ev_pos, ums[ums.acc.isin(ev_prot).values], any_lab, kcls, res_of, N_NEG, rng, MIN_POS)
    keepk = set(E.kinase_acc)
    Trows = sample_rows(tr[tr.kinase_acc.isin(keepk) & tr.site_id.isin(reps.pos)],
                        ums[ums.acc.isin(tr_prot).values], any_lab, kcls, res_of, N_NEG_TR, rng, MIN_POS)
    print(f"eval rows {len(E)} ({E.y.sum()} pos, {E.kinase_acc.nunique()} kinases); train rows {len(Trows)}",
          flush=True)
    timing = {}
    # ---- knn through the driver
    t = time.time()
    ftr, fev = crossfit.crossfit(lambda: knn.KnnSim(spaces=SPACES), tr, Trows[["kinase_acc", "site_id"]],
                                 E[["kinase_acc", "site_id"]], seed=SEED, outer_fold=0, clusters=cl)
    timing["knn_crossfit_5inner_rows_s"] = round(time.time() - t, 1)
    # ---- learned PSSM through the driver (same rows / seed / folds)
    t = time.time()
    kp = {a: i for i, a in enumerate(K.uniprot_acc)}
    sp = {s: j for j, s in enumerate(S.site_id)}
    E2 = E.assign(kinase_idx=E.kinase_acc.map(kp), site_idx=E.site_id.map(sp))
    T2 = Trows.assign(kinase_idx=Trows.kinase_acc.map(kp), site_idx=Trows.site_id.map(sp))
    ptr, pev = crossfit.crossfit(crossfit.REGISTRY["f2a_pssm"], tr, T2, E2, seed=SEED, outer_fold=0, clusters=cl)
    timing["pssm_crossfit_s"] = round(time.time() - t, 1)
    # ---- Cantley (label-free reference)
    C = load_matrix("motif/cantley_unprimed_pct_published")
    cant = np.asarray(C[E2.kinase_idx.to_numpy(), E2.site_idx.to_numpy()], np.float32) if hasattr(C, "shape") \
        else np.asarray(C[0][E2.kinase_idx.to_numpy(), E2.site_idx.to_numpy()], np.float32)
    ev = pd.concat([E.reset_index(drop=True), fev, pev[["f2a_pssm", "f2a_pssm_ph"]]], axis=1)
    ev["cantley_pct_published"] = cant
    kcols = [c for c in knn.columns(SPACES) if c != "knn_log_ntrain"]
    cols = kcols + ["f2a_pssm", "f2a_pssm_ph", "cantley_pct_published"]
    # simple complementarity probe: per-kinase rank average of knn lr (best space) and PSSM
    pk = per_kinase_auc(ev, cols)
    # ---- real-data leak demonstration on training rows
    t = time.time()
    fam_all = knn.KnnSim(spaces=SPACES).fit(tr, crossfit.FitContext(cl))
    f_rob = fam_all.transform(Trows[["kinase_acc", "site_id"]], crossfit.FitContext(cl))
    f_leak = knn.KnnSim(spaces=SPACES, exclude_same_cluster=False).fit(tr, crossfit.FitContext(cl)) \
        .transform(Trows[["kinase_acc", "site_id"]], crossfit.FitContext(cl))
    timing["leak_demo_s"] = round(time.time() - t, 1)
    leak = {}
    for name, f in (("leaky_no_exclusion", f_leak), ("robust_transform_no_driver", f_rob), ("driver", ftr)):
        d = pd.concat([Trows.reset_index(drop=True), f.reset_index(drop=True)], axis=1)
        a = per_kinase_auc(d, kcols)
        leak[name] = {c: round(float(a[c].mean()), 4) for c in kcols}
    # ---- timing: one fit (all fold-0 training labels, 503 kinases) + dense over all U_ms sites
    torch.cuda.synchronize() if torch.cuda.is_available() else None
    t = time.time()
    fam = knn.KnnSim(spaces=SPACES).fit(tr, crossfit.FitContext(cl))
    D = fam.dense(ums.site_id.to_numpy(), list(K.uniprot_acc), crossfit.FitContext(cl))
    torch.cuda.synchronize() if torch.cuda.is_available() else None
    timing["dense_503k_x_all_ums_s"] = round(time.time() - t, 1)
    timing["dense_shape"] = [len(K), int(len(ums))]
    timing["dense_nan_frac_max_first_space"] = round(float(np.isnan(D[kcols[0]]).mean()), 4)
    del D
    # ---- summaries
    ntr = tr.groupby("kinase_acc").size()
    pk["n_train"] = pk.kinase_acc.map(ntr).fillna(0).astype(int)
    pk["bin"] = pd.cut(pk.n_train, [-1, 19, 99, 10 ** 6], labels=["<20", "20-99", ">=100"]).astype(str)
    pk["residue_class"] = pk.kinase_acc.map(kcls)
    pk["gene"] = pk.kinase_acc.map(dict(zip(K.uniprot_acc, K.gene)))
    pk.to_csv(OUT / "sanity_fold0_per_kinase.tsv", sep="\t", index=False, float_format="%.4f")
    rows = []
    for c in cols:
        r = {"feature": c, "macro_auroc": pk[c].mean(), "median": pk[c].median(), "n_kinases": int(pk[c].notna().sum()),
             "frac_kinases_beat_pssm": float((pk[c] > pk["f2a_pssm"]).mean())}
        for b in ("<20", "20-99", ">=100"):
            r[f"macro_{b}"] = pk[pk.bin == b][c].mean()
        r["macro_Yclass"] = pk[pk.residue_class == "Y"][c].mean()
        rows.append(r)
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "sanity_fold0_summary.tsv", sep="\t", index=False, float_format="%.4f")
    # correlations with PSSM on eval rows (pooled within kinase, Spearman, mean over kinases)
    corr = {}
    for c in kcols:
        rs = []
        for k, idx in ev.groupby("kinase_acc").indices.items():
            a, b = ev[c].to_numpy()[idx], ev["f2a_pssm"].to_numpy()[idx]
            m = ~(np.isnan(a) | np.isnan(b))
            if m.sum() > 20:
                rs.append(spearmanr(a[m], b[m]).statistic)
        corr[c] = round(float(np.nanmean(rs)), 3)
    # rank-average probe (knn lr per space + PSSM), within kinase
    comb = {}
    for sp in SPACES:
        c = f"knn_{sp}_lr"
        ev["_r1"] = ev.groupby("kinase_acc")[c].rank(pct=True)
        ev["_r2"] = ev.groupby("kinase_acc")["f2a_pssm"].rank(pct=True)
        ev[f"rankavg_{sp}_lr+pssm"] = ev[["_r1", "_r2"]].mean(axis=1)
        comb[f"rankavg_{sp}_lr+pssm"] = round(float(per_kinase_auc(ev, [f"rankavg_{sp}_lr+pssm"])
                                                    [f"rankavg_{sp}_lr+pssm"].mean()), 4)
    meta = {"spaces": list(SPACES), "tau": fam.tau, "n_eval_rows": int(len(E)), "n_eval_pos": int(E.y.sum()),
            "n_kinases": int(E.kinase_acc.nunique()), "n_train_labels": int(len(tr)), "timing": timing,
            "leak_demo_macro_auroc_training_rows": leak, "spearman_with_pssm": corr,
            "rank_average_with_pssm_macro_auroc": comb, "runtime_s": round(time.time() - t00, 1),
            "note": "NON-OFFICIAL sanity: random in-class negatives, no paralog/Hamming masks, not ksp.eval"}
    (OUT / "sanity_fold0_summary.json").write_text(json.dumps(meta, indent=1, default=float))
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(summ.round(3).to_string(index=False))
    print(json.dumps({k: meta[k] for k in ("timing", "leak_demo_macro_auroc_training_rows", "spearman_with_pssm",
                                           "rank_average_with_pssm_macro_auroc")}, indent=1, default=float))


if __name__ == "__main__":
    main()
