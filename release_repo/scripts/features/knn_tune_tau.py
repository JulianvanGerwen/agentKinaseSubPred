"""Choose the kernel temperature tau of knn_sim ONCE, on CV fold 0 INNER data (outer fold-0 eval proteins never used).

Protocol (non-official; labels via ksp.crossfit.outer_split(0) = DEV labels of the 4 outer-training folds only):
  * inner folds = ksp.splits.inner_folds over the outer-training U_ms sites + label sites (cluster-grouped, seed 0),
    exactly as the crossfit driver assigns them; inner fold 0 is held out.
  * fit labels = outer-training labels outside inner fold 0; held-out positives = labels inside inner fold 0.
  * per kinase with >= 5 held-out in-class positives: 300 random held-out U_ms sites of the kinase's residue
    types (in class) that carry no DEV label of that kinase (load_dev('all'), outer-training proteins only).
  * per space, macro (per-kinase) AUROC of knn_<sp>_lr (primary) and knn_<sp>_kde over a tau grid;
    tau = argmax macro AUROC of lr. Writes results/features/knn/tau.json and tau_grid.tsv.
Usage: CUDA_VISIBLE_DEVICES=2 python -I scripts/features/knn_tune_tau.py
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

from ksp import crossfit, labels, splits  # noqa: E402
from ksp.eval import auroc  # noqa: E402
from ksp.features import knn  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402

torch.set_num_threads(6)
OUT = ROOT / "results" / "features" / "knn"
GRID = {"cos": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5], "blosum": [0.2, 0.5, 1.0, 2.0, 4.0, 8.0]}
SP = knn.ALL_SPACES       # tau for every space; the default space set (knn.SPACES) is chosen from the same table
N_NEG, MIN_POS, SEED = 300, 5, 0


def class_ok(res: np.ndarray, cls: np.ndarray) -> np.ndarray:
    isY = res == "Y"
    return np.where(isY, np.isin(cls, ["Y", "STY"]), np.isin(cls, ["ST", "STY"]))


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    S, K = load_sites(), load_kinases()
    reps = knn.get_reps(SP)
    cl = splits.load_clusters()
    tr, tr_prot, _ = crossfit.outer_split(0)
    tr = tr[["kinase_acc", "site_id"]].drop_duplicates()
    tr = tr[tr.kinase_acc.isin(set(K.uniprot_acc)) & tr.site_id.isin(reps.pos)].reset_index(drop=True)
    S_tr = S[S.in_U_ms.values & S.acc.isin(tr_prot).values].reset_index(drop=True)
    sites = pd.Series(pd.concat([S_tr.site_id, tr.site_id]).unique())
    f_of = pd.Series(splits.inner_folds(sites, 5, SEED, cl), index=sites.values)
    lab_f = tr.site_id.map(f_of).to_numpy()
    fit_lab, ho = tr[lab_f != 0].reset_index(drop=True), tr[lab_f == 0].reset_index(drop=True)
    ho_sites = S_tr[S_tr.site_id.map(f_of).to_numpy() == 0].reset_index(drop=True)
    allk = labels.load_dev("all")[["kinase_acc", "site_id"]]
    any_lab = allk.groupby("kinase_acc").site_id.apply(set).to_dict()
    kcls = dict(zip(K.uniprot_acc, K.residue_class))
    res_of = dict(zip(S.site_id, S.res))
    rows = []
    for k, g in ho.groupby("kinase_acc"):
        res = np.array([res_of[s] for s in g.site_id])
        pos = g.site_id.to_numpy()[class_ok(res, np.full(len(res), kcls[k]))]
        if len(pos) < MIN_POS:
            continue
        cand = ho_sites[class_ok(ho_sites.res.to_numpy(), np.full(len(ho_sites), kcls[k]))].site_id.to_numpy()
        cand = cand[~np.isin(cand, list(any_lab.get(k, ())))]
        neg = rng.choice(cand, min(N_NEG, len(cand)), replace=False)
        rows.append(pd.DataFrame({"kinase_acc": k, "site_id": np.r_[pos, neg], "y": np.r_[np.ones(len(pos)),
                                                                                            np.zeros(len(neg))]}))
    Q = pd.concat(rows, ignore_index=True)
    ctx = crossfit.FitContext(cl, "tau_inner0", 0, SEED)
    crossfit._check_disjoint(fit_lab, Q, cl, "tau_inner0")       # held-out block is cluster-disjoint from fit labels
    print(f"fit labels {len(fit_lab)}, held-out kinases {Q.kinase_acc.nunique()}, rows {len(Q)}", flush=True)
    tab = []
    for gi in range(len(GRID["cos"])):
        tau = {sp: GRID["blosum" if sp == "blosum" else "cos"][gi] for sp in SP}
        f = knn.KnnSim(spaces=SP, tau=tau).fit(fit_lab, ctx).transform(Q[["kinase_acc", "site_id"]], ctx)
        for sp in SP:
            for st in ("kde", "lr"):
                a = []
                for k, idx in Q.groupby("kinase_acc").indices.items():
                    v, y = f[f"knn_{sp}_{st}"].to_numpy()[idx], Q.y.to_numpy()[idx]
                    m = ~np.isnan(v)
                    if (y[m] == 1).sum() and (y[m] == 0).sum():
                        a.append(auroc(v[m][y[m] == 1], v[m][y[m] == 0]))
                tab.append({"space": sp, "stat": st, "tau": tau[sp], "macro_auroc": float(np.mean(a)),
                            "n_kinases": len(a)})
        print(pd.DataFrame(tab).tail(2 * len(SP)).to_string(index=False), flush=True)
    tab = pd.DataFrame(tab)
    tab.to_csv(OUT / "tau_grid.tsv", sep="\t", index=False, float_format="%.4f")
    best = {sp: float(tab[(tab.space == sp) & (tab.stat == "lr")].sort_values("macro_auroc").tau.iloc[-1])
            for sp in SP}
    meta = {"tau": best, "criterion": "max macro per-kinase AUROC of knn_<sp>_lr on CV fold 0, inner fold 0 held out",
            "grid": GRID, "n_fit_labels": int(len(fit_lab)), "n_kinases": int(Q.kinase_acc.nunique()),
            "n_rows": int(len(Q)), "neg_per_kinase": N_NEG, "min_pos": MIN_POS, "seed": SEED,
            "c": 1.0, "n_bg": 2000, "runtime_s": round(time.time() - t0, 1)}
    (OUT / "tau.json").write_text(json.dumps(meta, indent=1))
    print(meta)


if __name__ == "__main__":
    main()
