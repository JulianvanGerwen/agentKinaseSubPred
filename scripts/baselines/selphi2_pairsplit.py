"""NON-OFFICIAL reference: SELPHI2.0's own protocol (pair-level stratified CV) on our DEV data, one fold.

Purpose: size of the optimism of a pair-split evaluation relative to our substrate-cluster CV. NOT an endpoint.
Protocol (their learn.py, adapted):
  * rows = all DEV L2 positives on U_B4 + 10x uniform random masked negative pairs (DEV sites), one bag;
  * label-derived columns: ksp.crossfit over DEV clusters (inner 5-fold; no row's own cluster enters its feature) -
    i.e. features are NOT computed from all labels (SELPHI did), so this is a LOWER bound on their protocol's optimism;
  * StratifiedKFold(5, shuffle, seed 0) over PAIRS; fold 0 held out; RF with fold-0 grid winner of our run;
  * metric = their metric: pooled AUROC (and AP) of held-out positives vs held-out random negatives.
Our-split counterpart with the same metric: fold-0 OOF predictions of the faithful arm; held-out fold-0 L2 positives on
U_B4 vs 10x random masked fold-0 pairs.
Output: results/baselines/selphi2/pairsplit_reference.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

from ksp import crossfit, labels, splits  # noqa: E402
from ksp.baselines.selphi2 import features as F  # noqa: E402
from ksp.baselines.selphi2 import model as M  # noqa: E402
from ksp.index import load_kinases, load_sites  # noqa: E402


def main():
    t0 = time.time()
    res_dir = ROOT / "results" / "baselines" / "selphi2"
    best = json.loads((res_dir / "faithful" / "fold0.json").read_text())["best_params"]
    clusters = splits.load_clusters()
    lab, dev_prot, _ = crossfit.outer_split(None, "L2")
    S, K = load_sites(), load_kinases()
    kpos = pd.Series(np.arange(len(K)), index=K.uniprot_acc)
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    U = M.selphi_universe() & np.isin(S.acc.to_numpy(), list(dev_prot))
    pos = lab[lab.kinase_acc.isin(kpos.index)][["kinase_acc", "site_id"]].drop_duplicates()
    pos = pos.assign(kinase_idx=kpos[pos.kinase_acc].to_numpy(), site_idx=spos.reindex(pos.site_id).to_numpy())
    pos = pos[pos.site_idx.notna()].astype({"site_idx": np.int64})
    pos = pos[U[pos.site_idx.to_numpy()]].reset_index(drop=True)
    masks = labels.NegativeMasks(labels.load_dev("all"))
    rng = np.random.default_rng(12345)
    neg = M.sample_negatives(pos, U, masks, sorted(pos.kinase_acc.unique()), 10 * len(pos), rng)
    rows = pd.concat([pos, neg], ignore_index=True)
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    lab_tr, _ = M.label_features(lab, rows, None, None, clusters)
    X = M.Assembler().matrix(rows.kinase_idx.to_numpy(), rows.site_idx.to_numpy(), lab_tr)
    tr, te = next(StratifiedKFold(5, shuffle=True, random_state=0).split(X, y))
    rf = RandomForestClassifier(n_estimators=best["n_estimators"], max_depth=best["max_depth"],
                                min_samples_split=best["min_samples_split"], bootstrap=True, n_jobs=12, random_state=0)
    rf.fit(X[tr], y[tr])
    p = rf.predict_proba(X[te])[:, 1]
    # share of held-out positive pairs whose SITE also appears (with another kinase / as negative) in training rows
    tr_sites = set(rows.site_idx.to_numpy()[tr])
    te_pos_sites = rows.site_idx.to_numpy()[te][y[te] == 1]
    out = {"label": "NON-OFFICIAL reference (SELPHI pair-split protocol on our DEV data; one fold)",
           "best_params": best, "n_pos": int(len(pos)), "n_neg": int(len(neg)),
           "pairsplit": {"auroc": float(roc_auc_score(y[te], p)), "ap": float(average_precision_score(y[te], p)),
                         "n_test_pos": int(y[te].sum()), "n_test_neg": int((1 - y[te]).sum()),
                         "frac_test_pos_sites_seen_in_training": float(np.isin(te_pos_sites, list(tr_sites)).mean())}}
    # our-split counterpart, same metric
    oof = pd.read_parquet(F.CACHE / "faithful" / "oof_fold0.parquet")
    folds = splits.load_cv_folds()
    f0 = set(folds.index[folds == 0])
    ev = labels.load_dev("L2")
    ev = ev[splits.site_protein(ev.site_id) != None]  # noqa: E711
    ev = ev[pd.Series(splits.site_protein(ev.site_id)).isin(f0).to_numpy()]
    key = oof.kinase_acc + "|" + oof.site_id
    pk = set(ev.kinase_acc + "|" + ev.site_id)
    is_pos = key.isin(pk).to_numpy()
    oof = oof.assign(kinase_idx=kpos[oof.kinase_acc].to_numpy(), site_idx=spos[oof.site_id].to_numpy())
    cand = oof[~is_pos]
    keep = np.ones(len(cand), bool)
    for k, g in cand.groupby("kinase_acc"):
        keep[cand.index.get_indexer(g.index[np.isin(g.site_idx.to_numpy(), masks.excluded(k))])] = False
    cand = cand[keep]
    pp = oof[is_pos]
    nn = cand.sample(n=min(10 * len(pp), len(cand)), random_state=0)
    yy = np.r_[np.ones(len(pp)), np.zeros(len(nn))]
    ss = np.r_[pp.score.to_numpy(), nn.score.to_numpy()]
    out["our_cluster_split_fold0"] = {"auroc": float(roc_auc_score(yy, ss)), "ap": float(average_precision_score(yy, ss)),
                                      "n_test_pos": int(len(pp)), "n_test_neg": int(len(nn))}
    out["inflation_auroc"] = out["pairsplit"]["auroc"] - out["our_cluster_split_fold0"]["auroc"]
    out["seconds"] = round(time.time() - t0, 1)
    (res_dir / "pairsplit_reference.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
