"""Synthetic self-tests for ksp.downstream (formulas vs. reference implementations). Run: python scripts/downstream/selftest.py"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from ksp.downstream import evaluate, inference, libraries  # noqa: E402


def main():
    rng = np.random.default_rng(0)
    sites = [f"P{i:05d}_S{i}" for i in range(400)]
    X = pd.DataFrame(rng.normal(size=(400, 3)), index=sites, columns=["a", "b", "c"])
    X.iloc[rng.choice(400, 60, replace=False), 1] = np.nan
    lib = pd.DataFrame({"kinase_acc": ["K1"] * 20 + ["K2"] * 4 + ["K3"] * 30,
                        "site_id": list(rng.choice(sites, 20, replace=False)) + sites[:4] + sites[100:130]})
    for col in X.columns:
        x = X[col].dropna()
        for k in ["K1", "K3"]:
            t = x[x.index.isin(lib.site_id[lib.kinase_acc.eq(k)])]
            z_ref = t.sum() / (x.std(ddof=1) * np.sqrt(len(t)))
            w = x.index.isin(t.index).astype(float)
            ulm_ref = stats.linregress(w, x.to_numpy())
            assert np.isclose(inference.infer(X, lib, "zscore").loc[k, col], z_ref)
            assert np.isclose(inference.infer(X, lib, "ulm").loc[k, col], ulm_ref.slope / ulm_ref.stderr)
            assert np.isclose(inference.infer(X, lib, "wmean").loc[k, col], t.mean())
            assert np.isclose(inference.infer(X, lib, "ksea").loc[k, col],
                              (t.mean() - x.mean()) * np.sqrt(len(t)) / x.std(ddof=1))
    assert inference.infer(X, lib, "zscore").loc["K2"].isna().all()       # < 5 targets
    p, q = rng.normal(1, 1, 50), rng.normal(0, 1, 300)
    q[:20] = p[:20]                                                        # ties
    assert np.isclose(evaluate.auroc(p, q), roc_auc_score(np.r_[np.ones(50), np.zeros(300)], np.r_[p, q]))
    from ksp.downstream import paired
    w = rng.integers(0, 3, 350).astype(float)
    sc, yy = np.r_[p, q], np.r_[np.ones(50), np.zeros(300)]
    assert np.isclose(paired._wauc(paired._prep_auc(sc, yy), w), roc_auc_score(yy, sc, sample_weight=w))
    sh = libraries.shuffled(lib, 1)
    assert libraries.sizes(sh).equals(libraries.sizes(lib)) and set(sh.site_id) == set(lib.site_id)
    print("downstream selftest OK")


if __name__ == "__main__":
    main()
