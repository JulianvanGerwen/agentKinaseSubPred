"""Cross-fitting (DESIGN §7, §16 A1): planted-signal toy - a leaky label-derived feature gives AUROC ~1 when fitted on
all labels, ~0.5 under the cross-fitting driver; the driver raises LeakageError on overlapping blocks."""
import numpy as np
import pandas as pd
import pytest

from ksp import crossfit, splits
from ksp.eval import auroc


def _toy(n_prot=300, sites_per=4, frac=0.3, seed=0):
    rng = np.random.default_rng(seed)
    prots = [f"P{i:04d}" for i in range(n_prot)]
    # clusters: pairs of proteins are homologous (same cluster)
    clusters = pd.Series([i // 2 for i in range(n_prot)], index=prots)
    sites = [f"{p}_S{j + 1}" for p in prots for j in range(sites_per)]
    lab = rng.random(len(sites)) < frac
    rows = pd.DataFrame({"kinase_acc": "K1", "site_id": sites, "y": lab.astype(int)})
    return rows, rows[rows.y == 1][["kinase_acc", "site_id"]], clusters


def test_leaky_feature_without_crossfit_is_perfect():
    rows, labs, cl = _toy()
    f = crossfit.LeakyLabelFlag().fit(labs, None).transform(rows, None)
    assert auroc(f.is_labelled[rows.y == 1], f.is_labelled[rows.y == 0]) == 1.0


def test_leaky_feature_under_crossfit_is_chance():
    rows, labs, cl = _toy()
    tr, _ = crossfit.crossfit(crossfit.REGISTRY["toy_is_labelled"], labs, rows, None, clusters=cl)
    a = auroc(tr.is_labelled[rows.y == 1].values, tr.is_labelled[rows.y == 0].values)
    assert abs(a - 0.5) < 0.05
    assert tr.is_labelled.sum() == 0          # own cluster never fitted


def test_eval_rows_inner_mode_and_all_mode():
    rows, labs, cl = _toy(seed=1)
    # outer split: proteins of clusters >= 120 are eval
    ev_mask = cl[splits.site_protein(rows.site_id)].values >= 120
    tr_rows, ev_rows = rows[~ev_mask].reset_index(drop=True), rows[ev_mask].reset_index(drop=True)
    tr_labs = labs[~cl[splits.site_protein(labs.site_id)].values.__ge__(120)]
    a, e = crossfit.crossfit(crossfit.LeakyLabelFlag, tr_labs, tr_rows, ev_rows, clusters=cl)
    assert len(e) == len(ev_rows) and e.is_labelled.sum() == 0
    a2, e2 = crossfit.crossfit(crossfit.LeakyLabelFlag, tr_labs, tr_rows, ev_rows, clusters=cl, eval_mode="all")
    assert e2.is_labelled.sum() == 0


def test_driver_raises_on_leak():
    rows, labs, cl = _toy()
    with pytest.raises(crossfit.LeakageError):
        crossfit.crossfit(crossfit.LeakyLabelFlag, labs, rows, rows, clusters=cl, eval_mode="all")


class _CountFamily(crossfit.LabelDerivedFeature):
    """Legit protein-level signal: number of labelled sites of k on proteins of the same cluster *excluded* ->
    uses only other clusters; carries a real signal when labels concentrate on a protein class."""
    name = "toy_count"

    def fit(self, train_pairs, ctx):
        self.n = len(train_pairs)
        return self

    def transform(self, query, ctx):
        return pd.DataFrame({"n_train": np.full(len(query), float(self.n))})


def test_eval_inner_mode_matches_training_coverage():
    rows, labs, cl = _toy(seed=2)
    ev_mask = cl[splits.site_protein(rows.site_id)].values >= 120
    tr_rows, ev_rows = rows[~ev_mask].reset_index(drop=True), rows[ev_mask].reset_index(drop=True)
    tr_labs = labs[cl[splits.site_protein(labs.site_id)].values < 120]
    tr, ev = crossfit.crossfit(_CountFamily, tr_labs, tr_rows, ev_rows, clusters=cl)
    # eval rows see the same ~80% label coverage as training rows (A1), not 100%
    assert ev.n_train.max() < len(tr_labs)
    assert abs(ev.n_train.mean() - tr.n_train.mean()) / len(tr_labs) < 0.1
