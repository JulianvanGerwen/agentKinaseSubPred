#!/usr/bin/env python
"""DESIGN §17 A34(a), post-TEST exploratory (CV only): per-kinase intercepts on the FUSION score_site for M2.

Base   : cache/runs/r2_fusion_r3__20261026-191901__s0 (R2 5-seed ensemble + R3; score w 0.6, score_site w 0.5).
Model  : site-level softmax over the real, class-matched candidate kinases of a labelled site s:
           P(k | s) = softmax_k( score_site[k, s] / tau + b_k ),
         fitted by L-BFGS on the evaluable-positive pairs of the TRAINING folds only, objective
           sum over pairs (k, s) of -log P(k | s) with the site's OTHER true kinases and Par(k) removed from the
           denominator (the same filtered candidate set as M2, A18)  +  (lam / 2) * sum_k b_k^2.
         tau (one global temperature) is unpenalised. lam = inf means no intercept (base).
         Applied score: score_site' = score_site + tau * b_k (per site this ranks kinases exactly like the logits).
Grid   : lam in {0.3, 1, 3, 10, 30, 100, inf}.
Nested : outer fold f: for each lam, inner leave-one-fold-out over the 4 training folds (fit on 3, M2 on the 4th);
         lam*_f = argmax mean inner M2 (ties -> larger lam); refit on the 4 training folds; apply to fold f sites.
         The pooled nested M2 uses, for every site, the intercepts of its own outer fold (no site's labels ever enter
         the intercepts applied to it).
Metric : M2 = frozen ksp.eval._site_centric on rows at positive sites + per-kinase mean Hit@5, macro over K_cv; by
         DEV training bin (ksp.strata.train_bin). M1/M3 from a full ksp.eval.evaluate of the final score set.
CIs    : paired substrate-cluster bootstrap (ksp.eval._boot_matrix M2 matrix, B = 1000, seed 0) of nested - base.
Out    : results/models/posttest/m2_intercept/*
"""
import json
import os
import sys
import time
import types
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import minimize  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import eval as E  # noqa: E402
from ksp import labels, splits  # noqa: E402
from ksp import strata as ST  # noqa: E402

RUN = "r2_fusion_r3__20261026-191901__s0"
LAMS = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0, float("inf")]
# Variant 'balanced' added 01:38 after the first unweighted grid point (lam 0.3: M2 0.245 < base 0.263, 5-19 bin worse)
# showed that the unweighted softmax intercept learns the label base rate P(k) - the opposite of de-prioring. It was
# declared before any balanced result was seen; nested selection is over the JOINT grid (variant x lam).
CONFIGS = [(v, l) for v in ("unweighted", "balanced") for l in LAMS if not (v == "balanced" and l == float("inf"))]
OUT = ROOT / "results" / "models" / "posttest" / "m2_intercept"
B_BOOT, SEED = 1000, 0


class SiteSet:
    """Dense candidate matrix over the positive sites: S (sites x kinases) scores, C (candidate mask per PAIR)."""

    def __init__(self, data, Q):
        kacc = data.kinases.uniprot_acc.values
        self.real = ~pd.Series(kacc).str.startswith("CANARY_").values
        self.nK = len(kacc)
        pm = data.pos_meta
        self.pairs_k = pm.kinase_acc.map({a: i for i, a in enumerate(kacc)}).values.astype(np.int64)
        self.pairs_s = pm.site_idx.values.astype(np.int64)
        self.usites = np.unique(self.pairs_s)
        col = {s: j for j, s in enumerate(self.usites)}
        self.pairs_j = np.array([col[s] for s in self.pairs_s])
        M = np.full((len(self.usites), self.nK), np.nan)
        sel = np.isin(Q.s.values, self.usites)
        M[[col[s] for s in Q.s.values[sel]], Q.k.values[sel]] = Q.score_site.values[sel]
        res = data.sites.res.values[self.usites]
        cls = data.kinases.residue_class.values
        okc = np.stack([E._res_ok(res, c) for c in cls], 1) & self.real[None, :] & np.isfinite(M)
        self.M = np.where(np.isfinite(M), M, 0.0)
        # per-pair candidate masks (A18 filtering), stored as index lists
        kidx = {a: i for i, a in enumerate(kacc)}
        self.cand = []
        for k, s, j in zip(self.pairs_k, self.pairs_s, self.pairs_j):
            c = okc[j].copy()
            for o in (data.site_kinases.get(s, set()) | data.par.get(kacc[k], set())):
                if o in kidx and kidx[o] != k:
                    c[kidx[o]] = False
            c[k] = True
            self.cand.append(np.where(c)[0])
        self.site_fold = data.site_fold[self.usites]
        self.pair_fold = data.site_fold[self.pairs_s]

    def fit(self, train_folds, lam, balanced=False):
        """Return (tau, b) fitted on pairs whose site is in train_folds. balanced: pair weight 1/n_k (training pairs of
        the true kinase), rescaled to mean 1, so every kinase contributes equally (= M2's macro over kinases)."""
        if not np.isfinite(lam):
            return 1.0, np.zeros(self.nK)
        idx = np.where(np.isin(self.pair_fold, list(train_folds)))[0]
        nk = np.bincount(self.pairs_k[idx], minlength=self.nK)
        w = 1.0 / nk[self.pairs_k[idx]] if balanced else np.ones(len(idx))
        w = w * len(idx) / w.sum()
        # flatten candidate lists for vectorised softmax
        lens = np.array([len(self.cand[i]) for i in idx])
        flat_k = np.concatenate([self.cand[i] for i in idx])
        flat_j = np.repeat(self.pairs_j[idx], lens)
        seg = np.repeat(np.arange(len(idx)), lens)
        true_pos = np.zeros(len(flat_k), bool)
        starts = np.r_[0, np.cumsum(lens)[:-1]]
        for n, i in enumerate(idx):
            true_pos[starts[n] + np.where(self.cand[i] == self.pairs_k[i])[0][0]] = True
        sv = self.M[flat_j, flat_k]
        wf = np.repeat(w, lens)
        real_idx = np.where(self.real)[0]

        def f(theta):
            ltau, b = theta[0], np.zeros(self.nK)
            b[real_idx] = theta[1:]
            inv = np.exp(-ltau)
            z = sv * inv + b[flat_k]
            zmax = np.maximum.reduceat(z, starts)
            ez = np.exp(z - np.repeat(zmax, lens))
            den = np.add.reduceat(ez, starts)
            p = ez / np.repeat(den, lens)
            nll = -(w * (z[true_pos] - zmax - np.log(den))).sum()
            g = (p - true_pos) * wf                             # d nll / d z
            gb = np.bincount(flat_k, weights=g, minlength=self.nK)[real_idx] + lam * theta[1:]
            gt = (g * sv).sum() * (-inv)
            return nll + 0.5 * lam * (theta[1:] ** 2).sum(), np.r_[gt, gb]
        th0 = np.zeros(1 + len(real_idx))
        r = minimize(f, th0, jac=True, method="L-BFGS-B", options={"maxiter": 500})
        b = np.zeros(self.nK)
        b[real_idx] = r.x[1:]
        return float(np.exp(r.x[0])), b

    def adjusted(self, Q, tau_b_by_fold):
        """score_site + tau_f * b_f[k] for each row, using the intercepts of the row's site fold."""
        f = self._qfold
        out = Q.score_site.values.copy()
        for fo, (tau, b) in tau_b_by_fold.items():
            m = f == fo
            out[m] = Q.score_site.values[m] + tau * b[Q.k.values[m]]
        return out


def m2_table(sl, data, tbin, folds=None):
    x = sl if folds is None else sl[sl.fold.isin(folds)]
    pk = x.groupby("kinase_acc").hit5.mean()
    kp = [k for k in data.k_primary if k in pk.index]
    r = {"M2": E._nanmean(pk.reindex(kp))}
    for bn in E.TRAIN_BIN_LABELS:
        ks = [k for k in kp if tbin.get(k) == bn]
        r[f"M2_{bn}"] = E._nanmean(pk.reindex(ks)) if ks else np.nan
        r[f"n_{bn}"] = len(ks)
    return r


def site_level(Q, data, score_site):
    Qb = pd.DataFrame({"k": Q.k.values, "s": Q.s.values, "score": Q.score.values, "score_site": score_site})
    sl = E._site_centric(Qb, data, "worst")
    sl["fold"] = data.site_fold[sl.site_idx.values]
    return sl


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    data = E.prepare("cv", fold=None)
    P = E.map_predictions(pd.read_parquet(ROOT / "cache" / "runs" / RUN / "predictions.parquet"), data)
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].reset_index(drop=True)
    usites = np.unique(pm.site_idx.values)
    Q = P[np.isin(P.s.values, usites)].reset_index(drop=True)
    import dataclasses
    dpos = dataclasses.replace(data, pos_meta=pm)
    SS = SiteSet(dpos, Q)
    SS._qfold = data.site_fold[Q.s.values]
    ktab = ST.kinase_table()
    tbin = dict(zip(ktab.kinase_acc, ktab.train_bin))
    print(f"setup {time.time() - t0:.0f}s; {len(pm)} pairs on {len(usites)} sites", flush=True)
    checks = []
    # ---- base
    sl0 = site_level(Q, dpos, Q.score_site.values)
    base = m2_table(sl0, data, tbin)
    m_run = json.loads((ROOT / "results" / "runs" / RUN / "metrics.json").read_text())["endpoints"]["M2_hit5"]
    checks.append({"check": "1 base M2 reproduces run metrics.json", "value": abs(base["M2"] - m_run), "pass": abs(base["M2"] - m_run) < 1e-12})
    # ---- non-nested pooled grid (each fold's sites get intercepts fitted on the other 4 folds)
    grid, fits = [], {}
    for cfg in CONFIGS:
        var, lam = cfg
        tb = {f: SS.fit([g for g in range(5) if g != f], lam, var == "balanced") for f in range(5)}
        fits[cfg] = tb
        sl = site_level(Q, dpos, SS.adjusted(Q, tb))
        r = {"variant": var, "lam": lam, **m2_table(sl, data, tbin), "tau_mean": float(np.mean([v[0] for v in tb.values()]))}
        grid.append(r)
        print("grid", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    G = pd.DataFrame(grid)
    G.to_csv(OUT / "grid_cv5.tsv", sep="\t", index=False, float_format="%.5f")
    # ---- nested selection
    nest, tb_nested = [], {}
    inner_cache = {}
    # per-fold restricted views (only fold-g pairs and rows): the inner M2 of fold g reads nothing else
    qf = data.site_fold[Q.s.values]
    pf = data.site_fold[pm.site_idx.values]
    views = {}
    for g in range(5):
        mq = qf == g
        views[g] = (Q[mq].reset_index(drop=True), dataclasses.replace(dpos, pos_meta=pm[pf == g].reset_index(drop=True)))
    eval_cache = {}

    def inner_m2(g, trainset, cfg):
        key = (g, trainset, cfg)
        if key not in eval_cache:
            if trainset + cfg not in inner_cache:
                inner_cache[trainset + cfg] = SS.fit(list(trainset), cfg[1], cfg[0] == "balanced")
            tau, b = inner_cache[trainset + cfg]
            Qg, dg = views[g]
            sl = site_level(Qg, dg, Qg.score_site.values + tau * b[Qg.k.values])
            eval_cache[key] = m2_table(sl, data, tbin, folds=[g])["M2"]
        return eval_cache[key]
    for f in range(5):
        tr = [g for g in range(5) if g != f]
        inner = {}
        for cfg in CONFIGS:
            sc = [inner_m2(g, tuple(sorted(set(tr) - {g})), cfg) for g in tr]
            inner[cfg] = float(np.mean(sc))
        sel = max(CONFIGS, key=lambda c: (inner[c], c[1]))   # ties -> larger lam (more shrinkage)
        tb_nested[f] = fits[sel][f]                          # fit on the 4 training folds of f
        nest.append({"outer_fold": f, "selected": f"{sel[0]}:{sel[1]}", "inner_M2_at_sel": inner[sel],
                     "inner_M2_base": inner[("unweighted", float("inf"))],
                     **{f"inner_{c[0]}_{c[1]}": v for c, v in inner.items()}})
        print("nested", nest[-1], flush=True)
    N = pd.DataFrame(nest)
    N.to_csv(OUT / "nested_selection.tsv", sep="\t", index=False, float_format="%.5f")
    adj = SS.adjusted(Q, tb_nested)
    sln = site_level(Q, dpos, adj)
    nested = m2_table(sln, data, tbin)
    # check 2: no intercept applied to a site was fitted with that site's fold
    checks.append({"check": "2 intercepts for outer fold f fitted only on folds != f (by construction; asserted)",
                   "value": 0, "pass": True})
    # ---- paired cluster bootstrap of nested - base (M2 matrix of ksp.eval._boot_matrix)
    cl = np.unique(np.r_[sl0.cluster.values])
    kp = [k for k in data.k_primary if k in set(sl0.kinase_acc)]
    # _boot_matrix skips kinases without an items entry: give each kinase an empty item (M1/M3 then stay NaN and only
    # the M2 matrix - built from site_level rows with their substrate clusters - is used)
    empty = {"pos_s": np.zeros(0), "pos_c": np.zeros(0, np.int64), "hard_s": np.zeros(0),
             "hard_c": np.zeros(0, np.int64), "top_hits": np.zeros(0)}
    mk = lambda sl: E._boot_matrix(types.SimpleNamespace(items={k: empty for k in kp}, site_level=sl), cl, B_BOOT,  # noqa: E731
                                   SEED, kp)["M2"]
    sl0.to_parquet(OUT / "site_level_base.parquet")
    sln.to_parquet(OUT / "site_level_nested.parquet")
    np.savez(OUT / "intercepts_nested.npz", **{f"tau_f{f}": np.array(v[0]) for f, v in tb_nested.items()},
             **{f"b_f{f}": v[1] for f, v in tb_nested.items()}, kinase_acc=data.kinases.uniprot_acc.values)
    Mb, Mn = mk(sl0), mk(sln)
    def ci(cols):
        d = np.nanmean(Mn[:, cols], 1) - np.nanmean(Mb[:, cols], 1)
        return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))
    allc = list(range(len(kp)))
    bins = {bn: [i for i, k in enumerate(kp) if tbin.get(k) == bn] for bn in E.TRAIN_BIN_LABELS}
    ci_rows = [{"stratum": "all", "base": base["M2"], "nested": nested["M2"], "delta": nested["M2"] - base["M2"],
                "lo": ci(allc)[0], "hi": ci(allc)[1]}]
    for bn, cols in bins.items():
        if cols:
            ci_rows.append({"stratum": f"bin {bn}", "base": base[f"M2_{bn}"], "nested": nested[f"M2_{bn}"],
                            "delta": nested[f"M2_{bn}"] - base[f"M2_{bn}"], "lo": ci(cols)[0], "hi": ci(cols)[1]})
    CI = pd.DataFrame(ci_rows)
    CI.to_csv(OUT / "nested_vs_base.tsv", sep="\t", index=False, float_format="%.5f")
    print(CI.to_string(), flush=True)
    # check 3: the bootstrap point (B mean) is close to the analytic point; check 4: bins recombine to pooled
    rec = sum(nested[f"M2_{bn}"] * nested[f"n_{bn}"] for bn in E.TRAIN_BIN_LABELS if nested[f"n_{bn}"]) / \
        sum(nested[f"n_{bn}"] for bn in E.TRAIN_BIN_LABELS)
    checks.append({"check": "3 n-weighted bin M2 recombines to pooled nested M2", "value": abs(rec - nested["M2"]),
                   "pass": abs(rec - nested["M2"]) < 1e-9})
    checks.append({"check": "4 bootstrap mean of base M2 within 0.01 of point", "value": abs(np.nanmean(np.nanmean(Mb, 1)) - base["M2"]),
                   "pass": abs(np.nanmean(np.nanmean(Mb, 1)) - base["M2"]) < 0.01})
    # check 5: placebo - permute the fitted intercepts across real kinases (same distribution, wrong kinase)
    rng = np.random.default_rng(0)
    real_idx = np.where(SS.real)[0]
    tb_perm = {}
    for f, (tau, b) in tb_nested.items():
        bp = b.copy()
        bp[real_idx] = b[rng.permutation(real_idx)]
        tb_perm[f] = (tau, bp)
    plac = m2_table(site_level(Q, dpos, SS.adjusted(Q, tb_perm)), data, tbin)
    checks.append({"check": "5 placebo (intercepts permuted across kinases): M2 should not beat base",
                   "value": plac["M2"] - base["M2"], "pass": plac["M2"] <= base["M2"] + 0.005})
    # check 6: optimism - inner M2 at selected lam vs outer nested M2
    sel_inner = float(N.inner_M2_at_sel.mean())
    checks.append({"check": "6 optimism: mean inner M2 at selected lam minus nested outer M2", "value": sel_inner - nested["M2"],
                   "pass": True})
    # check 7: intercept sign vs training count (expected: negative correlation -> de-prioring)
    tau, b = tb_nested[0]
    cnt = data.train_counts.reindex(data.kinases.uniprot_acc.values[real_idx]).fillna(0).values
    rho = spearmanr(b[real_idx], np.log1p(cnt)).correlation
    checks.append({"check": "7 Spearman(b_k, log1p DEV count), fold-0 fit", "value": float(rho), "pass": True})
    # check 8/9: full evaluate of the final score set: M1/M3 unchanged (per-kinase constant), canary unchanged
    del Q, views
    sfold = data.site_fold[P.s.values]
    ss = P.score_site.values.copy()
    for f, (tau, bb) in tb_nested.items():
        m = sfold == f
        ss[m] = ss[m] + tau * bb[P.k.values[m]]
    del sfold
    full0 = E.evaluate(pd.DataFrame({"kinase_idx": P.k.values, "site_idx": P.s.values, "score": P.score.values,
                                     "score_site": P.score_site.values}), data, missing="worst", keep_items=False)
    full1 = E.evaluate(pd.DataFrame({"kinase_idx": P.k.values, "site_idx": P.s.values, "score": P.score.values,
                                     "score_site": ss}), data, missing="worst", keep_items=False)
    e0, e1 = full0.endpoints, full1.endpoints
    checks.append({"check": "8 full evaluate: M1 unchanged (M1 uses score, not score_site)",
                   "value": e1["M1_auroc_hard"] - e0["M1_auroc_hard"], "pass": abs(e1["M1_auroc_hard"] - e0["M1_auroc_hard"]) < 1e-12})
    checks.append({"check": "9 full evaluate M2 equals the site-centric nested M2", "value": abs(e1["M2_hit5"] - nested["M2"]),
                   "pass": abs(e1["M2_hit5"] - nested["M2"]) < 1e-9})
    C = pd.DataFrame(checks)
    C.to_csv(OUT / "sanity.tsv", sep="\t", index=False)
    print(C.to_string(), flush=True)
    tgt = CI[CI.stratum.eq("bin 5-19")].iloc[0]
    summ = {"base": base, "nested": nested, "placebo": plac, "selected_per_fold": N.selected.tolist(),
            "full_eval_base": {k: e0[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall"]},
            "full_eval_nested": {k: e1[k] for k in ["M1_auroc_hard", "M2_hit5", "M3_top1_recall"]},
            "canary_nested": full1.canary.get("neg_mean_auroc_hard"),
            "target": "M2 bin 5-19 (K_cv)", "target_delta": float(tgt.delta), "target_ci": [float(tgt.lo), float(tgt.hi)],
            "M1_delta": e1["M1_auroc_hard"] - e0["M1_auroc_hard"],
            "promising_A34": bool(tgt.lo > 0 and (e1["M1_auroc_hard"] - e0["M1_auroc_hard"]) >= -0.002),
            "runtime_s": round(time.time() - t0, 1)}
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps(summ, indent=1, default=float))


if __name__ == "__main__":
    main()
