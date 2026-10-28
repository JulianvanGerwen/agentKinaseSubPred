#!/usr/bin/env python
"""DESIGN §17 A32(b): prior-debiased site score for M2 (post hoc on the candidate final OOF fusion; no retraining).

Input : cache/runs/r2_fusion_r3__20261026-191901__s0/predictions.parquet  (R2 5-seed ensemble + R3; score = w 0.6
        blend for M1, score_site = w 0.5 blend for M2; z = per-fold SD-standardised logits; docs/model_r2.md §6).
Score : score_site' = score_site - beta * prior(k, fold), beta in {0, 0.1, 0.25, 0.5, 0.75, 1}.
        prior(k, f) = log(1 + n_train(k, f)), n_train = DEV L2 (U_ms) pairs of k on proteins outside CV fold f
        (= the fold model's training labels; +1 because some kinases have 0 training labels in a fold).
SCALE RULE (pre-declared, applied by the code BEFORE any M2 is computed; output scale_check.json):
        s_site  = median over positive sites of the SD of score_site across that site's real candidate kinases;
        s_prior = SD of prior(k, f) across real kinases (pooled over folds).
        If 0.25 <= s_prior / s_site <= 4 the log-count prior is used as is; otherwise the prior is rescaled to unit
        variance ((prior - mean) / SD; the constant does not change within-site ranks).
M2    : exactly ksp.eval: E._site_centric (frozen; only rows at positive sites are read) then per-kinase mean Hit@5,
        macro over K_cv (= _endpoints). The beta = 0 value is checked against the run's metrics.json.
        Selected beta is re-evaluated with the full E.evaluate (all endpoints + canaries).
Selection: beta* = argmax pooled M2 (A32b). Nested leave-one-fold-out: for each fold f, beta chosen on the other
        4 folds' positive site pairs, M2 of fold f's pairs at that beta (macro over K_cv kinases with pairs in f).
Adoption (A32b literal): adopt beta* > 0 only if M2(beta*) >= M2(0) - 0.003 AND 5-19-bin M2(beta*) > M2(0)-bin.
Out   : results/models/priming/m2_debias/{scale_check.json, grid.tsv, nested.tsv, decision.json}
"""
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import labels, splits  # noqa: E402

RUN = "r2_fusion_r3__20261026-191901__s0"
BETAS = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0]
OUT = ROOT / "results" / "models" / "priming" / "m2_debias"


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    data = E.prepare("cv", fold=None)
    P = pd.read_parquet(ROOT / "cache" / "runs" / RUN / "predictions.parquet")
    P = E.map_predictions(P, data)                       # k, s, score, score_site (validated)
    fold = data.site_fold[P.s.values]
    kacc = data.kinases.uniprot_acc.values
    nK = len(kacc)
    real = ~pd.Series(kacc).str.startswith("CANARY_").values
    # per-(kinase, fold) training counts: DEV L2 (U_ms) pairs on proteins outside fold f
    L2 = labels.load_dev("L2")
    L2 = L2[L2.kinase_acc.isin(set(kacc))]
    folds = splits.load_cv_folds()
    lf = L2.acc.map(folds).values
    kpos = {a: i for i, a in enumerate(kacc)}
    ki = L2.kinase_acc.map(kpos).values
    tot = np.bincount(ki, minlength=nK)
    ntr = np.stack([tot - np.bincount(ki[lf == f], minlength=nK) for f in range(5)], 1)   # nK x 5
    prior_kf = np.log1p(ntr.astype(float))
    # ---- scale check (before any M2)
    usites = np.unique(data.pos_meta.site_idx.values)
    at = np.isin(P.s.values, usites)
    Q = P[at].copy()
    Q["fold"] = fold[at]
    res_arr = data.sites.res.values
    cls = data.kinases.residue_class.values
    rq, cq = res_arr[Q.s.values], cls[Q.k.values]
    ok = ((cq == "STY") | ((cq == "ST") & np.isin(rq, ["S", "T"])) | ((cq == "Y") & (rq == "Y"))) & real[Q.k.values]
    s_site = float(Q[ok].groupby("s").score_site.std().median())
    s_prior = float(prior_kf[real].std())
    ratio = s_prior / s_site
    rescale = not (0.25 <= ratio <= 4)
    mu, sd = float(prior_kf[real].mean()), float(prior_kf[real].std())
    prior_used = (prior_kf - mu) / sd if rescale else prior_kf
    sc = {"s_site_median_within_site_sd": s_site, "s_prior_sd_log1p_count": s_prior, "ratio": ratio,
          "rule": "use log1p count if 0.25 <= ratio <= 4 else rescale prior to unit variance",
          "decision": "rescale_unit_variance" if rescale else "log1p_count_as_is",
          "score_site_sd_overall": float(P.score_site.std()), "prior_mean": mu, "prior_sd": sd,
          "decided_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "note": "written before any M2 was computed"}
    (OUT / "scale_check.json").write_text(json.dumps(sc, indent=1))
    print("scale check:", sc, flush=True)
    # ---- grid (M2 only needs rows at positive sites; _site_centric reads nothing else)
    pr = prior_used[Q.k.values, np.clip(Q.fold.values, 0, 4)]
    kp = data.k_primary
    tb = pd.Series(data.train_counts).reindex(kacc).fillna(0)
    tbin = pd.Series(pd.cut(tb, E.TRAIN_BINS, labels=E.TRAIN_BIN_LABELS).astype(str).values, index=kacc)
    rows, levels = [], {}
    for b in BETAS:
        Qb = pd.DataFrame({"k": Q.k.values, "s": Q.s.values, "score": Q.score.values,
                           "score_site": Q.score_site.values - b * pr})
        sl = E._site_centric(Qb, data, "worst")
        sl["fold"] = data.site_fold[sl.site_idx.values]
        levels[b] = sl
        pk = sl.groupby("kinase_acc").hit5.mean()
        r = {"beta": b, "M2": E._nanmean(pk.reindex([k for k in kp if k in pk.index]))}
        for bn in E.TRAIN_BIN_LABELS:
            ks = [k for k in kp if k in pk.index and tbin[k] == bn]
            r[f"M2_bin_{bn}"] = E._nanmean(pk.reindex(ks)) if ks else np.nan
            r[f"n_bin_{bn}"] = len(ks)
        r["MRR"] = E._nanmean(sl.groupby("kinase_acc").rr.mean().reindex([k for k in kp if k in pk.index]))
        r["hit1"] = E._nanmean(sl.groupby("kinase_acc").hit1.mean().reindex([k for k in kp if k in pk.index]))
        rows.append(r)
        print(r, flush=True)
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "grid.tsv", sep="\t", index=False, float_format="%.5f")
    # ---- nested leave-one-fold-out

    def m2_sub(sl, fs):
        x = sl[sl.fold.isin(fs)]
        pk = x.groupby("kinase_acc").hit5.mean()
        return E._nanmean(pk.reindex([k for k in kp if k in pk.index]))
    nest = []
    for f in range(5):
        inner = [g for g in range(5) if g != f]
        sc_in = {b: m2_sub(levels[b], inner) for b in BETAS}
        bsel = max(BETAS, key=lambda b: (sc_in[b], -b))
        nest.append({"held_out_fold": f, "beta_selected": bsel, "M2_inner_at_sel": sc_in[bsel],
                     "M2_heldout_at_sel": m2_sub(levels[bsel], [f]), "M2_heldout_beta0": m2_sub(levels[0.0], [f])})
    N = pd.DataFrame(nest)
    N.to_csv(OUT / "nested.tsv", sep="\t", index=False, float_format="%.5f")
    print(N.to_string(), flush=True)
    # ---- selection + A32(b) rule
    g0 = G[G.beta.eq(0.0)].iloc[0]
    best = G.loc[G.M2.idxmax()]
    bstar = float(best.beta)
    adopt = bool(bstar > 0 and best.M2 >= g0.M2 - 0.003 and best["M2_bin_5-19"] > g0["M2_bin_5-19"])
    # informational: largest beta within the 0.003 M2 tolerance that improves the 5-19 bin (NOT the A32b selection)
    tol = G[(G.M2 >= g0.M2 - 0.003) & (G["M2_bin_5-19"] > g0["M2_bin_5-19"])]
    dec = {"beta_star_argmax_pooled_M2": bstar, "M2_beta0": float(g0.M2), "M2_beta_star": float(best.M2),
           "M2_5_19_beta0": float(g0["M2_bin_5-19"]), "M2_5_19_beta_star": float(best["M2_bin_5-19"]),
           "adopt_beta_gt_0": adopt,
           "check_beta0_vs_run_metrics_M2": None,
           "informational_largest_beta_within_tolerance": float(tol.beta.max()) if len(tol) else None,
           "nested_mean_heldout_M2_selected": float(N.M2_heldout_at_sel.mean()),
           "nested_mean_heldout_M2_beta0": float(N.M2_heldout_beta0.mean())}
    m = json.loads((ROOT / "results" / "runs" / RUN / "metrics.json").read_text())["endpoints"]["M2_hit5"]
    dec["check_beta0_vs_run_metrics_M2"] = float(abs(m - g0.M2))
    # full evaluate at the selected beta (all endpoints + canary), and at beta 0 for reference
    for b in sorted({0.0, bstar}):
        pr_all = prior_used[P.k.values, np.clip(fold, 0, 4)]
        Pb = pd.DataFrame({"kinase_idx": P.k.values, "site_idx": P.s.values, "score": P.score.values,
                           "score_site": P.score_site.values - b * pr_all})
        res = E.evaluate(Pb, data, missing="worst", keep_items=False)
        e = res.endpoints
        dec[f"full_eval_beta_{b}"] = {"M1": e["M1_auroc_hard"], "M2": e["M2_hit5"], "M3": e["M3_top1_recall"],
                                      "M4": e["M4_auroc_hard_Kany"], "neg_canary_hard": res.canary.get("neg_mean_auroc_hard"),
                                      "pos_canary_hard": res.canary.get("pos_auroc_hard")}
        del Pb, pr_all
    dec["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "decision.json").write_text(json.dumps(dec, indent=1, default=float))
    print(json.dumps(dec, indent=1, default=float))


if __name__ == "__main__":
    main()
