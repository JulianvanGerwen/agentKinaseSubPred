#!/usr/bin/env python
"""DESIGN §17 A34(b), post-TEST exploratory (CV only): kinase-gated blend with the priming-aware Cantley arm.

For the curated priming kinases (ksp.strata.priming_curated: GSK3A/B, CK1 family, CK2 alpha, PLK1) only:
    score' = (1 - gamma) * z(fusion score) + gamma * z(Cantley primed-hard score)
z = standardisation within (CV fold, kinase) over that kinase's scored candidate rows. Fusion = final OOF fusion
(r2_fusion_r3__20261026-191901__s0, score channel, w 0.6). Cantley primed-hard = B1 arm simple_b1_primed_hard_reeval
(kinase-library pp=True with U_ms flanks lower-cased; label-free). Rows without a Cantley score keep z(fusion).
All other kinases, and score_site (M2), are unchanged.
Grid    : gamma in {0, 0.1, ..., 1}.
Nested  : outer fold f: gamma*_f = argmax over the grid of the mean, over the 4 other folds g, of the per-fold
          priming-kinase M1 (macro over K_cv priming kinases of ksp.eval per-fold E-hard AUROC 'auroc_hard_f{g}',
          A24-normalised pooled evaluation); ties -> smaller gamma. Fold f rows use gamma*_f.
Report  : overall M1/M2/M3 (full ksp.eval.evaluate, A24), priming-kinase M1, primed / unprimed positive strata
          (ksp.strata 'primed', curated offsets) with paired cluster-bootstrap CIs (B = 1000, seed 0) vs the fusion.
A34 rule: promising if the primed-stratum nested gain has CI > 0 and overall M1 loses <= 0.002.
Out     : results/models/posttest/priming_gate/*
"""
import dataclasses
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
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import strata as ST  # noqa: E402
import strata_cv_report as R  # noqa: E402

FUSION = "r2_fusion_r3__20261026-191901__s0"
CANT = "simple_b1_primed_hard_reeval__20261026-155750__s0"
GAMMAS = [round(0.1 * i, 1) for i in range(11)]
OUT = ROOT / "results" / "models" / "posttest" / "priming_gate"


def zscore(v, key):
    s = pd.Series(v)
    g = s.groupby(key)
    sd = g.transform("std").values
    return np.where(sd > 0, (v - g.transform("mean").values) / np.where(sd > 0, sd, 1), 0.0)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    data = E.prepare("cv", fold=None)
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].reset_index(drop=True)
    data = dataclasses.replace(data, pos_meta=pm,
                               positives={k: v for k, v in data.positives.items() if not k.startswith("CANARY_")})
    kacc = data.kinases.uniprot_acc.values
    real = ~pd.Series(kacc).str.startswith("CANARY_").values
    P = E.map_predictions(pd.read_parquet(ROOT / "cache" / "runs" / FUSION / "predictions.parquet"), data)
    P = P[real[P.k.values]].reset_index(drop=True)
    cur = ST.priming_curated(data.kinases[real])
    pk_idx = np.array(sorted(i for i, a in enumerate(kacc) if a in cur))
    gene = dict(zip(kacc, data.kinases.gene.values))
    m_pr = np.isin(P.k.values, pk_idx)
    Pp = P[m_pr].reset_index(drop=True)
    C = E.map_predictions(pd.read_parquet(ROOT / "cache" / "runs" / CANT / "predictions.parquet"), data)
    C = C[np.isin(C.k.values, pk_idx)]
    key_p = Pp.k.values.astype(np.int64) * len(data.sites) + Pp.s.values
    key_c = C.k.values.astype(np.int64) * len(data.sites) + C.s.values
    cs = pd.Series(C.score.values, index=key_c).reindex(key_p).values
    fold = data.site_fold[Pp.s.values]
    gk = Pp.k.values * 8 + fold
    zf = zscore(Pp.score.values, gk)
    has_c = np.isfinite(cs)
    zc = np.where(has_c, zscore(np.where(has_c, cs, 0.0), np.where(has_c, gk, -1)), zf)
    cov = {gene[kacc[i]]: float(has_c[Pp.k.values == i].mean()) for i in pk_idx}
    print("Cantley coverage per priming kinase:", cov, flush=True)
    del C
    kprim_cv = [k for k in data.k_primary if k in cur]
    # ---- grid: per-fold priming-kinase M1 (restricted evaluation: priming kinases only, their pairs only)
    dpr = dataclasses.replace(data, pos_meta=pm[pm.kinase_acc.isin(cur)].reset_index(drop=True),
                              positives={k: v for k, v in data.positives.items() if k in cur})
    rows, perk = [], {}
    for gm in GAMMAS:
        sc = (1 - gm) * zf + gm * zc
        res = E.evaluate(pd.DataFrame({"kinase_idx": Pp.k.values, "site_idx": Pp.s.values, "score": sc}), dpr,
                         missing="worst", kinases=sorted(dpr.positives), keep_items=False)
        pk = res.per_kinase.set_index("kinase_acc").loc[kprim_cv]
        r = {"gamma": gm, "M1_priming_pooled": float(pk.auroc_hard.mean())}
        for f in range(5):
            r[f"M1_priming_f{f}"] = float(pk[f"auroc_hard_f{f}"].mean())
        rows.append(r)
        perk[gm] = pk.auroc_hard
        print("grid", {k: round(v, 4) for k, v in r.items()}, flush=True)
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "grid.tsv", sep="\t", index=False, float_format="%.5f")
    pd.DataFrame(perk).rename(index=gene).to_csv(OUT / "grid_per_kinase.tsv", sep="\t", float_format="%.4f")
    nest = []
    gsel = {}
    for f in range(5):
        inner = {gm: float(np.mean([G.loc[G.gamma.eq(gm), f"M1_priming_f{g}"].iloc[0] for g in range(5) if g != f]))
                 for gm in GAMMAS}
        gsel[f] = max(GAMMAS, key=lambda x: (inner[x], -x))
        nest.append({"outer_fold": f, "gamma_selected": gsel[f], "inner_M1_at_sel": inner[gsel[f]],
                     "inner_M1_gamma0": inner[0.0]})
    N = pd.DataFrame(nest)
    N.to_csv(OUT / "nested_selection.tsv", sep="\t", index=False, float_format="%.5f")
    print(N.to_string(), flush=True)
    gam_row = np.array([gsel[f] if f >= 0 else 0.0 for f in fold])
    sc_nested = (1 - gam_row) * zf + gam_row * zc
    # ---- full evaluations: base fusion vs nested gated blend (all kinases)
    res = {}
    for tag in ["base", "nested"]:
        s = P.score.values.copy()
        if tag == "nested":
            s[np.where(m_pr)[0]] = sc_nested
        df = pd.DataFrame({"kinase_idx": P.k.values, "site_idx": P.s.values, "score": s,
                           "score_site": P.score_site.values})
        res[tag] = E.evaluate(df, data, missing="worst", keep_items=True)
        del df, s
    out = {"fusion": FUSION, "cantley": CANT, "gamma_selected": gsel, "cantley_coverage": cov}
    for t, r in res.items():
        e = r.endpoints
        out[t] = {"M1": e["M1_auroc_hard"], "M1_fold_mean": e.get("M1_fold_mean"), "M2": e["M2_hit5"],
                  "M3": e["M3_top1_recall"], "M4": e["M4_auroc_hard_Kany"]}
    out["paired_overall"] = E.paired_bootstrap(res["nested"], res["base"], B=1000, seed=0).to_dict(orient="records")
    # strata (positive subsets) with paired bootstrap on identical cluster weights
    ktab = ST.kinase_table()
    site_attr, _, _ = R.site_attributes(data, ktab)
    cl = E._universe_clusters(res["base"])
    srows = []
    specs = [("priming kinases (all positives)", None, None)] + [(f"{s}={v}", s, v) for s, v in
                                                                [("primed", "primed"), ("primed", "unprimed"),
                                                                 ("primed_cantley", "primed"),
                                                                 ("primed_cantley", "unprimed")]]
    for name, s, v in specs:
        subs, mats = {}, {}
        for t in res:
            if s is None:
                keep = {k: np.ones(len(data.positives[k]), bool) for k in cur if k in data.positives}
            else:
                kin_ok, attr = site_attr[(s, v)]
                keep = {k: attr[k] == v for k in kin_ok if k in data.positives}
            sub = R.subset_result(res[t], data, keep)
            spk = sub.per_kinase.set_index("kinase_acc")
            kp = [k for k in data.k_primary if k in spk.index]
            subs[t] = (spk, kp)
            mats[t] = E._boot_matrix(sub, cl, R.B_BOOT, R.SEED, kp)
        kp = subs["base"][1]
        pt = {t: E._nanmean(subs[t][0].loc[kp, "auroc_hard"]) for t in res}
        d = np.nanmean(mats["nested"]["M1"], 1) - np.nanmean(mats["base"]["M1"], 1)
        srows.append({"stratum": name, "n_kinases": len(kp), "n_pos": int(subs["base"][0].loc[kp, "n_pos"].sum()),
                      "M1_base": pt["base"], "M1_nested": pt["nested"], "delta": pt["nested"] - pt["base"],
                      "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
                      "per_kinase": {gene[k]: [round(float(subs["base"][0].loc[k, "auroc_hard"]), 4),
                                               round(float(subs["nested"][0].loc[k, "auroc_hard"]), 4)] for k in kp}})
    out["strata"] = srows
    pd.DataFrame([{k: v for k, v in r.items() if k != "per_kinase"} for r in srows]).to_csv(
        OUT / "strata_nested_vs_base.tsv", sep="\t", index=False, float_format="%.5f")
    # sanity: non-priming kinases' per-kinase M1 identical; base reproduces the run's metrics
    pb, pn = res["base"].per_kinase.set_index("kinase_acc"), res["nested"].per_kinase.set_index("kinase_acc")
    other = [k for k in pb.index if k not in cur]
    m_run = json.loads((ROOT / "results" / "runs" / FUSION / "metrics.json").read_text())["endpoints"]
    out["sanity"] = {
        "max_abs_diff_nonpriming_kinase_auroc": float(np.nanmax(np.abs(pb.loc[other, "auroc_hard"] - pn.loc[other, "auroc_hard"]))),
        "base_M1_vs_run_metrics": abs(out["base"]["M1"] - m_run["M1_auroc_hard"]),
        "base_M2_vs_run_metrics": abs(out["base"]["M2"] - m_run["M2_hit5"]),
        "M2_unchanged": abs(out["base"]["M2"] - out["nested"]["M2"]),
        "gamma0_grid_equals_base_priming_M1": abs(G.loc[G.gamma.eq(0.0), "M1_priming_pooled"].iloc[0] -
                                                    E._nanmean(pb.loc[kprim_cv, "auroc_hard"]))}
    tgt = srows[1]
    dM1 = out["nested"]["M1"] - out["base"]["M1"]
    out.update({"target": "primed=primed M1 (K_cv)", "target_delta": tgt["delta"], "target_ci": [tgt["lo"], tgt["hi"]],
                "M1_delta": dM1, "promising_A34": bool(tgt["lo"] > 0 and dM1 >= -0.002),
                "runtime_s": round(time.time() - t0, 1)})
    (OUT / "summary.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: v for k, v in out.items() if k != "strata"}, indent=1, default=float))
    for r in srows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
