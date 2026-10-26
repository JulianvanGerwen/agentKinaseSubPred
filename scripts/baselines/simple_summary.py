"""Compile results/baselines/simple/summary.md from the per-method json files written by simple_run.py.
All numbers are read from ksp.eval outputs (endpoints, bootstrap, strata); nothing is recomputed here."""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "baselines" / "simple"

ORDER = [("b0_prior", "B0-prior (kinase training count)"),
         ("b0_pop", "B0-pop (site popularity; diagnostic)"),
         ("b0_miss", "B0-miss (LightGBM on missingness only)"),
         ("b1_published", "B1 (i) Cantley as-published pct"),
         ("b1_universe", "B1 (ii) Cantley universe-recalibrated pct"),
         ("b1_primed_hard", "B1 (iii) Cantley priming-aware hard (universe pct)"),
         ("b1_primed_soft", "B1 (iii) Cantley priming-aware soft (universe pct)"),
         ("b1_primed_hard_published", "B1 suppl.: primed hard, published pct"),
         ("b1_primed_soft_published", "B1 suppl.: primed soft, published pct"),
         ("b1_imputed_universe", "B1 suppl.: universe pct + nearest-atlas imputed"),
         ("b1_imputed_primed_soft", "B1 suppl.: primed soft + nearest-atlas imputed"),
         ("b1_imputed_primed_hard", "B1 suppl.: primed hard + nearest-atlas imputed"),
         ("b2_pssm", "B2 learned PSSM (F2a, family prior)"),
         ("b2_pssm_ph", "B2 learned PSSM + phospho channel"),
         ("b7_lr", "B7 simple strong (LR)")]


def md(df: pd.DataFrame, index: bool = False) -> str:
    """Minimal markdown table (no tabulate dependency)."""
    d = df.reset_index() if index else df
    cols = [str(c) for c in d.columns]
    fmt = lambda v: "" if (isinstance(v, float) and pd.isna(v)) else str(v)  # noqa: E731
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    lines += ["| " + " | ".join(fmt(v) for v in r) + " |" for r in d.itertuples(index=False)]
    return "\n".join(lines)


def ci(boot, e):
    r = [b for b in boot if b["endpoint"] == e][0]
    return f"{r['point']:.3f} [{r['lo95']:.3f}, {r['hi95']:.3f}]"


def main():
    L = [f"# Simple baselines on the 5-fold DEV CV (B0 / B1 / B2 / B7)\n",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `scripts/baselines/simple_summary.py` from "
         "`results/baselines/simple/<method>.json` (each = one `ksp.manifest` run; metrics only from `ksp.eval`).\n",
         "Pooled out-of-fold predictions over **all** evaluation candidates (every kinase with >=1 DEV positive x "
         "every DEV U_ms site of its class, plus all real kinases x labelled sites for M2). "
         "Missing pairs = 'worst' (ties at the bottom). 95% CI: cluster bootstrap over substrate clusters, B=1000 "
         "(`ksp.eval.bootstrap`). M1 = M5 on CV (K_cv, E-hard); M4 = CV analogue over K_any_cv (>=1 DEV positive). "
         "M2 = filtered site-centric Hit@5 (A18), M3 = recall in top 1% of E-all.\n"]
    rows, cov = [], []
    for m, lab in ORDER:
        f = OUT / f"{m}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        b = d["bootstrap"]
        E = d["endpoints"]
        c = d["coverage"]
        rows.append({"method": f"`{m}`", "description": lab, "M1 (=M5) AUROC E-hard": ci(b, "M1"),
                     "M2 Hit@5": ci(b, "M2"), "M3 top-1% recall": ci(b, "M3"), "M4-CV AUROC (K_any)": ci(b, "M4"),
                     "MRR": f"{E['sec_mrr']:.3f}", "AUROC E-all": f"{E['sec_auroc_all']:.3f}",
                     "fit+predict s": f"{d['runtime_fit_predict_s']:.0f}", "total s": f"{d['runtime_total_s']:.0f}"})
        r = {"method": f"`{m}`", "kinases scored / with DEV pos": f"{c['kinases_scored']}/{c['kinases_with_positives']}",
             "K_cv scored": f"{c['K_cv_scored']}/{c['K_cv']}",
             "required pairs scored": f"{c['frac_required_pairs_scored']:.3f}"}
        if "endpoints_covered_only" in d:
            ec = d["endpoints_covered_only"]
            r.update({"covered-only M1": f"{ec['M1_auroc_hard']:.3f} (n={ec['n_K_primary_evaluated']})",
                      "covered-only M2": f"{ec['M2_hit5']:.3f}", "covered-only M3": f"{ec['M3_top1_recall']:.3f}"})
        cov.append(r)
    L.append("## Endpoints (CV, pooled OOF; point [95% CI])\n")
    L.append(md(pd.DataFrame(rows)) + "\n")
    L.append("## Coverage\n")
    L.append("Covered-only = `missing='drop'` restricted to kinases with >=1 scored pair (not a matched-universe "
             "comparison; head-to-heads must use `ksp.eval.matched_universe`).\n")
    L.append(md(pd.DataFrame(cov).fillna("")) + "\n")
    bf = OUT / "b1_best_arm.json"
    if bf.exists():
        bb = json.loads(bf.read_text())
        L.append("## Official B1 = Cantley best arm (A19; chosen on CV M1)\n")
        L.append(f"Rule: {bb['rule']}. **Best: `{bb['best_method']}`** (`{bb['store']}`; run `{bb['run_id']}`).\n")
        L.append(md(pd.DataFrame(bb["candidates"])[["method", "M1", "M2", "M3", "M4"]].round(4))
                 + "\n")
        pb = pd.DataFrame(bb["paired_bootstrap_best_vs_other"])
        L.append("Paired cluster bootstrap (B=1000), best minus other:\n")
        L.append(md(pb[["comparison", "endpoint", "diff", "lo95", "hi95", "p_boot"]].round(4))
                 + "\n")
    can = []
    for m in ["b2_pssm", "b2_pssm_ph", "b7_lr", "b0_miss"]:
        f = OUT / f"{m}_canary.json"
        if f.exists():
            c = json.loads(f.read_text())["canary"]
            can.append({"method": f"`{m}`", "neg mean AUROC E-hard": round(c["neg_mean_auroc_hard"], 3),
                        "neg mean AUROC E-all": round(c["neg_mean_auroc_all"], 3),
                        "pos AUROC E-hard": round(c["pos_auroc_hard"], 3),
                        "pos AUROC E-all": round(c["pos_auroc_all"], 3), "neg pass (<=0.57)": c["neg_pass"],
                        "pos pass (>=0.90)": c["pos_pass"]})
    if can:
        L.append("## Canaries (trained with `load_dev(with_canaries=True)`; `ksp.eval.canary_check`)\n")
        L.append(md(pd.DataFrame(can)) + "\n")
    # strata: M1 per training-label bin and residue class over K_any
    st = []
    for m, _ in ORDER:
        f = OUT / f"{m}.json"
        if not f.exists():
            continue
        s = pd.DataFrame(json.loads(f.read_text())["strata"])
        s = s[s.stratum.isin(["train_bin", "residue_class", "priming", "pikk"])]
        for r in s.itertuples():
            st.append({"method": m, "stratum": f"{r.stratum}={r.value}", "auroc_hard": r.auroc_hard,
                       "n": r.n_kinases})
    if st:
        S = pd.DataFrame(st)
        piv = S.pivot(index="method", columns="stratum", values="auroc_hard").reindex([m for m, _ in ORDER
                                                                                     if m in set(S.method)])
        n = S.drop_duplicates("stratum").set_index("stratum").n
        piv.columns = [f"{c} (n={n[c]})" for c in piv.columns]
        L.append("## Strata: macro AUROC (E-hard) over K_any kinases\n")
        L.append(md(piv.round(3), index=True) + "\n")
        L.append("Full strata (incl. manning group, M2/M3 per stratum, pair-level evidence/temporal/external strata) "
                 "are in each `<method>.json` (`strata`, `strata_K_cv`, `pair_strata`).\n")
    nf = OUT / "notes.md"
    if nf.exists():
        L.append(nf.read_text())
    (OUT / "summary.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
