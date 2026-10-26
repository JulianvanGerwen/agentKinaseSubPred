"""Collect pooled-CV summaries of every regime / variant run into results/models/regime/summary.tsv and the
per-evidence-class pair strata into results/models/regime/evidence_strata.tsv.
Optional --paired: paired cluster bootstrap (ksp.eval.paired_bootstrap, B=1000) of each run vs the L2 baseline."""
import argparse, json
import numpy as np, pandas as pd
from ksp.index import ROOT

ap = argparse.ArgumentParser()
ap.add_argument("--paired", action="store_true")
a = ap.parse_args()
RUNS = ROOT / "results" / "runs"
OUTD = ROOT / "results" / "models" / "regime"
rows = []
for st in sorted((OUTD / "state").glob("*.json")):
    s = json.loads(st.read_text())
    if "pooled" in s:
        rows.append({"config": st.stem, "run_id": s["pooled"], **s["summary"],
                     "train_wall_s": sum(v.get("t_fold_s", 0) for v in s.get("logs", {}).values())})
for p in sorted(RUNS.glob("regime_*_cv5__*")):
    if any(r["run_id"] == p.name for r in rows):
        continue
    sj = p / "summary.json"
    if sj.exists():
        rows.append({"config": p.name.split("__")[0], "run_id": p.name, **json.loads(sj.read_text())})
T = pd.DataFrame(rows)
T.to_csv(OUTD / "summary.tsv", sep="\t", index=False)
cols = ["config", "M1", "M1_fold_mean", "M2", "M3", "M4"] + [c for c in T.columns if c.startswith("M5bin_")] + \
       ["M1_ST", "M1_Y", "canary_neg_hard", "canary_pos", "canary_gate"]
print(T[[c for c in cols if c in T]].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
ev = []
for _, r in T.iterrows():
    ps = pd.read_csv(RUNS / r.run_id / "pair_strata.tsv", sep="\t")
    ps = ps[ps.stratum.str.startswith("evidence=") | ps.stratum.eq("external_curated")]
    for _, x in ps.iterrows():
        ev.append({"config": r.config, "stratum": x.stratum, "kinase_set": x.kinase_set, "n_kinases": x.n_kinases,
                   "n_pos": x.n_pos, "auroc_hard": x.auroc_hard, "top1_recall": x.top1_recall})
E = pd.DataFrame(ev)
E.to_csv(OUTD / "evidence_strata.tsv", sep="\t", index=False)
print(E[E.kinase_set == "K_primary"].pivot(index="config", columns="stratum", values="auroc_hard").to_string(
    float_format=lambda x: f"{x:.4f}"))

# markdown table (M1 / fold-mean / M2 / M3 / M4 / M5 by bin / ST / Y / canary) with deltas vs the matched L2 baseline
base = T[T.config == "regime_L2_L2fixed_s0"].iloc[0]
order = ["regime_L2_L2fixed_s0", "regime_L1_refit_s0", "regime_L1_L2fixed_s0", "regime_L2_vitro_refit_s0",
         "regime_L2_vitro_L2fixed_s0", "regime_L2_vivo_refit_s0", "regime_L2ext_refit_s0",
         "regime_perkinase_lgb_N20_qmap_cv5", "regime_perkinase_lr_N20_qmap_cv5", "regime_ppool_lam300_cv5",
         "regime_resboost_lam300_cv5", "regime_pergroup_raw_cv5"]
T2 = T.set_index("config")
lines = ["| config | M1 | ΔM1 | fold-mean | M2 | ΔM2 | M3 | M4 | M5 1-4 | 5-19 | 20-99 | ≥100 | M1 ST | M1 Y | canary neg/pos | adopt? |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
f = lambda x: f"{x:.4f}" if pd.notna(x) else ""
for c in order + [c for c in T2.index if c not in order]:
    if c not in T2.index:
        continue
    r = T2.loc[c]
    sd = c.rsplit("_s", 1)[-1] if c.rsplit("_s", 1)[-1] in ("1", "2", "3", "4") else "0"
    if sd == "0":
        bm1, bm2 = base.M1, base.M2
    else:            # matched-seed R2 baseline (frozen final config, same seed; R2's own run - not retrained)
        rid = [p.name for p in RUNS.glob(f"r2_final_s{sd}_cv5__*")][0]
        e = json.loads((RUNS / rid / "metrics.json").read_text())["endpoints"]
        bm1, bm2 = e["M1_auroc_hard"], e["M2_hit5"]
    dm1, dm2 = r.M1 - bm1, r.M2 - bm2
    adopt = "baseline" if c == "regime_L2_L2fixed_s0" else (
        "YES" if (dm1 >= 0.005 and dm2 >= -0.01 and bool(r.canary_gate)) else "no")
    lines.append(f"| {c} | {f(r.M1)} | {dm1:+.4f} | {f(r.M1_fold_mean)} | {f(r.M2)} | {dm2:+.4f} | {f(r.M3)} | {f(r.M4)} | "
                 f"{f(r['M5bin_1-4'])} | {f(r['M5bin_5-19'])} | {f(r['M5bin_20-99'])} | {f(r['M5bin_>=100'])} | "
                 f"{f(r.M1_ST)} | {f(r.M1_Y)} | {r.canary_neg_hard:.3f}/{r.canary_pos:.3f} {'PASS' if r.canary_gate else 'FAIL'} | {adopt} |")
(OUTD / "summary_table.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
