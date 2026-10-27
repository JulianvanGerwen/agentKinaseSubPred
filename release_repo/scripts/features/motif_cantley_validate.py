"""Reproduction check of the Cantley scoring engine (ksp.features.motif_cantley) against published outputs.

1. kinase-library 1.8.0 scored Ochoa reference parquets (82,735 S/T x 311; 5,431 Y x 93): log2 scores.
2. Johnson 2023 Suppl. Table 3 (89,784 S/T sites x 303 kinases): percentiles.
3. Yaron-Barir 2024 Suppl. Table 3, non-canonical sheet (7,315 Y sites x 93 kinases): percentiles.
4. Johnson 2023 Suppl. Table 4 (ATM percentiles of IR-induced sites).
Writes results/features/motif/reproduction_{summary.json,per_kinase.tsv}.
Supp tables are first extracted by scripts/features/motif_extract_supp.py into cache/features/motif/validation/.
"""
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.features import motif_cantley as mc  # noqa: E402

RES = ROOT / "results/features/motif"
VAL = ROOT / "cache/features/motif/validation"
summary, per_k = {}, []
t0 = time.time()

# ---- 1. scores vs kinase-library scored reference ----
for t in mc.KIN_TYPES:
    A = mc.load_atlas(t)
    R = mc.load_reference(t)
    w = list(R.index)
    S = mc.score(A, mc.encode(w), mc.center_residue(w))
    raw = mc.score(A, mc.encode(w), mc.center_residue(w), published_rounding=False)
    ref = R[A.names].to_numpy(float)
    d = np.abs(S - ref)
    summary[f"kl_reference_scores_{t}"] = {
        "n_sites": len(w), "n_kinases": len(A.names), "n_values": int(d.size),
        "exact_match_3dp": float((d < 1e-9).mean()), "within_0.001": float((d <= 0.001 + 1e-9).mean()),
        "max_abs_diff": float(d.max()), "max_abs_diff_unrounded": float(np.abs(raw - ref).max()),
        "pearson": float(np.corrcoef(S.ravel(), ref.ravel())[0, 1])}
    P_ref_self = mc.percentile(S, ref)
    P_pub_self = mc.percentile(ref, ref)
    dp = np.abs(P_ref_self - P_pub_self)
    summary[f"kl_reference_scores_{t}"]["percentile_exact"] = float((dp < 1e-9).mean())
    summary[f"kl_reference_scores_{t}"]["percentile_max_abs_diff"] = float(dp.max())


def compare_pct(tab, t, label):
    A = mc.load_atlas(t)
    ref = mc.load_reference(t)[A.names].to_numpy(float)
    w = tab["SITE_+/-7_AA"].astype(str).tolist()
    ok = np.array([len(x) == 15 and x[7].upper() in mc.TYPE_RES[t] for x in w])
    tab = tab[ok].reset_index(drop=True)
    w = [x for x, o in zip(w, ok) if o]
    lower_flanks = sum(any(c in "sty" for i, c in enumerate(x) if i != 7) for x in w)
    codes = mc.encode(w, keep_lower_flanks=True)   # tables give unprimed windows; keep as-is if any
    P = mc.percentile(mc.score(A, codes, mc.center_residue(w)), ref)
    kcols = [c[:-len("_percentile")] for c in tab.columns if c.endswith("_percentile") and c != "median_percentile"]
    common = [k for k in kcols if k.upper() in A.names]
    missing = [k for k in kcols if k.upper() not in A.names]
    pub = tab[[k + "_percentile" for k in common]].to_numpy(float)
    mine = P[:, [A.names.index(k.upper()) for k in common]]
    d = np.abs(pub - mine)
    for j, k in enumerate(common):
        rho = spearmanr(pub[:, j], mine[:, j]).correlation
        per_k.append(dict(table=label, matrix=k, n=len(w), exact=float((d[:, j] < 1e-6).mean()),
                          within_0p1=float((d[:, j] <= 0.1 + 1e-9).mean()), within_1=float((d[:, j] <= 1 + 1e-9).mean()),
                          max_abs=float(d[:, j].max()), spearman=float(rho)))
    pk = pd.DataFrame([r for r in per_k if r["table"] == label])
    good = pk.spearman > 0.999
    summary[label] = {"n_sites": len(w), "n_kinases_compared": len(common), "kinases_not_in_1.8.0": missing,
                      "windows_with_lowercase_flanks": lower_flanks, "n_values": int(d.size),
                      "exact_match_2dp": float((d < 1e-6).mean()), "within_0.1": float((d <= 0.1 + 1e-9).mean()),
                      "within_1.0": float((d <= 1 + 1e-9).mean()),
                      "pearson_all": float(np.corrcoef(pub.ravel(), mine.ravel())[0, 1]),
                      "median_per_kinase_spearman": float(pk.spearman.median()),
                      "n_kinases_spearman_gt_0.999": int(good.sum()),
                      "kinases_spearman_lt_0.999 (matrix revised after publication)": pk[~good].sort_values("spearman")[["matrix", "spearman"]].round(4).values.tolist(),
                      "exact_match_2dp_unrevised_kinases": float(np.mean([(d[:, j] < 1e-6).mean() for j in np.where(good)[0]])),
                      "within_0.1_unrevised_kinases": float(np.mean([(d[:, j] <= 0.1 + 1e-9).mean() for j in np.where(good)[0]]))}


compare_pct(pd.read_parquet(VAL / "johnson2023_st3.parquet"), "ser_thr", "johnson2023_ST3")
compare_pct(pd.read_parquet(VAL / "yaronbarir2024_st3_noncanonical.parquet"), "tyrosine", "yaronbarir2024_ST3")

# ---- 4. Johnson ST4 (ATM) ----
try:
    import openpyxl
    wb = openpyxl.load_workbook(ROOT / "external/cantley_johnson2023/41586_2022_5575_MOESM6_ESM.xlsx", read_only=True)
    rows = list(wb["Supplementary Table 4 Data"].iter_rows(values_only=True))
    t4 = pd.DataFrame(rows[1:], columns=rows[0]).dropna(subset=["SITE_+/-7_AA"])
    A = mc.load_atlas("ser_thr")
    ref = mc.load_reference("ser_thr")[A.names].to_numpy(float)
    w = t4["SITE_+/-7_AA"].tolist()
    P = mc.percentile(mc.score(A, mc.encode(w), mc.center_residue(w)), ref)[:, A.names.index("ATM")]
    pub = t4["ATM percentile score"].astype(float).to_numpy()
    summary["johnson2023_ST4_ATM"] = {"n": len(w), "max_abs_diff": float(np.abs(P - pub).max()),
                                     "within_0.1": float((np.abs(P - pub) <= 0.1).mean()),
                                     "note": "ST4 percentiles given to 3 dp; ours 2 dp"}
except Exception as e:  # pragma: no cover
    summary["johnson2023_ST4_ATM"] = {"error": repr(e)}

summary["runtime_s"] = round(time.time() - t0, 1)
pd.DataFrame(per_k).to_csv(RES / "reproduction_per_kinase.tsv", sep="\t", index=False, float_format="%.5g")
(RES / "reproduction_summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
