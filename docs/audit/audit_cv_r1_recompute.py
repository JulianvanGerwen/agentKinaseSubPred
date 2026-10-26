"""Independent audit recomputation of CV headline metrics from saved predictions with ksp.eval.
Usage: python -I recompute.py <spec_name> <out_dir>
Read-only on /workspace; writes only to out_dir (scratchpad)."""
import json
import sys
import time

sys.path.insert(0, "/workspace/src")
import numpy as np
import pandas as pd

from ksp import eval as kev, splits, index, canary

C = "/workspace/cache/runs/"
R = "/workspace/results/runs/"

R3F = ["r3_reg1x5_w7_nce05_cpt_aux_fold0__20261026-155231__s0", "r3_reg1x5_w7_nce05_cpt_aux_fold1__20261026-160057__s0",
       "r3_reg1x5_w7_nce05_cpt_aux_fold2__20261026-160712__s0", "r3_reg1x5_w7_nce05_cpt_aux_fold3__20261026-161523__s0",
       "r3_reg1x5_w7_nce05_cpt_aux_fold4__20261026-161009__s0"]
R3B = ['r3_reg1x5_w7_nce05_fold0__20261026-154250__s0', 'r3_reg1x5_w7_nce05_fold1__20261026-154825__s0',
       'r3_reg1x5_w7_nce05_fold2__20261026-155158__s0', 'r3_reg1x5_w7_nce05_fold3__20261026-155528__s0',
       'r3_reg1x5_w7_nce05_fold4__20261026-155809__s0']
R3L = ["r3_reg1x5_w7_nce05_cpt_aux_lko0__20261026-162151__s0", "r3_reg1x5_w7_nce05_cpt_aux_lko1__20261026-162828__s0",
       "r3_reg1x5_w7_nce05_cpt_aux_lko2__20261026-161641__s0", "r3_reg1x5_w7_nce05_cpt_aux_lko3__20261026-162328__s0",
       "r3_reg1x5_w7_nce05_cpt_aux_lko4__20261026-163141__s0"]
R2FULL = ['r2_a_l127_lr03_v2__20261026-153808__s0', 'r2_full_f1__20261026-160659__s0', 'r2_full_f2__20261026-161233__s0',
          'r2_full_f3__20261026-161712__s0', 'r2_full_f4__20261026-162043__s0']
R2G2 = ['r2_g2_main_f0__20261026-155427__s0', 'r2_g2_main_f1__20261026-155642__s0', 'r2_g2_main_f2__20261026-160129__s0',
        'r2_g2_main_f3__20261026-160725__s0', 'r2_g2_main_f4__20261026-161638__s0']

SIMPLE = {m: r for m, r in [
    ("b0_prior", "simple_b0_prior__20261026-150925__s0"), ("b0_pop", "simple_b0_pop__20261026-145746__s0"),
    ("b0_miss", "simple_b0_miss__20261026-150214__s0"), ("b1_published", "simple_b1_published__20261026-150119__s0"),
    ("b1_universe", "simple_b1_universe__20261026-150125__s0"),
    ("b1_imputed_universe", "simple_b1_imputed_universe__20261026-150833__s0"),
    ("b2_pssm", "simple_b2_pssm__20261026-150923__s0"), ("b7_lr", "simple_b7_lr__20261026-153330__s0")]}

SPECS = {
    "r3_final_cv": dict(fold_runs=R3F, mode="cv", reported=R + "r3_w7_nce05_cpt_aux_cv_pooled__20261026-162205__s0/metrics.json"),
    "r3_base_cv": dict(fold_runs=R3B, mode="cv", reported=R + "r3_w7_nce05_cv_pooled__20261026-160122__s0/metrics.json"),
    "r3_final_lko": dict(fold_runs=R3L, mode="lko", reported=R + "r3_w7_nce05_cpt_aux_lko_pooled__20261026-164851__s0/metrics.json"),
    "r2_full": dict(fold_runs=R2FULL, mode="cv", reported=R + "r2_full_cv5__20261026-162445__s0/metrics.json",
                    pooled=C + "r2_full_cv5__20261026-162445__s0/predictions.parquet"),
    "r2_g2_main": dict(fold_runs=R2G2, mode="cv", reported=R + "r2_g2_main_cv5__20261026-162249__s0/metrics.json",
                       pooled=C + "r2_g2_main_cv5__20261026-162249__s0/predictions.parquet"),
}
for m, r in SIMPLE.items():
    SPECS["simple_" + m] = dict(pooled=C + r + "/predictions.parquet", mode="simple",
                                reported=f"/workspace/results/baselines/simple/{m}.json")

GENES = ["PRKACA", "CDK1", "CSNK2A1", "MAPK1", "ATM", "SRC", "GSK3B", "CDK2", "PRKCA", "AKT1", "EGFR", "CSNK1A1",
         "PLK1", "ATR", "ABL1", "MAPK3", "CHEK1", "AURKA", "CAMK2A", "PRKDC"]


def main(name, out):
    sp = SPECS[name]
    t0 = time.time()
    log = {"spec": name}
    S = index.load_sites()
    folds = splits.load_cv_folds()
    site_fold = S.acc.map(folds).fillna(-1).astype(int).values
    K = canary.kinases_with_canaries()
    integ = {}
    if "fold_runs" in sp:
        parts = []
        for f, rid in enumerate(sp["fold_runs"]):
            p = pd.read_parquet(C + rid + "/predictions.parquet", columns=["kinase_idx", "site_idx", "score"])
            p["kinase_idx"] = p.kinase_idx.astype(np.int64)
            p["site_idx"] = p.site_idx.astype(np.int64)
            sf = site_fold[p.site_idx.values]
            d = {"rows": len(p), "n_sites": int(p.site_idx.nunique()), "n_kinases": int(p.kinase_idx.nunique()),
                 "site_fold_counts": {int(k): int(v) for k, v in zip(*np.unique(sf, return_counts=True))},
                 "nonfinite": int((~np.isfinite(p.score.values)).sum())}
            if sp["mode"] == "cv":
                # OOF integrity: every predicted site must be on a protein of eval fold f
                d["oof_ok"] = bool((sf == f).all())
                # coverage of the fold's U_ms sites
                fs = np.where((site_fold == f) & S.in_U_ms.values)[0]
                d["fold_ums_sites"] = len(fs)
                d["fold_ums_sites_scored"] = int(np.isin(fs, p.site_idx.unique()).sum())
            parts.append(p)
            integ[rid] = d
        P = pd.concat(parts, ignore_index=True)
        if "pooled" in sp:
            Q = pd.read_parquet(sp["pooled"], columns=["kinase_idx", "site_idx", "score"])
            a = P.sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
            b = Q.sort_values(["kinase_idx", "site_idx"]).reset_index(drop=True)
            integ["pooled_equals_concat"] = bool(len(a) == len(b) and (a.kinase_idx.values == b.kinase_idx.values).all()
                                                 and (a.site_idx.values == b.site_idx.values).all()
                                                 and np.allclose(a.score.values, b.score.values))
            del Q, a, b
    else:
        P = pd.read_parquet(sp["pooled"])
        P["kinase_idx"] = P.kinase_idx.astype(np.int64)
        P["site_idx"] = P.site_idx.astype(np.int64)
        sf = site_fold[P.site_idx.values]
        integ["pooled"] = {"rows": len(P), "n_kinases": int(P.kinase_idx.nunique()),
                           "site_fold_counts": {int(k): int(v) for k, v in zip(*np.unique(sf, return_counts=True))}}
    log["integrity"] = integ
    data = kev.prepare("cv", fold=None)
    if sp["mode"] == "lko":
        held = sorted(set(K.uniprot_acc.values[np.unique(P.kinase_idx.values)]) & set(data.positives))
        res = kev.evaluate(P, data, kinases=held, missing="drop")
        # also how many required pairs missing (coverage) under 'worst'
    elif sp["mode"] == "simple":
        rep = json.load(open(sp["reported"]))
        real_k = sorted(k for k in data.positives if not k.startswith("CANARY_"))
        res = kev.evaluate(P, data, kinases=real_k, missing="worst")
    else:
        res = kev.evaluate(P, data, missing="worst")
    E = res.endpoints
    log["endpoints"] = {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v) for k, v in E.items()}
    log["canary"] = res.canary
    log["missing"] = {k: v for k, v in res.missing.items() if k != "per_kinase"}
    log["missing_per_kinase_n"] = len(res.missing.get("per_kinase", {}))
    pk = res.per_kinase
    pk.to_csv(f"{out}/{name}_per_kinase.tsv", sep="\t", index=False)
    sel = pk[pk.gene.isin(GENES)][["gene", "n_pos", "n_neg_hard", "n_neg_all", "auroc_hard", "auroc_all",
                                    "top1_recall"] + (["hit5"] if "hit5" in pk else [])]
    log["named_kinases"] = sel.to_dict(orient="records")
    # reported
    try:
        rep = json.load(open(sp["reported"]))
        log["reported_endpoints"] = rep.get("endpoints", rep.get("metrics", {}).get("endpoints"))
        log["reported_canary"] = rep.get("canary")
    except Exception as e:
        log["reported_error"] = repr(e)
    # per-fold M1 for cv specs (fold-specific prepare)
    if sp["mode"] == "cv":
        pf = {}
        for f in range(5):
            dfold = kev.prepare("cv", fold=f)
            Pf = P[site_fold[P.site_idx.values] == f]
            rf = kev.evaluate(Pf, dfold, missing="worst", keep_items=False)
            pf[f] = {"M1": rf.endpoints["M1_auroc_hard"], "M2": rf.endpoints["M2_hit5"],
                     "M3": rf.endpoints["M3_top1_recall"], "neg_canary_hard": rf.canary.get("neg_mean_auroc_hard")}
        log["per_fold"] = pf
    log["t_s"] = round(time.time() - t0, 1)
    with open(f"{out}/{name}.json", "w") as fh:
        json.dump(log, fh, indent=1, default=kev._json_default)
    print(json.dumps({k: log[k] for k in ["spec", "endpoints", "canary", "missing", "t_s"]}, indent=1,
                     default=kev._json_default))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
