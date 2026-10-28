#!/usr/bin/env python
"""A32(a) comparison: priming-feature R2 (pooled 5-fold, seed 0) vs r2_final_s0_cv5, with the pre-declared rule:
adopt only if overall M1 delta >= -0.002, M2 delta >= -0.005, primed-positives-stratum M1 delta >= +0.02 with a
cluster-bootstrap CI above 0, and the canary passes (pooled metrics of the new run).
Primed stratum = ksp.strata 'primed' (PRIMARY curated priming kinases/offsets; docs/strata_cv.md §1), positives
subset against unchanged E-hard negatives, macro over K_cv. Paired CIs: ksp.eval.paired_bootstrap (overall) and
ksp.eval._boot_matrix on identical cluster weights for the stratum (B = 1000, seed 0). A24 as in strata_cv_report.
Usage: python -I scripts/models/priming_compare.py --new <pooled run_id>
"""
import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import strata as ST  # noqa: E402
import strata_cv_report as R  # noqa: E402

BASE = "r2_final_s0_cv5__20261026-181435__s0"
OUT = ROOT / "results" / "models" / "priming"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True)
    a = ap.parse_args()
    data = E.prepare("cv", fold=None)
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].reset_index(drop=True)
    data = dataclasses.replace(data, pos_meta=pm,
                               positives={k: v for k, v in data.positives.items() if not k.startswith("CANARY_")})
    ktab = ST.kinase_table()
    site_attr, _, _ = R.site_attributes(data, ktab)
    res = {}
    for tag, rid in [("base", BASE), ("new", a.new)]:
        Pn = R.load_normalised(rid, data, None)
        res[tag] = E.evaluate(Pn, data, missing="worst", fold_normalise=False, keep_items=True)
        del Pn
    out = {"base": BASE, "new": a.new}
    for t in res:
        e = res[t].endpoints
        out[t] = {"M1": e["M1_auroc_hard"], "M2": e["M2_hit5"], "M3": e["M3_top1_recall"], "M4": e["M4_auroc_hard_Kany"]}
    pb = E.paired_bootstrap(res["new"], res["base"], B=1000, seed=0)
    out["paired_overall"] = pb.to_dict(orient="records")
    cl = E._universe_clusters(res["base"])
    rows = []
    for (s, v) in [("primed", "primed"), ("primed", "unprimed"), ("primed_cantley", "primed"),
                   ("primed_cantley", "unprimed")]:
        kin_ok, attr = site_attr[(s, v)]
        subs, mats = {}, {}
        for t in res:
            pk = res[t].per_kinase.set_index("kinase_acc")
            keep = {k: attr[k] == v for k in pk.index if k in kin_ok and k in data.positives}
            sub = R.subset_result(res[t], data, keep)
            spk = sub.per_kinase.set_index("kinase_acc")
            kp = [k for k in data.k_primary if k in spk.index]
            subs[t] = (spk, kp)
            mats[t] = E._boot_matrix(sub, cl, R.B_BOOT, R.SEED, kp)
        kp = subs["base"][1]
        assert kp == subs["new"][1]
        pt = {t: E._nanmean(subs[t][0].loc[kp, "auroc_hard"]) for t in res}
        d = np.nanmean(mats["new"]["M1"], 1) - np.nanmean(mats["base"]["M1"], 1)
        rows.append({"stratum": f"{s}={v}", "n_kinases": len(kp), "n_pos": int(subs["base"][0].loc[kp, "n_pos"].sum()),
                     "M1_base": pt["base"], "M1_new": pt["new"], "delta": pt["new"] - pt["base"],
                     "delta_lo": float(np.percentile(d, 2.5)), "delta_hi": float(np.percentile(d, 97.5)),
                     "per_kinase": {subs['base'][0].loc[k, 'gene']: [round(subs['base'][0].loc[k, 'auroc_hard'], 4),
                                                                      round(subs['new'][0].loc[k, 'auroc_hard'], 4)]
                                    for k in kp}})
    out["strata"] = rows
    pool = json.loads((ROOT / "results" / "runs" / a.new / "metrics.json").read_text())
    can = pool.get("canary", {})
    out["canary_new"] = can
    dM1 = out["new"]["M1"] - out["base"]["M1"]
    dM2 = out["new"]["M2"] - out["base"]["M2"]
    pr = rows[0]
    crit = {"M1_delta>=-0.002": dM1 >= -0.002, "M2_delta>=-0.005": dM2 >= -0.005,
            "primed_M1_delta>=+0.02": pr["delta"] >= 0.02, "primed_CI_lo>0": pr["delta_lo"] > 0,
            "canary_pass": bool(can.get("neg_pass") and can.get("pos_pass"))}
    out.update({"dM1": dM1, "dM2": dM2, "criteria": crit, "adopt": all(crit.values())})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "priming_compare.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
