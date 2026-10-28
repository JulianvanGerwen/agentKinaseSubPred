"""Validate the benchmarKIN re-implementation with reference libraries (PSP, Cantley, shuffled) and write
results/downstream/{validation_summary.tsv, validation_per_kinase.tsv, library_stats.tsv, evaluable_kinases.tsv,
benchmark_mapping.json}. Usage: python scripts/downstream/run_validation.py [--n-iter 1000] [--methods zscore,ksea,ulm,wmean]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")

import pandas as pd  # noqa: E402

from ksp.downstream import api, bench_data  # noqa: E402

OUT = ROOT / "results" / "downstream"
SHUF = [f"PSP_all_shuffled_s{s}" for s in api.SEEDS]
MATCHED = {"psp": ["PSP_all", "PSP_dev"] + SHUF,
           "cantley_measured": ["PSP_all", "Cantley_topPSP_measured", "Cantley_top50_measured",
                                "Cantley_top100_measured", "KL15site_measured", "PSP_all+KL15site"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-iter", type=int, default=1000)
    ap.add_argument("--methods", default="zscore,ksea,ulm,wmean")
    ap.add_argument("--refresh", action="store_true", help="rebuild the mapped benchmark caches")
    a = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pert = bench_data.load_perturbation(refresh=a.refresh)
    tum = bench_data.load_tumour(refresh=a.refresh)
    meas = bench_data.measured_sites(pert, tum)
    libs = api.standard_libraries(universe_measured=meas)
    api.library_stats(libs, pert, tum).to_csv(OUT / "library_stats.tsv", sep="\t", index=False)
    res = api.run_benchmarks(libs, methods=tuple(a.methods.split(",")), rm_bg_modes=(True, False),
                             n_iter=a.n_iter, matched=MATCHED)
    s = res["summary"]
    s.to_csv(OUT / "validation_summary.tsv", sep="\t", index=False, float_format="%.4f")
    pk = res["per_kinase"]
    pk.to_csv(OUT / "validation_per_kinase.tsv", sep="\t", index=False, float_format="%.4f")
    # evaluable kinases (z-score, native): perturbation = perturbed kinases ranked in >=1 experiment;
    # tumour = kinases with >=30 scored tumours and >=1 GS+ and GS- pair in >=1 cohort
    ev = pk[pk.method.eq("zscore") & pk.dataset.isin(["all", "cptac7"])]
    ev = ev.groupby(["benchmark", "library"]).agg(n=("gene", "nunique"),
                                                   kinases=("gene", lambda x: ";".join(sorted(set(x))))).reset_index()
    ev.to_csv(OUT / "evaluable_kinases.tsv", sep="\t", index=False)
    (OUT / "benchmark_mapping.json").write_text(json.dumps({"perturbation": pert.mapping, "tumour": tum.mapping,
                                                            "n_experiments": pert.experiments.groupby("dataset")
                                                            .size().to_dict()}, indent=1, default=str))
    cols = ["benchmark", "dataset", "library", "method", "mode", "rm_bg", "auroc_mean", "scaled_rank_mean",
            "phit10", "n_perturbed_kinases_evaluable", "n_experiments_auroc", "n_kinases"]
    with pd.option_context("display.width", 250, "display.max_rows", 500, "display.max_columns", 20):
        q = s[s["mode"].eq("native") & s.method.eq("zscore") & (s.rm_bg.ne(False))]
        print(q[cols].round(3).to_string(index=False))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
