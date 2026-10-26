"""Evaluate a model's predictions as kinase-substrate libraries on benchmarKIN (DESIGN §6.7, A21).

  python scripts/downstream/run_predictions.py --pred PRED.parquet --tag R2clean [--universe U_obs|U_ms|measured]
         [--methods zscore,ulm] [--n-iter 1000] [--experiments-file exps.txt] [--no-reference]

PRED.parquet: long table with kinase_acc (or kinase_idx) + site_id (or site_idx) + score (or prob/p/pred), over
U_obs x kinases. It MUST come from the clean model variant (phospho features without the benchmark studies and
without CPTAC; results/features/phospho/exclude_lists.json 'clean_bench'); this script cannot check that, so the
caller records the model run id in --tag. Arms: {tag}_top{PSP,50,100}, ...+PSP (union with PSP_all) and ...-PSP
(PSP_all pairs removed before ranking); reference arms PSP_all, PSP_dev, Cantley_top{PSP,50,100}, 5 shuffles.
Both native (each library's own evaluable set, paper protocol) and matched (cells/kinases scored by all arms) rows
are written to results/downstream/{tag}_summary.tsv and {tag}_per_kinase.tsv.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")

import pandas as pd  # noqa: E402

from ksp.downstream import api, bench_data  # noqa: E402
from ksp.downstream import libraries as L  # noqa: E402

OUT = ROOT / "results" / "downstream"


def matched_groups(libs: dict, tag: str) -> dict:
    """Per N: the model arms vs PSP_all / PSP_dev / Cantley top-N, on the cells all of them score."""
    g = {}
    for n in ("PSP", 50, 100):
        names = [f"{tag}_top{n}", f"{tag}_top{n}+PSP", f"{tag}_top{n}-PSP", "PSP_all", "PSP_dev",
                 f"Cantley_top{n}_Uobs", f"Cantley_top{n}_measured"]
        g[f"N{n}"] = [x for x in names if x in libs]
        g[f"N{n}_vsPSPdev"] = [x for x in (f"{tag}_top{n}", "PSP_dev") if x in libs]
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--universe", default="U_obs", choices=["U_obs", "U_ms", "measured"])
    ap.add_argument("--methods", default="zscore,ulm")
    ap.add_argument("--n-iter", type=int, default=1000)
    ap.add_argument("--experiments-file", default=None, help="restrict the perturbation benchmark (leave-drug-out)")
    ap.add_argument("--no-reference", action="store_true")
    ap.add_argument("--benchmarks", default="perturbation,tumour")
    a = ap.parse_args()
    pert, tum = bench_data.load_perturbation(), bench_data.load_tumour()
    meas = bench_data.measured_sites(pert, tum)
    uni = {"U_obs": None, "U_ms": "U_ms", "measured": meas}[a.universe]
    psp = L.psp_all()
    libs = {} if a.no_reference else api.standard_libraries(universe_measured=meas)
    libs.update(api.prediction_arms(a.pred, a.tag, universe=uni, psp=psp))
    exps = None
    if a.experiments_file:
        exps = [x.strip() for x in Path(a.experiments_file).read_text().split() if x.strip()]
    OUT.mkdir(parents=True, exist_ok=True)
    api.library_stats(libs, pert, tum).to_csv(OUT / f"{a.tag}_library_stats.tsv", sep="\t", index=False)
    res = api.run_benchmarks(libs, methods=tuple(a.methods.split(",")), benchmarks=tuple(a.benchmarks.split(",")),
                             experiments=exps, n_iter=a.n_iter, matched=matched_groups(libs, a.tag))
    res["summary"].assign(pred=a.pred, universe=a.universe).to_csv(OUT / f"{a.tag}_summary.tsv", sep="\t",
                                                                   index=False, float_format="%.4f")
    res["per_kinase"].to_csv(OUT / f"{a.tag}_per_kinase.tsv", sep="\t", index=False, float_format="%.4f")
    s = res["summary"]
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(s[s.dataset.isin(["all", "cptac7"])][["benchmark", "library", "method", "mode", "auroc_mean",
                                                   "scaled_rank_mean", "n_perturbed_kinases_evaluable",
                                                   "n_kinases"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
