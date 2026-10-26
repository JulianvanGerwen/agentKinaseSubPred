"""Final pre-registered downstream analysis (DESIGN §6.7, §16 A21) for one predictions run.

  python scripts/downstream/run_final.py --tag r2_clean [--methods zscore,ulm] [--B 1000] [--n-iter 1000]

Inputs: cache/downstream/<tag>/{uobs_superset,measured_topk}.parquet from scripts/downstream/pred_topk.py (exhaustive
predictions streamed; truncation verified in results/downstream/<tag>_topk_check.json).
Arms (U = Uobs [pre-registered] | measured [sensitivity]):
  ours_top{PSP,50,100}_U, ours_top{N}+PSP_U (∪ PSP_all), ours_top{N}-PSP_U (PSP_all pairs removed before ranking),
  Cantley_top{N}_U, PSP_all, PSP_dev, shuffled PSP_all / ours_top50_Uobs / ours_topPSP_Uobs (seeds 0-4).
Perturbation subsets: all 230 experiments, per dataset, and the leave-drug-out strata (ksp.downstream.ldo):
  drug_unseen / drug_seen / no_drug.
Outputs (results/downstream/final_<tag>_*.tsv): native (each library alone, paper protocol), paired (matched kinase
sets, paired bootstrap), coverage, ldo_strata, per_kinase.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp.downstream import bench_data, evaluate, inference, ldo, paired  # noqa: E402
from ksp.downstream import libraries as L  # noqa: E402

OUT = ROOT / "results" / "downstream"
NS = ("PSP", 50, 100)
SEEDS = (0, 1, 2, 3, 4)


def build_libs(tag: str, meas: set) -> dict:
    psp, dev = L.psp_all(), L.psp_dev()
    libs = {"PSP_all": psp, "PSP_dev": dev}
    src = {"Uobs": pd.read_parquet(ROOT / "cache" / "downstream" / tag / "uobs_superset.parquet"),
           "measured": pd.read_parquet(ROOT / "cache" / "downstream" / tag / "measured_topk.parquet")}
    for u, d in src.items():
        d = d[["kinase_acc", "site_id", "score"]]
        for n in NS:
            nn = L.sizes(psp) if n == "PSP" else n
            top = L.topn_from_long(d, nn, lib_name=f"ours_top{n}_{u}")
            libs[f"ours_top{n}_{u}"] = top
            libs[f"ours_top{n}+PSP_{u}"] = L.union(top, psp, f"ours_top{n}+PSP_{u}")
            libs[f"ours_top{n}-PSP_{u}"] = L.topn_from_long(d, nn, exclude=psp, lib_name=f"ours_top{n}-PSP_{u}")
            libs[f"Cantley_top{n}_{u}"] = L.cantley_topn(nn, universe=None if u == "Uobs" else meas)
    for s in SEEDS:
        libs[f"PSP_all_shuf_s{s}"] = L.shuffled(psp, s)
        libs[f"ours_top50_Uobs_shuf_s{s}"] = L.shuffled(libs["ours_top50_Uobs"], s)
        libs[f"ours_topPSP_Uobs_shuf_s{s}"] = L.shuffled(libs["ours_topPSP_Uobs"], s)
    return libs


def pairs() -> list:
    out = []
    for u in ("Uobs", "measured"):
        for n in NS:
            o = f"ours_top{n}_{u}"
            out += [("PSP_all", o), ("PSP_dev", o), ("PSP_all", f"ours_top{n}+PSP_{u}"),
                    ("PSP_all", f"ours_top{n}-PSP_{u}"), (f"Cantley_top{n}_{u}", o),
                    (f"Cantley_top{n}_{u}", f"ours_top{n}-PSP_{u}")]
    for n in ("PSP", 50):
        out += [(f"ours_top{n}_Uobs_shuf_s0", f"ours_top{n}_Uobs")]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="r2_clean")
    ap.add_argument("--methods", default="zscore,ulm")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--n-iter", type=int, default=1000)
    a = ap.parse_args()
    t0 = time.time()
    pert, tum = bench_data.load_perturbation(), bench_data.load_tumour()
    meas = bench_data.measured_sites(pert, tum)
    strata = ldo.drug_strata(pert.experiments)
    strata.to_csv(OUT / f"final_{a.tag}_{a.methods.replace(",", "+")}_ldo_strata.tsv", sep="\t", index=False)
    subsets = {"all": list(pert.experiments.experiment)}
    for d in bench_data.PERTURBATION_DATASETS:
        subsets[d] = list(pert.experiments.experiment[pert.experiments.dataset.eq(d)])
    for s in ("drug_unseen", "drug_seen", "no_drug"):
        subsets[s] = list(strata.experiment[strata.stratum.eq(s)])
    libs = build_libs(a.tag, meas)
    print(f"{len(libs)} libraries built in {time.time() - t0:.0f}s", flush=True)
    methods = a.methods.split(",")
    native, per_k, prs = [], [], []
    for meth in methods:
        pa = {n: inference.infer(pert.data, lib, meth) for n, lib in libs.items()}
        ta = {n: {c: inference.infer(d, lib, meth, exclude_autophos=True) for c, d in tum.data.items()}
              for n, lib in libs.items()}
        for n in libs:
            for sname, exps in subsets.items():
                b = pert.subset(experiments=exps)
                o, pk = evaluate.perturbation_eval(pa[n], b, rm_bg=True, n_iter=a.n_iter)
                native.append({"benchmark": "perturbation", "subset": sname, "library": n, "method": meth, **o})
                if sname == "all":
                    per_k.append(pk.assign(benchmark="perturbation", library=n, method=meth))
            o, pk = evaluate.tumour_eval(ta[n], tum.gs, n_iter=a.n_iter)
            native.append({"benchmark": "tumour", "subset": "cptac7", "library": n, "method": meth, **o})
            per_k.append(pk.assign(benchmark="tumour", library=n, method=meth))
        print(f"{meth}: native done {time.time() - t0:.0f}s", flush=True)
        for ra, rb in pairs():
            for sname in ("all", "drug_unseen", "drug_seen", "no_drug", "hernandez", "hijazi", "tyrosine"):
                b = pert.subset(experiments=subsets[sname])
                o = paired.perturbation_paired(pa[ra], pa[rb], b, B=a.B)
                prs.append({"benchmark": "perturbation", "subset": sname, "method": meth, "ref": ra, "lib": rb, **o})
            o = paired.tumour_paired(ta[ra], ta[rb], tum.gs, B=a.B)
            prs.append({"benchmark": "tumour", "subset": "cptac7", "method": meth, "ref": ra, "lib": rb, **o})
        print(f"{meth}: paired done {time.time() - t0:.0f}s", flush=True)
    nat = pd.DataFrame(native)
    nat.to_csv(OUT / f"final_{a.tag}_{a.methods.replace(",", "+")}_native.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(prs).to_csv(OUT / f"final_{a.tag}_{a.methods.replace(",", "+")}_paired.tsv", sep="\t", index=False, float_format="%.4f")
    pk = pd.concat(per_k, ignore_index=True)
    pk.to_csv(OUT / f"final_{a.tag}_{a.methods.replace(",", "+")}_per_kinase.tsv", sep="\t", index=False, float_format="%.4f")
    # coverage
    ps, ts = set(pert.data.index), set().union(*[set(d.index) for d in tum.data.values()])
    cov = []
    z = nat[nat.method.eq(methods[0])]
    for n, lib in libs.items():
        sz = L.sizes(lib)
        mp = lib[lib.site_id.isin(ps)].groupby("kinase_acc").site_id.nunique()
        mt = lib[lib.site_id.isin(ts)].groupby("kinase_acc").site_id.nunique()
        r1 = z[(z.library == n) & (z.benchmark == "perturbation") & (z.subset == "all")].iloc[0]
        r2 = z[(z.library == n) & (z.benchmark == "tumour")].iloc[0]
        cov.append({"library": n, "n_pairs": len(lib), "n_kinases": len(sz), "n_sites": lib.site_id.nunique(),
                    "median_size": float(sz.median()), "frac_pairs_in_PSP_all": float(
                        pd.Index(list(zip(lib.kinase_acc, lib.site_id))).isin(
                            list(zip(libs["PSP_all"].kinase_acc, libs["PSP_all"].site_id))).mean()),
                    "kinases_ge5_pert_sites": int((mp >= 5).sum()), "kinases_ge5_tumour_sites": int((mt >= 5).sum()),
                    "pert_evaluable_perturbed_kinases": int(r1.n_perturbed_kinases_evaluable),
                    "pert_experiments": int(r1.n_experiments_auroc), "tumour_evaluable_kinases": int(r2.n_kinases)})
    pd.DataFrame(cov).to_csv(OUT / f"final_{a.tag}_{a.methods.replace(",", "+")}_coverage.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
