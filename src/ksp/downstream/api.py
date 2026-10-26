"""Entry point: run the benchmarKIN perturbation and tumour benchmarks for a set of libraries.

    from ksp.downstream import api, libraries as L
    libs = api.standard_libraries()                                  # PSP_all, PSP_dev, Cantley arms, shuffled
    libs.update(api.prediction_arms("cache/.../clean_preds.parquet", "R2clean"))   # top-N, ∪PSP, minus PSP
    res = api.run_benchmarks(libs, methods=("zscore", "ulm"))
    res["summary"]  # one row per (benchmark, dataset, library, method, mode)

Leakage rule (DESIGN §6.7): predictions evaluated here must come from the 'clean' model variant whose phospho
features exclude the benchmark studies (results/features/phospho/exclude_lists.json: clean_bench; tumour benchmark
also excludes all CPTAC-derived features). For leave-drug-out (A21) pass `experiments=` (the experiments of the
held-out drugs) together with predictions from the matching leave-drug-out model.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from . import bench_data, evaluate, inference
from . import libraries as L

SEEDS = (0, 1, 2, 3, 4)


def standard_libraries(universe_measured: set | None = None, shuffle_seeds=SEEDS, cantley_ns=("PSP", 50, 100),
                       include_dev: bool = True) -> dict:
    """Reference arms: PSP_all, PSP_dev, Cantley top-N per kinase (U_obs and, if given, the benchmark-measured
    universe), the paper's Kinase-Library arm (top-15 kinases per measured site, alone and ∪ PSP_all), and
    degree-preserving shuffles of PSP_all."""
    psp = L.psp_all()
    libs = {"PSP_all": psp}
    if include_dev:
        libs["PSP_dev"] = L.psp_dev()
    for n in cantley_ns:
        nn = L.sizes(psp) if n == "PSP" else n
        libs[f"Cantley_top{n}_Uobs"] = L.cantley_topn(nn)
        if universe_measured is not None:
            libs[f"Cantley_top{n}_measured"] = L.cantley_topn(nn, universe=universe_measured)
    if universe_measured is not None:          # the paper's Kinase Library arm: top-15 kinases per measured site
        kl = L.cantley_top_kinases_per_site(15, universe_measured)
        libs["KL15site_measured"] = kl
        libs["PSP_all+KL15site"] = L.union(psp, kl, "PSP_all+KL15site")
    for s in shuffle_seeds:
        libs[f"PSP_all_shuffled_s{s}"] = L.shuffled(psp, s)
    return libs


def prediction_arms(source: str, tag: str, ns=("PSP", 50, 100), universe=None, psp: pd.DataFrame | None = None) -> dict:
    """A21 arms from a predictions parquet over U_obs x kinases: top-N, top-N ∪ PSP_all, top-N of non-PSP pairs
    ('minus PSP': PSP_all pairs removed before ranking)."""
    psp = L.psp_all() if psp is None else psp
    pred = L.read_predictions(source)
    out = {}
    for n in ns:
        nn = L.sizes(psp) if n == "PSP" else n
        top = L.topn_from_long(pred, nn, universe, lib_name=f"{tag}_top{n}")
        out[f"{tag}_top{n}"] = top
        out[f"{tag}_top{n}+PSP"] = L.union(top, psp, f"{tag}_top{n}+PSP")
        out[f"{tag}_top{n}-PSP"] = L.topn_from_long(pred, nn, universe, exclude=psp, lib_name=f"{tag}_top{n}-PSP")
    return out


def library_stats(libs: dict, pert=None, tum=None, min_targets: int = 5) -> pd.DataFrame:
    rows = []
    ps = set(pert.data.index) if pert is not None else set()
    for name, lib in libs.items():
        sz = L.sizes(lib)
        r = {"library": name, "n_pairs": len(lib), "n_kinases": int(len(sz)), "median_size": float(sz.median())
             if len(sz) else 0, "n_sites": int(lib.site_id.nunique())}
        if pert is not None:
            m = lib[lib.site_id.isin(ps)].groupby("kinase_acc").site_id.nunique()
            r["pert_kinases_ge5_measured"] = int((m >= min_targets).sum())
            r["pert_frac_pairs_measured"] = float(lib.site_id.isin(ps).mean())
        if tum is not None:
            ts = set().union(*[set(d.index) for d in tum.data.values()])
            m = lib[lib.site_id.isin(ts)].groupby("kinase_acc").site_id.nunique()
            r["tum_kinases_ge5_measured"] = int((m >= min_targets).sum())
        rows.append(r)
    return pd.DataFrame(rows)


def _match_pert(acts: dict) -> dict:
    """Common (kinase, experiment) cells across libraries (NaN elsewhere)."""
    keys = list(acts)
    idx = sorted(set.intersection(*[set(a.index) for a in acts.values()]))
    cols = list(acts[keys[0]].columns)
    mask = np.ones((len(idx), len(cols)), bool)
    for a in acts.values():
        mask &= a.reindex(index=idx, columns=cols).notna().to_numpy()
    return {k: a.reindex(index=idx, columns=cols).where(mask) for k, a in acts.items()}


def _groups(libs: dict, matched) -> list:
    if not matched:
        return []
    if matched is True:
        return [("matched", list(libs))]
    return [(f"matched:{g}", [n for n in names if n in libs]) for g, names in dict(matched).items()]


def run_benchmarks(libs: dict, methods=("zscore",), benchmarks=("perturbation", "tumour"),
                   perturbation_datasets=(("all", ("hernandez", "hijazi", "tyrosine")), ("hernandez", ("hernandez",)),
                                          ("hijazi", ("hijazi",)), ("tyrosine", ("tyrosine",)),
                                          ("hernandez+hijazi", ("hernandez", "hijazi"))),
                   experiments=None, rm_bg_modes=(True,), n_iter: int = 1000, min_targets: int = 5,
                   matched=False, verbose: bool = True) -> dict:
    """Returns {'summary': DataFrame, 'per_kinase': DataFrame}.
    mode 'native' = each library on its own evaluable set (paper protocol). `matched` (DESIGN matched-universe rule):
    True = one group of all libraries, or {group_name: [library names]}; each group is additionally evaluated on the
    (experiment, kinase) cells (perturbation) / per-cohort kinases (tumour) scored by ALL its members
    (mode 'matched:<group>')."""
    groups = _groups(libs, matched)
    summ, perk = [], []
    if "perturbation" in benchmarks:
        pert = bench_data.load_perturbation()
        if experiments is not None:
            pert = pert.subset(experiments=experiments)
        for meth in methods:
            t0 = time.time()
            acts = {n: inference.infer(pert.data, lib, meth, min_targets) for n, lib in libs.items()}
            modes = [("native", acts)] + [(g, _match_pert({n: acts[n] for n in names})) for g, names in groups]
            for mode, aa in modes:
                for dname, dsets in perturbation_datasets:
                    b = pert.subset(datasets=dsets)
                    if b.experiments.empty:
                        continue
                    for rm_bg in rm_bg_modes:
                        for n, a in aa.items():
                            o, pk = evaluate.perturbation_eval(a, b, rm_bg=rm_bg, n_iter=n_iter)
                            summ.append({"benchmark": "perturbation", "dataset": dname, "library": n, "method": meth,
                                         "mode": mode, "rm_bg": rm_bg, **o})
                            if mode == "native" and rm_bg == rm_bg_modes[0]:
                                perk.append(pk.assign(benchmark="perturbation", dataset=dname, library=n,
                                                      method=meth))
            if verbose:
                print(f"perturbation {meth}: {time.time() - t0:.1f}s", flush=True)
    if "tumour" in benchmarks:
        tum = bench_data.load_tumour()
        for meth in methods:
            t0 = time.time()
            acts = {n: {c: inference.infer(d, lib, meth, min_targets, exclude_autophos=True)
                        for c, d in tum.data.items()} for n, lib in libs.items()}
            modes = [("native", list(libs), None)]
            for g, names in groups:
                common = {c: set.intersection(*[set(acts[n][c].index[acts[n][c].notna().sum(1) >= 30])
                                                for n in names]) for c in tum.data}
                modes.append((g, names, common))
            for mode, names, kins in modes:
                for n in names:
                    o, pk = evaluate.tumour_eval(acts[n], tum.gs, kinases=kins, n_iter=n_iter)
                    summ.append({"benchmark": "tumour", "dataset": "cptac7", "library": n, "method": meth,
                                 "mode": mode, "rm_bg": np.nan, **o})
                    if mode == "native":
                        perk.append(pk.assign(benchmark="tumour", dataset="cptac7", library=n, method=meth))
            if verbose:
                print(f"tumour {meth}: {time.time() - t0:.1f}s", flush=True)
    return {"summary": pd.DataFrame(summ), "per_kinase": pd.concat(perk, ignore_index=True) if perk else
            pd.DataFrame()}
