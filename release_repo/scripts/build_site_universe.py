#!/usr/bin/env python
"""Build the observed-phosphosite universe U_obs on Swiss-Prot canonical human sequences.

Outputs
  cache/sites/universe.parquet         one row per canonical site (site_id, acc, pos, res, gene, window15, window31,
                                       src_* flags, n_sources, per-source metadata — metadata are NOT features)
  cache/sites/mapping_log.parquet      every source site with its mapping result (method / failure reason)
  results/sites/universe_summary.tsv   per-source mapping counts and failure rates; universe totals
  results/sites/universe_overlap.tsv   source-combination counts (UpSet-style) and pairwise overlaps
Rules: src/ksp/ids.py, src/ksp/sites.py, docs/labels_universe.md.
Usage: python -I scripts/build_site_universe.py
"""
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import ids, kinases, sites  # noqa: E402

SOURCES = ["psp", "psp_ks", "compendium", "cptac", "ochoa"]
U_OBS = ["psp", "compendium", "cptac", "ochoa"]  # DESIGN §2 U_obs definition (psp_ks kept as a flag)


def first_nonnull(s):
    s = s.dropna()
    return s.iloc[0] if len(s) else None


def main():
    t0 = time.time()
    out_c = ROOT / "cache" / "sites"
    out_r = ROOT / "results" / "sites"
    out_c.mkdir(parents=True, exist_ok=True)
    out_r.mkdir(parents=True, exist_ok=True)
    st = sites.default_store()
    hg = kinases.load_hgnc()
    e2a = sites.ensg_candidates(st, hg)

    loaders = {"psp": sites.load_psp_sites_human, "psp_ks": sites.load_psp_ks_human_sites,
               "compendium": sites.load_compendium_sites, "cptac": sites.load_cptac_sites,
               "ochoa": sites.load_ochoa_sites}
    mapped, summ, logs = {}, [], []
    for name in SOURCES:
        t = time.time()
        src = loaders[name]()
        m = sites.map_source(st, src, name, ensg2acc=e2a)
        mapped[name] = m
        summ.append(sites.mapping_summary(m, name))
        logcols = ["source", "src_id", "acc", "pos", "res", "window", "map_method", "map_reason", "map_n_hits",
                   "site_id"]
        if name == "cptac":
            m = m.assign(acc=m.ensg)
        logs.append(m[logcols].astype({"pos": "Int64"}))
        print(f"{name}: {len(src)} source sites, mapped {m.site_id.notna().mean():.4f}, "
              f"{time.time() - t:.1f}s", file=sys.stderr)

    # ---- aggregate to canonical sites
    agg = {}
    m = mapped["psp"].dropna(subset=["site_id"])
    agg["psp"] = m.groupby("site_id").agg(
        psp_lt_lit=("psp_lt_lit", "max"), psp_ms_lit=("psp_ms_lit", "max"), psp_ms_cst=("psp_ms_cst", "max"),
        psp_site_grp_id=("psp_site_grp_id", "min"), psp_ambiguous=("psp_ambiguous", "any"),
        psp_n_rows=("src_id", "size"),
        psp_map_method=("map_method", lambda s: "exact_pos" if (s == "exact_pos").any() else s.iloc[0]))
    m = mapped["psp_ks"].dropna(subset=["site_id"])
    agg["psp_ks"] = m.groupby("site_id").agg(psp_ks_n_rows=("src_id", "size"))
    m = mapped["compendium"].dropna(subset=["site_id"])
    agg["compendium"] = m.groupby("site_id").agg(
        comp_uid=("src_id", "first"), comp_n_uids=("src_id", "size"), comp_n_conditions=("comp_n_conditions", "max"),
        comp_map_method=("map_method", lambda s: "exact_pos" if (s == "exact_pos").any() else s.iloc[0]))
    m = mapped["cptac"].dropna(subset=["site_id"])
    agg["cptac"] = m.groupby("site_id").agg(
        cptac_id=("src_id", "first"),
        cptac_cohorts=("cptac_cohorts", lambda s: ",".join(sorted(set(",".join(s).split(","))))),
        cptac_map_method=("map_method", lambda s: "exact_pos" if (s == "exact_pos").any() else s.iloc[0]))
    agg["cptac"]["cptac_n_cohorts"] = agg["cptac"].cptac_cohorts.str.count(",") + 1
    m = mapped["ochoa"].dropna(subset=["site_id"])
    agg["ochoa"] = m.groupby("site_id").agg(
        ochoa_functional_score=("ochoa_functional_score", "max"),
        ochoa_biological_samples=("ochoa_biological_samples", "max"),
        ochoa_best_loc_prob=("ochoa_best_loc_prob", "max"))

    all_ids = sorted(set().union(*[set(a.index) for a in agg.values()]))
    u = pd.DataFrame({"site_id": all_ids})
    p = u.site_id.str.extract(ids.SITE_RE.pattern)
    u["acc"], u["res"], u["pos"] = p[0], p[1], p[2].astype(int)
    u["gene"] = u.acc.map(st.acc_to_gene).fillna("")
    for name in SOURCES:
        u[f"src_{name}"] = u.site_id.isin(agg[name].index)
    u["n_sources"] = u[[f"src_{s}" for s in U_OBS]].sum(axis=1)
    u["in_U_obs"] = u.n_sources > 0
    for name in SOURCES:
        u = u.merge(agg[name], left_on="site_id", right_index=True, how="left")
    # sanity: canonical residue
    chk = [st.residue(a, p) for a, p in zip(u.acc, u.pos)]
    assert all(c == r for c, r in zip(chk, u.res)), "canonical residue mismatch"
    u["window15"] = [st.window(a, p, 7) for a, p in zip(u.acc, u.pos)]
    u["window31"] = [st.window(a, p, 15) for a, p in zip(u.acc, u.pos)]
    for c in ["psp_lt_lit", "psp_ms_lit", "psp_ms_cst", "psp_n_rows", "psp_ks_n_rows", "comp_n_uids",
              "comp_n_conditions", "cptac_n_cohorts", "ochoa_biological_samples", "psp_site_grp_id"]:
        u[c] = u[c].astype("Int64")
    u["psp_ambiguous"] = u.psp_ambiguous.astype("boolean")
    u.to_parquet(out_c / "universe.parquet", index=False)
    pd.concat(logs, ignore_index=True).to_parquet(out_c / "mapping_log.parquet", index=False)

    # ---- summaries
    s = pd.concat(summ, ignore_index=True)
    tot = [("universe", "sites_all_sources", len(u), np.nan),
           ("universe", "sites_U_obs", int(u.in_U_obs.sum()), np.nan),
           ("universe", "sites_only_psp_ks (labels outside U_obs)", int((~u.in_U_obs).sum()), np.nan),
           ("universe", "proteins_U_obs", int(u.loc[u.in_U_obs, "acc"].nunique()), np.nan)]
    for r, n in u[u.in_U_obs].res.value_counts().items():
        tot.append(("universe", f"U_obs_res:{r}", int(n), n / u.in_U_obs.sum()))
    ev = (u.psp_lt_lit.fillna(0) >= 1) | ((u.psp_ms_lit.fillna(0) + u.psp_ms_cst.fillna(0)) >= 2)
    tot.append(("universe", "psp_LT>=1_or_MS>=2", int(ev.sum()), np.nan))
    tot.append(("universe", "U_obs_supported_by>=2_sources", int((u.n_sources >= 2).sum()), np.nan))
    for name in SOURCES:
        only = u[f"src_{name}"] & (u[[f"src_{x}" for x in SOURCES if x != name]].sum(axis=1) == 0)
        tot.append((name, "sites_unique_to_source", int(only.sum()), np.nan))
    s = pd.concat([s, pd.DataFrame(tot, columns=s.columns)], ignore_index=True)
    s.to_csv(out_r / "universe_summary.tsv", sep="\t", index=False, float_format="%.4f")

    flags = [f"src_{x}" for x in SOURCES]
    combo = u.groupby(flags).size().rename("n_sites").reset_index().sort_values("n_sites", ascending=False)
    combo.insert(0, "kind", "combination")
    pw = []
    for a_, b_ in itertools.product(SOURCES, SOURCES):
        n_ab = int((u[f"src_{a_}"] & u[f"src_{b_}"]).sum())
        pw.append({"kind": "pairwise", "source_a": a_, "source_b": b_, "n_both": n_ab,
                   "frac_of_a": n_ab / max(int(u[f"src_{a_}"].sum()), 1)})
    pd.concat([combo, pd.DataFrame(pw)], ignore_index=True).to_csv(out_r / "universe_overlap.tsv", sep="\t",
                                                                   index=False, float_format="%.4f")
    print(s.to_string(index=False))
    print(f"total {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
