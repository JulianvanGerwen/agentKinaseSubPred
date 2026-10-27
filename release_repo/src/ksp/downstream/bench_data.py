"""benchmarKIN data (Zenodo 12566560; Müller-Dott et al. 2025 Nat Commun 16:4771) mapped to canonical IDs.

Perturbation benchmark (one column per experiment; values = log2FC as distributed):
  hernandez  benchmark_data_hernandez.csv  rows 'GENE|S103|ENSG|ENSP' (no flanking sequence) -> candidates from the
             ENSG (HGNC uniprot_ids) else the gene symbol; the stated position is accepted on the Swiss-Prot canonical
             when the residue matches and exactly one candidate matches (ASSUMPTION: same as the paper, which maps
             these by symbol + residue + position). Metadata: benchmark_metadata_hernandez.csv (one row per
             experiment x target; sign +1 = activation expected, -1 = inhibition).
  hijazi     benchmark_data_hijazi.csv rows 'GENE_S12|GENE|S12' -> gene symbol candidates, same position rule.
             Metadata targets ';'-separated (the paper's manually curated TTD-based target list).
  tyrosine   18 tyrosine-kinase experiments (Dasatinib/Imatinib/Bosutinib/Nilotinib, EGF, HRG, TCR stimulation;
             6 PMIDs, none in our compendium) from the benchmarKIN GitHub repo data/*.rda (external/benchmarkin_github),
             same position rule.
Tumour benchmark (CPTAC BCM v1.2, 7 cohorts, site median-centred, >=30 samples):
  rows 'ENSG.v|ENSP.v|S1053|15-mer|k' -> 15-mer window mapping on the ENSG candidates (ksp.ids map_window), i.e. the
  same rule used for the CPTAC part of the site universe. Gold standard: benchmarking_gold_standard_top5per.Rds
  (per cohort GS_pos_pairs / GS_neg_pairs: kinase symbol -> sample ids; top/bottom 5% kinase protein level), read
  with the pure-Python `rdata` parser (no R code is executed).
Kinase names are resolved with ksp.kinases.resolve_kinase_names(source='benchmarkin'); unresolved names (e.g. the
lipid kinases PIK3CB/CD/G, metabolic kinases in the tumour GS) are reported and dropped.
Several source rows mapping to one site_id are averaged (nan-mean).
"""
from __future__ import annotations

import json
import re
import warnings
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
BK = ROOT / "external" / "benchmarkin"
BK_GH = ROOT / "external" / "benchmarkin_github"     # tyrosine set (not in the Zenodo record)
CACHE = ROOT / "cache" / "downstream"
RESULTS = ROOT / "results" / "downstream"
TUMOUR_COHORTS = ("BRCA", "CCRCC", "GBM", "HNSCC", "LSCC", "LUAD", "UCEC")


@dataclass
class PerturbationBenchmark:
    data: pd.DataFrame            # site_id x experiment (float)
    meta: pd.DataFrame            # one row per (experiment, target kinase): experiment, dataset, kinase_acc, gene, sign, ...
    experiments: pd.DataFrame     # one row per experiment: experiment, dataset, description, pmid, sign, n_sites
    mapping: dict = field(default_factory=dict)

    def subset(self, datasets=None, experiments=None) -> "PerturbationBenchmark":
        e = self.experiments
        if datasets is not None:
            e = e[e.dataset.isin(list(datasets))]
        if experiments is not None:
            e = e[e.experiment.isin(list(experiments))]
        keep = list(e.experiment)
        return PerturbationBenchmark(self.data[keep], self.meta[self.meta.experiment.isin(keep)].reset_index(drop=True),
                                     e.reset_index(drop=True), self.mapping)


@dataclass
class TumourBenchmark:
    data: dict                    # cohort -> DataFrame site_id x sample
    gs: pd.DataFrame              # cohort, kinase_acc, gene, sample, label (1 = GS+, 0 = GS-)
    mapping: dict = field(default_factory=dict)


# ------------------------------------------------------------------------------------------------ mapping helpers
@lru_cache(maxsize=1)
def _store():
    from ksp import sites
    return sites.default_store(verbose=False)


@lru_cache(maxsize=1)
def _ensg2acc() -> dict:
    from ksp import sites
    from ksp.kinases import load_hgnc
    return sites.ensg_candidates(_store(), load_hgnc())


def map_by_position(gene: list, res: list, pos: list, ensg: list | None = None) -> pd.DataFrame:
    """Map (gene, residue, position[, ENSG]) without a window: candidates = ENSG->acc (HGNC) else gene symbol;
    accept iff exactly one candidate canonical has `res` at `pos`. Returns site_id, method, reason."""
    st = _store()
    e2a = _ensg2acc()
    ensg = ensg if ensg is not None else [None] * len(gene)
    out = []
    for g, r, p, e in zip(gene, res, pos, ensg):
        cands = list(e2a.get(e, [])) if e else []
        how = "ensg"
        if not cands:
            cands = list(st.gene_to_accs.get(str(g).upper(), [])) if isinstance(g, str) else []
            how = "symbol"
        if not cands:
            out.append((None, "unmapped", "no_candidate"))
            continue
        hits = [a for a in cands if 1 <= p <= len(st.seqs[a]) and st.seqs[a][p - 1] == r]
        if len(hits) == 1:
            out.append((f"{hits[0]}_{r}{p}", f"pos_residue_{how}", ""))
        elif len(hits) > 1:
            out.append((None, "unmapped", "multiple_candidates"))
        else:
            out.append((None, "unmapped", "residue_mismatch"))
    return pd.DataFrame(out, columns=["site_id", "map_method", "map_reason"])


def _collapse(values: pd.DataFrame, site_id: pd.Series) -> pd.DataFrame:
    v = values.copy()
    v.index = site_id.values
    v = v[v.index.notna()]
    return v.groupby(level=0, sort=True).mean()          # nan-mean over duplicate source rows


def _summary(m: pd.DataFrame, name: str) -> dict:
    return {"source_rows": int(len(m)), "mapped_rows": int(m.site_id.notna().sum()),
            "mapped_frac": round(float(m.site_id.notna().mean()), 4),
            "unique_sites": int(m.site_id.dropna().nunique()),
            "methods": m.map_method.value_counts().to_dict(),
            "unmapped_reasons": m.loc[m.site_id.isna(), "map_reason"].value_counts().to_dict(), "name": name}


def resolve_targets(names: list[str]) -> tuple[pd.DataFrame, dict]:
    from ksp.kinases import resolve_kinase_names
    uniq = sorted({n for n in names if isinstance(n, str) and n})
    r, info = resolve_kinase_names(uniq, source="benchmarkin")
    return r, info


# ------------------------------------------------------------------------------------------------ perturbation
def _load_hernandez() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    d = pd.read_csv(BK / "benchmark_data_hernandez.csv")
    p = d.ID.str.split("|", expand=True)
    rp = p[1].str.extract(r"^([STY])(\d+)$")
    ok = rp[0].notna()
    m = pd.DataFrame({"site_id": None, "map_method": "unmapped", "map_reason": "bad_site"}, index=d.index)
    mm = map_by_position(p.loc[ok, 0].tolist(), rp.loc[ok, 0].tolist(), rp.loc[ok, 1].astype(int).tolist(),
                         p.loc[ok, 2].tolist())
    m.loc[ok, ["site_id", "map_method", "map_reason"]] = mm.values
    data = _collapse(d.drop(columns="ID"), m.site_id)
    meta = pd.read_csv(BK / "benchmark_metadata_hernandez.csv")
    meta = meta.rename(columns={"id": "experiment", "target": "gene_in"})
    meta["dataset"] = "hernandez"
    meta["drug"] = meta.Description
    return data, meta, _summary(m, "hernandez")


def _load_hijazi() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    d = pd.read_csv(BK / "benchmark_data_hijazi.csv")
    p = d.ID.str.split("|", expand=True)
    rp = p[2].str.extract(r"^([STY])(\d+)$")
    ok = rp[0].notna()
    m = pd.DataFrame({"site_id": None, "map_method": "unmapped", "map_reason": "bad_site"}, index=d.index)
    mm = map_by_position(p.loc[ok, 1].tolist(), rp.loc[ok, 0].tolist(), rp.loc[ok, 1].astype(int).tolist())
    m.loc[ok, ["site_id", "map_method", "map_reason"]] = mm.values
    data = _collapse(d.drop(columns="ID"), m.site_id)
    meta = pd.read_csv(BK / "benchmark_metadata_hijazi.csv")
    meta["gene_in"] = meta.target.str.split(";")
    meta = meta.explode("gene_in").rename(columns={"id": "experiment"}).drop(columns="target")
    meta["dataset"] = "hijazi"
    meta["Description"] = meta.drug
    meta["PMID"] = 31959955
    meta["Time_min"] = 60
    return data, meta, _summary(m, "hijazi")


def _load_tyrosine() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """18 tyrosine-kinase experiments (paper 'Tyrosine' set; benchmarKIN GitHub data/tyrosine{Data,Meta}.rda, read with
    the pure-Python rdata parser). Rows 'GENE_Y115|GENE|Y115' -> gene symbol + position rule."""
    import rdata
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = rdata.read_rda(str(BK_GH / "tyrosineData.rda"))["tyrosineData"].reset_index(drop=True)
        meta = rdata.read_rda(str(BK_GH / "tyrosineMeta.rda"))["tyrosineMeta"].reset_index(drop=True)
    p = d.ID.str.split("|", expand=True)
    rp = p[2].str.extract(r"^([STY])(\d+)$")
    ok = rp[0].notna()
    m = pd.DataFrame({"site_id": None, "map_method": "unmapped", "map_reason": "bad_site"}, index=d.index)
    mm = map_by_position(p.loc[ok, 1].tolist(), rp.loc[ok, 0].tolist(), rp.loc[ok, 1].astype(int).tolist())
    m.loc[ok, ["site_id", "map_method", "map_reason"]] = mm.values
    data = _collapse(d.drop(columns="ID").astype(float), m.site_id)
    meta = meta.rename(columns={"id": "experiment", "target": "gene_in"})
    meta["dataset"] = "tyrosine"
    meta["Description"] = meta.treatment
    meta["drug"] = meta.treatment
    meta["PMID"] = pd.to_numeric(meta.PMID)
    meta["Time_min"] = np.nan
    meta["sign"] = meta.sign.astype(int)
    return data, meta, _summary(m, "tyrosine")


PERTURBATION_DATASETS = ("hernandez", "hijazi", "tyrosine")


def load_perturbation(refresh: bool = False) -> PerturbationBenchmark:
    """The three perturbation datasets (Hernandez-Armenta 2017 compilation, Hijazi 2020, tyrosine set) mapped to site_id and kinase
    accessions. Experiment ids are prefixed by dataset ('hernandez:102_117', 'hijazi:HL60.MK2206'). Cached."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f_data, f_meta, f_map = CACHE / "perturbation_data.parquet", CACHE / "perturbation_meta.parquet", \
        CACHE / "perturbation_mapping.json"
    if not refresh and f_data.exists():
        data = pd.read_parquet(f_data)
        meta = pd.read_parquet(f_meta)
        mapping = json.loads(f_map.read_text())
    else:
        parts, metas, mapping = [], [], {}
        for name, fn in [("hernandez", _load_hernandez), ("hijazi", _load_hijazi), ("tyrosine", _load_tyrosine)]:
            dd, mt, info = fn()
            dd.columns = [f"{name}:{c}" for c in dd.columns]
            mt["experiment"] = name + ":" + mt.experiment.astype(str)
            parts.append(dd)
            metas.append(mt)
            mapping[name] = info
        data = pd.concat(parts, axis=1, join="outer").sort_index()
        meta = pd.concat(metas, ignore_index=True)
        r, info = resolve_targets(meta.gene_in.tolist())
        rmap = dict(zip(r.name, zip(r.gene, r.uniprot_acc, r.how)))
        meta["gene"] = meta.gene_in.map(lambda n: rmap.get(n, (None,))[0])
        meta["kinase_acc"] = meta.gene_in.map(lambda n: rmap.get(n, (None, None))[1])
        meta["resolve_how"] = meta.gene_in.map(lambda n: rmap.get(n, (None, None, "unknown"))[2])
        mapping["kinase_resolution"] = {"counts": info["counts"], "unresolved": info["unresolved"],
                                        "dropped_target_rows": int(meta.kinase_acc.isna().sum())}
        meta = meta[["experiment", "dataset", "Description", "drug", "PMID", "Time_min", "sign", "gene_in", "gene",
                     "kinase_acc", "resolve_how"]].astype({"PMID": "Int64"})
        meta["Time_min"] = pd.to_numeric(meta.Time_min, errors="coerce")
        data.to_parquet(f_data)
        meta.to_parquet(f_meta, index=False)
        f_map.write_text(json.dumps(mapping, indent=1, default=str))
    # experiments: only those with metadata and data; sign must be constant within an experiment
    sg = meta.groupby("experiment").sign.nunique()
    if (sg > 1).any():
        raise ValueError(f"experiments with mixed signs: {list(sg[sg > 1].index)}")
    e = meta.groupby("experiment", sort=False).agg(dataset=("dataset", "first"), description=("Description", "first"),
                                                   pmid=("PMID", "first"), sign=("sign", "first"),
                                                   targets=("gene_in", lambda s: ";".join(s))).reset_index()
    e = e[e.experiment.isin(data.columns)].reset_index(drop=True)
    e["n_sites"] = data[e.experiment].notna().sum().values
    m2 = meta[meta.kinase_acc.notna() & meta.experiment.isin(e.experiment)].drop_duplicates(
        ["experiment", "kinase_acc"]).reset_index(drop=True)
    return PerturbationBenchmark(data[e.experiment], m2, e, mapping)


# ------------------------------------------------------------------------------------------------ tumour
def _cohort_file(c: str) -> Path:
    return BK / f"{c.lower()}_original_medcent_30plus.tsv"


def _map_cptac_ids(ids: pd.Series) -> pd.DataFrame:
    p = ids.str.split("|", expand=True)
    sid = pd.Series([None] * len(ids), index=ids.index, dtype=object)
    meth = pd.Series("unmapped", index=ids.index)
    reason = pd.Series("", index=ids.index)
    need = sid.isna()
    if need.any():
        st, e2a = _store(), _ensg2acc()
        rp = p.loc[need, 2].str.extract(r"^([STY])(\d+)$")
        for i, e, r, pos, w in zip(p.index[need], p.loc[need, 0].str.split(".").str[0], rp[0], rp[1],
                                   p.loc[need, 3]):
            if not isinstance(r, str):
                reason[i] = "bad_site"
                continue
            mr = st.map_window(None, int(pos), str(w).upper(), r, candidates=e2a.get(e, []))
            if mr.site_id:
                sid[i], meth[i] = mr.site_id, "window_" + mr.method
            else:
                reason[i] = mr.reason
    return pd.DataFrame({"site_id": sid.values, "map_method": meth.values, "map_reason": reason.values})


def read_gold_standard(path: Path = BK / "benchmarking_gold_standard_top5per.Rds") -> dict:
    import rdata
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rdata.read_rds(str(path))


def load_tumour(cohorts=TUMOUR_COHORTS, refresh: bool = False) -> TumourBenchmark:
    CACHE.mkdir(parents=True, exist_ok=True)
    f_gs, f_map = CACHE / "tumour_gs.parquet", CACHE / "tumour_mapping.json"
    data, mapping = {}, {}
    if not refresh and f_gs.exists() and all((CACHE / f"tumour_{c}.parquet").exists() for c in cohorts):
        for c in cohorts:
            data[c] = pd.read_parquet(CACHE / f"tumour_{c}.parquet")
        gs = pd.read_parquet(f_gs)
        mapping = json.loads(f_map.read_text())
        return TumourBenchmark(data, gs[gs.cohort.isin(cohorts)].reset_index(drop=True), mapping)
    for c in cohorts:
        d = pd.read_csv(_cohort_file(c), sep="\t")
        m = _map_cptac_ids(d.site)
        data[c] = _collapse(d.drop(columns="site"), m.site_id)
        data[c].to_parquet(CACHE / f"tumour_{c}.parquet")
        mapping[c] = _summary(m, c)
    o = read_gold_standard()
    rows = []
    for c in cohorts:
        for lab, key in [(1, "GS_pos_pairs"), (0, "GS_neg_pairs")]:
            for k, samples in o[c][key].items():
                for s in np.atleast_1d(samples):
                    rows.append((c, str(k), str(s), lab))
    gs = pd.DataFrame(rows, columns=["cohort", "gene_in", "sample", "label"])
    r, info = resolve_targets(gs.gene_in.tolist())
    rmap = dict(zip(r.name, zip(r.gene, r.uniprot_acc)))
    gs["gene"] = gs.gene_in.map(lambda n: rmap.get(n, (None, None))[0])
    gs["kinase_acc"] = gs.gene_in.map(lambda n: rmap.get(n, (None, None))[1])
    mapping["kinase_resolution"] = {"counts": info["counts"], "n_unresolved_names": len(info["unresolved"]),
                                    "unresolved": info["unresolved"],
                                    "gs_pairs_total": int(len(gs)), "gs_pairs_dropped": int(gs.kinase_acc.isna().sum())}
    gs = gs[gs.kinase_acc.notna()].reset_index(drop=True)
    gs.to_parquet(f_gs, index=False)
    f_map.write_text(json.dumps(mapping, indent=1, default=str))
    return TumourBenchmark(data, gs, mapping)


def measured_sites(pert: PerturbationBenchmark | None = None, tum: TumourBenchmark | None = None) -> set:
    s = set()
    if pert is not None:
        s |= set(pert.data.index[pert.data.notna().any(axis=1)])
    if tum is not None:
        for d in tum.data.values():
            s |= set(d.index)
    return s
