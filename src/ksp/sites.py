"""Loaders for phosphosite sources of the site universe and their mapping to canonical site_ids (src/ksp/ids.py).

Each loader returns one row per *source site* with columns:
  src_id (source identifier), acc (source accession), pos (stated position), res, window (15-mer, or None when the
  source has no window), gene (source gene symbol or ''), plus source-specific metadata columns (prefixed).
`map_source` adds the canonical mapping columns (map_acc, map_pos, map_res, map_method, map_reason, site_id).

Sources:
  psp      PhosphoSitePlus Phosphorylation_site_dataset, human, S/T/Y, MOD_RSD '-p' (io_psp.load_psp_sites)
  psp_ks   PhosphoSitePlus Kinase_Substrate_Dataset, human substrate rows, S/T/Y (io_psp.load_psp)
  compendium  data/02_phosphosite_information.csv (peptide_UID 'ACC_POS_M1', 15-mer); n conditions from
           cache/compendium_long.parquet (metadata only)
  cptac    CPTAC pan-cancer BCM-harmonised freeze v1.2 (LinkedOmics mirror), all cohorts, tumour + normal
           phospho_site files; idx 'ENSG|ENSP|S1053|15-mer flank|k' -> mapped by flank via ENSG -> HGNC -> UniProt
  ochoa    Ochoa et al. 2020 annotated phosphoproteome (MOESM4 'annotated_phosphoproteome' + MOESM5 functional
           score); NO window in the source -> position + residue on the current canonical, accepted only when the
           canonical length equals Ochoa's prot_length (sequence-version proxy).
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import io_psp
from .ids import CanonicalStore, map_table

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
EXT = ROOT / "external"
CACHE = ROOT / "cache"
COLS = ["src_id", "acc", "pos", "res", "window", "gene"]


def default_store(verbose: bool = True) -> CanonicalStore:
    """Swiss-Prot canonical store (UniProt 2026_03) with the gene-symbol fallback index extended by UniProt gene
    synonyms and HGNC previous symbols (never overriding a primary symbol)."""
    import sys
    from .kinases import load_hgnc
    st = CanonicalStore.from_uniprot()
    ptm = pd.read_csv(EXT / "uniprot/ptm_tsv/UP000005640_reviewed_ptm.tsv.gz", sep="\t", dtype=str,
                      keep_default_na=False)
    pairs = [(s, a) for a, syn in zip(ptm.Entry, ptm["Gene Names (synonym)"]) for s in syn.split()]
    hg = load_hgnc()
    for u, prev in zip(hg.uniprot_ids, hg.prev_symbol):
        for a in u.split("|"):
            pairs += [(p, a) for p in prev.split("|") if p]
    n = st.add_gene_synonyms(pairs)
    if verbose:
        print(f"store: {len(st.seqs)} canonical sequences; +{n} synonym gene->acc links", file=sys.stderr)
    return st


# ------------------------------------------------------------------------------------------------- loaders
def load_psp_sites_human() -> pd.DataFrame:
    s = io_psp.load_psp_sites()
    s = s[s.organism.eq("human") & s.is_sty & s.mod_type.eq("p")]
    return pd.DataFrame({
        "src_id": s.acc_full + "_" + s.mod_rsd, "acc": s.acc_full, "pos": s.position, "res": s.residue,
        "window": s.window, "gene": s.gene, "psp_lt_lit": s.lt_lit, "psp_ms_lit": s.ms_lit, "psp_ms_cst": s.ms_cst,
        "psp_site_grp_id": s.site_grp_id, "psp_ambiguous": s.ambiguous}).reset_index(drop=True)


def load_psp_ks_human_sites(d: pd.DataFrame | None = None) -> pd.DataFrame:
    d = io_psp.load_psp() if d is None else d
    h = d[d.sub_organism.eq("human") & d.is_sty]
    h = h.drop_duplicates(["sub_acc_full", "position", "window"])
    return pd.DataFrame({"src_id": h.sub_acc_full + "_" + h.mod_rsd, "acc": h.sub_acc_full, "pos": h.position,
                         "res": h.residue, "window": h.window, "gene": h.sub_gene,
                         "psp_site_grp_id": h.site_grp_id}).reset_index(drop=True)


def load_compendium_sites(with_counts: bool = True) -> pd.DataFrame:
    c = pd.read_csv(DATA / "02_phosphosite_information.csv", usecols=["peptide_UID", "sequence_window"], dtype=str)
    p = c.peptide_UID.str.extract(r"^([^_]+)_(\d+)_M(\d+)$")
    out = pd.DataFrame({"src_id": c.peptide_UID, "acc": p[0], "pos": pd.to_numeric(p[1]).astype("Int64"),
                        "res": c.sequence_window.str[7].str.upper(), "window": c.sequence_window.str.upper(),
                        "gene": ""})
    if with_counts and (CACHE / "compendium_long.parquet").exists():
        lg = pd.read_parquet(CACHE / "compendium_long.parquet", columns=["peptide_UID"])
        n = lg.peptide_UID.value_counts()
        out["comp_n_conditions"] = out.src_id.map(n).fillna(0).astype(int)
    return out


def cptac_files(base=EXT / "cptac_pancan_v1.2") -> list[str]:
    return sorted(glob.glob(str(Path(base) / "*" / "*_phospho_site_abundance_*.txt")))


def load_cptac_sites(base=EXT / "cptac_pancan_v1.2") -> pd.DataFrame:
    recs = []
    for f in cptac_files(base):
        cohort = Path(f).parent.name
        tissue = "Normal" if f.endswith("_Normal.txt") else "Tumor"
        ids = pd.read_csv(f, sep="\t", usecols=[0], dtype=str).iloc[:, 0]
        recs.append(pd.DataFrame({"idx": ids, "cohort": cohort, "tissue": tissue}))
    a = pd.concat(recs, ignore_index=True)
    p = a.idx.str.split("|", expand=True)
    a["ensg"] = p[0].str.split(".").str[0]
    a["ensp"] = p[1]
    m = p[2].str.extract(r"^([STY])(\d+)$")
    a["res"], a["pos"] = m[0], pd.to_numeric(m[1]).astype("Int64")
    a["window"] = p[3].str.upper()
    key = ["ensg", "ensp", "res", "pos", "window"]
    g = a.groupby(key, dropna=False)
    out = g.agg(cptac_n_cohorts=("cohort", "nunique"),
                cptac_cohorts=("cohort", lambda s: ",".join(sorted(set(s)))),
                cptac_n_files=("idx", "size"), src_id=("idx", "first")).reset_index()
    out["acc"] = ""
    out["gene"] = ""
    return out


def load_ochoa_sites(cache_path=CACHE / "sites" / "ochoa_raw.parquet") -> pd.DataFrame:
    cache_path = Path(cache_path)
    if cache_path.exists():
        return pd.read_parquet(cache_path)
    f4 = EXT / "ochoa2020" / "41587_2019_344_MOESM4_ESM.xlsx"
    f5 = EXT / "ochoa2020" / "41587_2019_344_MOESM5_ESM.xlsx"
    a = pd.read_excel(f4, sheet_name="annotated_phosphoproteome",
                      usecols=["uniprot", "position", "residue", "biological_samples", "spectralcounts",
                               "prot_length", "best_localization_prob"])
    fs = pd.read_excel(f5, sheet_name="functional_score")
    a = a.merge(fs, on=["uniprot", "position"], how="outer")
    out = pd.DataFrame({"src_id": a.uniprot.astype(str) + "_" + a.position.astype("Int64").astype(str),
                        "acc": a.uniprot.astype(str), "pos": a.position.astype("Int64"),
                        "res": a.residue.astype("string").str.upper(), "window": None, "gene": "",
                        "ochoa_prot_length": a.prot_length.astype("Int64"),
                        "ochoa_biological_samples": a.biological_samples.astype("Int64"),
                        "ochoa_spectralcounts": a.spectralcounts.astype("Int64"),
                        "ochoa_best_loc_prob": a.best_localization_prob.astype(float),
                        "ochoa_functional_score": a.functional_score.astype(float)})
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(cache_path, index=False)
    return out


# ------------------------------------------------------------------------------------------------- mapping
def ensg_candidates(store: CanonicalStore, hgnc: pd.DataFrame) -> dict:
    """ENSG (unversioned) -> list of Swiss-Prot canonical accessions (HGNC uniprot_ids, else HGNC symbol)."""
    out = {}
    for e, u, s in zip(hgnc.ensembl_gene_id, hgnc.uniprot_ids, hgnc.symbol):
        if not e:
            continue
        c = [x for x in u.split("|") if x in store.seqs] if u else []
        if not c:
            c = [a for a in store.gene_to_accs.get(s.upper(), [])]
        if c:
            out.setdefault(e, [])
            out[e] += [x for x in c if x not in out[e]]
    return out


def map_source(store: CanonicalStore, src: pd.DataFrame, name: str, ensg2acc: dict | None = None) -> pd.DataFrame:
    """Map a loaded source table to canonical sites. Returns src with mapping columns and `source`=name."""
    src = src.reset_index(drop=True)
    if name == "ochoa":
        recs = []
        for a, p, r, L in zip(src.acc, src.pos, src.res, src.ochoa_prot_length):
            m = store.map_position(a, p, r, expected_length=L)
            recs.append((m.acc, m.pos, m.res, m.method, m.reason, m.n_hits, m.site_id))
        o = pd.DataFrame(recs, columns=["map_acc", "map_pos", "map_res", "map_method", "map_reason", "map_n_hits",
                                        "site_id"])
        o["map_pos"] = o.map_pos.astype("Int64")
    elif name == "cptac":
        recs = []
        cache = {}
        for e, p, w, r in zip(src.ensg, src.pos, src.window, src.res):
            k = (e, p, w)
            m = cache.get(k)
            if m is None:
                c = (ensg2acc or {}).get(e, [])
                m = store.map_window(None, p, w, r, candidates=c)
                cache[k] = m
            recs.append((m.acc, m.pos, m.res, m.method, m.reason, m.n_hits, m.site_id))
        o = pd.DataFrame(recs, columns=["map_acc", "map_pos", "map_res", "map_method", "map_reason", "map_n_hits",
                                        "site_id"])
        o["map_pos"] = o.map_pos.astype("Int64")
    else:
        o = map_table(store, src, "acc", "pos", "window", "res", "gene")
        o = o.reset_index(drop=True)
    out = pd.concat([src, o], axis=1)
    # residue must be S/T/Y on the canonical
    bad = out.site_id.isna() & out.map_acc.notna()
    out.loc[bad, ["map_method", "map_reason"]] = ["unmapped", "non_sty"]
    out["source"] = name
    return out


def mapping_summary(mapped: pd.DataFrame, name: str) -> pd.DataFrame:
    n = len(mapped)
    rows = [(name, "source_sites", n, 1.0)]
    for k, v in mapped.map_method.value_counts().items():
        rows.append((name, f"method:{k}", int(v), v / n))
    for k, v in mapped.loc[mapped.map_method.eq("unmapped"), "map_reason"].value_counts().items():
        rows.append((name, f"unmapped:{k}", int(v), v / n))
    ok = mapped.site_id.notna()
    rows.append((name, "mapped_source_sites", int(ok.sum()), ok.mean()))
    rows.append((name, "unique_canonical_sites", int(mapped.loc[ok, "site_id"].nunique()), np.nan))
    return pd.DataFrame(rows, columns=["source", "metric", "count", "fraction"])
