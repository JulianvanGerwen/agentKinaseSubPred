"""Kinase–substrate label construction (PSP Kinase_Substrate_Dataset + external curated pairs).

Output schema (cache/labels/ks_pairs.parquet): one row per (kinase, site_id, provenance), provenance in
{'PSP', 'OmniPath', 'UniProt'}. Rules are documented in docs/labels_universe.md; summary:

PSP rows (S/T/Y only)
  kinase   : human rows -> kinase list by UniProt accession, else symbol/synonym; non-human rows -> human kinase by
             (upper-case) symbol/synonym (ASSUMPTION: same symbol => ortholog); empty GENE -> PSP kinase name.
  substrate: human rows -> canonical site by 15-mer (ids.CanonicalStore.map_window: exact at stated position, else
             unique elsewhere in the same gene's canonical); non-human rows -> SITE_GRP_ID -> human sites of that
             group (PSP Phosphorylation_site_dataset human rows + KS human rows, each mapped to canonical by 15-mer);
             else ortholog sequence fallback (io_psp.map_sites_by_sequence: same symbol, ungapped 15-mer,
             >= min_identity/15 identical, same residue) on the Swiss-Prot canonical.
  flags    : native_human (>=1 row with kinase AND substrate organism human), sub_native_human, species_origin,
             in_vivo / in_vitro (any row) and their _native / _hmr versions, autophos (kinase accession == substrate
             accession), kinase_included, residue_in_class, in_universe.
  variants : L1 = native rows; L2 = kinase & substrate organism in {human, mouse, rat}; L3 = any organism;
             L2_vitro = L2 restricted to HMR rows with in vitro evidence; L2_noauto = L2 without autophos.
             Every variant additionally requires kinase_included & residue_in_class (& a mapped site).
External (never PSP variants): OmniPath enzsub rows with >= 1 curated non-PSP source; UniProt MOD_RES
  'Phospho...; by <KINASE>' notes with ECO:0000269 evidence. Variant column EXT = kinase_included & residue_in_class.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import io_psp
from .ids import CanonicalStore, map_table, site_ids
from .kinases import KinaseResolver

HMR = ("human", "mouse", "rat")
OMNIPATH_CURATED = ("SIGNOR", "HPRD", "phosphoELM", "NCI-PID_ProtMapper", "Reactome_ProtMapper", "SIGNOR_ProtMapper")
OMNIPATH_OTHER_AGGREGATORS = ("dbPTM", "Li2012")  # neither PSP/text-mining nor primary curation; off by default
PSP_TOKENS = ("PhosphoSite", "PhosphoSite_MIMP", "PhosphoSite_ProtMapper")
_MOD = re.compile(r'MOD_RES (\d+); /note="(Phospho(serine|threonine|tyrosine)[^"]*)"(?:; /evidence="([^"]*)")?')
_RES = {"serine": "S", "threonine": "T", "tyrosine": "Y"}


def residue_in_class(res: pd.Series, cls: pd.Series) -> pd.Series:
    return pd.Series([(isinstance(c, str) and r in c) if c != "ST" else r in ("S", "T")
                      for r, c in zip(res, cls)], index=res.index)


# ---------------------------------------------------------------------------------------------- PSP
def resolve_psp_kinases(d: pd.DataFrame, res: KinaseResolver) -> pd.DataFrame:
    """Add kinase (kinase-list gene or None) and kin_map_how."""
    hum_name2gene = (d[d.kin_organism.eq("human") & d.kin_gene_uc.ne("")]
                     .groupby(d.kin_name.str.upper()).kin_gene_uc.unique())
    name2gene = {k: v[0] for k, v in hum_name2gene.items() if len(v) == 1}
    out_g, out_h = [], []
    for org, g, acc, name in zip(d.kin_organism, d.kin_gene, d.kin_acc, d.kin_name):
        if org == "human":
            k, how = res.resolve(g, acc=acc)
        else:
            k, how = res.resolve(g)
            if k is None and not g:
                gg = name2gene.get(str(name).upper())
                if gg:
                    k, how = res.resolve(gg)
                    how = f"{how}(via_psp_name)"
            how = f"ortholog_{how}" if k is not None else how
        out_g.append(k)
        out_h.append(how)
    d = d.copy()
    d["kinase"], d["kin_map_how"] = out_g, out_h
    return d


def sgid_to_human_sites(store: CanonicalStore, psp_sites_human: pd.DataFrame, ks: pd.DataFrame) -> dict:
    """SITE_GRP_ID -> set of canonical human site_ids (PSP site dataset human rows + KS human rows, mapped by 15-mer)."""
    a = psp_sites_human[["acc_full", "position", "residue", "window", "gene", "site_grp_id"]].rename(
        columns={"acc_full": "acc", "position": "pos", "residue": "res"})
    h = ks[ks.sub_organism.eq("human") & ks.is_sty][["sub_acc_full", "position", "residue", "window", "sub_gene",
                                                      "site_grp_id"]].rename(
        columns={"sub_acc_full": "acc", "position": "pos", "residue": "res", "sub_gene": "gene"})
    x = pd.concat([a, h], ignore_index=True).dropna(subset=["site_grp_id"]).drop_duplicates()
    m = map_table(store, x, "acc", "pos", "window", "res", "gene")
    x = pd.concat([x.reset_index(drop=True), m.reset_index(drop=True)], axis=1).dropna(subset=["site_id"])
    return x.groupby("site_grp_id").site_id.agg(lambda s: sorted(set(s))).to_dict()


def map_psp_substrates(d: pd.DataFrame, store: CanonicalStore, sgid_map: dict, min_identity: int = 11) -> pd.DataFrame:
    """Return rows exploded over mapped canonical sites: original columns + site_id, sub_map_method, sub_map_reason,
    sgid_n_human_sites, seq_identity."""
    d = d.reset_index(drop=True).copy()
    d["_row"] = np.arange(len(d))
    is_h = d.sub_organism.eq("human")
    # human rows: window mapping
    hm = map_table(store, d[is_h].rename(columns={"sub_acc_full": "acc", "position": "pos", "window": "window"}),
                   "acc", "pos", "window", "residue", "sub_gene")
    hum = d[is_h].copy()
    hum["site_id"] = hm.site_id.values
    hum["sub_map_method"] = np.where(hm.site_id.notna(), "human_window_" + hm.map_method.values, "unmapped")
    hum["sub_map_reason"] = np.where(hm.site_id.notna(), "", "human_window_" + hm.map_reason.values)
    # non-human rows: SGID, then ortholog sequence fallback
    non = d[~is_h].copy()
    sg = non.site_grp_id.map(lambda g: sgid_map.get(g, []) if pd.notna(g) else [])
    non["_sg"] = sg.values
    has = non._sg.map(len) > 0
    a = non[has].explode("_sg").rename(columns={"_sg": "site_id"})
    a["sub_map_method"] = "sgid"
    a["sub_map_reason"] = ""
    a["sgid_n_human_sites"] = a.groupby("_row").site_id.transform("nunique")
    b = non[~has].drop(columns="_sg")
    if len(b):
        prot = [(g, acc, acc, store.seqs[acc]) for g, accs in store.gene_to_accs.items() for acc in accs]
        hp = pd.DataFrame(prot, columns=["gene_uc", "acc_full", "acc", "seq"])
        need = b[b.sub_gene_uc.isin(set(hp.gene_uc))]
        uniq = need.drop_duplicates(["sub_gene_uc", "window", "residue", "position"])[
            ["sub_gene_uc", "window", "residue", "position"]]
        hits = io_psp.map_sites_by_sequence(uniq, hp[hp.gene_uc.isin(set(uniq.sub_gene_uc))], min_identity)
        b = b.merge(hits, on=["sub_gene_uc", "window", "residue", "position"], how="left")
        ok = b.seq_h_sub_acc.notna()
        b["site_id"] = None
        b.loc[ok, "site_id"] = site_ids(b.loc[ok, "seq_h_sub_acc"], b.loc[ok, "residue"], b.loc[ok, "seq_h_position"])
        b["sub_map_method"] = np.where(ok, "ortholog_seq", "unmapped")
        b["sub_map_reason"] = np.where(ok, "", np.where(b.site_grp_id.notna(), "sgid_no_human_site;no_ortholog_hit",
                                                        "no_sgid;no_ortholog_hit"))
        b = b.drop(columns=["seq_h_sub_acc_full", "seq_h_sub_acc", "seq_h_position"])
    out = pd.concat([hum, a, b], ignore_index=True).sort_values("_row", kind="stable").reset_index(drop=True)
    out["sgid_n_human_sites"] = out.get("sgid_n_human_sites", pd.Series(np.nan, index=out.index))
    return out


def psp_pairs(d: pd.DataFrame, kl: pd.DataFrame, universe_ids: set | None = None) -> pd.DataFrame:
    """Aggregate mapped PSP rows (columns kinase, site_id, ...) to one row per (kinase, site_id)."""
    k = d.dropna(subset=["kinase", "site_id"]).copy()
    k["site_res"] = k.site_id.str.extract(r"_([STY])\d+$")[0]
    k["site_acc"] = k.site_id.str.split("_").str[0]
    acc = dict(zip(kl.gene, kl.uniprot_acc))
    k["kin_acc_h"] = k.kinase.map(acc)
    k["native"] = k.kin_organism.eq("human") & k.sub_organism.eq("human")
    k["hmr"] = k.kin_organism.isin(HMR) & k.sub_organism.isin(HMR)
    k["sub_native"] = k.sub_organism.eq("human")
    k["autophos_row"] = k.kin_acc_h.eq(k.site_acc)
    k["res_changed"] = k.residue.ne(k.site_res)
    k["origin"] = k.kin_organism.str.replace(" ", "_") + ">" + k.sub_organism.str.replace(" ", "_")
    for ev in ["in_vivo", "in_vitro"]:
        k[f"{ev}_native"] = k[ev] & k.native
        k[f"{ev}_hmr"] = k[ev] & k.hmr
    g = k.groupby(["kinase", "site_id"], sort=False)
    p = g.agg(native_human=("native", "any"), hmr=("hmr", "any"), sub_native_human=("sub_native", "any"),
              species_origin=("origin", lambda s: ";".join(sorted(set(s)))),
              in_vivo=("in_vivo", "any"), in_vitro=("in_vitro", "any"),
              in_vivo_native=("in_vivo_native", "any"), in_vitro_native=("in_vitro_native", "any"),
              in_vivo_hmr=("in_vivo_hmr", "any"), in_vitro_hmr=("in_vitro_hmr", "any"),
              autophos=("autophos_row", "any"), psp_autophos_flag=("autophos", "any"),
              sub_map_methods=("sub_map_method", lambda s: ";".join(sorted(set(s)))),
              kin_map_hows=("kin_map_how", lambda s: ";".join(sorted(set(s)))),
              res_changed=("res_changed", "any"), sgid_multi=("sgid_n_human_sites", lambda s: bool((s > 1).any())),
              n_rows=("_row", "nunique"), psp_mod_rsd=("mod_rsd", lambda s: ";".join(sorted(set(s)))),
              psp_sub_accs=("sub_acc_full", lambda s: ";".join(sorted(set(s)))),
              psp_sgids=("site_grp_id", lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x)})))).reset_index()
    p["provenance"] = "PSP"
    return p


def evidence_class(vivo: pd.Series, vitro: pd.Series) -> pd.Series:
    return pd.Series(np.select([vivo & vitro, vivo, vitro], ["both", "in_vivo_only", "in_vitro_only"], "none"),
                     index=vivo.index)


# ---------------------------------------------------------------------------------------------- OmniPath
def omnipath_pairs(path, store: CanonicalStore, res: KinaseResolver, curated=OMNIPATH_CURATED,
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Curated non-PSP OmniPath enzyme-substrate phosphorylation rows -> (pairs, unmapped log)."""
    o = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    o = o[o.modification.eq("phosphorylation")].copy()
    toks = o.sources.str.split(";")
    o["curated_sources"] = toks.map(lambda t: ";".join(x for x in t if x in curated))
    o["omnipath_lists_psp"] = toks.map(lambda t: any(x in PSP_TOKENS for x in t))
    log = []
    o_drop = o[o.curated_sources.eq("")]
    log.append(pd.DataFrame({"provenance": "OmniPath", "reason": "no_curated_non_psp_source",
                             "kinase_name": o_drop.enzyme_genesymbol, "substrate": o_drop.substrate,
                             "site": o_drop.residue_type + o_drop.residue_offset, "detail": o_drop.sources}))
    o = o[o.curated_sources.ne("")]
    cur_refs = []
    for refs, cs in zip(o.references, o.curated_sources):
        keep = set(c.replace("_ProtMapper", "") for c in cs.split(";")) | set(cs.split(";"))
        cur_refs.append(";".join(sorted({r.split(":", 1)[1] for r in refs.split(";") if ":" in r and
                                         r.split(":", 1)[0] in keep})))
    o["curated_refs"] = cur_refs
    kin, how, sid, reason = [], [], [], []
    for e, eg, s, rt, ro in zip(o.enzyme, o.enzyme_genesymbol, o.substrate, o.residue_type, o.residue_offset):
        k, h = res.resolve(eg, acc=e)
        m = store.map_position(s, ro, rt)
        kin.append(k), how.append(h), sid.append(m.site_id), reason.append(m.reason)
    o["kinase"], o["kin_map_how"], o["site_id"], o["sub_map_reason"] = kin, how, sid, reason
    bad = o.kinase.isna() | o.site_id.isna()
    b = o[bad]
    log.append(pd.DataFrame({"provenance": "OmniPath",
                             "reason": np.where(b.kinase.isna(), "kinase_" + b.kin_map_how,
                                                "site_" + b.sub_map_reason),
                             "kinase_name": b.enzyme_genesymbol, "substrate": b.substrate,
                             "site": b.residue_type + b.residue_offset, "detail": b.sources}))
    o = o[~bad]
    p = o.groupby(["kinase", "site_id"], as_index=False).agg(
        ext_sources=("curated_sources", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),
        ext_refs=("curated_refs", lambda s: ";".join(sorted(set(x for v in s for x in v.split(";") if x)))),
        omnipath_all_sources=("sources", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),
        omnipath_lists_psp=("omnipath_lists_psp", "any"), kin_map_hows=("kin_map_how", "first"))
    p["provenance"] = "OmniPath"
    return p, pd.concat(log, ignore_index=True)


# ---------------------------------------------------------------------------------------------- UniProt
def parse_uniprot_mod_res(path) -> pd.DataFrame:
    """Rows: acc, pos, res, note, evidence, kinase_name (one per 'by' name; 'autocatalysis' kept), eco269, pmids."""
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    recs = []
    for acc, gene, f in zip(t.Entry, t["Gene Names (primary)"], t["Modified residue"]):
        for pos, note, rname, ev in _MOD.findall(f):
            if "; by " not in note:
                continue
            by = note.split("; by ", 1)[1].split(";")[0].strip()
            pmids = ";".join(sorted(set(re.findall(r"ECO:0000269\|PubMed:(\d+)", ev or ""))))
            for name in [x.strip() for x in re.split(r",\s*|\s+and\s+|\s+or\s+", by, flags=re.I) if x.strip()]:
                recs.append((acc, gene, int(pos), _RES[rname], note, ev, name, "ECO:0000269" in (ev or ""), pmids))
    return pd.DataFrame(recs, columns=["acc", "gene", "pos", "res", "note", "evidence", "kinase_name", "eco269",
                                       "pmids"])


def uniprot_pairs(path, store: CanonicalStore, res: KinaseResolver) -> tuple[pd.DataFrame, pd.DataFrame]:
    u = parse_uniprot_mod_res(path)
    log = [pd.DataFrame({"provenance": "UniProt", "reason": "no_ECO:0000269", "kinase_name": u.kinase_name[~u.eco269],
                         "substrate": u.acc[~u.eco269], "site": u.res[~u.eco269] + u.pos[~u.eco269].astype(str),
                         "detail": u.note[~u.eco269]})]
    u = u[u.eco269].copy()
    from .kinases import resolve_kinase_names
    rr, _ = resolve_kinase_names(u.kinase_name, "uniprot", res)   # strict resolver + homonym rules (A4)
    kin, how = [], []
    for name, acc, g, h in zip(u.kinase_name, u.acc, rr.gene, rr.how):
        if name.lower() == "autocatalysis":
            k, h = res.resolve(acc=acc)
            h = "autocatalysis" if k is not None else "autocatalysis_not_kinase"
        else:
            k = g if isinstance(g, str) else None
            h = "family_name" if h == "family_name" else h
        kin.append(k), how.append(h)
    u["kinase"], u["kin_map_how"] = kin, how
    sid, reason = [], []
    for a, p, r in zip(u.acc, u.pos, u.res):
        m = store.map_position(a, p, r)
        sid.append(m.site_id), reason.append(m.reason)
    u["site_id"], u["sub_map_reason"] = sid, reason
    bad = u.kinase.isna() | u.site_id.isna()
    b = u[bad]
    log.append(pd.DataFrame({"provenance": "UniProt",
                             "reason": np.where(b.kinase.isna(), "kinase_" + b.kin_map_how, "site_" + b.sub_map_reason),
                             "kinase_name": b.kinase_name, "substrate": b.acc, "site": b.res + b.pos.astype(str),
                             "detail": b.note}))
    u = u[~bad]
    p = u.groupby(["kinase", "site_id"], as_index=False).agg(
        ext_sources=("kinase_name", lambda s: "UniProt:" + ";".join(sorted(set(s)))),
        ext_refs=("pmids", lambda s: ";".join(sorted(set(x for v in s for x in v.split(";") if x)))),
        uniprot_note=("note", "first"), kin_map_hows=("kin_map_how", lambda s: ";".join(sorted(set(s)))))
    p["provenance"] = "UniProt"
    return p, pd.concat(log, ignore_index=True)


# ---------------------------------------------------------------------------------------------- summaries
def variant_counts(df: pd.DataFrame, mask: pd.Series, name: str) -> dict:
    x = df[mask]
    per_k = x.groupby("kinase").site_id.nunique()
    return {"variant": name, "pairs": len(x.drop_duplicates(["kinase", "site_id"])), "sites": x.site_id.nunique(),
            "substrate_proteins": x.site_id.str.split("_").str[0].nunique(), "kinases": int((per_k > 0).sum()),
            "kinases_ge5": int((per_k >= 5).sum()), "kinases_ge20": int((per_k >= 20).sum()),
            "kinases_ge50": int((per_k >= 50).sum()),
            "top10_share": float(per_k.sort_values(ascending=False).head(10).sum() / max(per_k.sum(), 1))}


# ================================================================================================================
# Core-infra extensions (DESIGN §3-§5): DEV loader, known-positive sets KP(k), paralog sets Par(k), negative masks.
# Owner: core-infra. The DEV label table (cache/labels/dev_pairs.parquet) is written by scripts/build_splits.py and
# physically excludes TEST proteins; TEST labels are sealed (read only by ksp.eval_test).
# ================================================================================================================
from functools import lru_cache  # noqa: E402
from pathlib import Path  # noqa: E402

from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
FULL_PAIRS = ROOT / "cache" / "labels" / "ks_pairs.parquet"     # pre-split full table (removed after sealing)
DEV_PAIRS = ROOT / "cache" / "labels" / "dev_pairs.parquet"     # DEV proteins only (all provenances / flags)
KINCORE_MSA = ROOT / "external" / "kincore" / "alignment-files_Human-PK-alignment.fasta"
PSP2021 = (ROOT / "external" / "kinase_library" / "src" / "kinase_library-1.8.0" / "src" / "kinase_library" /
           "databases" / "substrates" / "Kinase_Substrate_Dataset_count_07_2021.txt")
VARIANTS = ("L1", "L2", "L3", "L2_vitro", "L2_noauto", "EXT", "EXT_not_in_psp")
PAR_IDENTITY = 0.60          # DESIGN §5 paralog threshold (kinase-domain identity)
HAMMING_MAX = 2              # DESIGN §5 Hamming mask (15-mer)


def evaluable_positive(df: pd.DataFrame) -> pd.Series:
    """Evaluation positives (DESIGN §3 'Test positives'): human-native PSP pairs, kinase included, residue in class,
    not autophosphorylation (L1 already implies included & in class). Used identically for TEST and CV."""
    return df.provenance.eq("PSP") & df.L1.astype(bool) & ~df.autophos.astype(bool)


def select_variant(df: pd.DataFrame, variant: str) -> pd.DataFrame:
    """One row per (kinase_acc, site_id) for a label variant. 'all' returns every row unchanged."""
    if variant == "all":
        return df
    if variant == "eval":
        x = df[evaluable_positive(df)]
    elif variant in VARIANTS:
        x = df[df[variant].astype(bool)]
    else:
        raise ValueError(f"unknown label variant {variant!r}; use one of {VARIANTS + ('eval', 'all')}")
    # EXT may have OmniPath + UniProt rows for one pair; keep the first by provenance order
    return x.sort_values(["kinase_acc", "site_id", "provenance"]).drop_duplicates(["kinase_acc", "site_id"]) \
        .reset_index(drop=True)


def load_dev(variant: str = "L2", with_canaries: bool = False, path: Path | None = None,
             ums_only: bool | None = None) -> pd.DataFrame:
    """THE canonical DEV label loader (DESIGN §4.3). Returns DEV-protein labels only (TEST proteins are absent from the
    file). variant: L1/L2/L3/L2_vitro/L2_noauto/EXT/EXT_not_in_psp (one row per pair), 'eval' (evaluation positives),
    or 'all' (every provenance row; use for KP/negative masks). with_canaries=True appends the canary pseudo-kinase
    labels (ksp.canary; kinase_acc 'CANARY_*'; flag column is_canary) - they must never enter reported metrics.
    ums_only (default: True for every variant except 'all'): keep only sites in U_ms (A6). Use variant='all' (all
    sites, all rows) for KP / negative masks."""
    p = Path(path) if path else DEV_PAIRS
    if not p.exists():
        raise FileNotFoundError(f"{p} missing: run scripts/build_splits.py (splits are frozen at G1)")
    d = pd.read_parquet(p)
    d["is_canary"] = False
    d = select_variant(d, variant)
    if ums_only is None:
        ums_only = variant != "all"
    if ums_only:   # A6 (DESIGN §16): training/eval positives restricted to U_ms (flag written at split freeze)
        d = d[d.in_U_ms.astype(bool)].reset_index(drop=True)
    if with_canaries:
        from . import canary
        c = canary.canary_labels()
        if variant not in ("all",):
            c = select_variant(c, variant) if variant != "eval" else c[evaluable_positive(c)]
        d = pd.concat([d, c], ignore_index=True)
    return d


def training_counts(variant: str = "L2") -> pd.Series:
    """kinase_acc -> number of DEV label sites in `variant` (0 for kinases without labels)."""
    from . import index
    d = load_dev(variant)
    n = d.groupby("kinase_acc").site_id.nunique()
    return index.load_kinases().uniprot_acc.map(n).fillna(0).astype(int).set_axis(index.load_kinases().uniprot_acc)


# ------------------------------------------------------------------------------------------------ windows / Hamming
@lru_cache(maxsize=1)
def _canonical_seqs() -> dict:
    from . import uniprot
    sp = uniprot.read_fasta(uniprot.DEFAULT_DIR / "proteome" / "UP000005640_reviewed_isoforms.fasta.gz")
    can = sp[(sp.isoform == "") & (sp.db == "sp")]
    return dict(zip(can.acc, can.seq))


def site_windows(site_ids_: pd.Series | list, flank: int = 7) -> pd.Series:
    """15-mer ('_'-padded) for arbitrary canonical site_ids (index sites read from cache/index, others from sequence)."""
    from . import index
    from .ids import window
    s = pd.Series(list(site_ids_), dtype=object)
    sites_ = index.load_sites()[["site_id", "window15"]]
    w = s.map(dict(zip(sites_.site_id, sites_.window15)))
    miss = w.isna()
    if miss.any():
        seqs = _canonical_seqs()
        out = []
        for sid in s[miss]:
            acc, rp = sid.rsplit("_", 1)
            if acc not in seqs:
                raise KeyError(f"site {sid}: accession not in canonical store")
            out.append(window(seqs[acc], int(rp[1:]), flank))
        w[miss] = out
    return w


def _to_u8(windows) -> np.ndarray:
    a = np.frombuffer("".join(windows).encode("ascii"), dtype=np.uint8)
    return a.reshape(len(windows), -1) if len(windows) else np.zeros((0, 15), np.uint8)


def hamming_neighbors(query: list[str], ref: list[str], max_d: int = HAMMING_MAX) -> list[np.ndarray]:
    """For each query 15-mer, indices into `ref` with Hamming distance <= max_d (exact pigeonhole search over
    max_d+1 blocks: any pair within distance max_d shares >= 1 identical block). Distance 0 (identical) included."""
    if not query:
        return []
    L = len(query[0])
    if any(len(w) != L for w in query) or any(len(w) != L for w in ref):
        raise ValueError("all windows must have equal length")
    nb = max_d + 1
    edges = np.linspace(0, L, nb + 1).astype(int)
    R = _to_u8(ref)
    Q = _to_u8(query)
    buckets = []
    for b in range(nb):
        lo, hi = edges[b], edges[b + 1]
        key = pd.Series([w[lo:hi] for w in ref])
        buckets.append(key.groupby(key).indices)
    out = []
    for qi, w in enumerate(query):
        cand = [buckets[b].get(w[edges[b]:edges[b + 1]]) for b in range(nb)]
        cand = [c for c in cand if c is not None]
        if not cand:
            out.append(np.zeros(0, dtype=np.int64))
            continue
        c = np.unique(np.concatenate(cand))
        d = (R[c] != Q[qi]).sum(1)
        out.append(c[d <= max_d])
    return out


# ------------------------------------------------------------------------------------------------ KP(k)
def known_positive_sets(pairs: pd.DataFrame, expand_identical: bool = True, include_sensitivity: bool = False) -> dict:
    """KP(k) (DESIGN §3): kinase_acc -> set of site_ids labelled to k in ANY row of `pairs` (all sources PSP/OmniPath/
    UniProt, all species, any flag incl. out-of-class / not-in-variant), plus (expand_identical) every U_obs site whose
    15-mer is identical to a labelled site of k. Pass the unfiltered table (load_dev('all') or DEV+TEST at TEST time)."""
    from . import index
    p = pairs.dropna(subset=["kinase_acc"])
    if not include_sensitivity and "mask_sensitivity_only" in p:
        p = p[~p.mask_sensitivity_only.fillna(False).astype(bool)]
    kp = p.groupby("kinase_acc").site_id.agg(lambda s: set(s)).to_dict()
    if expand_identical:
        sites_ = index.load_sites()[["site_id", "window15"]]
        by_w = sites_.groupby("window15").site_id.agg(list).to_dict()
        ls = pd.Series(sorted(set(p.site_id)))
        w = dict(zip(ls, site_windows(ls)))
        for k, ss in kp.items():
            extra = set()
            for s in ss:
                extra.update(by_w.get(w[s], ()))
            ss |= extra
    return kp


# ------------------------------------------------------------------------------------------------ Par(k)
def _read_kincore_msa(path: Path = KINCORE_MSA) -> pd.DataFrame:
    recs, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            h = line[1:].split()
            cur = [h[0], h[-1], []]
            recs.append(cur)
        elif cur is not None:
            cur[2].append(line.strip())
    df = pd.DataFrame([(a, b, "".join(c)) for a, b, c in recs], columns=["name", "acc", "aln"])
    return df[df.acc.ne("P00000")].reset_index(drop=True)   # drop the ANNOTATION row


def kinase_domain_identity(cache: Path | None = ROOT / "cache" / "splits" / "kinase_domain_identity.parquet"
                           ) -> pd.DataFrame:
    """Pairwise kinase-domain identity from the KinCore human PK multiple alignment (DESIGN §4.5/§5):
    identity = identical aligned residues / columns where both have a residue (gaps '-'/'.' ignored), requiring the
    shared columns to cover >= 50% of the shorter domain (else identity 0). Multi-domain kinases (e.g. JAK JH1/JH2,
    RSK N/C): max over domain pairs. Returns long table (acc_a, acc_b, identity) for acc_a < acc_b, identity > 0.3,
    restricted to kinases in cache/index/kinases.parquet. Kinases absent from KinCore (28 atypicals: PIKK, PDHK,
    Alpha, ABC1, RIO, FAM20, FAST) get no identity edges and rely on the Manning family rule."""
    if cache is not None and Path(cache).exists():
        return pd.read_parquet(cache)
    from . import index
    ka = set(index.load_kinases().uniprot_acc)
    m = _read_kincore_msa()
    m = m[m.acc.isin(ka)].reset_index(drop=True)
    A = np.array([list(s.upper()) for s in m.aln])
    res = ~np.isin(A, ["-", "."])
    aa = sorted(set(A[res].ravel()))
    n, L = A.shape
    X = np.zeros((n, L * len(aa)), dtype=np.float32)
    for i, a in enumerate(aa):
        X[:, i::len(aa)] = (A == a)
    match = X @ X.T
    R = res.astype(np.float32)
    both = R @ R.T
    lens = R.sum(1)
    short = np.minimum.outer(lens, lens)
    ident = np.where((both >= 0.5 * short) & (both > 0), match / np.maximum(both, 1), 0.0)
    ii, jj = np.triu_indices(n, k=1)
    long = pd.DataFrame({"acc_a": m.acc.values[ii], "acc_b": m.acc.values[jj], "identity": ident[ii, jj]})
    sw = long.acc_a > long.acc_b
    long.loc[sw, ["acc_a", "acc_b"]] = long.loc[sw, ["acc_b", "acc_a"]].values
    long = long[long.acc_a < long.acc_b].groupby(["acc_a", "acc_b"], as_index=False).identity.max()
    long = long[long.identity > 0.3].reset_index(drop=True)
    long["method"] = "kincore_msa"
    extra = _atypical_identity(index.load_kinases(), set(m.acc))
    long = pd.concat([long, extra], ignore_index=True)
    if cache is not None:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        long.to_parquet(cache, index=False)
    return long


def _atypical_identity(k: pd.DataFrame, in_kincore: set) -> pd.DataFrame:
    """Kinases absent from KinCore: pairwise global alignment (Biopython PairwiseAligner, BLOSUM62, open -10,
    extend -0.5) of UniProt kinase domains (full sequence when no domain feature) against every kinase of the same
    Manning family. identity = identical / aligned (both non-gap) columns, 0 if those cover < 50% of the shorter."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.mode = "global"
    seqs = _canonical_seqs()

    def dom(r):
        sq = seqs[r.uniprot_acc]
        if pd.notna(r.kinase_domain_start) and pd.notna(r.kinase_domain_end):
            return sq[int(r.kinase_domain_start) - 1:int(r.kinase_domain_end)]
        return sq
    rows = []
    fam = k.manning_family.map(_real_family)
    for f, sub in k.assign(f=fam).dropna(subset=["f"]).groupby("f"):
        recs = list(sub.itertuples(index=False))
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                a, b = recs[i], recs[j]
                if a.uniprot_acc in in_kincore and b.uniprot_acc in in_kincore:
                    continue
                sa, sb = dom(a), dom(b)
                aln = al.align(sa, sb)[0]
                ca = aln.counts()
                both = ca.identities + ca.mismatches
                ident = ca.identities / both if both >= 0.5 * min(len(sa), len(sb)) else 0.0
                x, y = sorted([a.uniprot_acc, b.uniprot_acc])
                rows.append((x, y, float(ident), "pairwise_uniprot_domain"))
    return pd.DataFrame(rows, columns=["acc_a", "acc_b", "identity", "method"])


def _real_family(f) -> str | None:
    if not isinstance(f, str) or not f or f.endswith("-Unique") or f == "Unique":
        return None
    return f


def paralog_sets(identity: float = PAR_IDENTITY, kinases: pd.DataFrame | None = None,
                 fam_min_identity: float | None = 0.50) -> dict:
    """Par(k) (DESIGN §5): kinase_acc -> set of OTHER kinase accessions with kinase-domain identity >= `identity`
    (KinCore MSA) OR the same Manning subfamily; for kinases without a subfamily, the same Manning family
    ('*-Unique' pseudo-families never group) AND (A8, DESIGN §16) >= fam_min_identity (50%) KinCore domain identity
    (kinases absent from KinCore - atypicals - therefore get no family-fallback paralogs unless fam_min_identity=None).
    Symmetrised (j in Par(k) <=> k in Par(j))."""
    from . import index
    k = kinases if kinases is not None else index.load_kinases()
    par = {a: set() for a in k.uniprot_acc}
    ide = kinase_domain_identity()
    ide_d = dict(zip(zip(ide.acc_a, ide.acc_b), ide.identity))
    for a, b in zip(*[ide[ide.identity >= identity][c] for c in ("acc_a", "acc_b")]):
        if a in par and b in par:
            par[a].add(b), par[b].add(a)
    fam = pd.Series([_real_family(f) for f in k.manning_family], index=k.index, dtype=object)
    sub = pd.Series([s if isinstance(s, str) and s else None for s in k.manning_subfamily], index=k.index,
                    dtype=object)
    by_sub = k.assign(f=fam, s=sub).dropna(subset=["s"]).groupby(["f", "s"]).uniprot_acc.agg(list)
    by_fam = k.assign(f=fam).dropna(subset=["f"]).groupby("f").uniprot_acc.agg(list)
    for a, f, s in zip(k.uniprot_acc, fam, sub):
        if s is not None and f is not None:  # noqa: E501
            mem = by_sub.get((f, s), [])
        elif f is not None:
            mem = by_fam.get(f, [])
        else:
            mem = []
        for b in mem:
            if b == a:
                continue
            if s is None and f is not None and fam_min_identity is not None:
                # A8 (DESIGN §16): the family fallback requires >= fam_min_identity domain identity
                if ide_d.get((min(a, b), max(a, b)), 0.0) < fam_min_identity:
                    continue
            par[a].add(b), par[b].add(a)
    return par


# ------------------------------------------------------------------------------------------------ negative masks
class NegativeMasks:
    """Per-kinase negative exclusion over the canonical site index (DESIGN §5).

    excluded(k) = KP(k) (incl. identical 15-mers) U {U_obs sites within Hamming <= 2 of a KP(k) site's 15-mer}
                  U {sites labelled (any source) to a kinase in Par(k)}
    hard_pool    : per site, the set of kinases with a human-native PSP non-autophosphorylation label ('eval' labels;
                   PI r2): a site is an E-hard negative for k if it has such a label to a kinase outside {k} U Par(k).
    `pairs` must be the unfiltered label table relevant for the evaluation (CV: load_dev('all'); TEST: DEV U TEST).
    Arrays are site-index positions (row order of cache/index/sites.parquet)."""

    def __init__(self, pairs: pd.DataFrame, paralogs: dict | None = None, max_d: int = HAMMING_MAX,
                 include_sensitivity: bool = False):
        from . import index
        self.sites = index.load_sites()
        self.spos = pd.Series(np.arange(len(self.sites)), index=self.sites.site_id)
        self.par = paralogs if paralogs is not None else paralog_sets()
        p = pairs.dropna(subset=["kinase_acc"])
        if not include_sensitivity and "mask_sensitivity_only" in p:   # A7: text-mining/KEA rows: sensitivity only
            p = p[~p.mask_sensitivity_only.fillna(False).astype(bool)]
        self.kp = known_positive_sets(p, expand_identical=False)
        lab_sites = sorted(set(p.site_id))
        w = site_windows(lab_sites)
        nb = hamming_neighbors(list(w), list(self.sites.window15), max_d)
        self._nb = dict(zip(lab_sites, nb))
        self._own = {s: self.spos.get(s) for s in lab_sites}
        self._excl_cache: dict = {}
        hp = p[evaluable_positive(p)]          # PI r2: E-hard membership from human-native PSP labels
        hp = hp[hp.site_id.isin(self.spos.index)]
        self.hard_labelers = hp.groupby("site_id").kinase_acc.agg(lambda s: frozenset(s))
        self.hard_idx = self.spos[self.hard_labelers.index].values

    def _site_block(self, sites_) -> np.ndarray:
        parts = [self._nb[s] for s in sites_]
        parts += [np.array([self._own[s]]) for s in sites_ if self._own[s] is not None]
        return np.unique(np.concatenate(parts)).astype(np.int64) if parts else np.zeros(0, np.int64)

    def kp_block(self, k: str) -> np.ndarray:
        """Index positions masked because of KP(k): KP sites, identical and Hamming<=2 15-mers."""
        return self._site_block(self.kp.get(k, ()))

    def excluded(self, k: str) -> np.ndarray:
        if k not in self._excl_cache and str(k).startswith("CANARY_"):
            # canary v2: KP mask = the canary's own sites only (no Hamming / paralog expansion)
            self._excl_cache[k] = np.unique(self.spos.reindex(sorted(self.kp.get(k, ()))).dropna()
                                            .values.astype(np.int64))
        if k not in self._excl_cache:
            sites_ = set(self.kp.get(k, ()))
            par_sites = set()
            for j in self.par.get(k, ()):
                par_sites |= self.kp.get(j, set())
            a = self.kp_block(k)
            b = self.spos.reindex(sorted(par_sites)).dropna().values.astype(np.int64)
            self._excl_cache[k] = np.union1d(a, b)
        return self._excl_cache[k]

    def hard_candidates(self, k: str) -> np.ndarray:
        """Index positions with >= 1 human-native non-autophos label to a kinase outside {k} U Par(k)."""
        own = {k} | set(self.par.get(k, ()))
        ok = np.fromiter((len(ls - own) > 0 for ls in self.hard_labelers.values), bool, len(self.hard_labelers))
        return self.hard_idx[ok]


# ------------------------------------------------------------------------------------------------ temporal / external
def psp2021_keys(path: Path = PSP2021) -> set:
    """Keys of the PSP Kinase_Substrate_Dataset bundled with kinase-library 1.8.0 (07/2021), ALL species (A10):
    {(KIN_ACC, SUB_ACC, WINDOW15_upper)}, {(KIN_ACC, SUB_ACC, residue+pos)}, and {('G', KINASE_GENE_upper, SGID)}
    (kinase symbol upper-cased = ortholog assumption, as in the label build) for cross-species matching."""
    d = pd.read_csv(path, sep="\t", skiprows=3, dtype=str, encoding="latin1", keep_default_na=False)
    ka = d.KIN_ACC_ID.str.split("-").str[0]
    sa = d.SUB_ACC_ID.str.split("-").str[0]
    keys = set(zip(ka, sa, d["SITE_+/-7_AA"].str.upper()))
    keys |= set(zip(ka, sa, d.SUB_MOD_RSD.str.upper()))
    keys |= {("G", g, sg) for g, sg in zip(d.GENE.str.upper(), d.SITE_GRP_ID) if sg}
    return keys


def new_since_2021(pairs: pd.DataFrame, keys: set | None = None) -> pd.Series:
    """True for PSP human-native pairs absent from the 2021 PSP release in ALL species (A10): no 2021 row with the
    same kinase accession + substrate accession + (15-mer or residue/position), and no 2021 row with the same
    kinase symbol (upper) + SITE_GRP_ID (any organism). Non-PSP / non-native -> False.
    ASSUMPTION: 2021 KIN_ACC_ID / symbols are comparable with today's (kinase symbol = ortholog, as in the build)."""
    keys = keys if keys is not None else psp2021_keys()
    w = site_windows(pairs.site_id).str.upper().values
    rp = pairs.site_id.str.rsplit("_", n=1).str[1].values
    subs = pairs.psp_sub_accs.fillna("").str.split(";").values
    sg = (pairs.psp_sgids.fillna("") if "psp_sgids" in pairs else pd.Series("", index=pairs.index)).str.split(";")
    out = []
    for ka, kg, acc, ww, r, ss, sgs, prov, nat in zip(pairs.kinase_acc, pairs.kinase.astype(str).str.upper(),
                                                     pairs.acc, w, rp, subs, sg.values, pairs.provenance,
                                                     pairs.native_human):
        if prov != "PSP" or not nat:
            out.append(False)
            continue
        cand = {acc} | {s.split("-")[0] for s in ss if s}
        old = any((ka, a, ww) in keys or (ka, a, r) in keys for a in cand)
        old = old or any(("G", kg, x.split(".")[0]) in keys for x in sgs if x)
        out.append(not old)
    return pd.Series(out, index=pairs.index, dtype=bool)


def external_test_filter(pairs: pd.DataFrame, paralogs: dict | None = None, window: int = 5) -> pd.Series:
    """External curated stratum (DESIGN §3, DL-5): EXT_not_in_psp pairs (OmniPath curated / UniProt ECO:0000269,
    in class, included kinase, single-kinase names) minus near-copies of PSP: (i) the same kinase + protein has a PSP
    pair within +-`window` residues (renumbering), (ii) the same site is a PSP pair of a paralog (Par) of the kinase.
    (iii) A9: OmniPath rows listing PhosphoSite among their sources. Then ext_test_tiers() applies the per-PMID cap and
    tier A/B. Returns a bool mask over `pairs` (False for non-EXT rows)."""
    par = paralogs if paralogs is not None else paralog_sets()
    psp = pairs[pairs.provenance.eq("PSP")]
    pos_by = psp.groupby(["kinase_acc", "acc"]).pos.agg(lambda s: np.array(sorted(set(s)))).to_dict()
    lists_psp = pairs.omnipath_lists_psp.fillna(False).astype(bool) if "omnipath_lists_psp" in pairs else \
        pd.Series(False, index=pairs.index)
    kin_by_site = psp.groupby("site_id").kinase_acc.agg(set).to_dict()
    out = []
    for ka, acc, pos, sid, ext, lp in zip(pairs.kinase_acc, pairs.acc, pairs.pos, pairs.site_id,
                                          pairs.EXT_not_in_psp.astype(bool), lists_psp):
        if not ext or lp:                      # A9: drop OmniPath rows that list PhosphoSite among sources
            out.append(False)
            continue
        near = pos_by.get((ka, acc))
        if near is not None and np.any((np.abs(near - pos) <= window) & (near != pos)):
            out.append(False)
            continue
        if kin_by_site.get(sid, set()) & par.get(ka, set()):
            out.append(False)
            continue
        out.append(True)
    return pd.Series(out, index=pairs.index, dtype=bool)


def ext_test_tiers(pairs: pd.DataFrame, keep: pd.Series, cap_per_pmid: int = 20, seed: int = 0) -> pd.Series:
    """A9 (DESIGN §16): among rows with keep=True, cap 20 pairs per PMID (pairs supported by any PMID with <= cap
    pairs are kept; pairs supported only by larger PMIDs: at most `cap` per PMID, chosen by a seeded shuffle), then
    tier 'A' = SIGNOR / HPRD (OmniPath curated) or UniProt ECO:0000269 (single named kinase + named residue - all
    EXT rows have both), tier 'B' = other curated sources (phosphoELM, ProtMapper-derived). Returns '' / 'A' / 'B'."""
    x = pairs[keep].copy()
    refs = x.ext_refs.fillna("").str.split(";").map(lambda r: [t for t in r if t])
    cnt = pd.Series([t for r in refs for t in r]).value_counts()
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(x))
    used: dict = {}
    ok = np.zeros(len(x), bool)
    for i in order:
        r = refs.iloc[i]
        if not r or any(cnt[t] <= cap_per_pmid for t in r):
            ok[i] = True
            continue
        big = [t for t in r if used.get(t, 0) < cap_per_pmid]
        if big:
            ok[i] = True
            for t in r:
                used[t] = used.get(t, 0) + 1
    src = x.ext_sources.fillna("")
    tierA = src.str.contains(r"(?:^|;)(?:SIGNOR|HPRD)(?:;|$)") | x.provenance.eq("UniProt")
    tier = pd.Series("", index=pairs.index, dtype=object)
    tier.loc[x.index[ok]] = np.where(tierA.values[ok], "A", "B")
    return tier


MASK_PROVENANCES = ("MASK:PSP_subunit", "MASK:EXT_subunit", "MASK:OmniPath_subunit", "MASK:UniProt_family",
                    "MASK:OmniPath_textmining")


def mask_only_pairs(allp: pd.DataFrame, olog: pd.DataFrame, ulog: pd.DataFrame, store, res,
                    textmining: pd.DataFrame | None = None) -> pd.DataFrame:
    """A7 (DESIGN §16) mask-only rows (never positives, never negatives): kinase_acc, site_id, provenance, flags.
    (1) PSP/EXT pairs whose kinase is a regulatory/complex subunit (config/complex_subunits.tsv) -> catalytic kinases;
    (2) OmniPath rows whose enzyme is such a subunit (unresolved enzymes) -> catalytic kinases;
    (3) UniProt ECO:0000269 family-level notes ('by PKA', 'by CK2') -> all family members;
    (4) sensitivity-only: OmniPath rows without a curated non-PSP source (text mining, KEA, ...) resolved by symbol,
        flagged mask_sensitivity_only=True (used only when NegativeMasks(include_sensitivity=True))."""
    from .kinases import family_members, load_complex_subunits, resolve_kinase_names
    cx = load_complex_subunits()
    fam = family_members(res)
    acc_of = dict(zip(res.kl.gene, res.kl.uniprot_acc))
    rec = []
    sub = allp[allp.kinase.astype(str).str.upper().isin(cx)]
    for g, sid, prov in zip(sub.kinase.str.upper(), sub.site_id, sub.provenance):
        for c in cx[g]:
            rec.append((c, sid, "MASK:PSP_subunit" if prov == "PSP" else "MASK:EXT_subunit", False, g))
    for lg, kind in [(olog, "omni"), (ulog, "uni")]:
        if lg is None or not len(lg):
            continue
        if kind == "omni":
            x = lg[lg.reason.str.startswith("kinase_") & lg.kinase_name.str.upper().isin(cx)]
            for g, a, st in zip(x.kinase_name.str.upper(), x.substrate, x.site):
                m = store.map_position(a, int(st[1:]), st[0])
                if m.site_id:
                    for c in cx[g]:
                        rec.append((c, m.site_id, "MASK:OmniPath_subunit", False, g))
        else:
            x = lg[lg.reason.eq("kinase_family_name")]
            rr, _ = resolve_kinase_names(x.kinase_name, "uniprot", res)
            for mem, a, st, nm in zip(rr.family_members, x.substrate, x.site, x.kinase_name):
                if not mem:
                    mem = " ".join(fam.get(str(nm).upper(), []))
                if not mem:
                    continue
                m = store.map_position(a, int(st[1:]), st[0])
                if m.site_id:
                    for c in mem.split():
                        rec.append((c, m.site_id, "MASK:UniProt_family", False, nm))
    if textmining is not None and len(textmining):
        rr, _ = resolve_kinase_names(textmining.kinase_name, "omnipath", res)
        for g, a, st in zip(rr.gene, textmining.substrate, textmining.site):
            if not isinstance(g, str):
                continue
            try:
                m = store.map_position(a, int(st[1:]), st[0])
            except (ValueError, KeyError):
                continue
            if m.site_id:
                rec.append((g, m.site_id, "MASK:OmniPath_textmining", True, g))
    d = pd.DataFrame(rec, columns=["kinase", "site_id", "provenance", "mask_sensitivity_only", "mask_source_name"])
    d = d[d.kinase.isin(acc_of)].drop_duplicates(["kinase", "site_id", "provenance"]).reset_index(drop=True)
    d["mask_only"] = True
    return d
