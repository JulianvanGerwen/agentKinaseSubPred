"""Loader for the PhosphoSitePlus Kinase_Substrate_Dataset.

File format (release header "051925"): 3 preamble lines (release date, licence text, blank),
then a tab-separated table with header; encoding latin1.

Columns in the raw file -> parsed columns
  GENE          -> kin_gene        kinase gene symbol (organism-specific capitalisation, e.g. 'Dyrk2')
  KINASE        -> kin_name        PSP kinase protein name (e.g. 'PKACA')
  KIN_ACC_ID    -> kin_acc_full, kin_acc (base accession), kin_isoform (suffix '' / '2' ...)
  KIN_ORGANISM  -> kin_organism
  SUBSTRATE     -> sub_name
  SUB_GENE_ID   -> sub_gene_id     (Entrez gene id, str)
  SUB_ACC_ID    -> sub_acc_full, sub_acc (base accession), sub_isoform (suffix, '' for canonical)
  SUB_GENE      -> sub_gene
  SUB_ORGANISM  -> sub_organism
  SUB_MOD_RSD   -> residue (one letter), position (int), mod_rsd (raw, e.g. 'S336')
  SITE_GRP_ID   -> site_grp_id     (Int64; PSP site group = orthologous sites across species/isoforms)
  SITE_+/-7_AA  -> window_raw (15-mer as in file), window (upper-case, '_' padding kept),
                   center_lower (centre char is lower-case = recorded as modified),
                   n_other_phos_in_window (lower-case residues in window other than centre),
                   other_phos_offsets (tuple of offsets -7..7 of other lower-case residues)
  DOMAIN        -> domain
  IN_VIVO_RXN   -> in_vivo  (bool, 'X' in file)
  IN_VITRO_RXN  -> in_vitro (bool)
  CST_CAT#      -> cst_cat

Derived:
  is_sty         residue in {S,T,Y}
  evidence       'both' / 'in_vivo_only' / 'in_vitro_only' / 'none'
  kin_gene_uc / sub_gene_uc : upper-cased gene symbols (used for cross-species symbol mapping)
  autophos       kinase and substrate are the same gene (upper-cased symbol) in the same organism

No rows are dropped by default; use `filter_psp` for the usual cleaning steps.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_PATH = Path(os.environ.get("KSP_DATA", Path(__file__).resolve().parents[2] / "data")) / "Kinase_Substrate_Dataset"
N_PREAMBLE = 3
_RSD = re.compile(r"^([A-Za-z])(\d+)$")


def read_release(path: str | Path = DEFAULT_PATH) -> str:
    """Return the release tag on the first line of the file (e.g. '051925' = 2025-05-19)."""
    with open(path, encoding="latin1") as f:
        return f.readline().strip()


def _split_acc(s: pd.Series) -> tuple[pd.Series, pd.Series]:
    parts = s.str.split("-", n=1, expand=True)
    base = parts[0]
    iso = parts[1] if parts.shape[1] > 1 else pd.Series("", index=s.index)
    return base, iso.fillna("")


def load_psp(path: str | Path = DEFAULT_PATH) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", skiprows=N_PREAMBLE, encoding="latin1", dtype=str,
                      keep_default_na=False)
    raw.columns = [c.strip() for c in raw.columns]
    expected = ["GENE", "KINASE", "KIN_ACC_ID", "KIN_ORGANISM", "SUBSTRATE", "SUB_GENE_ID", "SUB_ACC_ID",
                "SUB_GENE", "SUB_ORGANISM", "SUB_MOD_RSD", "SITE_GRP_ID", "SITE_+/-7_AA", "DOMAIN",
                "IN_VIVO_RXN", "IN_VITRO_RXN", "CST_CAT#"]
    missing = [c for c in expected if c not in raw.columns]
    if missing:
        raise ValueError(f"PSP file missing columns {missing}; got {raw.columns.tolist()}")
    raw = raw.apply(lambda c: c.str.strip())
    d = pd.DataFrame(index=raw.index)
    d["kin_gene"] = raw["GENE"]
    d["kin_name"] = raw["KINASE"]
    d["kin_acc_full"] = raw["KIN_ACC_ID"]
    d["kin_acc"], d["kin_isoform"] = _split_acc(raw["KIN_ACC_ID"])
    d["kin_organism"] = raw["KIN_ORGANISM"]
    d["sub_name"] = raw["SUBSTRATE"]
    d["sub_gene_id"] = raw["SUB_GENE_ID"]
    d["sub_acc_full"] = raw["SUB_ACC_ID"]
    d["sub_acc"], d["sub_isoform"] = _split_acc(raw["SUB_ACC_ID"])
    d["sub_gene"] = raw["SUB_GENE"]
    d["sub_organism"] = raw["SUB_ORGANISM"]
    d["mod_rsd"] = raw["SUB_MOD_RSD"]
    m = raw["SUB_MOD_RSD"].str.extract(_RSD)
    d["residue"] = m[0].str.upper()
    d["position"] = pd.to_numeric(m[1], errors="coerce").astype("Int64")
    d["site_grp_id"] = pd.to_numeric(raw["SITE_GRP_ID"], errors="coerce").astype("Int64")
    w = raw["SITE_+/-7_AA"]
    d["window_raw"] = w
    d["window"] = w.str.upper()
    d["center_lower"] = w.str[7].str.islower()
    offs = w.map(lambda s: tuple(i - 7 for i, ch in enumerate(s) if ch.islower() and i != 7) if len(s) == 15 else ())
    d["other_phos_offsets"] = offs
    d["n_other_phos_in_window"] = offs.map(len).astype(int)
    d["domain"] = raw["DOMAIN"]
    d["in_vivo"] = raw["IN_VIVO_RXN"].str.upper().eq("X")
    d["in_vitro"] = raw["IN_VITRO_RXN"].str.upper().eq("X")
    d["cst_cat"] = raw["CST_CAT#"]
    d["is_sty"] = d["residue"].isin(["S", "T", "Y"])
    d["evidence"] = np.select([d.in_vivo & d.in_vitro, d.in_vivo, d.in_vitro],
                              ["both", "in_vivo_only", "in_vitro_only"], "none")
    d["kin_gene_uc"] = d["kin_gene"].str.upper()
    d["sub_gene_uc"] = d["sub_gene"].str.upper()
    d["autophos"] = (d.kin_gene_uc == d.sub_gene_uc) & (d.kin_organism == d.sub_organism)
    d["window_ok"] = (w.str.len() == 15) & (d["window"].str[7] == d["residue"])
    return d


def filter_psp(d: pd.DataFrame, sty_only: bool = True, kin_organisms=None, sub_organisms=None,
               evidence=None, canonical_sub_only: bool = False) -> pd.DataFrame:
    """Common filters. evidence: None (any), or a subset of {'both','in_vivo_only','in_vitro_only'}."""
    m = pd.Series(True, index=d.index)
    if sty_only:
        m &= d.is_sty
    if kin_organisms is not None:
        m &= d.kin_organism.isin(list(kin_organisms))
    if sub_organisms is not None:
        m &= d.sub_organism.isin(list(sub_organisms))
    if evidence is not None:
        m &= d.evidence.isin(list(evidence))
    if canonical_sub_only:
        m &= d.sub_isoform.eq("")
    return d[m].copy()


def human_kinase_symbols(d: pd.DataFrame) -> set[str]:
    """Upper-case gene symbols of kinases with at least one human-kinase row in PSP."""
    return set(d.loc[d.kin_organism == "human", "kin_gene_uc"])


def map_kinases_to_human(d: pd.DataFrame, human_symbols: set[str] | None = None,
                         symbol_to_acc: dict | None = None) -> pd.DataFrame:
    """Add kin_human_gene / kin_human_mapped. Non-human kinases are mapped to a human kinase when the
    upper-cased GENE symbol matches a human kinase symbol (orthology by symbol; ASSUMPTION: same
    symbol => 1:1 ortholog; not checked by sequence). `human_symbols` defaults to the human kinase
    symbols present in PSP; pass a UniProt-derived set to extend. `symbol_to_acc` optionally maps the
    human symbol to a human UniProt accession (human rows keep their own accession)."""
    hs = human_kinase_symbols(d) if human_symbols is None else set(human_symbols)
    out = d.copy()
    # rows with empty GENE (KIN_ORGANISM 'no information'): fall back to the PSP kinase name, mapped to the
    # gene symbol it carries in human rows when that name is unambiguous.
    hn = d[d.kin_organism.eq("human") & d.kin_gene_uc.ne("")].groupby(d.kin_name.str.upper()).kin_gene_uc.unique()
    name2gene = {k: v[0] for k, v in hn.items() if len(v) == 1}
    sym = out.kin_gene_uc.where(out.kin_gene_uc.ne(""), out.kin_name.str.upper().map(name2gene)).fillna("")
    out["kin_gene_from_name"] = out.kin_gene_uc.eq("") & sym.ne("")
    is_h = out.kin_organism.eq("human")
    mapped = is_h | (sym.isin(hs) & sym.ne(""))
    out["kin_human_gene"] = np.where(mapped, sym, None)
    out["kin_human_mapped"] = mapped
    if symbol_to_acc is not None:
        acc = out.kin_human_gene.map(symbol_to_acc)
        out["kin_human_acc"] = np.where(is_h, out.kin_acc, acc)
    return out


def human_site_table(d: pd.DataFrame) -> pd.DataFrame:
    """Human substrate sites in PSP keyed by site group: site_grp_id -> human (sub_acc_full, residue, position, window).
    A site group can contain >1 human site (isoforms); all are returned."""
    h = d[(d.sub_organism == "human") & d.site_grp_id.notna()]
    return (h[["site_grp_id", "sub_acc_full", "sub_acc", "sub_isoform", "sub_gene", "residue", "position", "window"]]
            .drop_duplicates(["site_grp_id", "sub_acc_full", "position"]))


def map_substrates_to_human(d: pd.DataFrame, extra_human_sites: pd.DataFrame | None = None) -> pd.DataFrame:
    """Map every row's substrate site to human site(s) via SITE_GRP_ID, using human sites seen anywhere in
    PSP (optionally plus `extra_human_sites` with columns site_grp_id, sub_acc_full, sub_acc, residue,
    position). Human rows map to themselves. Returns rows exploded over human sites, with columns
    h_sub_acc_full, h_sub_acc, h_residue, h_position, sub_human_mapped (unmapped rows kept with NaN)."""
    hs = human_site_table(d)
    if extra_human_sites is not None:
        hs = pd.concat([hs, extra_human_sites], ignore_index=True).drop_duplicates(
            ["site_grp_id", "sub_acc_full", "position"])
    hs = hs.rename(columns={"sub_acc_full": "h_sub_acc_full", "sub_acc": "h_sub_acc", "residue": "h_residue",
                            "position": "h_position"})[["site_grp_id", "h_sub_acc_full", "h_sub_acc", "h_residue",
                                                        "h_position"]]
    out = d.copy()
    out["_row"] = np.arange(len(out))
    is_h = out.sub_organism.eq("human")
    hum = out[is_h].assign(h_sub_acc_full=out.sub_acc_full, h_sub_acc=out.sub_acc, h_residue=out.residue,
                           h_position=out.position)
    non = out[~is_h].merge(hs, on="site_grp_id", how="left")
    res = pd.concat([hum, non], ignore_index=True).sort_values("_row", kind="stable").drop(columns="_row")
    res["sub_human_mapped"] = res.h_sub_acc_full.notna()
    return res.reset_index(drop=True)


def _best_window_hit(window: str, residue: str, seq: str, prefer_pos: int | None = None):
    """Ungapped scan of a 15-mer (centre = residue, '_' = padding/wildcard) along `seq`.
    Returns (best_pos_1based, n_identical, n_compared) or None. Ties are broken by closeness to prefer_pos."""
    L = len(seq)
    if L == 0:
        return None
    s = np.frombuffer(("_" * 7 + seq + "_" * 7).encode(), dtype=np.uint8)
    w = np.frombuffer(window.encode(), dtype=np.uint8)
    centres = np.flatnonzero(np.frombuffer(seq.encode(), dtype=np.uint8) == ord(residue))
    if len(centres) == 0:
        return None
    idx = centres[:, None] + np.arange(15)[None, :]
    sub = s[idx]
    cmp_mask = w != ord("_")
    ident = ((sub == w[None, :]) & cmp_mask[None, :]).sum(1)
    best = ident.max()
    cand = centres[ident == best]
    if prefer_pos is not None and len(cand) > 1:
        cand = cand[np.argsort(np.abs(cand + 1 - prefer_pos), kind="stable")]
    return int(cand[0]) + 1, int(best), int(cmp_mask.sum())


def map_sites_by_sequence(sites: pd.DataFrame, human_prot: pd.DataFrame, min_identity: int = 11) -> pd.DataFrame:
    """Map (non-human) substrate sites to human by sequence, as a fallback when the site group has no
    human member in the KS dataset. Candidate human proteins = `human_prot` rows (columns gene_uc, acc_full,
    acc, seq) whose upper-case gene symbol equals the site's sub_gene_uc (ASSUMPTION: same symbol =>
    ortholog). The 15-mer is scanned ungapped; the hit needs the same centre residue and >= min_identity
    identical residues among non-padding window positions (of 15).
    Input `sites` needs columns sub_gene_uc, window, residue, position. Returns `sites` with columns
    seq_h_sub_acc_full, seq_h_sub_acc, seq_h_position, seq_identity (NaN when no hit)."""
    by_gene = {g: grp for g, grp in human_prot.groupby("gene_uc")}
    out = []
    for r in sites.itertuples(index=False):
        best = None
        grp = by_gene.get(r.sub_gene_uc)
        if grp is not None:
            for p in grp.itertuples(index=False):
                hit = _best_window_hit(r.window, r.residue, p.seq, int(r.position) if pd.notna(r.position) else None)
                if hit and hit[1] >= min_identity and (best is None or hit[1] > best[3]):
                    best = (p.acc_full, p.acc, hit[0], hit[1])
        out.append(best if best else (None, None, None, None))
    res = sites.copy()
    o = pd.DataFrame(out, columns=["seq_h_sub_acc_full", "seq_h_sub_acc", "seq_h_position", "seq_identity"],
                     index=sites.index)
    o["seq_h_position"] = o["seq_h_position"].astype("Int64")
    return pd.concat([res, o], axis=1)


def check_windows_against_sequences(sites: pd.DataFrame, seqs: dict, acc_col: str = "sub_acc_full") -> pd.DataFrame:
    """Compare each site's 15-mer with the sequence `seqs[acc]` (dict acc_full -> sequence) at `position`.
    Adds columns: seq_found, pos_in_range, centre_match, window_match (exact over non-padding positions),
    window_identity (n identical / n compared), and status in
    {'no_sequence','pos_out_of_range','exact','centre_only','centre_mismatch'};
    and for non-exact rows, shifted_pos = position of the exact 15-mer elsewhere in the sequence (if unique/any)."""
    rec = []
    for acc, pos, win in zip(sites[acc_col], sites["position"], sites["window"]):
        seq = seqs.get(acc)
        if seq is None:
            rec.append((False, False, False, False, np.nan, "no_sequence", pd.NA))
            continue
        pos = int(pos)
        if pos < 1 or pos > len(seq):
            status, cm, wm, ident = "pos_out_of_range", False, False, np.nan
        else:
            padded = "_" * 7 + seq + "_" * 7
            ref = padded[pos - 1: pos + 14]
            comp = [(a, b) for a, b in zip(win, ref) if a != "_"]
            ident = sum(a == b for a, b in comp) / max(len(comp), 1)
            cm = win[7] == seq[pos - 1]
            wm = ident == 1.0
            status = "exact" if wm else ("centre_only" if cm else "centre_mismatch")
        shifted = pd.NA
        if status != "exact":
            core = win.strip("_")
            j = seq.find(core)
            if j >= 0:
                lead = len(win) - len(win.lstrip("_"))
                shifted = j + (7 - lead) + 1
        rec.append((True, status != "pos_out_of_range", cm, wm, ident, status, shifted))
    o = pd.DataFrame(rec, index=sites.index, columns=["seq_found", "pos_in_range", "centre_match", "window_match",
                                                      "window_identity", "status", "shifted_pos"])
    return pd.concat([sites, o], axis=1)
