#!/usr/bin/env python
"""Build F6 data blocks (deliverable 1): compendium study blocks + CPTAC residualised tumour cohorts.

Outputs
  cache/features/phospho/conditions.parquet        condition metadata (PMID, block_id, perturbation type, minutes, scale)
  cache/features/phospho/blocks/*.npy|parquet      one block per study / cohort (+ manifest.parquet)
  cache/features/phospho/cptac_kinase_protein.npz  per cohort: kinase idx + confounder-residualised kinase protein
  results/features/phospho/blocks_summary.tsv, cptac_residualisation.tsv, conditions_by_block.tsv

Usage: CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_build_blocks.py
       [--skip-compendium] [--skip-cptac] [--min-block-cols 8]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.features import phospho as ph  # noqa: E402
from ksp.features import phospho_cptac as pc  # noqa: E402
from profile_compendium import classify_condition  # noqa: E402

DECRYPTM = 36926954


def build_conditions(min_block_cols: int) -> pd.DataFrame:
    c = pd.read_csv(ROOT / "data" / "01_phosphoproteomics_conditions.csv")
    c = c.drop(columns=[x for x in c.columns if x.startswith("Unnamed") or x == ""], errors="ignore")
    c["PMID"] = pd.to_numeric(c.PMID, errors="coerce").astype("Int64")
    c["ptype"] = c.Condition.map(classify_condition)
    t = [ph.parse_minutes(tp, cond) for tp, cond in zip(c.Timepoint, c.Condition)]
    c["minutes"] = [x[0] for x in t]
    c["minutes_how"] = [x[1] for x in t]
    # blocks: PMID; decryptM split per drug series; small studies pooled
    blk = "pmid:" + c.PMID.astype(str)
    dm = c.PMID == DECRYPTM
    drug = c.Condition.astype(str).str.split("_").str[0]
    blk[dm] = "pmid:36926954:" + drug[dm]
    n = blk.map(blk.value_counts())
    small_dm = dm & (n < min_block_cols)
    blk[small_dm] = "pmid:36926954:other"
    n = blk.map(blk.value_counts())
    blk[(~dm) & (n < min_block_cols)] = "pooled_small"
    c["block_id"] = blk
    return c


def build_compendium(c: pd.DataFrame, sites: pd.DataFrame, threads: int = 10) -> list[dict]:
    t0 = time.time()
    tab = pq.read_table(ROOT / "cache" / "compendium_wide.parquet", use_threads=True)
    uids = tab.column("peptide_UID").to_pylist()
    cols = [x for x in tab.column_names if x != "peptide_UID"]
    assert set(cols) == set(c.conditionID), "condition mismatch"
    X = np.column_stack([tab.column(k).to_numpy(zero_copy_only=False).astype(np.float32) for k in c.conditionID])
    del tab
    print(f"read wide {X.shape} in {time.time()-t0:.1f}s", flush=True)
    m = pd.read_parquet(ROOT / "cache" / "sites" / "mapping_log.parquet", columns=["source", "src_id", "site_id"])
    m = m[(m.source == "compendium") & m.site_id.notna()].set_index("src_id").site_id
    spos = {s: j for j, s in enumerate(sites.site_id)}
    sidx = pd.Series(uids).map(m).map(spos)
    ok = sidx.notna().to_numpy()
    print(f"uids mapped to index sites: {ok.sum()} / {len(uids)}", flush=True)
    X = X[ok]
    sidx = sidx[ok].astype(np.int64).to_numpy()
    # robust per-condition scale on the raw values (all mapped uids)
    med = np.nanmedian(X, 0)
    mad = np.nanmedian(np.abs(X - med), 0) * 1.4826
    c["scale"] = mad
    c["n_sites_raw"] = np.isfinite(X).sum(0)
    # aggregate duplicate uids per site (nan-mean)
    uniq, inv = np.unique(sidx, return_inverse=True)
    S = np.zeros((len(uniq), X.shape[1]), np.float32)
    N = np.zeros((len(uniq), X.shape[1]), np.float32)
    fin = np.isfinite(X)
    np.add.at(S, inv, np.where(fin, X, 0))
    np.add.at(N, inv, fin.astype(np.float32))
    with np.errstate(invalid="ignore"):
        A = np.where(N > 0, S / np.maximum(N, 1), np.nan).astype(np.float32)
    A /= mad[None, :]
    print(f"aggregated to {len(uniq)} sites ({(N > 1).any(1).sum()} with merged uids) {time.time()-t0:.1f}s",
          flush=True)
    rep = []
    for bid, g in c.groupby("block_id", sort=True):
        j = g.index.to_numpy()
        Xb = A[:, j]
        rowok = np.isfinite(Xb).any(1)
        b = ph.Block(bid, "compendium", uniq[rowok].astype(np.int32), Xb[rowok], g.conditionID.tolist(),
                     tuple(sorted(int(p) for p in g.PMID.dropna().unique())), None)
        r = ph.save_block(b)
        r["median_obs_per_site"] = float(np.median(np.isfinite(b.X).sum(1)))
        r["n_sites_ge8"] = int((np.isfinite(b.X).sum(1) >= 8).sum())
        rep.append(r)
    return rep


def build_cptac(sites: pd.DataFrame, kin: pd.DataFrame, device: str) -> tuple[list[dict], list[dict]]:
    from ksp.kinases import load_hgnc
    spos = {s: j for j, s in enumerate(sites.site_id)}
    idmap = pc.cohort_site_map()
    hg = load_hgnc()
    acc2ensg = {}
    for u, e in zip(hg.uniprot_ids, hg.ensembl_gene_id):
        for a in u.split("|"):
            if a and e:
                acc2ensg.setdefault(a, e)
    kin_ensg = [acc2ensg.get(a, "") for a in kin.uniprot_acc]
    reps, man = [], []
    kp = {}
    for coh in pc.COHORTS:
        t0 = time.time()
        blk, rep, pr, C, samples = pc.build_cohort_block(coh, spos, idmap, device=device)
        r = ph.save_block(blk)
        man.append(r)
        # kinase protein abundance (residualised on confounders only)
        pr_base = pr.copy()
        pr_base.index = [i.split(".")[0] for i in pr_base.index]
        pr_base = pr_base[~pr_base.index.duplicated()]
        ki = np.array([i for i, e in enumerate(kin_ensg) if e and e in pr_base.index], dtype=np.int64)
        K = pr_base.reindex([kin_ensg[i] for i in ki]).to_numpy(np.float32)
        Kr, _ = pc.residualise(K, None, C, min_obs=20, device=device)
        good = np.isfinite(Kr).sum(1) >= 20
        kp[f"{coh}_kidx"] = ki[good]
        kp[f"{coh}_X"] = Kr[good]
        rep["n_kinases_with_protein"] = int(good.sum())
        rep["seconds"] = round(time.time() - t0, 1)
        reps.append(rep)
        print(json.dumps(rep), flush=True)
    np.savez(ph.CACHE / "cptac_kinase_protein.npz", **kp)
    return reps, man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-block-cols", type=int, default=8)
    ap.add_argument("--skip-compendium", action="store_true")
    ap.add_argument("--skip-cptac", action="store_true")
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    ph.CACHE.mkdir(parents=True, exist_ok=True)
    ph.RESULTS.mkdir(parents=True, exist_ok=True)
    sites = ph.site_index()
    kin = ph.kinase_index()
    man_path = ph.BLOCK_DIR / "manifest.parquet"
    old = pd.read_parquet(man_path) if man_path.exists() else pd.DataFrame(columns=["block_id", "kind"])
    man = []
    if not a.skip_compendium:
        c = build_conditions(a.min_block_cols)
        rep = build_compendium(c, sites)
        c.to_parquet(ph.CACHE / "conditions.parquet")
        man += rep
        cb = c.groupby("block_id").agg(n_conditions=("conditionID", "size"), pmids=("PMID", lambda x: ";".join(
            map(str, sorted(set(x.dropna()))))), ptypes=("ptype", lambda x: "; ".join(
                f"{k}:{v}" for k, v in x.value_counts().items())), samples=("Sample", "nunique")).reset_index()
        cb.to_csv(ph.RESULTS / "conditions_by_block.tsv", sep="\t", index=False)
        c[["conditionID", "PMID", "block_id", "ptype", "Timepoint", "minutes", "minutes_how", "scale"]].to_csv(
            ph.RESULTS / "conditions_meta.tsv", sep="\t", index=False)
    else:
        man += old[old.kind == "compendium"].to_dict("records")
    if not a.skip_cptac:
        reps, m2 = build_cptac(sites, kin, a.device)
        man += m2
        pd.DataFrame(reps).to_csv(ph.RESULTS / "cptac_residualisation.tsv", sep="\t", index=False)
    else:
        man += old[old.kind == "cptac"].to_dict("records")
    man = pd.DataFrame(man)
    man.to_parquet(man_path)
    man.to_csv(ph.RESULTS / "blocks_summary.tsv", sep="\t", index=False)
    print(man.groupby("kind").agg(n=("block_id", "size"), cols=("n_cols", "sum"), vals=("n_values", "sum")))


if __name__ == "__main__":
    main()
