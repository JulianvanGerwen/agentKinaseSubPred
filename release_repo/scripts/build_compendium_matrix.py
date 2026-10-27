#!/usr/bin/env python
"""Convert the long-format compendium quantifications (03_phosphosite_quantifications.csv)
into a wide float32 matrix: rows = peptide_UID (sites), columns = conditionID.

Outputs (cache/ is gitignored):
  cache/compendium_wide.parquet   index column 'peptide_UID' + one float32 column per condition
                                  (row order = 02_phosphosite_information.csv order; column order =
                                  01_phosphoproteomics_conditions.csv order; NaN = not quantified)
  cache/compendium_long.parquet   the long table as parquet (peptide_UID, conditionID dictionary-encoded)

Duplicate (peptide_UID, conditionID) pairs are checked; if present they are averaged and reported.
Usage: python scripts/build_compendium_matrix.py [--data data] [--out cache] [--threads 16]
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="cache")
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    pa.set_cpu_count(a.threads)
    pa.set_io_thread_count(a.threads)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    cond = pd.read_csv(os.path.join(a.data, "01_phosphoproteomics_conditions.csv"), usecols=["conditionID"])
    site = pd.read_csv(os.path.join(a.data, "02_phosphosite_information.csv"), usecols=["peptide_UID"])
    tab = pacsv.read_csv(
        os.path.join(a.data, "03_phosphosite_quantifications.csv"),
        read_options=pacsv.ReadOptions(block_size=64 << 20, use_threads=True),
        convert_options=pacsv.ConvertOptions(
            include_columns=["peptide_UID", "conditionID", "log2FC"],
            column_types={"peptide_UID": pa.string(), "conditionID": pa.string(), "log2FC": pa.float64()},
            null_values=["NA", "", "NaN"], strings_can_be_null=False),
    )
    t_read = time.time() - t0
    n_long = tab.num_rows
    pq.write_table(tab.combine_chunks().cast(pa.schema([("peptide_UID", pa.dictionary(pa.int32(), pa.string())),
                                                        ("conditionID", pa.dictionary(pa.int32(), pa.string())),
                                                        ("log2FC", pa.float32())])),
                   os.path.join(a.out, "compendium_long.parquet"))
    df = tab.to_pandas()
    del tab
    sites = pd.Index(site["peptide_UID"].astype(str))
    conds = pd.Index(cond["conditionID"].astype(str))
    extra_sites = pd.Index(df["peptide_UID"].unique()).difference(sites)
    extra_conds = pd.Index(df["conditionID"].unique()).difference(conds)
    if len(extra_sites):
        sites = sites.append(extra_sites)
    if len(extra_conds):
        conds = conds.append(extra_conds)
    ri = sites.get_indexer(df["peptide_UID"])
    ci = conds.get_indexer(df["conditionID"])
    v = df["log2FC"].to_numpy(np.float64)
    nan_vals = int(np.isnan(v).sum())
    key = ri.astype(np.int64) * len(conds) + ci
    uniq, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    n_dup_keys = int((cnt > 1).sum())
    M = np.full((len(sites), len(conds)), np.nan, dtype=np.float32)
    if n_dup_keys:
        ok = ~np.isnan(v)
        s = np.bincount(inv[ok], weights=v[ok], minlength=len(uniq))
        c = np.bincount(inv[ok], minlength=len(uniq))
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = (s / c).astype(np.float32)
        M.flat[uniq] = mean
    else:
        M[ri, ci] = v.astype(np.float32)
    t_pivot = time.time() - t0 - t_read
    wide = pd.DataFrame(M, index=pd.Index(sites, name="peptide_UID"), columns=list(conds))
    wide.reset_index().to_parquet(os.path.join(a.out, "compendium_wide.parquet"), index=False)
    t_total = time.time() - t0
    info = {
        "n_long_rows": n_long, "n_nan_log2FC_in_long": nan_vals,
        "n_sites": len(sites), "n_conditions": len(conds),
        "sites_not_in_02": len(extra_sites), "conditions_not_in_01": len(extra_conds),
        "n_duplicate_site_condition_pairs": n_dup_keys,
        "n_nonmissing_wide": int(np.isfinite(M).sum()),
        "density": float(np.isfinite(M).mean()),
        "seconds_read_csv": round(t_read, 1), "seconds_pivot": round(t_pivot, 1),
        "seconds_total": round(t_total, 1), "threads": a.threads,
    }
    json.dump(info, open(os.path.join(a.out, "compendium_wide.info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
