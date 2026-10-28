#!/usr/bin/env python
"""Download AlphaFold DB predicted-aligned-error (PAE) JSON files for a list of UniProt accessions (F1 models).

URL pattern (from the AFDB API field `paeDocUrl`): https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_v{ver}.json
(the server returns 403 without a browser-like User-Agent). Files are saved gzip-compressed to
external/alphafold/pae/AF-{acc}-F1-predicted_aligned_error_v{ver}.json.gz; a manifest TSV records acc, http status,
bytes and SHA-256 of the uncompressed JSON.

Usage: python -I scripts/download_afdb_pae.py --accs accs.txt [--threads 2] [--ver 6]
(AFDB returns HTTP 429 at 8 parallel connections; 2 threads with 429 back-off worked.)
"""
import argparse
import gzip
import hashlib
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
UA = {"User-Agent": "Mozilla/5.0 (ksp research pipeline)"}


def fetch(acc, out, ver):
    f = out / f"AF-{acc}-F1-predicted_aligned_error_v{ver}.json.gz"
    if f.exists():
        return acc, "cached", f.stat().st_size, ""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_v{ver}.json"
    for attempt in range(5):
        try:
            r = requests.get(url, headers=UA, timeout=120)
            if r.status_code == 200:
                tmp = f.with_suffix(".tmp")
                with gzip.open(tmp, "wb", compresslevel=3) as g:
                    g.write(r.content)
                tmp.rename(f)
                return acc, 200, len(r.content), hashlib.sha256(r.content).hexdigest()
            if r.status_code == 404:
                return acc, 404, 0, ""
            if r.status_code == 429:  # rate limited: back off politely
                time.sleep(60 * (attempt + 1))
                continue
        except requests.RequestException:
            pass
        time.sleep(2 * (attempt + 1))
    return acc, "fail", 0, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accs", required=True)
    ap.add_argument("--out", default=str(ROOT / "external/alphafold/pae"))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--ver", type=int, default=6)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    accs = Path(a.accs).read_text().split()
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex, open(out / "MANIFEST.tsv", "a") as man:
        for i, (acc, st, nb, sha) in enumerate(ex.map(lambda x: fetch(x, out, a.ver), accs)):
            man.write(f"{acc}\t{st}\t{nb}\t{sha}\t{time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
            if i % 100 == 0:
                print(f"{i}/{len(accs)} {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(accs)} in {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
