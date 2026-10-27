#!/usr/bin/env python
"""Download UniProt human reference proteome + kinase lists into external/uniprot/.

Files (each in its own subdirectory, as downloaded files are untrusted):
  proteome/UP000005640_all_isoforms.fasta.gz   all entries of reference proteome UP000005640
                                               (Swiss-Prot + TrEMBL), canonical + isoforms
  proteome/UP000005640_reviewed_isoforms.fasta.gz  Swiss-Prot only, canonical + isoforms
  human_tsv/UP000005640_annotations.tsv.gz     per-entry table (gene names, keywords, domains, families)
  kinase_kw0418/human_kw0418_reviewed.tsv.gz   reviewed human entries with keyword KW-0418 (Kinase)
  pkinfam/pkinfam.txt                          UniProt "Protein kinase families" doc (all species)
  ptm_tsv/UP000005640_reviewed_ptm.tsv.gz      reviewed human entries: gene synonyms, HGNC xref, MOD_RES features
                                               (incl. "Phosphoserine; by PKA" notes + ECO evidence codes)
  kinase_detail/human_kw0418_detail.tsv.gz     reviewed KW-0418 entries: domains, active site, binding, function,
                                               caution text (used for kinase list / pseudokinase flags)

A manifest (external/uniprot/MANIFEST.json) records URL, UniProt release, release date,
download timestamp and SHA-256 for every file.

Usage: python -I scripts/download_uniprot.py [--out external/uniprot] [--force]
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import sys

import requests

REST = "https://rest.uniprot.org/uniprotkb/stream"
TSV_FIELDS = ("accession,id,reviewed,gene_primary,gene_names,protein_name,organism_id,length,"
              "sequence_version,date_sequence_modified,keywordid,keyword,ft_domain,protein_families,ec")

FILES = [
    ("proteome/UP000005640_all_isoforms.fasta.gz", REST,
     {"query": "proteome:UP000005640", "format": "fasta", "includeIsoform": "true", "compressed": "true"}),
    ("proteome/UP000005640_reviewed_isoforms.fasta.gz", REST,
     {"query": "proteome:UP000005640 AND reviewed:true", "format": "fasta", "includeIsoform": "true",
      "compressed": "true"}),
    ("human_tsv/UP000005640_annotations.tsv.gz", REST,
     {"query": "proteome:UP000005640", "format": "tsv", "fields": TSV_FIELDS, "compressed": "true"}),
    ("kinase_kw0418/human_kw0418_reviewed.tsv.gz", REST,
     {"query": "organism_id:9606 AND reviewed:true AND keyword:KW-0418", "format": "tsv",
      "fields": TSV_FIELDS, "compressed": "true"}),
    ("ptm_tsv/UP000005640_reviewed_ptm.tsv.gz", REST,
     {"query": "proteome:UP000005640 AND reviewed:true", "format": "tsv",
      "fields": "accession,id,gene_primary,gene_synonym,xref_hgnc,sequence_version,length,ft_mod_res",
      "compressed": "true"}),
    ("kinase_detail/human_kw0418_detail.tsv.gz", REST,
     {"query": "organism_id:9606 AND reviewed:true AND keyword:KW-0418", "format": "tsv",
      "fields": "accession,id,gene_primary,gene_synonym,protein_name,length,ft_domain,ft_act_site,ft_binding,"
                "cc_function,cc_caution,ec,protein_families,keyword",
      "compressed": "true"}),
    ("pkinfam/pkinfam.txt", "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/pkinfam.txt",
     None),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="external/uniprot")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    man_path = os.path.join(a.out, "MANIFEST.json")
    manifest = json.load(open(man_path)) if os.path.exists(man_path) else {}
    for rel, url, params in FILES:
        path = os.path.join(a.out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path) and rel in manifest and not a.force:
            print("skip (exists)", rel)
            continue
        print("GET", url, params, file=sys.stderr)
        with requests.get(url, params=params, stream=True, timeout=600) as r:
            r.raise_for_status()
            tmp = path + ".part"
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
            os.replace(tmp, path)
            manifest[rel] = {
                "url": r.url,
                "uniprot_release": r.headers.get("X-UniProt-Release"),
                "uniprot_release_date": r.headers.get("X-UniProt-Release-Date"),
                "last_modified": r.headers.get("Last-Modified"),
                "downloaded_utc": dt.datetime.utcnow().isoformat(timespec="seconds"),
                "bytes": os.path.getsize(path),
                "sha256": sha256(path),
            }
        print(rel, manifest[rel])
        json.dump(manifest, open(man_path, "w"), indent=1)


if __name__ == "__main__":
    main()
