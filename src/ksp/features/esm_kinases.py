"""ESM-2 kinase descriptors (DESIGN §7 F7): kinase-domain and full-protein mean-pooled embeddings.

Kinase set: results/labels/kinase_list.tsv (column with UniProt accession: 'uniprot'/'acc'/'uniprot_acc'/
'UniProt') if it exists; otherwise all reviewed human UniProt entries with keyword KW-0418 that carry a
kinase-catalytic domain feature (fallback; stated in the output column `kinase_source`).
Kinase-catalytic domain = UniProt DOMAIN feature whose note matches KINASE_DOMAIN_RE ('Protein kinase',
'Protein kinase 1/2', 'PI3K/PI4K catalytic', 'Alpha-type protein kinase', 'Histidine kinase' (PDHK/BCKDK)).

One row per (kinase, kinase domain); kinases without such a domain get one row with domain fields NA.
Embedding columns (float32 lists, D=1280 for 650M):
  emb_full        mean over all residues of the full-protein embedding (windows for >1022 aa, overlaps averaged)
  emb_domain_ctx  mean over domain residues, taken from the full-protein embedding (domain in context)
  emb_domain_iso  mean over residues of the domain sequence embedded on its own
Usage: CUDA_VISIBLE_DEVICES=2 python -I src/ksp/features/esm_kinases.py --out cache/esm2_650m/kinases.parquet
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from ksp.features.esm_embed import canonical_swissprot, embed_sequences, load_model  # noqa: E402
from ksp.uniprot import read_annotations  # noqa: E402

KINASE_DOMAIN_RE = re.compile(r"^(Protein kinase( \d+)?|PI3K/PI4K catalytic|Alpha-type protein kinase|Histidine kinase)$")
DOM_RE = re.compile(r'DOMAIN (\d+)\.\.(\d+); /note="([^"]+)"')


def kinase_domains(ft: str) -> list[tuple[str, int, int]]:
    return [(n, int(a), int(b)) for a, b, n in DOM_RE.findall(ft or "") if KINASE_DOMAIN_RE.match(n)]


def kinase_set() -> tuple[pd.DataFrame, str]:
    kl = ROOT / "results/labels/kinase_list.tsv"
    ann = read_annotations()
    ann = ann[ann.Reviewed == "reviewed"]
    if kl.exists():
        k = pd.read_csv(kl, sep="\t", dtype=str)
        col = next(c for c in ["uniprot", "acc", "uniprot_acc", "UniProt", "kinase_acc", "uniprot_id"] if c in k.columns)
        gcol = next((c for c in ["gene", "hgnc_symbol", "symbol", "kinase"] if c in k.columns), None)
        k = pd.DataFrame({"acc": k[col].str.split("-").str[0], "gene": k[gcol] if gcol else ""}).drop_duplicates("acc")
        src = str(kl.relative_to(ROOT))
    else:
        kw = pd.read_csv(ROOT / "external/uniprot/kinase_kw0418/human_kw0418_reviewed.tsv.gz", sep="\t", dtype=str,
                         keep_default_na=False)
        kw = kw[kw["Domain [FT]"].map(lambda s: len(kinase_domains(s)) > 0)]
        k = pd.DataFrame({"acc": kw.Entry, "gene": kw["Gene Names (primary)"]})
        src = "UniProt KW-0418 with kinase-catalytic DOMAIN (fallback)"
    k = k.merge(ann[["Entry", "Domain [FT]"]].rename(columns={"Entry": "acc", "Domain [FT]": "ft"}), on="acc", how="left")
    return k, src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--out", default=str(ROOT / "cache/esm2_650m/kinases.parquet"))
    a = ap.parse_args()
    k, src = kinase_set()
    seqs = canonical_swissprot().set_index("acc").seq
    missing = sorted(set(k.acc) - set(seqs.index))
    if missing:
        print(f"WARNING {len(missing)} kinase accessions not in Swiss-Prot canonical: {missing[:10]}", flush=True)
    k = k[k.acc.isin(seqs.index)].reset_index(drop=True)
    tok, model = load_model(a.model, "fp16")
    t0 = time.time()
    full = embed_sequences([seqs[x] for x in k.acc], tok, model)
    rows, iso_seqs = [], []
    for (acc, gene, ft), E in zip(k[["acc", "gene", "ft"]].itertuples(index=False), full):
        doms = kinase_domains(ft if isinstance(ft, str) else "")
        base = dict(acc=acc, gene=gene, length=len(seqs[acc]), n_kinase_domains=len(doms), kinase_source=src,
                    emb_full=E.mean(0))
        if not doms:
            rows.append({**base, "domain": None, "dom_start": np.nan, "dom_end": np.nan, "emb_domain_ctx": None})
            iso_seqs.append(None)
        for i, (n, s, e) in enumerate(doms):
            rows.append({**base, "domain_idx": i, "domain": n, "dom_start": s, "dom_end": e,
                         "emb_domain_ctx": E[s - 1:e].mean(0)})
            iso_seqs.append(seqs[acc][s - 1:e])
    todo = [i for i, s in enumerate(iso_seqs) if s]
    iso = embed_sequences([iso_seqs[i] for i in todo], tok, model)
    for i, E in zip(todo, iso):
        rows[i]["emb_domain_iso"] = E.mean(0)
    df = pd.DataFrame(rows)
    for c in ["emb_full", "emb_domain_ctx", "emb_domain_iso"]:
        df[c] = df[c].map(lambda v: None if v is None or (isinstance(v, float)) else np.asarray(v, np.float32).tolist())
    df["model"] = a.model
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out, index=False)
    print(f"{len(df)} rows, {df.acc.nunique()} kinases ({(df.n_kinase_domains == 0).sum()} without kinase domain), "
          f"source={src}, {time.time() - t0:.0f}s -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
