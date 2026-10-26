"""ESM-2 kinase descriptors (DESIGN §7 F7): kinase-domain and full-protein mean-pooled embeddings.

Kinase set: results/labels/kinase_list.tsv (all candidate rows with a UniProt accession, `include` flag carried;
domains from its `kinase_domains` column, `is_primary` = its kinase_domain_start/end) if it exists; otherwise all reviewed human UniProt entries with keyword KW-0418 that carry a
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


def parse_kl_domains(s) -> list[tuple[str, int, int]]:
    """kinase_list.tsv `kinase_domains`: '545-809:Protein kinase 1;849-1124:Protein kinase 2'."""
    out = []
    if isinstance(s, str) and s:
        for part in s.split(";"):
            m = re.match(r"\s*(\d+)-(\d+):(.*)$", part)
            if m:
                out.append((m.group(3).strip(), int(m.group(1)), int(m.group(2))))
            elif out:  # a ';' inside a domain name, e.g. 'Protein kinase; truncated'
                n, a, b = out[-1]
                out[-1] = (n + ";" + part, a, b)
    return out


def kinase_set() -> tuple[pd.DataFrame, str]:
    """Return DataFrame acc, gene, include, domains [(name, start, end)], primary (start, end) or None."""
    kl = ROOT / "results/labels/kinase_list.tsv"
    if kl.exists():
        k = pd.read_csv(kl, sep="\t", dtype=str, keep_default_na=False)
        k = k[k.uniprot_acc != ""]
        k = pd.DataFrame({"acc": k.uniprot_acc.str.split("-").str[0], "gene": k.gene, "include": k.include == "True",
                          "domains": k.kinase_domains.map(parse_kl_domains),
                          "primary": [(int(a), int(b)) if a and b else None
                                      for a, b in zip(k.kinase_domain_start, k.kinase_domain_end)]})
        src = str(kl.relative_to(ROOT))
    else:
        ann = read_annotations()
        ann = ann[ann.Reviewed == "reviewed"].set_index("Entry")["Domain [FT]"]
        kw = pd.read_csv(ROOT / "external/uniprot/kinase_kw0418/human_kw0418_reviewed.tsv.gz", sep="\t", dtype=str,
                         keep_default_na=False)
        kw = kw[kw["Domain [FT]"].map(lambda s: len(kinase_domains(s)) > 0)]
        doms = [kinase_domains(ann.get(a, "")) for a in kw.Entry]
        k = pd.DataFrame({"acc": kw.Entry.values, "gene": kw["Gene Names (primary)"].values, "include": True,
                          "domains": doms, "primary": [(d[0][1], d[0][2]) if d else None for d in doms]})
        src = "UniProt KW-0418 with kinase-catalytic DOMAIN (fallback)"
    return k.drop_duplicates("acc").reset_index(drop=True), src


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
    for (acc, gene, inc, doms, prim), E in zip(k[["acc", "gene", "include", "domains", "primary"]].itertuples(
            index=False), full):
        L = len(seqs[acc])
        doms = [(n, s, e) for n, s, e in doms if 1 <= s <= e <= L]
        base = dict(acc=acc, gene=gene, include=inc, length=L, n_kinase_domains=len(doms), kinase_source=src,
                    emb_full=E.mean(0))
        if not doms:
            rows.append({**base, "domain": None, "dom_start": np.nan, "dom_end": np.nan, "is_primary": False,
                         "emb_domain_ctx": None})
            iso_seqs.append(None)
        for i, (n, s, e) in enumerate(doms):
            rows.append({**base, "domain_idx": i, "domain": n, "dom_start": s, "dom_end": e,
                         "is_primary": prim is not None and (s, e) == tuple(prim),
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
