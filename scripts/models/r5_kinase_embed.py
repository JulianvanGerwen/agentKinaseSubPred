"""R5: cache per-token ESM-2 650M embeddings of the 503 human kinase catalytic domains.

The R5 kinase tower does late interaction over kinase-domain *tokens*, so it needs per-residue (not pooled)
embeddings. They are computed once and cached; the tower on top of them is trained.

Domain = primary kinase domain of results/labels/kinase_list.tsv (cache/index/kinases.parquet columns
kinase_domain_start/end, 1-based inclusive) on the canonical Swiss-Prot sequence. The domain sequence is embedded
ON ITS OWN (a catalytic domain is a folding unit); this is the `emb_domain_iso` convention of
src/ksp/features/esm_kinases.py. Domains longer than L_KIN keep the N-terminal and C-terminal halves of the budget
(kinase inserts sit in the middle: SRPK1/2 spacer, MASTL, CDC7). Kinases with no annotated domain (5 of 503) fall
back to the first L_KIN residues of the protein and are flagged.

Outputs in cache/r5/:
  kin_tok.npy        float16 (503, L_KIN, 1280) per-token last-layer embeddings, 0 where padded
  kin_mask.npy       bool    (503, L_KIN)
  kin_meta.json      model id, L_KIN, per-kinase length / truncation / fallback flags

  CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_kinase_embed.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import torch  # noqa: E402

from ksp import index, uniprot  # noqa: E402

L_KIN = 320
MODEL = "facebook/esm2_t33_650M_UR50D"
OUT = ROOT / "cache" / "r5"


def domain_sequences() -> tuple[list[str], list[dict]]:
    K = index.load_kinases()
    fa = uniprot.read_fasta(ROOT / "external/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz")
    seqs = dict(zip(fa[fa.isoform == ""].acc, fa[fa.isoform == ""].seq))
    out, meta = [], []
    for acc, a, b in zip(K.uniprot_acc, K.kinase_domain_start, K.kinase_domain_end):
        s = seqs.get(acc)
        if s is None:
            raise KeyError(f"no canonical sequence for kinase {acc}")
        fallback = not (np.isfinite(a) and np.isfinite(b))
        d = s[: L_KIN] if fallback else s[int(a) - 1: int(b)]
        full = len(d)
        trunc = full > L_KIN
        if trunc:                       # keep both lobes, drop the middle insert
            h = L_KIN // 2
            d = d[:h] + d[full - (L_KIN - h):]
        out.append(d)
        meta.append({"acc": acc, "domain_len": full, "truncated": bool(trunc), "no_domain_annotation": bool(fallback),
                     "used_len": len(d)})
    return out, meta


def main():
    from transformers import AutoTokenizer, EsmModel
    OUT.mkdir(parents=True, exist_ok=True)
    seqs, meta = domain_sequences()
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = EsmModel.from_pretrained(MODEL, torch_dtype=torch.float16).cuda().eval()
    E = np.zeros((len(seqs), L_KIN, model.config.hidden_size), np.float16)
    M = np.zeros((len(seqs), L_KIN), bool)
    t0 = time.time()
    B = 16
    for a in range(0, len(seqs), B):
        batch = seqs[a:a + B]
        enc = tok(batch, return_tensors="pt", padding=True, add_special_tokens=True)
        enc = {k: v.cuda() for k, v in enc.items()}
        with torch.no_grad():
            h = model(**enc).last_hidden_state                       # (b, 2+L, D), BOS/EOS included
        for i, s in enumerate(batch):
            n = len(s)
            E[a + i, :n] = h[i, 1:1 + n].float().cpu().numpy().astype(np.float16)
            M[a + i, :n] = True
    np.save(OUT / "kin_tok.npy", E)
    np.save(OUT / "kin_mask.npy", M)
    (OUT / "kin_meta.json").write_text(json.dumps(
        {"model": MODEL, "L_KIN": L_KIN, "n_kinases": len(seqs), "d": int(E.shape[2]),
         "n_truncated": int(sum(m["truncated"] for m in meta)),
         "n_no_domain": int(sum(m["no_domain_annotation"] for m in meta)),
         "seconds": round(time.time() - t0, 1), "kinases": meta}, indent=1))
    print(f"kin_tok {E.shape} {E.nbytes/2**30:.2f} GiB in {time.time()-t0:.0f}s; "
          f"truncated {sum(m['truncated'] for m in meta)}, no-domain {sum(m['no_domain_annotation'] for m in meta)}")


if __name__ == "__main__":
    main()
