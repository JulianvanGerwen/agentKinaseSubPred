"""B3-pub on DEV CV (pooled): released fine-tuned Phosformer, zero-shot (contaminated reference).

Scores set A = (all DEV-universe positive sites U E-hard sites) x every residue-matched kinase, which makes M1/M2/M4/M5
complete; E-all metrics (M3, AUROC_all) are NOT complete (missing='drop') unless --eall-frac is used, which adds a
uniform random fraction of the remaining E-all sites (reported as an estimate only).
CAVEAT: the released model was trained on PSP/Phospho.ELM/PhosphoNetworks (13 organisms) -> its training positives
overlap our DEV and TEST labels; never compared on equal footing.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
import torch

from ksp import eval as E
from ksp import index
from ksp.baselines.phosformer.inputs import kinase_domains, window15_to_11mer
from ksp.baselines.phosformer.model import PhosformerClassifier
from ksp.baselines.phosformer.score import score_pairs
from ksp.manifest import RunManifest

ap = argparse.ArgumentParser()
ap.add_argument("--eall-frac", type=float, default=0.0)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
torch.set_num_threads(4)
REL = ROOT / "external/phosformer/hf_phosformer"

K = index.load_kinases()
S = index.load_sites()
kd, kfb = kinase_domains(K)
dom = dict(zip(kd.uniprot_acc, kd.domain_seq))
pep = np.array([window15_to_11mer(w) for w in S.window15.values], dtype=object)
res = S.res.values
data = E.prepare("cv", fold=None, with_canaries=False)
uni = data.universe
hard_u = data.masks.hard_idx[uni[data.masks.hard_idx]]
pos_u = np.unique(np.concatenate(list(data.positives.values())))
setA = np.zeros(len(S), bool)
setA[np.union1d(hard_u, pos_u)] = True


def cm(c):
    return np.ones(len(S), bool) if c == "STY" else (np.isin(res, ["S", "T"]) if c == "ST" else res == "Y")


rows = [pd.DataFrame({"kinase_acc": a, "site_idx": np.where(setA & cm(c))[0]})
        for a, c in zip(K.uniprot_acc, K.residue_class)]
if args.eall_frac > 0:
    r = np.random.default_rng(args.seed)
    extra = uni & ~setA & (r.random(len(S)) < args.eall_frac)
    kc = dict(zip(K.uniprot_acc, K.residue_class))
    rows += [pd.DataFrame({"kinase_acc": a, "site_idx": np.where(extra & cm(kc[a]))[0]}) for a in data.positives]
pairs = pd.concat(rows, ignore_index=True)
pairs["peptide"] = pep[pairs.site_idx.values]
print(f"scoring {len(pairs)} pairs ({int(setA.sum())} set-A sites)", flush=True)
model = PhosformerClassifier.from_released(REL).cuda().eval()
model.roberta = torch.compile(model.roberta, dynamic=True)
t0 = time.time()
pairs["score"] = score_pairs(model, dom, pairs, batch_size=1024, log_every=500)
el = time.time() - t0
preds = pd.DataFrame({"kinase_acc": pairs.kinase_acc.values, "site_id": S.site_id.values[pairs.site_idx.values],
                      "score": pairs.score.values.astype(np.float64)})
ev = E.evaluate(preds, data, missing="drop")
cfg = {"arm": "B3-pub", "weights": str(REL.relative_to(ROOT)), "set": "A=positives+E-hard x all kinases",
       "eall_frac": args.eall_frac, "n_pairs": len(pairs), "score_s": el, "pairs_per_s": len(pairs) / el,
       "contaminated": True}
with RunManifest("b3pub_phosformer_cv_pooled", config=cfg, seed=args.seed,
                 inputs=[REL / "pytorch_model.bin", ROOT / "results/splits/cv_folds.tsv",
                         ROOT / "cache/index/sites.parquet", ROOT / "cache/index/kinases.parquet"],
                 split="cv_pooled", notes="released fine-tuned weights, zero-shot; CONTAMINATED (trained on PSP)") as run:
    run.count_fallback("kinase_domain_full_seq", len(kfb))
    run.save_predictions(preds)
    run.save_metrics(ev)
    run.set("endpoints", ev.endpoints)
print(json.dumps({"run_id": run.run_id, "endpoints": ev.endpoints}, default=str, indent=1))
