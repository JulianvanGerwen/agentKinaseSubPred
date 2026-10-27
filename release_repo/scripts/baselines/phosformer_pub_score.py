"""B3-pub: released fine-tuned Phosformer weights, inference only (as-published arm).

Steps (Phase A):
 1. tokenizer check: our encoder vs transformers.RobertaTokenizer loaded from the released vocab;
 2. sanity: reproduce the paper's test-set AUROC on their released test_pos_neg.csv (paper Fig.4A: 0.885);
 3. build kinase-domain inputs for our 503 kinases (cache/index/kinases.parquet + UniProt);
 4. throughput benchmark (fp32 vs bf16, batch sizes) and score a 1k-site sample x all residue-matched kinases.
Outputs: results/baselines/phosformer/pub_phaseA.json, cache/phosformer/pub/sample1k_scores.parquet,
         cache/phosformer/kinase_domains.parquet.
CAVEAT: the released model was trained on PSP/Phospho.ELM/PhosphoNetworks positives (13 organisms) that overlap our
test labels -> optimistic, contaminated reference, never an equal-footing comparator.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from ksp.baselines.phosformer.inputs import kinase_domains, window15_to_11mer
from ksp.baselines.phosformer.model import PhosformerClassifier
from ksp.baselines.phosformer.score import score_pairs
from ksp.baselines.phosformer.tokens import collate_pairs, encode_pair
from ksp.index import load_sites

REL = ROOT / "external/phosformer/hf_phosformer"
SPLITS = ROOT / "external/phosformer/hf_dataset_curated/splits"
OUT = ROOT / "results/baselines/phosformer"
CACHE = ROOT / "cache/phosformer"

ap = argparse.ArgumentParser()
ap.add_argument("--n-sites", type=int, default=1000)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
OUT.mkdir(parents=True, exist_ok=True)
(CACHE / "pub").mkdir(parents=True, exist_ok=True)
torch.set_num_threads(4)
dev = "cuda"
res = {"fallbacks": []}

# 1. tokenizer check ------------------------------------------------------------------------------------------
from transformers import RobertaTokenizer  # library tokenizer with the released vocab files

tk = RobertaTokenizer.from_pretrained(REL)
test = pd.read_csv(SPLITS / "test_pos_neg.csv")
ex = test.sample(50, random_state=0)
n_bad = 0
for k, p in zip(ex.kinase, ex.substrate):
    ref = tk(p.replace("-", "<pad>"), k, truncation=False)["input_ids"]
    n_bad += int(list(encode_pair(p, k)) != ref)
pads = test[test.substrate.str.contains("-")].head(20)
for k, p in zip(pads.kinase, pads.substrate):
    ref = tk(p.replace("-", "<pad>"), k, truncation=False)["input_ids"]
    n_bad += int(list(encode_pair(p, k)) != ref)
res["tokenizer_mismatches"] = n_bad
print("tokenizer mismatches vs HF RobertaTokenizer:", n_bad, flush=True)
if n_bad:
    raise SystemExit("tokenizer mismatch")

# 2. sanity on their test split --------------------------------------------------------------------------------
model = PhosformerClassifier.from_released(REL).to(dev).eval()
t = test.rename(columns={"uniprot": "kinase_acc", "substrate": "peptide"})
kseq = {f"r{i}": s for i, s in enumerate(t.kinase.unique())}
inv = {s: k for k, s in kseq.items()}
t["kinase_acc"] = t.kinase.map(inv)
for amp in (False, True):
    pr = score_pairs(model, kseq, t, batch_size=256, amp=amp)
    res[f"their_test_pos_neg_{'bf16' if amp else 'fp32'}"] = {
        "n": len(t), "auroc": float(roc_auc_score(t.label, pr)), "auprc": float(average_precision_score(t.label, pr))}
    if not amp:
        pr32 = pr
    else:
        res["bf16_vs_fp32_max_abs_diff"] = float(np.abs(pr - pr32).max())
        res["bf16_vs_fp32_spearman"] = float(pd.Series(pr).corr(pd.Series(pr32), method="spearman"))
print(json.dumps({k: v for k, v in res.items() if k.startswith("their") or k.startswith("bf16")}, indent=1), flush=True)
neg = pd.read_csv(SPLITS / "test_neg_only.csv").rename(columns={"substrate": "peptide"})
kseq2 = {f"n{i}": s for i, s in enumerate(neg.kinase.unique())}
inv2 = {s: k for k, s in kseq2.items()}
neg["kinase_acc"] = neg.kinase.map(inv2)
prn = score_pairs(model, kseq2, neg, batch_size=256, amp=True)
res["their_test_neg_only_FPR@0.5"] = float((prn >= 0.5).mean())

# 3. our kinase domains -----------------------------------------------------------------------------------------
kd, fb = kinase_domains()
res["fallbacks"] += fb
kd.to_parquet(CACHE / "kinase_domains.parquet", index=False)
res["kinase_domains"] = {"n": len(kd), "n_fallback_full_seq": len(fb),
                         "len_median": float(kd.domain_seq.str.len().median()),
                         "len_max": int(kd.domain_seq.str.len().max())}
k_index = pd.read_parquet(ROOT / "cache/index/kinases.parquet")[["uniprot_acc", "residue_class"]]
kd = kd.merge(k_index, on="uniprot_acc")
dom = dict(zip(kd.uniprot_acc, kd.domain_seq))

# 4. benchmark + 1k sample ----------------------------------------------------------------------------------------
sites = load_sites()
smp = sites.sample(args.n_sites, random_state=args.seed)
smp = smp.assign(peptide=smp.window15.map(window15_to_11mer))


def matching(res_, cls):
    return (cls == "STY") or (res_ in "ST" and cls == "ST") or (res_ == "Y" and cls == "Y")


rows = [(k, c, s, r, p) for k, c in zip(kd.uniprot_acc, kd.residue_class)
        for s, r, p in zip(smp.site_id, smp.res, smp.peptide) if matching(r, c)]
pairs = pd.DataFrame(rows, columns=["kinase_acc", "kinase_class", "site_id", "res", "peptide"])
res["sample"] = {"n_sites": len(smp), "n_pairs_residue_matched": len(pairs)}
bench = {}
sub = pairs.iloc[:40000]
for amp in (False, True):
    for bs in (256, 1024):
        torch.cuda.synchronize()
        t0 = time.time()
        score_pairs(model, dom, sub, batch_size=bs, amp=amp)
        torch.cuda.synchronize()
        bench[f"{'bf16' if amp else 'fp32'}_bs{bs}"] = len(sub) / (time.time() - t0)
        print("bench", amp, bs, bench, flush=True)
res["throughput_pairs_per_s"] = bench
best = max(bench.values())
t0 = time.time()
pairs["score"] = score_pairs(model, dom, pairs, batch_size=1024, amp=True)
res["sample_scoring_s"] = time.time() - t0
pairs[["kinase_acc", "site_id", "score"]].to_parquet(CACHE / "pub/sample1k_scores.parquet", index=False)
# full-universe extrapolation (residue-class matched)
n_st = int(sites.res.isin(["S", "T"]).sum())
n_y = int((sites.res == "Y").sum())
cls = k_index.residue_class.value_counts().to_dict()
n_full = (cls.get("ST", 0) + cls.get("STY", 0)) * n_st + (cls.get("Y", 0) + cls.get("STY", 0)) * n_y
res["full_universe"] = {"n_pairs_residue_matched": n_full, "gpu_hours_at_best": n_full / best / 3600}
# Sanity on our sample: motif checks for a few canonical kinases (score distribution, top peptides)
motif = {}
for gene in ["PRKACA", "CDK1", "CSNK2A1", "ATM", "SRC", "MAPK1"]:
    acc = kd.loc[kd.gene == gene, "uniprot_acc"]
    if len(acc) == 0:
        continue
    q = pairs[pairs.kinase_acc == acc.iloc[0]].sort_values("score", ascending=False)
    motif[gene] = {"n": len(q), "mean": float(q.score.mean()), "top5": q.peptide.head(5).tolist()}
res["motif_peek"] = motif
(OUT / "pub_phaseA.json").write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
