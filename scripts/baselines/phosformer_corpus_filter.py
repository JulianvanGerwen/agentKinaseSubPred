"""Remove TEST proteins (and homologs) from the Phosformer MLM corpus (PI AIFM2-8), then tokenise it.

Inputs: external/phosformer/hf_dataset_curated/corpus/{kinases,substrates}.fa (their released corpus),
        results/splits/test_proteins.tsv + substrate_clusters.tsv (frozen split), canonical UniProt sequences.
TEST set T = TEST proteins ∪ every protein in a substrate cluster that contains a TEST protein.

Rules (stated, not tuned):
- substrates.fa (5,338 full-length proteins, 29 organisms): drop if the UniProt accession (header field 2, isoform
  suffix stripped) is in T, or if mmseqs2 finds a hit to any T protein with identity >= 0.40 and query AND target
  coverage >= 0.50 (= the substrate-cluster edge rule of ksp.splits, so non-human orthologs of TEST proteins go too).
- kinases.fa (295,320 NCBI-nr kinase DOMAINS): drop if it matches a T protein (TEST substrate proteins that are
  themselves kinases) with identity >= 0.90 over >= 0.80 of the domain (the protein itself and near-identical
  orthologs/alleles). The 40%/50% cluster rule is NOT applied to kinase domains: all eukaryotic protein kinase
  domains are ~25-50% identical to each other, so it would delete most of the corpus and change the method; counts
  at 0.40 and 0.50 are reported for transparency.
Outputs: cache/phosformer/corpus_filtered/{kinases.fa,substrates.fa,removed.tsv,tokens.npy,offsets.npy,source.npy,
         meta.json}; results/baselines/phosformer/corpus_filter.json.
Usage: python -I scripts/baselines/phosformer_corpus_filter.py [--test-proteins FILE] [--out DIR] [--threads 6]
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd

from ksp.baselines.phosformer.corpus import CORPUS_SRC, build_corpus, read_fasta
from ksp.baselines.phosformer.inputs import canonical_seqs

MMSEQS = "/home/tools/mmseqs/bin/mmseqs"
ap = argparse.ArgumentParser()
ap.add_argument("--test-proteins", default=str(ROOT / "results/splits/test_proteins.tsv"))
ap.add_argument("--clusters", default=str(ROOT / "results/splits/substrate_clusters.tsv"))
ap.add_argument("--out", default=str(ROOT / "cache/phosformer/corpus_filtered"))
ap.add_argument("--report", default=str(ROOT / "results/baselines/phosformer/corpus_filter.json"))
ap.add_argument("--threads", type=int, default=6)
args = ap.parse_args()
out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
tmp = out / "tmp"
tmp.mkdir(exist_ok=True)

tp = pd.read_csv(args.test_proteins, sep="\t", comment="#")
T = set(tp.acc)
n_test_direct = len(T)
if Path(args.clusters).exists():
    cl = pd.read_csv(args.clusters, sep="\t", comment="#")
    test_cl = set(cl.loc[cl.acc.isin(T), "cluster_id"])
    T |= set(cl.loc[cl.cluster_id.isin(test_cl), "acc"])
seqs = canonical_seqs(T)
missing = sorted(T - set(seqs))
if missing:
    raise KeyError(f"{len(missing)} TEST proteins without canonical sequence: {missing[:5]}")
with open(tmp / "test.fa", "w") as f:
    for a in sorted(T):
        f.write(f">{a}\n{seqs[a]}\n")

kin = read_fasta(CORPUS_SRC / "kinases.fa")
sub = read_fasta(CORPUS_SRC / "substrates.fa")
for name, recs, pre in [("kin", kin, "k"), ("sub", sub, "s")]:
    with open(tmp / f"{name}.fa", "w") as f:
        for i, (_, s) in enumerate(recs):
            f.write(f">{pre}{i}\n{s}\n")


def search(q: str, min_id: float, cov: float, cov_mode: int, tag: str) -> pd.DataFrame:
    m8 = tmp / f"{tag}.m8"
    wd = tmp / f"wd_{tag}"
    if wd.exists():
        shutil.rmtree(wd)
    cmd = [MMSEQS, "easy-search", str(tmp / f"{q}.fa"), str(tmp / "test.fa"), str(m8), str(wd),
           "--threads", str(args.threads), "-s", "7.5", "--max-seqs", "300", "-e", "1e-3",
           "--min-seq-id", str(min_id), "-c", str(cov), "--cov-mode", str(cov_mode),
           "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return pd.read_csv(m8, sep="\t", header=None,
                       names=["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"])


removed = []
# substrates: accession match + 40%/50% (both coverages, cov-mode 0)
sub_acc = []
for h, _ in sub:
    parts = h.split("|")
    sub_acc.append(parts[1].split("-")[0] if len(parts) > 2 else h.split()[0])
for i, a in enumerate(sub_acc):
    if a in T:
        removed.append(("substrates.fa", i, sub[i][0][:80], "accession_in_TEST", a, 1.0))
hs = search("sub", 0.40, 0.50, 0, "sub")
hs = hs[(hs.fident >= 0.40) & (hs.qcov >= 0.50) & (hs.tcov >= 0.50)]
done = {r[1] for r in removed}
for q, g in hs.groupby("query"):
    i = int(q[1:])
    if i not in done:
        b = g.sort_values("fident", ascending=False).iloc[0]
        removed.append(("substrates.fa", i, sub[i][0][:80], "homolog_40id_50cov", b.target, float(b.fident)))
# kinase domains: query-coverage mode (cov-mode 2 = coverage of query)
hk = search("kin", 0.40, 0.50, 2, "kin")
best = hk.sort_values("fident", ascending=False).drop_duplicates("query")
info = {"kin_domains_with_hit_id>=0.40_qcov>=0.50": int(len(best)),
        "kin_domains_with_hit_id>=0.50_qcov>=0.50": int((best.fident >= 0.50).sum()),
        "kin_domains_with_hit_id>=0.90_qcov>=0.80": int(((best.fident >= 0.90) & (best.qcov >= 0.80)).sum())}
kb = best[(best.fident >= 0.90) & (best.qcov >= 0.80)]
for b in kb.itertuples():
    i = int(b.query[1:])
    removed.append(("kinases.fa", i, kin[i][0][:80], "near_identical_90id_80qcov", b.target, float(b.fident)))

rm = pd.DataFrame(removed, columns=["file", "index", "header", "reason", "test_acc", "fident"])
rm.to_csv(out / "removed.tsv", sep="\t", index=False)
drop = {(f, i) for f, i in zip(rm.file, rm["index"])}
for fname, recs in [("kinases.fa", kin), ("substrates.fa", sub)]:
    with open(out / fname, "w") as f:
        for i, (h, s) in enumerate(recs):
            if (fname, i) not in drop:
                f.write(f">{h}\n{s}\n")
meta = build_corpus(out, out)
rep = {"n_test_proteins_direct": n_test_direct, "n_T_with_cluster_members": len(T),
       "removed": rm.groupby(["file", "reason"]).size().to_dict(),
       "removed_residues": {fn: int(sum(len(r[1]) for j, r in enumerate(recs) if (fn, j) in drop))
                            for fn, recs in [("kinases.fa", kin), ("substrates.fa", sub)]},
       "kinase_domain_hit_counts": info, "filtered_corpus": meta["counts"],
       "original": {"kinases.fa": len(kin), "substrates.fa": len(sub)}}
rep["removed"] = {f"{k[0]}|{k[1]}": int(v) for k, v in rep["removed"].items()}
Path(args.report).parent.mkdir(parents=True, exist_ok=True)
Path(args.report).write_text(json.dumps(rep, indent=1))
shutil.rmtree(tmp)
print(json.dumps(rep, indent=1))
