#!/bin/bash
# Unattended, resumable B3a MLM pretraining on GPU 0 on the TEST-filtered corpus (cache/phosformer/corpus_filtered,
# built by phosformer_corpus_filter.py). Retries up to 5x; each retry resumes from ckpt_last.pt.
cd /workspace
export CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=4 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
for i in 1 2 3 4 5; do
  .venv/bin/python -I scripts/baselines/phosformer_pretrain.py --run-dir cache/phosformer/pretrain \
     --corpus-dir cache/phosformer/corpus_filtered --epochs 50 --compile --token-budget 49152 --threads 4 && exit 0
  echo "[wrapper] attempt $i failed, resuming in 30s" >&2; sleep 30
done
exit 1
