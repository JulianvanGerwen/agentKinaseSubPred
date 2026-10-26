#!/bin/bash
# Run R3 final TEST components sequentially on one GPU. usage: r3_test_queue.sh <gpu> <final|J> [<final|J> ...]
cd /workspace
gpu=$1; shift
for j in "$@"; do
  if [ "$j" = final ]; then extra=""; else extra="--lko $j"; fi
  CUDA_VISIBLE_DEVICES=$gpu OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 /workspace/.venv/bin/python -I scripts/models/r3_cv.py \
    test --config w7_nce05_cpt_aux $extra > logs/r3/test_${j}.log 2>&1
done
