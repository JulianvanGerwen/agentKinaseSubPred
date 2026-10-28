#!/bin/bash
# Run a queue of R3 train jobs on one GPU. usage: r3_queue.sh <gpu> "<args1>" "<args2>" ...
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
gpu=$1; shift
for job in "$@"; do
  tag=$(echo "$job" | tr ' ' '_' | tr -d '-')
  CUDA_VISIBLE_DEVICES=$gpu OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "$PY" -I "$REPO"/scripts/models/r3_cv.py train $job > logs/r3/q_${tag}.log 2>&1
done
