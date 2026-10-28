#!/bin/bash
# R3 base-hyper-parameter dev sweep on fold-0 TRAINING rows only (inner validation; eval fold untouched).
# usage: scripts/models/r3_dev_sweep.sh <gpu> <name> <json> [<name> <json> ...]
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
gpu=$1; shift
while [ $# -gt 0 ]; do
  CUDA_VISIBLE_DEVICES=$gpu OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "$PY" -I "$REPO"/scripts/models/r3_cv.py \
    train --fold 0 --config w7_nce05 --inner-only --set-name "$1" --set "$2" > logs/r3/dev_$1.log 2>&1
  shift 2
done
