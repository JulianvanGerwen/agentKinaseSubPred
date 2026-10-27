#!/bin/bash
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
[ -n "$WAIT_PID" ] && while kill -0 $WAIT_PID 2>/dev/null; do sleep 20; done
CUDA_VISIBLE_DEVICES="" PYTHONPATH="$REPO/src" OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 "$PY" "$REPO"/scripts/models/r2_xfit.py --family knn --folds 0 1 2 3 4
