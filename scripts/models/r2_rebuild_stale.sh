#!/bin/bash
# Rebuild xfit stores built from the intermediate canary table (audit 17:30), then start the final-seed queue.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
E="PYTHONPATH=$REPO/src OMP_NUM_THREADS=2 MKL_NUM_THREADS=2"
env PYTHONPATH="$REPO/src" OMP_NUM_THREADS=2 "$PY" "$REPO"/scripts/models/r2_xfit.py --overwrite --family f3nb f6a f6a_clean --folds 0 > logs/r2/xfit4_f0.log 2>&1 &
p1=$!
env PYTHONPATH="$REPO/src" OMP_NUM_THREADS=2 "$PY" "$REPO"/scripts/models/r2_xfit.py --overwrite --family f3nb f2a --folds 3 > logs/r2/xfit4_f3.log 2>&1 &
p2=$!
wait $p1 $p2
env PYTHONPATH="$REPO/src" OMP_NUM_THREADS=2 "$PY" "$REPO"/scripts/models/r2_xfit.py --overwrite --family f2a --folds 0 > logs/r2/xfit4_f0b.log 2>&1
echo "rebuild done $(date +%H:%M)"
"$REPO"/scripts/models/r2_cv5_queue.sh cache/r2/queues/final_a.q > logs/r2/cv5q_final_a.log 2>&1
