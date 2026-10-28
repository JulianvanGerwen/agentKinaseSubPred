#!/bin/bash
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
until grep -q "clean library done" logs/r2/chain_lko_clean.log; do sleep 30; done
RID=$(grep "^RUN" logs/r2/final_clean_uobs.log | awk '{print $2}')
echo "clean run $RID"
OMP_NUM_THREADS=8 "$PY" "$REPO"/scripts/downstream/run_predictions.py --pred cache/runs/$RID/uobs_predictions.parquet --tag r2_clean > logs/r2/downstream_r2_clean.log 2>&1
echo "downstream done rc=$? $(date +%H:%M)"
