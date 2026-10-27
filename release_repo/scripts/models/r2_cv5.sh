#!/bin/bash
# 5-fold CV of one R2 config: per-fold runs (not registered) + pooled evaluation (registered).
# Usage: r2_cv5.sh <name> <threads> <stage> "<set args>" [fold0_run_id_to_reuse]
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
NAME=$1; T=$2; STAGE=$3; SETS=$4; REUSE0=$5
L=logs/r2/cv5_${NAME}
mkdir -p $L
RUNS=()
for f in 0 1 2 3 4; do
  if [ "$f" == "0" ] && [ -n "$REUSE0" ]; then RUNS+=("$REUSE0"); continue; fi
  PYTHONPATH="$REPO/src" R2_THREADS=$T OMP_NUM_THREADS=$T "$PY" "$REPO"/scripts/models/r2_run.py --name "${NAME}_f$f" --folds $f \
     --stage "$STAGE" --no-register --set $SETS > $L/fold$f.log 2>&1
  rid=$(grep RESULT $L/fold$f.log | python3 -c "import sys,json; print(json.loads(sys.stdin.read().split('RESULT ',1)[1])['run_id'])")
  [ -z "$rid" ] && { echo "fold $f failed"; tail -5 $L/fold$f.log; exit 1; }
  RUNS+=("$rid")
  echo "$(date +%H:%M:%S) $NAME fold $f -> $rid: $(grep 'M1/M5' $L/fold$f.log | cut -c1-120)"
done
PYTHONPATH="$REPO/src" OMP_NUM_THREADS=4 "$PY" "$REPO"/scripts/models/r2_pool.py --name "${NAME}_cv5" --runs "${RUNS[@]}" --stage "$STAGE" > $L/pool.log 2>&1
echo "$(date +%H:%M:%S) $NAME pooled: $(grep -E 'M1/M5|Error' $L/pool.log | cut -c1-600)"
