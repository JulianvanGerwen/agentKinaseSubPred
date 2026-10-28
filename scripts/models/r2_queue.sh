#!/bin/bash
# Run a queue of R2 configs sequentially. Usage: r2_queue.sh <queue_file> <threads> <log_prefix>
# queue file lines: name|folds|stage|set-args...
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # code checkout
PY="${KSP_PYTHON:-python}"                                      # interpreter with the requirements installed
cd "${KSP_DATA_DIR:-$REPO}"                                     # project data/working dir (cache/, logs/, results/)
Q=$1; T=${2:-12}; L=${3:-logs/r2/queue}
while IFS='|' read -r name folds stage sets; do
  [ -z "$name" ] && continue
  [[ "$name" == \#* ]] && continue
  echo "=== $(date +%H:%M:%S) start $name folds=$folds stage=$stage sets=$sets"
  PYTHONPATH="$REPO/src" R2_THREADS=$T OMP_NUM_THREADS=$T "$PY" "$REPO"/scripts/models/r2_run.py --name "$name" --folds $folds --stage "$stage" --set $sets > "${L}_${name}.log" 2>&1
  echo "=== $(date +%H:%M:%S) end $name rc=$?"; grep -E "M1/M5|RESULT" "${L}_${name}.log" | head -2 | cut -c1-400
done < "$Q"
