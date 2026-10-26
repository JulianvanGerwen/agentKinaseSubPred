#!/bin/bash
# R5: run the remaining CV folds of one config over GPUs 1-3 (GPU 0 is Phosformer's), <= 12 cores total.
#   bash scripts/models/r5_queue.sh <config> "<folds>"      e.g. bash scripts/models/r5_queue.sh li35_cpt_aux "1 2 3 4"
set -u
cd "$(dirname "$0")/../.."
CFG=${1:?config}
FOLDS=${2:-"1 2 3 4"}
mkdir -p logs/r5
GPUS=(1 2 3)
i=0
for f in $FOLDS; do
  g=${GPUS[$((i % 3))]}
  CUDA_VISIBLE_DEVICES=$g OMP_NUM_THREADS=4 nohup .venv/bin/python -I scripts/models/r5_cv.py \
      train --fold "$f" --config "$CFG" --n-seeds 3 --note "5-fold CV, 3 seeds" \
      > "logs/r5/f${f}_${CFG}.log" 2>&1 &
  echo "fold $f -> GPU $g (pid $!)"
  i=$((i + 1))
  # stagger: the frozen-650M window store is 9.4 GiB per job, so do not load three at once
  sleep 45
  if [ $((i % 3)) -eq 0 ]; then wait; fi
done
wait
echo "done $CFG folds: $FOLDS"
