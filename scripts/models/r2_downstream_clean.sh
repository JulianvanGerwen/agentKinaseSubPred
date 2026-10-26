#!/bin/bash
cd /workspace
until grep -q "clean library done" logs/r2/chain_lko_clean.log; do sleep 30; done
RID=$(grep "^RUN" logs/r2/final_clean_uobs.log | awk '{print $2}')
echo "clean run $RID"
OMP_NUM_THREADS=8 .venv/bin/python scripts/downstream/run_predictions.py --pred cache/runs/$RID/uobs_predictions.parquet --tag r2_clean > logs/r2/downstream_r2_clean.log 2>&1
echo "downstream done rc=$? $(date +%H:%M)"
