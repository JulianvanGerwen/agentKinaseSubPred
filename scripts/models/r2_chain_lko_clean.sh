#!/bin/bash
# waits for LKO xfit / bench xfit / final xfit, then runs LKO T2, clean CV, clean U_obs library (+downstream)
cd /workspace
V3='num_leaves=127 learning_rate=0.03 f2a_family=f2a_sep desc_v=2 es_real_only=true'
n_json() { ls $1 2>/dev/null | wc -l; }
(
  until [ $(n_json "cache/r2/xfit/lko*/all/*.json") -ge 15 ]; do sleep 30; done
  echo "$(date +%H:%M) LKO start"
  PYTHONPATH=src R2_THREADS=8 OMP_NUM_THREADS=8 .venv/bin/python scripts/models/r2_lko.py --name r2_final_lko --set $V3 seed=0 > logs/r2/lko_final.log 2>&1
  echo "$(date +%H:%M) LKO done: $(grep LKO-T2 logs/r2/lko_final.log)"
) &
(
  until [ $(n_json "cache/r2/xfit/cv/fold*/f6a_bench.json") -ge 5 ]; do sleep 30; done
  echo "$(date +%H:%M) clean CV start"
  scripts/models/r2_cv5.sh r2_final_cleanbench_s0 8 clean_bench "$V3 seed=0 clean=\"bench\""
  until [ $(n_json "cache/r2/xfit/final/all/*.json") -ge 4 ]; do sleep 30; done
  echo "$(date +%H:%M) clean U_obs library start"
  PYTHONPATH=src R2_THREADS=8 OMP_NUM_THREADS=8 .venv/bin/python scripts/models/r2_final.py --name r2_final_cleanbench_uobs --target uobs --seeds 0 --set $V3 'clean="bench"' > logs/r2/final_clean_uobs.log 2>&1
  echo "$(date +%H:%M) clean library done: $(grep RUN logs/r2/final_clean_uobs.log)"
) &
wait
