cd /workspace
until grep -q "DONE\|Traceback" logs/regime/run_L2_vitro_L2fixed.log 2>/dev/null; do sleep 20; done
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2ext --feat refit --seed 2 > logs/regime/run_L2ext_refit_s2.log 2>&1
