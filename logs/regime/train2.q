PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2ext --feat refit --fold0-eval > logs/regime/run_L2ext_refit.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L1 --feat L2fixed --fold0-eval > logs/regime/run_L1_L2fixed.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2_vitro --feat L2fixed --fold0-eval > logs/regime/run_L2_vitro_L2fixed.log 2>&1
