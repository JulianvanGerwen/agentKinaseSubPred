PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L1 --feat refit --fold0-eval > logs/regime/run_L1_refit.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2_vitro --feat refit --fold0-eval > logs/regime/run_L2_vitro_refit.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2_vivo --feat refit --fold0-eval > logs/regime/run_L2_vivo_refit.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2ext --feat refit --fold0-eval > logs/regime/run_L2ext_refit.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_pergroup.py --folds 0 1 2 3 4 --threads 6 > logs/regime/pergroup.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L1 --feat L2fixed --fold0-eval > logs/regime/run_L1_L2fixed.log 2>&1
PYTHONPATH=src R2_THREADS=6 R2_GATHER_THREADS=4 OMP_NUM_THREADS=6 .venv/bin/python scripts/models/regime_run.py --variant L2_vitro --feat L2fixed --fold0-eval > logs/regime/run_L2_vitro_L2fixed.log 2>&1
