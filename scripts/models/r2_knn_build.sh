#!/bin/bash
cd /workspace
[ -n "$WAIT_PID" ] && while kill -0 $WAIT_PID 2>/dev/null; do sleep 20; done
CUDA_VISIBLE_DEVICES="" PYTHONPATH=src OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 .venv/bin/python scripts/models/r2_xfit.py --family knn --folds 0 1 2 3 4
