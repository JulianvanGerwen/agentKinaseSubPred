cd /workspace
export PYTHONPATH=src
.venv/bin/python scripts/models/regime_pergroup_eval.py > logs/regime/pergroup_eval.log 2>&1
.venv/bin/python scripts/models/regime_resboost_eval.py --stage cv5 > logs/regime/resboost_eval_cv5.log 2>&1
.venv/bin/python scripts/models/regime_perkinase_eval.py --folds 0 1 2 3 4 > logs/regime/perkinase_eval_cv5b.log 2>&1
