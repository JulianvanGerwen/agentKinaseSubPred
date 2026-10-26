sleep 60
while pgrep -f "[r]egime_partial_pool_eval.py|[r]egime_resboost_eval.py" > /dev/null; do sleep 15; done
PYTHONPATH=src .venv/bin/python scripts/models/regime_pergroup_eval.py > logs/regime/pergroup_eval.log 2>&1
