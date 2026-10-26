cd /workspace
export PYTHONPATH=src
while pgrep -f "[c]ompare_chain.sh" > /dev/null; do sleep 20; done
B=$(ls -d results/runs/regime_L2_L2fixed_s0_cv5__* | head -1 | xargs basename)
A=$(ls -d results/runs/regime_L2ext_refit_s0_cv5__* | head -1 | xargs basename)
.venv/bin/python scripts/models/r2_compare.py --a $A --b $B --B 1000 --out results/models/regime/compare__regime_L2ext_refit_s0_cv5__vs__L2.tsv > logs/regime/compare_L2ext_s0.log 2>&1
until grep -q "DONE\|Traceback" logs/regime/run_L2ext_refit_s1.log; do sleep 20; done
A=$(ls -d results/runs/regime_L2ext_refit_s1_cv5__* | head -1 | xargs basename)
.venv/bin/python scripts/models/r2_compare.py --a $A --b r2_final_s1_cv5__20261026-181401__s1 --B 1000 --out results/models/regime/compare__regime_L2ext_refit_s1_cv5__vs__r2_final_s1.tsv > logs/regime/compare_L2ext_s1.log 2>&1
