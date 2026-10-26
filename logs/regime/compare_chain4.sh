cd /workspace
export PYTHONPATH=src
B=$(ls -d results/runs/regime_L2_L2fixed_s0_cv5__* | head -1 | xargs basename)
for r in "$@"; do
  A=$(ls -d results/runs/${r}__* | head -1 | xargs basename)
  .venv/bin/python scripts/models/r2_compare.py --a $A --b $B --B 1000 --out results/models/regime/compare__${r}__vs__L2_s0.tsv > logs/regime/compare_${r}.log 2>&1
done
