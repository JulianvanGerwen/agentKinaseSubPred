cd /workspace
export PYTHONPATH=src
while pgrep -f "[e]val_chain.sh" > /dev/null; do sleep 20; done
B=$(ls -d results/runs/regime_L2_L2fixed_s0_cv5__* | head -1 | xargs basename)
for r in regime_L1_refit_s0_cv5 regime_L2_vitro_refit_s0_cv5 regime_perkinase_lgb_N20_qmap_cv5 regime_ppool_lam300_cv5 regime_L2_vivo_refit_s0_cv5; do
  A=$(ls -d results/runs/${r}__* | head -1 | xargs basename)
  .venv/bin/python scripts/models/r2_compare.py --a $A --b $B --B 1000 --out results/models/regime/compare__${r}__vs__L2.tsv > logs/regime/compare_${r}.log 2>&1
done
