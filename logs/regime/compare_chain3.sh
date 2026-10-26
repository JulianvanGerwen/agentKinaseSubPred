cd /workspace
export PYTHONPATH=src
B=$(ls -d results/runs/regime_L2_L2fixed_s0_cv5__* | head -1 | xargs basename)
for r in regime_L2ext_refit_s0_cv5 regime_L1_refit_s0_cv5 regime_L2_vitro_refit_s0_cv5 regime_perkinase_lgb_N20_qmap_cv5 regime_ppool_lam300_cv5; do
  A=$(ls -d results/runs/${r}__* | head -1 | xargs basename)
  .venv/bin/python scripts/models/r2_compare.py --a $A --b $B --B 1000 --out results/models/regime/compare__${r}__vs__L2_s0.tsv > logs/regime/compare_${r}.log 2>&1
  if [ "$r" == "regime_L2ext_refit_s0_cv5" ]; then
    A=$(ls -d results/runs/regime_L2ext_refit_s1_cv5__* | head -1 | xargs basename)
    .venv/bin/python scripts/models/r2_compare.py --a $A --b r2_final_s1_cv5__20261026-181401__s1 --B 1000 --out results/models/regime/compare__regime_L2ext_refit_s1_cv5__vs__r2_final_s1.tsv > logs/regime/compare_L2ext_s1.log 2>&1
  fi
done
