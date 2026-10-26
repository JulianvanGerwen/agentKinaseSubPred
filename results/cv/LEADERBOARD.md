# CV leaderboard (pooled OOF, re-evaluated with ksp.eval A24)

Generated 2026-10-26 17:07:53 CET by `scripts/reeval_cv_runs.py` from saved predictions (metrics_a24.json per run). Sorted by M1 (macro E-hard AUROC over K_cv, fold rank-normalised; LKO pooled runs not normalised). M5bin = macro E-hard AUROC over K_any by DEV training-label bin; M1_ST/M1_Y = K_any macro by residue class. Coverage = fraction of required (kinase, U_ms site) pairs scored (missing = ranked last). Canary gate (v2): mean negative E-hard <= 0.57 and positive >= 0.90; '(v1-trained)' = run trained before the v2 canary labels (gate not meaningful). Official baselines per A23 (latest run of each arm).

| # | run | family | official | M1 | M1 fold-mean | M2 Hit@5 | M3 top1% | M4 | M5 0 | 1-4 | 5-19 | 20-99 | >=100 | M1 ST | M1 Y | canary neg / pos | gate | coverage | runtime s | started |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | b3pub_phosformer_cv_pooled__20261026-151731__s0 | B3 Phosformer |  | 0.796 | 0.796 | 0.248 | 0.852 | 0.804 |  | 0.807 | 0.805 | 0.797 | 0.817 | 0.846 | 0.616 | 0.500 / 0.500 | FAIL (v1-trained) | 0.020 | 0.18 | 2026-10-26 15:17:37 CET |

Config summaries: see LEADERBOARD.tsv (column `config`).
