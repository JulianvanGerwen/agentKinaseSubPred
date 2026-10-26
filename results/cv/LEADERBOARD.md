# CV leaderboard (pooled OOF, re-evaluated with ksp.eval A24)

Generated 2026-10-26 17:22:14 CET by `scripts/reeval_cv_runs.py` from saved predictions (metrics_a24.json per run). Sorted by M1 (macro E-hard AUROC over K_cv, fold rank-normalised; LKO pooled runs not normalised). M5bin = macro E-hard AUROC over K_any by DEV training-label bin; M1_ST/M1_Y = K_any macro by residue class. Coverage = fraction of required (kinase, U_ms site) pairs scored (missing = ranked last); runs with coverage < 0.95 have M3 blanked (E-all metrics not comparable: compare them on a matched universe, ksp.eval.matched_universe). Canary gate (v2): mean negative E-hard <= 0.57 and positive >= 0.90; '(v1-trained)' = run trained before the v2 canary labels (gate not meaningful). Official baselines per A23 (latest run of each arm).

| # | run | family | official | M1 | M1 fold-mean | M2 Hit@5 | M3 top1% | M4 | M5 0 | 1-4 | 5-19 | 20-99 | >=100 | M1 ST | M1 Y | canary neg / pos | gate | coverage | runtime s | started |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | b3pub_phosformer_cv_pooled__20261026-151731__s0 | B3 Phosformer |  | 0.796 | 0.796 | 0.248 |  | 0.804 |  | 0.807 | 0.805 | 0.797 | 0.817 | 0.846 | 0.616 |  /  | canaries not scored | 0.020 | 0.18 | 2026-10-26 15:17:37 CET |
| 2 | r2_alarm_compknow_cv5__20261026-165100__s0 | R2 LightGBM |  | 0.776 | 0.773 | 0.204 | 0.320 | 0.764 |  | 0.740 | 0.754 | 0.783 | 0.812 | 0.800 | 0.599 | 0.531 / 0.998 | PASS | 1.000 | 274.54 | 2026-10-26 16:51:14 CET |
| 3 | r2_g2_main_cv5__20261026-162249__s0 | R2 LightGBM |  | 0.754 | 0.751 | 0.146 | 0.272 | 0.735 |  | 0.689 | 0.739 | 0.750 | 0.789 | 0.774 | 0.555 | 0.539 / 0.997 | PASS | 1.000 | 272.38 | 2026-10-26 16:22:53 CET |
| 4 | r2_abl_F1F0F2a_cv5__20261026-165405__s0 | R2 LightGBM |  | 0.738 | 0.733 | 0.094 | 0.239 | 0.720 |  | 0.683 | 0.720 | 0.733 | 0.775 | 0.756 | 0.554 | 0.532 / 0.793 | FAIL | 1.000 | 236.33 | 2026-10-26 16:54:08 CET |
| 5 | r2_abl_F1F0_cv5__20261026-163715__s0 | R2 LightGBM |  | 0.733 | 0.729 | 0.169 | 0.236 | 0.717 |  | 0.678 | 0.718 | 0.734 | 0.764 | 0.756 | 0.544 | 0.545 / 0.450 | FAIL | 1.000 | 262.83 | 2026-10-26 16:37:18 CET |
| 6 | simple_b1_imputed_universe_reeval__20261026-160413__s0 | B1 Cantley |  | 0.722 | 0.723 | 0.149 | 0.132 | 0.701 |  | 0.652 | 0.698 | 0.726 | 0.766 | 0.742 | 0.542 |  /  | canaries not scored | 0.966 | 196.64 | 2026-10-26 16:04:16 CET |
| 7 | simple_b1_imputed_published__20261026-151630__s0 | B1 Cantley |  | 0.721 | 0.723 | 0.134 | 0.132 | 0.700 |  | 0.651 | 0.696 | 0.725 | 0.766 | 0.742 | 0.542 |  /  | canaries not scored | 0.966 | 211.97 | 2026-10-26 15:16:30 CET |
| 8 | b4_selphi2_faithful__20261026-164243__s0 | B4 SELPHI2 |  | 0.587 | 0.590 | 0.112 |  | 0.582 |  | 0.573 | 0.591 | 0.572 | 0.596 | 0.582 | 0.545 |  /  | canaries not scored | 0.329 | 672.87 | 2026-10-26 16:42:43 CET |
| 9 | simple_b0_miss__20261026-150214__s0 | B0 floor |  | 0.539 | 0.544 | 0.040 | 0.053 | 0.541 |  | 0.561 | 0.529 | 0.546 | 0.526 | 0.524 | 0.576 |  /  | canaries not scored | 0.966 | 320.65 | 2026-10-26 15:02:14 CET |
| 10 | simple_b0_miss_reeval__20261026-155334__s0 | B0 floor | B0-miss floor | 0.539 | 0.544 | 0.040 | 0.053 | 0.541 |  | 0.561 | 0.529 | 0.546 | 0.526 | 0.524 | 0.576 |  /  | canaries not scored | 0.966 | 190.97 | 2026-10-26 15:53:38 CET |
| 11 | simple_b0_prior__20261026-150125__s0 | B0 floor | B0-prior floor | 0.500 | 0.500 | 0.064 | 0.010 | 0.500 |  | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |  /  | canaries not scored | 0.966 | 131.59 | 2026-10-26 15:01:25 CET |

Config summaries: see LEADERBOARD.tsv (column `config`).
