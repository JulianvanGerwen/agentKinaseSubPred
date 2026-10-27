# Pipeline scripts

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Each script states its inputs, outputs and usage in its docstring. Every script reads and writes under the project
directory `ksp.paths.ROOT`, which is `$KSP_DATA_DIR`, or the repository root if the variable is unset. The shell
wrappers change into that directory, and they call Python through `$KSP_PYTHON` (default: `python`). Run Python
scripts with `python -I` from any directory after `pip install -e .`. Without an install, the scripts add
`<repo>/src` to `sys.path` themselves.

None of these scripts can run without the data in `DATA_SOURCES.md`. Scripts marked GPU need a CUDA device.
Exact command lines for the frozen configuration are in `docs/REPORT_internal.md` §12, `docs/model_r2.md` §6–7,
`docs/model_r3.md`, `docs/baseline_*.md` and `docs/G4_RUNBOOK.md`.

## Order of stages

1. **Downloads** (public data only): `download_uniprot.py`, `download_afdb_pae.py`, `dl.sh`, `dl_w2.sh`,
   `dl_w2b.sh`, `dl_pdc_pancancer.py`. Licensed data (PhosphoSitePlus, Klaeger, the compendium) must be placed by
   hand.
2. **Foundations:** `build_kinase_list.py` → `build_site_universe.py` → `build_labels.py` → `build_splits.py`
   (needs MMseqs2; writes the frozen splits and seals the TEST labels) → `build_index.py`. Then
   `strata_definitions.py` and `build_compendium_matrix.py`.
3. **Site-intrinsic features:** `python -m ksp.features.structure` and `python -m ksp.features.esm_embed` and `python -m ksp.features.esm_kinases` (GPU), then
   `export_structure_esm_features.py` and `validate_structure_esm_features.py`; `models/r3_window_embed.py` (GPU).
4. **Feature families:**
   - motif: `features/motif_extract_supp.py`, `features/motif_cantley_build.py`, `features/motif_cantley_validate.py`,
     `features/motif_validate.py`;
   - context: `features/context_ppi.py`, `features/context_loc.py`, `features/context_coexpr.py`, with checks in
     `features/context_ppi_label.py` and `features/context_sanity.py`;
   - phosphoproteomics: `features/phospho_study_lists.py`, `features/phospho_build_blocks.py`, `features/phospho_f6b.py`,
     `features/phospho_f6c.py`, with checks in `features/phospho_selftest.py` and `features/phospho_sanity.py`;
   - similarity to training substrates: `features/knn_build_reps.py`, `features/knn_tune_tau.py`,
     `features/knn_xfit.py`, `features/knn_sanity.py`;
   - Tyr features (evaluated, not adopted): `features/tyr_build.py`, `features/tyr_eval.py`.
5. **Cross-fitted label-derived stores:** `models/r2_xfit.py` (CV folds; `--toy` runs the leak test),
   `models/r2_final_xfit.py`, `models/r2_lko_xfit.py`, `models/r2_t3_xfit.py`; wrapper `models/r2_knn_build.sh`.
6. **R2 on CV:** `models/r2_run.py` per fold and `models/r2_pool.py`, or the wrappers `models/r2_cv5.sh`,
   `models/r2_cv5_queue.sh` and `models/r2_queue.sh`; then `models/r2_seed_ens.py`, `models/r2_reeval.py`,
   `models/r2_combine_m2.py`, `models/r2_compare.py`, `models/r2_partial_pool.py` and `models/r2_lko.py` (LKO T2).
7. **R3 on CV (GPU):** `models/r3_cv.py train|assemble|test` with the wrappers `models/r3_queue.sh`,
   `models/r3_dev_sweep.sh` and `models/r3_test_queue.sh`; `models/r3_paired.py`.
8. **Fusion and calibration:** `models/r2_fusion.py`, `models/r2_fusion_lko.py`, `models/r2_calib.py`.
9. **Baselines:** `baselines/simple_*.py` (B0 controls, B1 Cantley arms, B2, B7); `baselines/phosformer_*.py` and
   `baselines/phosformer_pretrain_run.sh` (B3, GPU); `baselines/selphi2_*.py` (B4).
10. **CV reports:** `pool_fold_runs.py`, `reeval_cv_runs.py` (leaderboard), `strata_cv_report.py`.
11. **TEST deliverables and FREEZE:** `models/r2_final.py`, `models/r2_fuse_test.py`, `models/r2_assign_test.py`,
    `models/r2_freeze_deliverable.py`, `models/r2_topk_filter.py`, `baselines/simple_test.py`,
    `baselines/selphi2_final.py`, `baselines/phosformer_test.py`, `baselines/simple_freeze_entries.py`,
    `build_freeze.py`; then the single batch `python -I src/ksp/eval_test.py --real`. The batch for this study has
    been run, so it must not be repeated on the same split.
12. **Downstream benchmark:** `downstream/selftest.py`, `downstream/run_validation.py`, `downstream/pred_topk.py`,
    `downstream/run_predictions.py`, `downstream/run_final.py`.
13. **Deployment and release:** `models/r2_deploy.py --stage xfit|train`, then `models/r2_release.py`. Together they
    write the files in `predictions/`.
14. **Leakage audits:** `audit_ppi_ks_leakage.py`, `phosformer_split_leakage_check.py`; data profiling:
    `profile_compendium.py` (also imported by `features/phospho_build_blocks.py`).

## Not included

These scripts are in the project record but not in this MVP. They belong to exploratory or negative branches, to
one-off maintenance, or to data profiling of PhosphoSitePlus:
- the R5 cross-encoder (negative; `docs/model_r5.md`);
- the training-regime and label-variant sweep (`docs/regimes_labels.md`);
- the post-TEST experiments (`docs/posttest_cv.md`, `docs/posttest_tyr.md`);
- priming and M2-prior diagnostics, canary diagnostics and Tyr diagnostics;
- the PhosphoSitePlus profiling and mapping checks;
- the generators of the data register;
- run-queue chains that were tied to the original machine.

The modules `ksp.models.r5` and `ksp.models.regime` are left out for the same reason.
