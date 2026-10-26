# G4 runbook: FREEZE.json and the single TEST batch

Owner: core-infra. The PI commits FREEZE.json and runs one command. TEST labels are read only by `src/ksp/eval_test.py`, once.

## Deadline: 04:30 CET (G4 at ~05:20)

By 04:30, each model owner delivers the items below for **every frozen config and every baseline** that FREEZE.json will list: the final model, B1 arms (i) and (ii), B3-retrained, B4-retrained, B7, B0 floors, and any secondary configs.

### 1. TEST predictions: `cache/runs/<run_id>/test_predictions.parquet`

**Model:** trained on all of DEV with the frozen configuration (§8.6). Deployment models trained on DEV+TEST must not be used.

**Columns:**
- `kinase_idx` (int, row of `cache/index/kinases.parquet`), `site_idx` (int, row of `cache/index/sites.parquet`) and `score` (float, higher = more likely).
- `kinase_acc`/`site_id` can be used instead of the two index columns.
- Optional `score_site`: the cross-kinase assignment score used for M2 (A20). If absent, M2 uses `score`. It may also be supplied as a separate file (`score_site_predictions`).
- No canary kinases. No NaN or inf. No duplicate (kinase, site) rows.

**Rows must be exhaustive over the TEST candidate set:**
- The set is every U_ms site on a TEST protein (`results/splits/test_proteins.tsv`, 69,430 sites) × every kinase whose residue class covers the site. That gives 25,922,452 pairs over 503 kinases (`ksp.eval_test.test_candidate_pairs()`).
- **All-or-nothing per kinase:** a method either scores all of a kinase's candidates or leaves that kinase unscored.
- Never score a label-selected subset, such as positives + E-hard only. That inflates E-all metrics (M3), as seen on CV for B3-pub and for non-atlas B1 kinases.
- Unscored kinases count as missing. They are ranked last and reported. Head-to-head tests use the matched universe (PRIOR-1, A27).

**Pre-flight:** run this before you report the file. It reads no labels, and `"ok": true` is required.

    /workspace/.venv/bin/python -I src/ksp/eval_test.py --check-predictions cache/runs/<run_id>/test_predictions.parquet

### 2. T3 predictions (final model, and any baseline with an LKO variant): `cache/runs/<run_id>/test_t3_predictions.parquet`

- **Same format** as item 1.
- **Each kinase is scored by the model trained with that kinase's LKO family removed** (`results/splits/kinase_folds.tsv`, fold = `lko_fold`). Training data is DEV.
- **Rows:** TEST candidates of the held-out kinases of each LKO fold. All five folds are concatenated, so every kinase appears exactly once.

### 3. Run manifest: `results/runs/<run_id>/manifest.json` (ksp.manifest)

It must contain the git hash, config, seed, input checksums and fallbacks. It also needs the **wall-clock runtimes** that FREEZE.json requires: `train_wall_s` (training on DEV) and `inference_wall_s` (scoring the TEST candidates).

### 4. A FREEZE.json entry, sent to the PI

The schema is documented in `src/ksp/eval_test.py`.

    {"name": "final", "role": "final", "run_id": "<run_id>",
     "predictions": "cache/runs/<run_id>/test_predictions.parquet",
     "score_site_predictions": null, "t3_predictions": "cache/runs/<run_id>/test_t3_predictions.parquet",
     "runtime": {"train_wall_s": 1234.5, "inference_wall_s": 210.0}}

## Pre-declared comparisons (DESIGN §6.8, A19, A23, A27)

There are exactly **9 primary tests**, each with B = 10,000 paired cluster-bootstrap resamples on the matched universe. Holm correction covers those 9; secondary comparisons get Holm separately.

| a | b | endpoints |
|---|---|---|
| final | B1 Cantley arm (i) as published | M1, M3 |
| final | B1 Cantley arm (ii) universe-recalibrated | M2 |
| final | B3 Phosformer retrained | M1, M2, M3 |
| final | B4 SELPHI2.0 retrained | M1, M2, M3 |

"Beats" = Holm-adjusted p < 0.05 **and** paired 95% CI of the difference > 0. Everything else (B7, B0 floors, ablations, B3-pub contaminated) goes in `"primary": false`.

## What the PI does at G4

1. Assemble FREEZE.json from the owners' entries and the comparisons above, and set `"harness_commit"` (A28).
2. Validate it without touching sealed/:

       /workspace/.venv/bin/python -I -c "import sys,json; sys.path.insert(0,'src'); from ksp import eval_test as T; T.validate_freeze(json.load(open('FREEZE.json')))"

3. Commit FREEZE.json.
4. Run the batch:

       /workspace/.venv/bin/python -I src/ksp/eval_test.py

   It refuses to run if FREEZE.json is untracked or modified, if a seal hash mismatches, or if this FREEZE.json was already scored.
5. Outputs:
   - `results/test/<freeze_sha12>/`: per-config metrics, CIs, strata and T3; `superiority.tsv`; `secondary.tsv`; `summary.md`.
   - One block appended to `results/test/LEDGER.md`.

## Rehearsal (any time; never reads sealed/)

    /workspace/.venv/bin/python -I src/ksp/eval_test.py --dry-run --fold 0 --B 1000

This uses DEV CV fold 0 as a fake TEST, with existing fold-0 / pooled-OOF runs as configs, under `cache/eval_test_dryrun/`. `tests/test_eval_test_batch.py` runs the same path on synthetic scorers with known answers.
