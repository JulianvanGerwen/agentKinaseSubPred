> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# G4 runbook: FREEZE.json and the single TEST batch

Owner: core-infra. The PI commits FREEZE.json and runs one command. TEST labels are read only by `src/ksp/eval_test.py`, once.

## Deadline: 04:30 CET (G4 at ~05:20)

By 04:30, each model owner delivers the items below for **every frozen config and every baseline** that FREEZE.json will list: the final model, B1 arms (i) and (ii), B3-retrained, B4-retrained, B7, B0 floors, and any secondary configs.

### 1. TEST predictions: `cache/runs/<run_id>/test_predictions.parquet`

**Model:** trained on all of DEV with the frozen configuration (§8.6). Deployment models trained on DEV+TEST must not be used.

**Columns:**
- `kinase_idx` (int, row of `cache/index/kinases.parquet`), `site_idx` (int, row of `cache/index/sites.parquet`) and `score` (float, higher = more likely).
- `kinase_acc`/`site_id` can be used instead of the two index columns.
- Optional `score_site`: the cross-kinase score used for M2 (A20). It can be **any finite value**, e.g. R2's raw fused score chosen on CV. If absent, M2 uses `score`. It may also be supplied as a separate file (`score_site_predictions`).
- Optional `assignment_prob`: the A17 site-level assignment probability P(k | s), **in [0, 1]**. It is used only for the A17 reliability analysis (with ECE) and Coverage@τ, never for an endpoint. It may also be supplied as a separate file (`assignment_prob_predictions`).
- A separate file must use the same key columns as the main file, carry exactly one value column, and have no duplicate keys. A column cannot be given both inline and as a separate file.
- No canary kinases. No NaN or inf. No duplicate (kinase, site) rows.

**Rows must be exhaustive over the TEST candidate set:**
- The set is every U_ms site on a TEST protein (`results/splits/test_proteins.tsv`, 69,430 sites) × every kinase whose residue class covers the site. That gives 25,922,452 pairs over 503 kinases (`ksp.eval_test.test_candidate_pairs()`).
- **All-or-nothing per kinase:** a method either scores all of a kinase's candidates or leaves that kinase unscored.
- Never score a label-selected subset, such as positives + E-hard only. That inflates E-all metrics (M3), as seen on CV for B3-pub and for non-atlas B1 kinases.
- Unscored kinases count as missing. They are ranked last and reported. Head-to-head tests use the matched universe (PRIOR-1, A27).

**Pre-flight:** run this before you report the file. It reads no labels, and `"ok": true` is required.

    python -I src/ksp/eval_test.py --check-predictions cache/runs/<run_id>/test_predictions.parquet

### 2. T3 predictions (final model, and any baseline with an LKO variant): `cache/runs/<run_id>/test_t3_predictions.parquet`

- **Same format** as item 1.
- **Each kinase is scored by the model trained with that kinase's LKO family removed** (`results/splits/kinase_folds.tsv`, fold = `lko_fold`). Training data is DEV.
- **Rows:** TEST candidates of the held-out kinases of each LKO fold. All five folds are concatenated, so every kinase appears exactly once.

### 3. Run manifest: `results/runs/<run_id>/manifest.json` (ksp.manifest)

It must contain the git hash, config, seed, input checksums and fallbacks. It also needs the **wall-clock runtimes** that FREEZE.json requires: `train_wall_s` (training on DEV) and `inference_wall_s` (scoring the TEST candidates).

### 4. A FREEZE.json entry, sent to the PI

The schema is documented in `src/ksp/eval_test.py`. Every file is **pinned by SHA-256**: parquets are gitignored, so this pin is the only record of their contents. Compute the pins with `ksp.eval_test.pin_files(cfg)`.

    {"name": "final", "role": "final", "run_id": "<run_id>",
     "predictions": "cache/runs/<run_id>/test_predictions.parquet", "predictions_sha256": "<64 hex>",
     "score_site_predictions": null, "assignment_prob_predictions": null,
     "t3_predictions": "cache/runs/<run_id>/test_t3_predictions.parquet", "t3_predictions_sha256": "<64 hex>",
     "runtime": {"train_wall_s": 1234.5, "inference_wall_s": 210.0}, "missing": "worst"}

**Baselines in the primary family also carry `baseline_id`:** `B1_arm_i`, `B1_arm_ii`, `B3_retrained` or `B4_retrained`. Other baselines can use any id, e.g. `B7`.

**`score_site`**, if given, is any finite cross-kinase score used for M2. **`assignment_prob`**, if given, is the A17 probability in [0, 1], used for reliability and coverage only.

After you deliver, do not modify the files. Any byte change breaks the pin, and the batch refuses to run.

## Pre-declared comparisons (DESIGN §6.8, A19, A23, A27)

There are exactly **9 primary tests**, each with B = 10,000 paired cluster-bootstrap resamples on the matched universe. Holm correction covers those 9; secondary comparisons get Holm separately. `validate_freeze` checks:
- exactly one config has `role: "final"`, and it is `a` in every primary comparison;
- each primary `b` carries a `baseline_id` whose endpoints match the table below;
- there are no duplicate (baseline_id, endpoint) pairs and no M4.

| a | b (baseline_id) | endpoints |
|---|---|---|
| final | B1_arm_i: Cantley as published | M1, M3 |
| final | B1_arm_ii: Cantley universe-recalibrated | M2 |
| final | B3_retrained: Phosformer | M1, M2, M3 |
| final | B4_retrained: SELPHI2.0 | M1, M2, M3 |

**Contingency (PI decision):** if B3 or B4 predictions are not ready by FREEZE, list the baseline as `"unavailable_primary": [{"baseline_id": "B4_retrained", "reason": "..."}]`. The family stays at m = 9. The absent tests enter Holm with p = 1 and are reported as "NOT TESTED (baseline unavailable)".

"Beats" = Holm-adjusted p < 0.05 **and** an (unadjusted) paired 95% CI lower bound > 0. Everything else goes in `"primary": false`: B7, B0 floors, ablations, and the contaminated B3-pub.

## What the PI does at G4

1. Assemble FREEZE.json from the owners' entries and the comparisons above.
   - Set `"B": 10000` and `"harness_commit": "974d520"` (A28).
   - Optionally set `"B_secondary": 1000`, the B used for standalone, strata and sensitivity CIs.
2. Validate it without touching sealed/. This checks the schema, the primary family and the file pins:

       python -I -c "import sys,json; sys.path.insert(0,'src'); from ksp import eval_test as T; T.validate_freeze(json.load(open('FREEZE.json')))"

3. Commit FREEZE.json.
4. Run the batch:

       python -I src/ksp/eval_test.py --real

   Without `--real` it only prints usage. What the run does, in order:
   - **Refuses** if:
     - FREEZE.json is untracked or modified;
     - B ≠ 10000;
     - harness_commit ≠ 974d520, or the frozen harness files differ from it;
     - any prediction file fails pre-flight (exhaustiveness, finite scores, no canaries; T3 must cover all of K_any);
     - a seal hash mismatches.
   - **Ledger:** writes a STARTED block before reading sealed/ (FREEZE sha, input hashes, git HEAD, SHA-256 of every `src/ksp/*.py`), and a COMPLETED block at the end.
   - **Outputs:** written to `results/test/<sha12>.partial/`, then atomically renamed to `results/test/<sha12>/`.
5. **Crash:** rerun the same command unchanged. An identical FREEZE sha, input hashes and code hash is logged as attempt n; anything else is refused.
6. **After a COMPLETED TEST block, every further batch is refused**, whatever the FREEZE sha. The only escape is `--override-completed "<reason>"`, which is logged in the ledger.
7. Outputs in `results/test/<sha12>/`:
   - per config: metrics json, cluster CIs, kinase-resampling CIs (§6.8), per-kinase and site-level tables, strata, ext_tier A/B, A17 reliability + Coverage@τ (configs with assignment_prob), and T3 (+ missing-kinase count);
   - `superiority.tsv`, `secondary.tsv`, `summary.md`, `summary.json`, `preflight.json`.

Runtime: the B = 10,000 rehearsal on DEV fold 0 (8.5M required pairs) took about 21 min with a 5.9 GB peak. TEST has about 3× the pairs, so allow roughly 45–60 min.

## Rehearsal (any time; never reads sealed/)

    python -I src/ksp/eval_test.py --dry-run --fold 0 --B 1000

`--B` and `--fold` are allowed only with `--dry-run`. B3/B4 have no fold-0 runs, so the rehearsal exercises the unavailable_primary contingency. The rehearsal uses DEV CV fold 0 as a fake TEST, with existing fold-0 / pooled-OOF runs as configs, under `cache/eval_test_dryrun/`. `tests/test_eval_test_batch.py` runs the same path on synthetic scorers with known answers.
