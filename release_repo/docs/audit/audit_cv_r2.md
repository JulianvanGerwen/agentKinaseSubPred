> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Independent audit, CV round 2: R2 final candidate and R2+R3 late fusion (pre-FREEZE)

**Auditor:** an independent results-auditor agent. I produced none of the audited work.
**When:** 2026-10-26, 21:24–22:15 CET.
**Scope:** `docs/model_r2.md` §6 onwards, `results/models/r2/`, the runs listed there, and the code paths they use.

**Constraints I kept:**
- at most 3 evaluation workers, under 60 GB RAM;
- no git, no reads of `sealed/`, no edits to other agents' files.

**How numbers were checked:**
- Every number was recomputed from the saved `cache/runs/*/predictions.parquet` files with the frozen harness: `ksp.eval` with A24 on, `eval.py` dated 17:01.
- Bootstrap CIs were **not** recomputed; only point estimates were checked.

**Scripts:** copied with their JSON/TSV outputs to `docs/audit/audit_cv_r2_scripts/`. They import from the session scratchpad path, so adjust `sys.path` before re-running them.

| Script | What it does |
|---|---|
| `pooled.py` | Reproduces the pooled numbers and checks R2/R3 row alignment |
| `fold_sweep.py` | Per-fold sweep of the fusion weight |
| `nested.py` | Leave-one-fold-out (LOFO) weight selection |
| `lko_strat.py` | T2 diagnostic |
| `clean_eval.py` | Clean-variant recomputation |
| `oof_check.py` | Out-of-fold (OOF) integrity of the 25 seed-fold runs |
| `calib_check*.py` | Assignment-probability reliability and recalibration |

## 1. Verified numbers (pooled 5-fold DEV OOF, A24; K_cv = 203)

| Claim (doc) | Reported | Recomputed | Notes |
|---|---|---|---|
| R2 5-seed ensemble `r2_final_ens5` M1 / M2 / M3 / M4 | 0.776 / 0.200 / 0.312 / 0.767 | **0.7760 / 0.2000 / 0.3115 / 0.7671** ✓ | Fold-mean M1 0.7715 ✓ |
| Fusion `r2_fusion_r3__20261026-191901` M1 | 0.786 | **0.7864** ✓ | Per fold 0.767 / 0.788 / 0.782 / 0.794 / 0.773, mean 0.7807 ✓ |
| Fusion M2 (score_site, w = 0.5) | 0.263 | **0.2626** ✓ | |
| Fusion M3 / M4 | 0.361 / 0.774 | **0.3609 / 0.7744** ✓ | |
| Fusion negative / positive canary (pooled, E-hard) | 0.539 / 1.000 | **0.5391 / 0.9996** ✓, PASS | Per fold (w = 0.6): 0.551 / **0.575** / 0.481 / 0.520 / **0.676**; mean 0.561 |
| LKO T2, `r2_final_lko` M1 / M3 / M4 | 0.751 / 0.295 / 0.749 | **0.7506 / 0.2946 / 0.7490** ✓ | With `fold_normalise=False`: 0.7504 / 0.2908 / 0.7493 |
| Clean variant `r2_final_cleanbench_s0_cv5` M1 / M2 / M3 / M4 | 0.771 / 0.189 / 0.300 / 0.768 | **0.7714 / 0.1891 / 0.2995 / 0.7675** ✓ | Canary 0.551 / 0.999 (pass) |
| Assignment tier [0.8, 1] (A17) | claimed 0.87, observed 0.54 (n = 167) | **0.868 vs 0.545** ✓ | Wilson 95% CI [0.47, 0.62] |
| **Fusion, LOFO-selected weights (audit)** | — | **M1 0.7864, M2 0.2589, M3 0.3609**, M4 0.7744 | §3: optimism ≈ 0 for M1, ≈ 0.004 for M2 |

**Integrity checks:**
- **Seed-run OOF.** All 25 seed × fold R2 runs contain only sites of their own fold. Each covers 100% of the fold's U_ms sites × 516 kinases.
- **Seed-run configs.** The configs are identical across the 25 runs; only the seed differs.
- **Stale stores.** The xfit stores they use (`f2a_sep`, `f3nb`, `f6a`) carry v2 canary counts for every fold. The round-1 issue M-2 (stale stores) is fixed; the rebuild finished at 17:47, before the first final-seed fold run started.
- **R2/R3 row alignment, per fold.** The R2 ensemble rows and the R3 `FINAL_OOF` fold rows have identical key sets. On every fold, n_r2 = n_r3 = n_common (21.2–21.4 M rows). Every R3 row lies in its own CV fold.
- **Folds match.** The fold is a function of the site's protein, so R2 and R3 rows are in the same folds by construction. The check above also verifies this directly.
- **Rebuilt blend.** My blend, rebuilt from the ensemble and R3 fold files, equals the saved fusion `score` and `score_site` to within 2.4e-7.
- **No missing pairs.** There are 0 missing required pairs across all 71.3 M rows.

## 2. Issues ranked by severity

### Critical
None. I found no label leakage, no eval-fold information in any fitted component, and no TEST access.

### Major
**R2-M1. The assignment probability (A17) is miscalibrated at the top, and on Y sites it points the wrong way. Coverage@0.8 as reported (1,674) is not supportable.**
- **Confirmed.** The [0.8, 1] tier claims 0.868 but observes 0.545 (n = 167, CI [0.47, 0.62]).
- **Stratified by residue:**

  | Sites | [0.5, 0.8) claimed / observed | [0.8, 1] claimed / observed |
  |---|---|---|
  | ST | 0.63 / 0.63 | 0.85 / 0.72 (n = 126) |
  | **Y** | 0.59 / 0.43 | **0.91 / 0.00 (n = 41)** |

  - Y top calls (EGFR, ABL1, ABL2, SYK, FYN) are never the curated kinase.
  - A cross-fitted Platt fit on Y sites gets a **negative** slope in 4 of 5 folds: a higher p_top means a *lower* chance of being correct.
- **Not explained by label strictness.** Counting a top kinase labelled to the site in any DEV source, or a paralog of a true kinase, raises the top tier only to 0.61–0.72.
- Sources: `results/models/r2/calib/r2_fusion_r3__20261026-191901__s0__calibrated/assignment.json`, `docs/model_r2.md` §6.1.
- **Minimal fix, CV-only:** add a **residue-stratified (ST / Y), fold-cross-fitted isotonic recalibration of p_top**, fitted on labelled DEV sites. For TEST and deployment, use the isotonic map fitted on all DEV OOF. It only post-processes the top-kinase probability: **no model, temperature or M1–M3 number changes.** Audit result for this fix:

  | Stratum | [0.2, 0.5) claimed / observed | [0.5, 0.8) claimed / observed | [0.8, 1] |
  |---|---|---|---|
  | ST | 0.32 / 0.32 | 0.62 / 0.61 (n = 1,272) | 2 sites |
  | Y | 0.41 / 0.41 | none reach 0.5 | none |

  - Coverage@0.5 becomes **21,115** (19,841 unlabelled), all of them ST sites.
  - Coverage@0.8 becomes **9**.
  - In other words, the model has **no reliable ≥ 0.8 tier and no reliable Y-site assignment ≥ 0.5**. That is what should be reported.
- **Alternatives I tested and rejected:**
  - Global Platt: [0.8, 1] still observes only 0.40.
  - Residue-stratified Platt: the ST top tier is 0.85 vs 0.71, and the slope on Y sites is negative.
- **Caveat:** reliability is measured on labelled sites only (PU, selection on being labelled), as A17 prescribes.

**R2-M2. The canary gate passes pooled but not robustly per fold (carried over from round 1, M-3).**
- **Pooled** (the gated value): negative canary 0.539, which passes. The positive canary is 0.9996.
- **Per fold** at w = 0.6: two folds exceed 0.57 (fold 1: 0.575; fold 4: 0.676), and the mean is 0.561. The fold-4 excess comes from R3: R3 alone scores 0.667 on fold 4, R2 alone 0.648.
- Every arm sits above 0.5. The v2-sampler explanation from the round-1 audit still applies, and there is no evidence of leakage.
- **The PI must rule before FREEZE** whether the pooled gate stands as written. I recommend recording the per-fold mean (0.561) next to it in FREEZE.

### Minor
**R2-m1. The M2 component deviates from the letter of A20 (my answer to claim 2).**
- **Facts:**
  - A20 says "M2 uses the cross-kinase assignment probability (A17)".
  - The candidate instead uses the raw fused `score_site` at w = 0.5, which scores M2 0.263 against 0.205 for the calibrated assignment logit.
  - Both arms were evaluated on cross-fitted CV OOF only (`r2_calib.py`: calibrators and temperature fitted on the other folds). TEST was not involved.
- **Optimism.** The choice between 2 arms on CV is immaterial next to a 0.058 gap. The w_M2 grid adds about 0.004 (§3).
- **Status.** The choice is legitimate under A20's second sentence ("each component is selected on its own CV endpoint"), but it contradicts its first. Record it as a DESIGN §17 clarification before FREEZE.

**R2-m2. The LKO T2 headline needs qualifying; no leak was found (claim 3).**
- **Label hold-out, verified in code and stores.** The held-out kinase fold's labels are removed from all of the following:
  - Training positives: `sample.FoldData(exclude_kinases)`, with an assert in `r2_lko.py:35`.
  - Negative masks: `sample.py:71`.
  - Every cross-fitted store (`xfit.build(exclude_kinases)`, `xfit.py:285`). The stores `lko{0..4}/all/{f2a_sep,f3nb,f6a}` carry n_excluded_kinases of 102 / 100 / 97 / 103 / 101 and 2,304–2,569 removed pairs. In every fold, fitted plus removed equals 12,223.
- **Per-kinase families.** F2a, F3-label and F6a are fitted per kinase:
  - F3-label (`ppi_label.py:79`) and F6a use only kinase k's own training substrates.
  - F2a's family and group priors cannot draw on the held-out family, because LKO folds merge whole Manning families.
- **Label-free inputs.** The F7 descriptors, version 2 (`assemble.py:248–300`), are label-free: ESM, Cantley-profile PCA, Manning categories. F1, F3-raw, F4, F5 and F6b/c are label-free too.
- **Early stopping** uses training kinases only.
- **The fused LKO (0.749)** uses R3's LKO runs, which also exclude held-out labels (round 1). Their pretraining does include the held-out kinases' label-free Cantley percentiles.
- **Qualification.** T2 is "seen proteins", and 45% of held-out positives are training positives of some *other* kinase. Splitting M1 by that:

  | Positives (n = held-out positives) | T2 M1 | CV control (same split) |
  |---|---|---|
  | also a training positive of another kinase (45%) | 0.798 | 0.799 |
  | not a training positive of any other kinase (55%) | **0.703** | 0.757 |

  - Memorising seen sites brings no extra benefit (0.798 against 0.799).
  - The novelty cost is concentrated on sites that no other kinase labels. The doc's "loses only 0.02 M1" is true on average, but on those sites it is about 0.05.
- **Harness detail.** The R2 LKO was evaluated with A24 on (`r2_lko.py:51`), although `evaluate()`'s docstring says LKO should pass `fold_normalise=False`. The effect is immaterial (0.7506 against 0.7504), but the R3 LKO (0.686) was evaluated with A24 off.

**R2-m3. Clean variant: the claims hold, with residual channels to disclose (claim 4).**
- **PMID matching.**
  - `exclude_lists.json` `clean_bench` lists 22 PMIDs: the 21 benchmarKIN Hernandez-Armenta PMIDs present in the compendium, plus Hijazi 2020 (31959955).
  - The 3 Hernandez PMIDs not in the list (23592840, 23612710, 24204825) do not occur in the compendium, so there is nothing to exclude for them.
  - Matching is by PMID only. Neither the benchmark metadata nor the compendium carries PXDs, so PXD matching cannot be verified.
- **F6b.** In `f6b_conditions_clean_bench.tsv`, 0 of the 483 used conditions carry an excluded PMID (209 rows are flagged `exclude_pmids`).
- **CPTAC columns.** F6c is dropped (the `F6c_dropped_clean_bench(CPTAC)` fallback appears in all 5 fold manifests). `f6a_cptac_top` and `f6a_cptac_med` are all-NaN, with zero split count and zero gain in every fold. No CPTAC-derived column is active.
- **Residual site-existence channels.** These are label-free and do not come from fitting on benchmark data:
  - U_ms membership (it includes compendium and CPTAC observation);
  - the F1 primed-soft weights (MS-evidence channels include the compendium and CPTAC);
  - R3's phospho channel. This matters only if R3 is ever fused into the clean library; the current clean library is R2-only.
- **Bookkeeping.** The U_obs library manifest `r2_final_cleanbench_uobs` records no F6c fallback (`fallbacks: {}`), because `r2_final.py` does not count Assembler fallbacks. The code path does drop F6c (`assemble.py:138`).

**R2-m4. REGISTRY undercounts the search (claim 7).**
- Distinct config names with split = cv:

  | Model | Configs | Breakdown |
  |---|---|---|
  | R2 | 48 | about 18 fold-0 tuning configs (stages a–d); 12 ablations and alarms; 10 per-seed runs plus 1 ensemble; knn × 5 seeds; F9; clean; full and G2 main; the fusion |
  | REGIME | 22 | label variants; partial-pool λ × 7; residual boosting λ × 7; per-kinase models |
  | R3 | 10 | |
  | R5 | 5 | |

  - The fusion weight sweeps count as **one** row each (`r2_fusion.py:68–71`), although they covered 9 + 6 grid points.
  - The M2 arm choice (calibrated vs raw) is not registered.
  - The 51 per-fold R2 rows were back-filled at 21:07; several are duplicates (for example `knn_s2`).
- **Optimism statement for FREEZE:**
  - About 90 configurations were evaluated on CV (selection-eligible: about 55).
  - The final pick (fusion, 0.786) leads the next-best arm (R2 ensemble, 0.776) by +0.010 [0.005, 0.016] (paired bootstrap, `compare__r2_fusion_r3__vs__r2_final_ens5.tsv`).
  - Upstream tuning happened on fold 0 only, within a noise band of about ±0.005. Its carry-over into pooled 5-fold M1 is at most about 0.002.
  - The fusion weight adds no measurable M1 optimism (§3).
  - **Expected CV→TEST optimism for M1 is therefore ≲ 0.005**, before any split-shift effects.

**R2-m5. Smaller code points:**
- `r2_fusion_lko.py:21` uses an inner merge without asserting equal lengths. It prints the row report but does not save it in the manifest (`rows: None`), unlike `r2_fusion.py:30`, which asserts.
- At TEST, `r2_fuse_test.py:32–33` z-scores over the 503 real kinases, whereas the CV z-scores included the 13 canaries. The effect is negligible.

## 3. Fusion weight selection: leakage and optimism (claim 1)

**Mechanism** (`scripts/models/r2_fusion.py`):
- Logits are z-scored per CV fold over all of that fold's rows. This is label-free and transductive, and it is mirrored at TEST over TEST rows.
- One scalar w is chosen by argmax of pooled-CV M1 over the grid, and w_site by argmax of M2.
- No labels enter the blend. The only fold information used is the fold's own score distribution.

**LOFO estimate.**
- I swept w = 0.0 to 1.0 in steps of 0.1 separately on each fold with `prepare('cv', fold=f)`.
- For each fold f, w was chosen on the mean of the other four folds and then applied to f.

| | w_M1 chosen per held-out fold | w_M2 chosen | Nested pooled M1 | Nested pooled M2 | Nested M3 |
|---|---|---|---|---|---|
| LOFO | 0.6 in all 5 folds | 0.3 / 0.5 / 0.4 / 0.3 / 0.4 | **0.7864** (= reported) | **0.2589** (reported 0.2626) | 0.3609 (= reported) |

Fold-mean M1 by w is flat near the optimum: 0.778 / 0.780 / 0.781 / 0.780 / 0.778 for w = 0.4–0.8. **The optimism from selecting w is therefore 0.000 for M1 and M3, and about 0.004 for M2.**

## 4. Verdicts

| Result | Verdict |
|---|---|
| R2 5-seed ensemble, M1 0.776 (M2 0.200, M3 0.312, M4 0.767) | **TRUSTED** |
| R2+R3 fusion, M1 0.786 / M3 0.361 / M4 0.774 (w = 0.6) | **TRUSTED**: reproduced exactly, no selection optimism, rows aligned |
| Fusion M2 0.263 (raw score_site, w = 0.5) | **TRUSTED WITH CAVEATS**: nested value 0.259, and the A20 wording needs a PI ruling (R2-m1) |
| LKO T2, R2 M1 0.751 (fusion 0.749) | **TRUSTED WITH CAVEATS**: no hold-out leak; qualify with the seen/unseen split (R2-m2) |
| Clean variant, CV M1 0.771 | **TRUSTED WITH CAVEATS**: exclusion list correct, no active CPTAC columns; disclose the residual site-existence channels (R2-m3) |
| Calibration in the large / ECE (A16) | **TRUSTED** (cross-fitted, monotone per kinase; slopes b > 0 for all 503 kinases) |
| Assignment tiers / coverage@0.8 = 1,674 (A17) | **NOT TRUSTED**. Apply the stratified isotonic fix: coverage@0.5 = 21,115 (ST only), coverage@0.8 ≈ 0 |
| Canary gate on the final candidate (pooled) | **PASS, WITH CAVEATS** (per-fold mean 0.561; two folds above 0.57; R2-M2) |

**Pre-FREEZE actions, in priority order:**
1. Replace the A17 tiers with residue-stratified, cross-fitted isotonic recalibration of p_top. Report coverage@0.5 and @0.8 from it, and state that no Y-site or ≥ 0.8 tier is reliable.
2. Get PI rulings on the canary gate definition (pooled only, or also the per-fold mean) and on the A20 M2 component. Log both in DESIGN §17.
3. Put the nested fusion numbers (M2 0.259) and the optimism statement in FREEZE. Make REGISTRY count each fusion grid point and the M2-arm choice.
4. Add the seen/unseen qualification to the T2 claims, and make the LKO evaluation harness setting consistent between R2 and R3.
