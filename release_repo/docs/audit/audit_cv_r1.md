> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Independent audit, CV round 1: R3 neural, simple baselines, R2 (in progress)

Auditor: independent results-auditor agent. I produced none of these results. Date: 2026-10-26, 16:45–17:30 CET.
Scope: (A) R3 (`docs/model_r3.md`, `src/ksp/models/r3/`, `scripts/models/r3_*.py`, `results/models/r3/FINAL_OOF.json` and its runs); (B) simple baselines (`results/baselines/simple/summary.md`, `src/ksp/baselines/simple.py`, `scripts/baselines/simple_*.py`); (C) R2 as of 16:55 (`docs/model_r2.md`, `src/ksp/models/r2/`, runs `r2_full_cv5*` and `r2_g2_main_cv5*`).
Constraints I kept to: read-only on every agent's files, ≤8 cores, no GPU, no git, `sealed/` never touched.

Method:
- Metrics were recomputed from the saved `cache/runs/*/predictions.parquet` with `ksp.eval.prepare('cv')` + `evaluate`.
- For R2 and R3, I concatenated the per-fold prediction files myself instead of using the pooled file. Separately, I checked that the R2 pooled parquet equals that concatenation (it does).
- Scripts: `docs/audit/audit_cv_r1_recompute.py` (integrity checks, per-fold evaluation, A24-default evaluation) and `docs/audit/audit_cv_r1_recompute2.py` (evaluation with `fold_normalise=False` and `True`).
- Bootstrap CIs were **not** recomputed; only point estimates were checked.

**Harness change during the audit.** `src/ksp/eval.py` was modified at 17:01 and `src/ksp/canary.py` at 16:55, and DESIGN §17 (A23–A25, post-freeze) was added at 16:56. A24 makes pooled-CV M1/M3/M4 use scores **rank-normalised within (fold, kinase)**, and this is now the default for `prepare('cv')`. Every reported number in scope predates A24. I therefore give both views:
- `fold_normalise=False`: the harness the authors used. Reproduces their numbers.
- A24: the current default.

## 1. Verified numbers (pooled 5-fold DEV OOF; K_cv = 203, K_any_cv = 339)

"Reported" values come from the run's `metrics.json` or summary; the doc values are the same numbers rounded. "Recomputed (raw)" means `fold_normalise=False`. Per-fold mean M1 is the mean over folds of macro M1 from `prepare('cv', fold=f)`.

| Result | Endpoint | Reported | Recomputed (raw) | Recomputed (A24, current default) | Per-fold mean M1 |
|---|---|---|---|---|---|
| **R3 final** `cpt_aux` | M1 | 0.7396 | **0.7396** ✓ | 0.7540 | 0.7494 (0.741/0.756/0.750/0.748/0.752) |
| | M2 | 0.2229 | 0.2229 ✓ | 0.2229 | |
| | M3 | 0.2820 | 0.2820 ✓ | 0.2936 | |
| | M4 | 0.7333 | 0.7333 ✓ | 0.7306 | |
| **R3 base** `w7_nce05` (canary v2) | M1 / M2 / M3 / M4 | 0.7208 / 0.1828 / 0.2334 / 0.7225 | identical ✓ | 0.7346 / 0.1828 / 0.2473 / 0.7222 | 0.7275 |
| **R3 LKO T2** (final) | M1 / M3 / M4 | 0.686 / 0.226 / 0.676 | 0.6862 / 0.2264 / 0.6764 ✓ | n/a (one model per kinase; A24 not applicable) | — |
| **R2 full** `r2_full_cv5` | M1 | 0.7631 | **0.7631** ✓ | 0.7729 | 0.7678 (0.754/0.773/0.767/0.790/0.756) |
| | M2 | 0.2075 | 0.2075 ✓ | 0.2075 | |
| | M3 | 0.3081 | 0.3081 ✓ | 0.3117 | |
| | M4 | **0.7822** | 0.7822 ✓ | **0.7697** | |
| **R2 G2 main** `r2_g2_main_cv5` | M1 / M2 / M3 / M4 | 0.7456 / 0.1457 / 0.2667 / **0.7516** | identical ✓ | 0.7537 / 0.1457 / 0.2718 / **0.7351** | 0.7508 |
| B0-prior | M1 / M2 / M3 / M4 | 0.500 / 0.0641 / 0.010 / 0.500 | identical ✓ | identical | 0.500 |
| B0-pop | M1 / M2 / M3 / M4 | 0.5521 / 0.0201 / 0.1399 / 0.5629 | identical ✓ | 0.5518 / 0.0201 / 0.1421 / 0.5617 | 0.556 |
| B0-miss | M1 / M2 / M3 / M4 | 0.5380 / 0.0402 / 0.0392 / 0.5397 | identical ✓ | 0.5390 / 0.0402 / 0.0528 / 0.5408 | 0.544 |
| B1 (i) published (official) | M1 / M2 / M3 / M4 | 0.7166 / 0.1570 / 0.1326 / 0.6901 | identical ✓ | 0.7166 / 0.1570 / 0.1323 / 0.6900 | 0.721 |
| B1 (ii) universe | M1 / M2 / M3 / M4 | 0.7166 / 0.1706 / 0.1326 / 0.6901 | identical ✓ | ≈ identical | 0.721 |
| B1 imputed-universe | M1 / M2 / M3 / M4 | 0.7219 / 0.1489 / 0.1324 / 0.7011 | identical ✓ | ≈ identical | 0.723 |
| B2 PSSM | M1 / M2 / M3 / M4 | 0.6605 / 0.1183 / 0.1083 / 0.6392 | identical ✓ | 0.6602 / 0.1183 / 0.1093 / 0.6383 | 0.673 |
| B7 LR | M1 / M2 / M3 / M4 | 0.7260 / 0.1959 / 0.1767 / 0.7086 | identical ✓ | 0.7258 / 0.1959 / 0.1788 / 0.7080 | 0.731 |

**Integrity checks** (`audit_cv_r1_recompute.py`):
- **OOF.** Every R3 and R2 fold file contains only sites on proteins of its own CV fold (`oof_ok = True` for all 15 fold files).
- **Coverage.** Every fold file covers 100% of the fold's U_ms sites × all 516 kinases (503 real + 13 canaries) of matching class. `n_missing_pairs = 0` out of 71,331,065 required pairs, and there are no missing site-centric (M2) scores.
- **Scores.** No non-finite scores.
- **R2 pooled file.** It equals the concatenation of the 5 fold files.
- **R3 LKO.** 0 required pairs are missing under `missing='drop'`.
- **Missing predictions** are treated as the worst score (`missing='worst'`, `eval.py:407` in the 17:01 version). Only the B1 own-matrix arms have gaps: 6,023,595 of 68.9 M required pairs (about 8.7%; 21 of 339 kinases are unscored), and they are counted and set to −inf. Covered-only numbers are reported separately in summary.md. No silent dropping was found.

**Canary v2 gate.** The gate is the negative-canary mean on E-hard ≤ 0.57; the positive canary must be ≥ 0.90 on E-hard and E-all.

| Run | Pooled neg E-hard (reported = recomputed raw) | Pooled neg E-hard, A24 | Per-fold neg E-hard (f0…f4) | Per-fold mean | Pos E-hard / E-all | Gate (pooled) |
|---|---|---|---|---|---|---|
| R3 final | 0.501 | 0.524 | 0.498 / 0.522 / 0.456 / 0.530 / **0.667** | 0.535 | 0.929 / 0.974 | PASS / PASS |
| R3 base | 0.486 | 0.499 | 0.522 / 0.546 / 0.407 / 0.499 / **0.585** | 0.512 | 0.830 / 0.836 | PASS / **FAIL** (as reported) |
| R2 full | 0.525 | 0.544 | **0.635** / 0.564 / 0.443 / 0.495 / **0.613** | 0.550 | 0.998 / 0.999 | PASS / PASS |
| R2 G2 main | 0.530 | 0.540 | 0.557 / **0.582** / 0.439 / **0.601** / **0.615** | 0.559 | 0.997 / 0.999 | PASS / PASS |
| B2 / B2+ph / B7 (canary runs; reported only, not recomputed) | 0.541 / 0.547 / 0.530 | — | — | — | 0.953 / 0.958 / 0.999 (E-hard) | PASS / PASS |

**Per-kinase sanity** (pooled raw E-hard AUROC, as reported; my recomputation is identical; n = DEV evaluation positives):

| Kinase (n) | R3 final | R2 full | B1 (i) | B7 | Plausible? |
|---|---|---|---|---|---|
| PRKACA (365) | 0.855 | 0.865 | **0.875** | 0.847 | Yes. R-x-x-S/T is strong; Cantley is best here. |
| CDK1 (515) | 0.888 | 0.906 | 0.865 | 0.866 | Yes (S/T-P). |
| CSNK2A1 (365) | 0.885 | 0.885 | **0.909** | 0.883 | Yes. The acidic motif is Cantley-dominated. |
| MAPK1 (305) | 0.868 | 0.878 | 0.860 | 0.866 | Yes. |
| ATM (186) | 0.957 | 0.954 | 0.922 | 0.949 | Yes. S/T-Q is very distinctive; E-hard negatives are rarely SQ. |
| SRC (334) | 0.589 | 0.597 | 0.574 | 0.579 | Yes. Its E-hard negatives are mostly other TK substrates; all Y-class motif models are near 0.53–0.60. |
| GSK3B (230) | 0.751 | 0.770 | 0.578 | 0.701 | Yes. Priming is learned (+4 pS); unprimed Cantley fails, as expected. |

None of these looks too good. The pattern matches the biology: motif-defined kinases are high, TKs are low, and the learned models gain most for priming-dependent GSK3B. **One subgroup does look too good:** R2's raw M4 in the 1–4-label bin (§2, M-1).

## 2. Issues by severity

### Critical
None. I found no label leakage into eval-fold features, no use of eval-fold information in training, early stopping or pretraining, no TEST access, and no silent dropping of predictions. Details of the leakage review are in §3.

### Major
**M-1. Pooled-OOF offset artefact inflates R2's M4, especially for understudied kinases. The reported R2 M4 and every "1–4 label" claim are not trustworthy as stated.**
- Within-fold rank normalisation (A24) removes fold-varying per-kinase offsets. The changes in M4 are:
  - R2 full: 0.782 → **0.770**;
  - R2 G2 main: 0.752 → **0.735**;
  - R3 final: 0.733 → 0.731.
- By training bin (K_any; raw → A24), the 1–4 bin moves sharply:

  | Model | 1–4 (n = 79) | 5–19 (n = 137) | 20–99 (n = 92) | ≥100 (n = 31) |
  |---|---|---|---|---|
  | R2 full | **0.821 → 0.764** | 0.762 → 0.760 | 0.770 → 0.777 | 0.808 → 0.809 |
  | R2 G2 main | **0.766 → 0.689** | 0.738 → 0.739 | 0.747 → 0.750 | 0.788 → 0.789 |
  | R3 final | **0.696 → 0.647** | 0.743 → 0.755 | 0.735 → 0.748 | 0.784 → 0.786 |

- A raw 1–4-bin AUROC (0.821) that beats the ≥100 bin is the "too good" signal. The cause is that a kinase's few positives sit in the fold where it has 0–3 training labels. R2 and R3 score such low- or zero-count kinases at a different overall level, and pooling across folds turns that level into apparent discrimination. This is the same artefact as B0-prior's 0.40, but here it runs in the inflating direction.
- The artefact deflates M1 (R2 +0.010 and R3 +0.014 after A24) and the pooled negative-canary statistic.
- Locations: `docs/model_r2.md:214,226` (M4 0.752 / 0.782); `docs/model_r3.md:16,31` (M4 and the "1–4 bin" statement); `results/models/r3/paired_cptaux_vs_base_v2.tsv` (raw M4 Δ).
- **Fix:** restate M1, M3 and M4 under A24 and add the per-fold mean, in both docs. Make no understudied-kinase claim from raw pooled M4.

**M-2. R2 cross-fitted stores for fold 0 (all four families) and fold 3 (F2a, F3-nb) were built from an intermediate canary table, not canary v2.**
- Evidence (`cache/r2/xfit/cv/fold*/*.json`, file mtimes):
  - Fold-0 f2a/f3nb/f6a/f6a_clean, and fold-3 f2a/f3nb, record 331 / 329 training canary pairs.
  - The final `canary_dev.parquet` (written 15:34:35) gives 348 / 350 for those folds, and v1 gives 288 / 289.
  - Fold-3 f6a (built 15:49) has 350, and folds 1, 2 and 4 match v2 exactly (332, 337, 345).
  - The fold-0 builds started at about 15:32:48 (mtime minus runtime), before the v2 write.
- Real-kinase label counts are identical (fold 0: 9,433; fold 3: 9,436), so real-kinase features are unaffected. The only possible exception is F2a's shared α/β tuning, which includes canaries; the effect is tiny.
- No leak follows from this: the fitted canary labels lie on outer-training proteins only. But the canary gate for R2 folds 0 and 3 is **not end-to-end**, because the canaries' F2a and F3-nb features are fitted on labels other than the ones they are trained and evaluated on.
- `docs/model_r2.md:108` ("all xfit stores were rebuilt with v2") is inaccurate.
- **Fix:** rebuild `cv/fold0/*` and `cv/fold3/{f2a,f3nb}` and re-run the canary check, before R2 is frozen.

**M-3. The canary gate passes only in its pooled form, and the pooled statistic is biased toward passing.**
- Per fold, the negative-canary E-hard mean exceeds 0.57 in:
  - 3 of 5 folds for R2 G2 (mean 0.559);
  - 2 of 5 for R2 full (mean 0.550);
  - 1 of 5 for R3 (fold 4: 0.667).
- The pooled value is lower because of the M-1 artefact. Under A24 it rises to 0.540 / 0.544 / 0.524, which still passes.
- Every model, including the label-light B7 (0.530) and R2's label-free F0+F7 control (0.539 on fold 0), sits above 0.5. That fits a sampler-population difference (protein-first canary sampling against site-weighted E-hard negatives), as the R2 agent diagnosed in `docs/model_r2.md` §5.2, not leakage. R2 is still within 0.01–0.03 of the threshold.
- **Recommendation (PI):** gate on the A24 pooled value *and* report the per-fold mean. Consider a v3 negative canary sampled site-uniformly from the E-hard-eligible pool, so the null is exactly 0.5.

**M-4. Post-freeze metric change (A24); every number must state its harness version.**
- A24 is principled and is applied to all methods. It changes almost nothing for fixed-score baselines (B1 ΔM1 = 0.000) but adds +0.010 to +0.014 M1 for the learned models. Each learned-versus-Cantley margin therefore widens under A24: R2 full − B1 goes from +0.046 to +0.056.
- The brief forbids changing metrics to flatter results. This change is defensible (it is a documented artefact and TEST is unaffected), but it needs explicit PI sign-off, and the docs must give both values side by side.
- Every number in `model_r2.md`, `model_r3.md`, `summary.md` and REGISTRY is pre-A24. Pre-A24 and A24 numbers must not be compared with each other.

### Minor
- **m-1. M2 treats unscored competitor kinases as −inf** (`src/ksp/eval.py:544`). For a method with partial kinase coverage, missing competitors can never outrank the true kinase, so its M2 is inflated.
  - Evidence: on the B1-matched universe, R2 G2's own M2 rises from 0.146 (full) to 0.166 (`results/models/r2/compare__r2_g2_main_cv5__vs__simple_b1_universe.tsv`). B1 (i) and (ii) M2 (0.157 / 0.171) are therefore optimistic by roughly 0.02 relative to full-coverage methods; B1-imputed (full coverage) gets 0.149.
  - Head-to-heads already use `matched_universe` + `drop`, which is fair. The headline M2 tables are not.
- **m-2. R3 selection is within noise.** cpt_aux beats cpt by only 0.0017 on fold 0 (`docs/model_r3.md` §2), and the pre-declared rule was applied correctly. The +0.019 gain over base is confirmed out of selection: on folds 1–4 the per-fold differences are +0.011, +0.029, +0.010 and +0.014.
  - The round-2 configurations were designed after the round-1 pooled 5-fold results were seen. That is mild optimism over 10 eval-fold configurations, all logged correctly in REGISTRY.
  - Folds 1–4 of cpt_aux were launched (16:00:57) before the two 3B fold-0 results existed (16:03, 16:10). This pre-commits the choice; it is not an optimism problem.
- **m-3. R2 selection bookkeeping.**
  - Stage (a) was selected under canary v1 and the pre-15:31 harness. The chosen config re-run under v2 scores 0.7537, below several un-re-run v1 alternatives; all differences are within the ±0.005–0.01 noise floor the author states.
  - The fold-1–4 runs of `r2_full`/`r2_g2_main` are not REGISTRY rows; only the pooled rows are.
  - `n_configs_tried_model` counts rows, not distinct configs. There are about 20 fold-0 configs, plus 1 control and 7 pooled runs (main, full, ablations, alarm).
  - `r2_full` reuses the fold-0 run `r2_a_l127_lr03_v2`. Its config is identical to f1–f4 (verified from the manifests).
- **m-4. LKO T2 is not kinase-naive.**
  - Held-out kinases' label-free Cantley percentiles enter pretraining (disclosed, `docs/model_r3.md` §5.3).
  - Held-out kinases' E-hard negatives are training *positives* of other kinases, because the proteins are seen. This probably deflates T2, so it is not leakage.
  - The pooled LKO M2 (0.054) is meaningless, because only held-out kinases compete, and should not be used. It is not reported in the doc.
- **m-5. Naming.** "M4 (K_any, 339)" in CV is the CV analogue over kinases with ≥1 DEV positive (`eval.py` k_any for cv). It is not the K_any = 256 TEST set. Label it "M4-CV".
- **m-6. Stale REGISTRY values for the simple baselines.** REGISTRY holds the pre-15:31 harness values (for example b1 0.71875 against the re-evaluated 0.7166). summary.md documents this.
- **m-7. Canary results for B2/B7 were not recomputed in this audit.** Those are separate canary runs; their code paths were reviewed (§3).

## 3. Leakage and discipline review (checks 2–3)

**Grep.** No module in `src/ksp/features/` imports `ksp.labels` or reads label tables. Label reads in model and baseline code are only these:
- `r2/sample.py:54,68`: training-fold rows, filtered with `acc.isin(train_prot)`;
- `r3/data.py:210–223`: filtered to training proteins (and training kinases for LKO);
- `baselines/simple.py:140–145`: masks filtered to training proteins.

No `sealed/` access exists in scope. `tests/test_no_test_access.py` and `tests/test_crossfit.py` pass (8/8).

**Cross-fitting.**
- R2's F2a, F3-label and F6a go through `ksp.crossfit.crossfit` (`r2/xfit.py:278`; the file was edited by its owner at 17:21, after my review of the 16:5x version), with eval rows set to the A1 inner-fold fit.
- `DenseExpand` expands only over kinases at fixed query sites, and the driver's protein/cluster disjointness check covers every block.
- Outer training labels come from `crossfit.outer_split` (eval-fold proteins excluded).
- B7's B2 feature is cross-fitted (`simple.py:478`). B2-as-model is fit on fold-training labels only.
- No all-DEV labels feed features for eval-fold rows. Eval masks use all DEV labels by design (`eval.py` prepare), and that affects only metric definitions.

**R3.**
- Targets and masks come from training-protein labels only (`r3/data.py:209–215`).
- Early stopping uses only the inner validation split of training clusters (`r3/train.py:170, 307–318`).
- Distillation targets are the label-free `motif/cantley_*_pct_universe` (`r3/train.py:45–57`; universe percentile over U_ms, `docs/features_motif.md:59`).
- Pretraining samples training-fold rows only (`pre_rows=tr_idx`, `r3/train.py:201, 270`).
- The phospho channel uses U_ms membership of flanking residues, with the centre forced to 0 (`r3/train.py:99`).
- Kinase descriptors (ESM, Cantley profile, Manning group/family) are label-free, and families with no trained kinase are zeroed.

**R2 feature importance.** No single feature dominates. The top gains are f7_family 0.10, f1_unprimed_pct 0.06–0.07, f1_score 0.035 and f3_rwr_dc_pct 0.035; every other feature is ≤ 0.026. The COMPARTMENTS-knowledge alarm adds nothing (M1 0.761 against 0.763), which argues against curated-knowledge leakage.

## 4. Verdicts

| Result | Verdict | Conditions |
|---|---|---|
| R3 final (`cpt_aux`), CV M1/M2/M3 | **TRUSTED WITH CAVEATS** | Exactly reproduced and no leakage found. Restate under A24 (M1 0.754, M3 0.294; per-fold mean M1 0.749). Fold-4 canary is 0.667 (M-3). |
| R3 cpt_aux vs base gain (+0.019 M1, +0.040 M2, +0.049 M3) | **TRUSTED** | Holds out of the selection fold and under A24 (0.754 vs 0.735). |
| R3 M4 and the "1–4 label bin" statements | **NOT TRUSTED as stated** | Pooling artefact (M-1). Under A24 the 1–4 bin is 0.647. |
| R3 LKO T2 (M1 0.686) | **TRUSTED WITH CAVEATS** | Not kinase-naive (m-4). Do not use LKO M2. |
| Simple baselines B0/B1/B2/B7 (all endpoints) | **TRUSTED WITH CAVEATS** | Exactly reproduced and essentially A24-invariant. B1 M2 is optimistic by about 0.02 for partial coverage (m-1). The B1 arm choice is a tie on M1 (A23 resolves it). |
| R2 G2 gate (beats B1 on M1, +0.023 [0.009, 0.034]) | **TRUSTED WITH CAVEATS** | The point M1 is reproduced and the gate holds more strongly under A24. Per-fold canary failures (M-3) and stale fold-0 canary stores (M-2) need closing before freeze. |
| R2 full, M1 0.763 / M3 0.308 / M2 0.208 | **TRUSTED WITH CAVEATS** | Reproduced. A24: M1 0.773, M3 0.312. Same canary caveats (M-2, M-3). |
| R2 full / G2 M4 (0.782 / 0.752) | **NOT TRUSTED as stated** | Inflated by the pooling artefact. A24: 0.770 / 0.735, and the 1–4 bin falls from 0.821 to 0.764. |

**Requested actions before FREEZE:**
1. Rebuild R2 xfit fold 0 (all four families) and fold 3 (f2a, f3nb), then re-run the R2 canary check.
2. Re-issue every CV table under A24, with the per-fold mean and the harness version (eval.py mtime 17:01 or later).
3. The PI signs off A24 and decides whether the canary gate also requires the per-fold mean.
4. Report M2 for partial-coverage methods with the m-1 caveat, or fix `eval.py:544`, for example by restricting candidates to kinases the method scores.
