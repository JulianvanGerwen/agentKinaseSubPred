# Simple baselines on the 5-fold DEV CV (B0 / B1 / B2 / B7)

Generated 2026-10-26 16:43 by `scripts/baselines/simple_summary.py` from `results/baselines/simple/<method>.json` (each = one `ksp.manifest` run; metrics only from `ksp.eval`).

Pooled out-of-fold predictions over **all** evaluation candidates (every kinase with >=1 DEV positive x every DEV U_ms site of its class, plus all real kinases x labelled sites for M2). Missing pairs = 'worst' (ties at the bottom). 95% CI: cluster bootstrap over substrate clusters, B=1000 (`ksp.eval.bootstrap`). M1 = M5 on CV (K_cv, E-hard); M4 = CV analogue over K_any_cv (>=1 DEV positive). M2 = filtered site-centric Hit@5 (A18), M3 = recall in top 1% of E-all.

## Endpoints (CV, pooled OOF; point [95% CI])

| method | description | M1 (=M5) AUROC E-hard | M2 Hit@5 | M3 top-1% recall | M4-CV AUROC (K_any) | MRR | AUROC E-all | fit+predict s | total s |
|---|---|---|---|---|---|---|---|---|---|
| `b0_prior` | B0-prior (kinase training count) | 0.500 [0.500, 0.500] | 0.064 [0.062, 0.069] | 0.010 [0.010, 0.010] | 0.500 [0.500, 0.500] | 0.056 | 0.500 | 2 | 219 |
| `b0_pop` | B0-pop (site popularity; diagnostic) | 0.552 [0.539, 0.573] | 0.020 [0.020, 0.020] | 0.140 [0.112, 0.170] | 0.563 [0.545, 0.585] | 0.008 | 0.782 | 0 | 281 |
| `b0_miss` | B0-miss (LightGBM on missingness only) | 0.538 [0.523, 0.560] | 0.040 [0.029, 0.049] | 0.039 [0.028, 0.050] | 0.540 [0.522, 0.560] | 0.034 | 0.694 | 138 | 292 |
| `b1_published` | B1 (i) Cantley as-published pct | 0.717 [0.701, 0.731] | 0.157 [0.141, 0.174] | 0.133 [0.117, 0.151] | 0.690 [0.681, 0.713] | 0.112 | 0.757 | 1 | 439 |
| `b1_universe` | B1 (ii) Cantley universe-recalibrated pct | 0.717 [0.701, 0.731] | 0.171 [0.154, 0.189] | 0.133 [0.117, 0.151] | 0.690 [0.681, 0.713] | 0.124 | 0.757 | 1 | 415 |
| `b1_primed_hard` | B1 (iii) Cantley priming-aware hard (universe pct) | 0.714 [0.698, 0.729] | 0.157 [0.141, 0.175] | 0.116 [0.103, 0.132] | 0.688 [0.680, 0.711] | 0.115 | 0.750 | 7 | 378 |
| `b1_primed_soft` | B1 (iii) Cantley priming-aware soft (universe pct) | 0.711 [0.696, 0.727] | 0.160 [0.146, 0.179] | 0.124 [0.109, 0.142] | 0.686 [0.677, 0.710] | 0.120 | 0.749 | 2 | 399 |
| `b1_primed_hard_published` | B1 suppl.: primed hard, published pct | 0.713 [0.698, 0.729] | 0.122 [0.109, 0.137] | 0.116 [0.103, 0.132] | 0.688 [0.680, 0.711] | 0.089 | 0.750 | 7 | 412 |
| `b1_primed_soft_published` | B1 suppl.: primed soft, published pct | 0.711 [0.696, 0.727] | 0.137 [0.122, 0.154] | 0.123 [0.109, 0.141] | 0.686 [0.677, 0.710] | 0.100 | 0.749 | 8 | 379 |
| `b1_imputed_universe` | B1 suppl.: universe pct + nearest-atlas imputed | 0.722 [0.705, 0.737] | 0.149 [0.133, 0.168] | 0.132 [0.117, 0.150] | 0.701 [0.689, 0.722] | 0.106 | 0.763 | 5 | 274 |
| `b2_pssm` | B2 learned PSSM (F2a, family prior) | 0.660 [0.646, 0.679] | 0.118 [0.108, 0.131] | 0.108 [0.093, 0.125] | 0.639 [0.630, 0.663] | 0.094 | 0.735 | 146 | 324 |
| `b2_pssm_ph` | B2 learned PSSM + phospho channel | 0.654 [0.642, 0.672] | 0.128 [0.115, 0.140] | 0.129 [0.105, 0.151] | 0.633 [0.625, 0.658] | 0.098 | 0.751 | 142 | 289 |
| `b7_lr` | B7 simple strong (LR) | 0.726 [0.710, 0.742] | 0.196 [0.179, 0.213] | 0.177 [0.155, 0.199] | 0.709 [0.695, 0.731] | 0.146 | 0.799 | 925 | 290 |

## Coverage

Covered-only = `missing='drop'` restricted to kinases with >=1 scored pair (not a matched-universe comparison; head-to-heads must use `ksp.eval.matched_universe`).

| method | kinases scored / with DEV pos | K_cv scored | required pairs scored | covered-only M1 | covered-only M2 | covered-only M3 |
|---|---|---|---|---|---|---|
| `b0_prior` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b0_pop` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b0_miss` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b1_published` | 318/339 | 199/203 | 0.913 | 0.723 (n=199) | 0.161 | 0.135 |
| `b1_universe` | 318/339 | 199/203 | 0.913 | 0.723 (n=199) | 0.174 | 0.135 |
| `b1_primed_hard` | 318/339 | 199/203 | 0.913 | 0.716 (n=199) | 0.160 | 0.118 |
| `b1_primed_soft` | 318/339 | 199/203 | 0.913 | 0.715 (n=199) | 0.163 | 0.126 |
| `b1_primed_hard_published` | 318/339 | 199/203 | 0.913 | 0.716 (n=199) | 0.125 | 0.118 |
| `b1_primed_soft_published` | 318/339 | 199/203 | 0.913 | 0.715 (n=199) | 0.140 | 0.126 |
| `b1_imputed_universe` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b2_pssm` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b2_pssm_ph` | 339/339 | 203/203 | 1.000 |  |  |  |
| `b7_lr` | 339/339 | 203/203 | 1.000 |  |  |  |

## Official B1 = Cantley best arm (A19; chosen on CV M1)

Rule: max pooled-CV M1 (K_cv, E-hard, missing=worst) among DESIGN arms (i),(ii),(iii); tie-break M2. **Best: `b1_published`** (`motif/cantley_unprimed_pct_published`; run `simple_b1_published_reeval__20261026-155352__s0`).

| method | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| b1_published | 0.7166 | 0.157 | 0.1326 | 0.6901 |
| b1_universe | 0.7166 | 0.1706 | 0.1326 | 0.6901 |
| b1_primed_hard | 0.7135 | 0.157 | 0.1163 | 0.688 |
| b1_primed_soft | 0.7111 | 0.1602 | 0.124 | 0.6862 |

Paired cluster bootstrap (B=1000), best minus other:

| comparison | endpoint | diff | lo95 | hi95 | p_boot |
|---|---|---|---|---|---|
| b1_published - b1_universe | M1 | 0.0 | -0.0 | 0.0 | 0.6633 |
| b1_published - b1_universe | M2 | -0.0135 | -0.0201 | -0.007 | 0.002 |
| b1_published - b1_universe | M3 | 0.0 | -0.0002 | 0.0002 | 0.7433 |
| b1_published - b1_universe | M4 | -0.0 | -0.0 | 0.0 | 0.9331 |
| b1_published - b1_primed_hard | M1 | 0.0031 | -0.0024 | 0.0082 | 0.2897 |
| b1_published - b1_primed_hard | M2 | 0.0 | -0.0127 | 0.0111 | 0.9211 |
| b1_published - b1_primed_hard | M3 | 0.0164 | 0.0066 | 0.0275 | 0.006 |
| b1_published - b1_primed_hard | M4 | 0.002 | -0.004 | 0.0068 | 0.5994 |
| b1_published - b1_primed_soft | M1 | 0.0055 | 0.0011 | 0.0092 | 0.016 |
| b1_published - b1_primed_soft | M2 | -0.0031 | -0.0146 | 0.0061 | 0.4316 |
| b1_published - b1_primed_soft | M3 | 0.0087 | 0.0013 | 0.0177 | 0.024 |
| b1_published - b1_primed_soft | M4 | 0.0039 | -0.001 | 0.0075 | 0.1578 |

## Canaries (trained with `load_dev(with_canaries=True)`; `ksp.eval.canary_check`)

| file | what | run_id | neg mean AUROC E-hard (gate <=0.57) | neg mean AUROC E-all (not gated) | pos AUROC E-hard | pos AUROC E-all | neg pass | pos pass (>=0.90) |
|---|---|---|---|---|---|---|---|---|
| `b2_pssm_canary` | B2, canary v2 (official: canary-free tuning, isolated canary priors) | simple_b2_pssm_canary__20261026-154412__s0 | 0.541 | 0.6 | 0.953 | 0.981 | True | True |
| `b2_pssm_ph_canary` | B2+phospho, canary v2 (official) | simple_b2_pssm_ph_canary__20261026-154453__s0 | 0.547 | 0.62 | 0.958 | 0.98 | True | True |
| `b7_lr_canary` | B7, canary v2 (official) | simple_b7_lr_canary__20261026-154540__s0 | 0.53 | 0.583 | 0.999 | 1.0 | True | True |
| `b2_pssm_canary_pooledprior` | diag: B2 v2, canaries share group 'CANARY' prior (tuning canary-free) | simple_b2_pssm_canary__20261026-154216__s0 | 0.537 | 0.654 | 0.849 | 0.909 | True | False |
| `b2_pssm_ph_canary_pooledprior` | diag: B2+ph v2, shared canary prior | simple_b2_pssm_ph_canary__20261026-154256__s0 | 0.541 | 0.667 | 0.776 | 0.868 | True | False |
| `b2_pssm_canary_tunedwithcanaries` | diag: B2 v2, shared prior + canaries in alpha/beta tuning | simple_b2_pssm_canary__20261026-153933__s0 | 0.538 | 0.654 | 0.832 | 0.897 | True | False |
| `b2_pssm_ph_canary_tunedwithcanaries` | diag: B2+ph v2, shared prior + canaries in tuning | simple_b2_pssm_ph_canary__20261026-154014__s0 | 0.542 | 0.669 | 0.789 | 0.874 | True | False |
| `b2_pssm_canary_v1` | superseded: B2, canary v1 | simple_b2_pssm_canary__20261026-150138__s0 | 0.498 | 0.53 | 0.932 | 0.958 | True | True |
| `b7_lr_canary_v1` | superseded: B7, canary v1 | simple_b7_lr_canary__20261026-153214__s0 | 0.476 | 0.587 | 0.997 | 0.999 | True | True |

## Strata: macro AUROC (E-hard) over K_any kinases

| method | pikk=False (n=334) | pikk=True (n=5) | priming=False (n=327) | priming=True (n=12) | residue_class=ST (n=265) | residue_class=STY (n=9) | residue_class=Y (n=65) | train_bin=1-4 (n=79) | train_bin=20-99 (n=92) | train_bin=5-19 (n=137) | train_bin=>=100 (n=31) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| b0_prior | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| b0_pop | 0.563 | 0.532 | 0.566 | 0.468 | 0.549 | 0.812 | 0.583 | 0.568 | 0.537 | 0.582 | 0.542 |
| b0_miss | 0.54 | 0.495 | 0.541 | 0.504 | 0.522 | 0.779 | 0.577 | 0.562 | 0.543 | 0.528 | 0.526 |
| b1_published | 0.687 | 0.866 | 0.689 | 0.706 | 0.732 | 0.573 | 0.534 | 0.624 | 0.725 | 0.688 | 0.766 |
| b1_universe | 0.687 | 0.866 | 0.689 | 0.706 | 0.732 | 0.573 | 0.534 | 0.624 | 0.725 | 0.688 | 0.766 |
| b1_primed_hard | 0.685 | 0.857 | 0.684 | 0.794 | 0.73 | 0.556 | 0.537 | 0.628 | 0.723 | 0.681 | 0.769 |
| b1_primed_soft | 0.684 | 0.861 | 0.683 | 0.784 | 0.728 | 0.555 | 0.534 | 0.625 | 0.719 | 0.681 | 0.77 |
| b1_primed_hard_published | 0.685 | 0.857 | 0.684 | 0.792 | 0.729 | 0.556 | 0.537 | 0.628 | 0.723 | 0.681 | 0.768 |
| b1_primed_soft_published | 0.684 | 0.861 | 0.683 | 0.783 | 0.728 | 0.555 | 0.534 | 0.625 | 0.719 | 0.681 | 0.77 |
| b1_imputed_universe | 0.699 | 0.866 | 0.701 | 0.706 | 0.742 | 0.64 | 0.542 | 0.652 | 0.726 | 0.698 | 0.766 |
| b2_pssm | 0.636 | 0.843 | 0.635 | 0.762 | 0.66 | 0.82 | 0.53 | 0.578 | 0.67 | 0.629 | 0.748 |
| b2_pssm_ph | 0.63 | 0.833 | 0.629 | 0.752 | 0.648 | 0.843 | 0.544 | 0.576 | 0.663 | 0.624 | 0.737 |
| b7_lr | 0.706 | 0.899 | 0.707 | 0.761 | 0.744 | 0.851 | 0.543 | 0.654 | 0.727 | 0.712 | 0.776 |

Full strata (incl. manning group, M2/M3 per stratum, pair-level evidence/temporal/external strata) are in each `<method>.json` (`strata`, `strata_K_cv`, `pair_strata`).

## Notes, assumptions and caveats

- **Candidate set.** `ksp.baselines.simple.candidate_pairs`: 339 real kinases with >=1 evaluable DEV positive x every
  DEV U_ms site of their class (= positives + E-all), plus all 503 real kinases x every labelled positive site of a
  matching class (M2). About 71M pairs. Each site is scored by the model of its protein's CV fold. Evaluation is
  `ksp.eval.prepare('cv')` (default: canary labels are in the mask table) with `evaluate(..., kinases=<real kinases with positives>, missing='worst')`.
- **Pooled-OOF offset artefact (important for every per-kinase-offset model).** The first B0-prior run
  (`simple_b0_prior__20261026-150125__s0`, superseded) used the *fold-training* label count as the within-kinase score.
  It scored M1 = 0.400 and M4 = 0.307 instead of 0.5. The cause: a kinase's held-out positives sit in the folds where
  its training count is lowest. Pooled OOF AUROC therefore penalises (or rewards) any score whose per-kinase level
  depends on the fold's training labels; examples are `f2a_pssm_log_ntrain`, per-kinase intercepts and
  count-type features. Fix used here: B0-prior `score` = constant DEV count (M1 = 0.5 exactly), and
  `score_site` = fold-training count (OOF), which `ksp.eval` uses for M2 only. Likewise, B7 `score` has no per-kinase
  intercept (M1/M3/M4), while `score_site` adds the intercept (M2), as A20 prescribes. **Recommendation to PI / R2:**
  report per-fold M1 (`prepare('cv', fold=f)`) next to pooled M1, or keep fold-dependent per-kinase offsets out of the
  within-kinase score.
- **B0-pop** = PSP LT_LIT + MS_LIT + MS_CST + compendium `comp_n_conditions` (NA -> 0). It is constant across kinases,
  so M2 is ties (expected-rank tie handling in `ksp.eval`). Diagnostic only.
- **B0-miss** uses 57 binary indicators:
  - pair level: NaN or not for F6b ridge_z / contrast and F6c prot / aloop (kinase x site), and for PPI RWR, the 6
    localisation cosines and the 2 co-expression measures (kinase x protein);
  - kinase level and protein/site level: "has any value" for the same stores, PPI source membership, structure and
    topology NaN, and whether the site is in the compendium or CPTAC.

  PSP and Ochoa site presence are excluded because they measure popularity (B0-pop), not feature missingness. F6a
  (label-derived, GPU) missingness is not included. LightGBM uses lr 0.05, 31 leaves, min 100 per leaf, and early
  stopping on a 20% cluster-grouped inner split. It was not tuned (1 configuration). Kinase identity is not a feature.
- **Training rows (B0-miss, B7).** Per kinase: fold-training L2 positives (U_ms, in class), plus 50 negatives per
  positive (50% E-hard-type, filled per A15). Negatives come from training-protein U_ms sites, with masks built from
  training-protein labels only. In about 140 of ~1,700 kinase-folds the hard pool was short and was filled from E-all.
- **B1.** Arms (i) and (ii) are monotone in the same unprimed score for each kinase. Their M1/M3/M4 are therefore
  identical up to ties: the published percentile is rounded to 2 dp against an 82k-site reference. They differ only in
  M2 (cross-kinase comparability). Coverage is the same for all own-matrix arms. The nearest-atlas imputed variants are
  reported separately and are not candidates for "best arm".
- **B2** = `motif_learned.fit_pssm` (family prior; alpha and beta tuned internally on fold-training pairs grouped by
  substrate cluster), fit on all fold-training L2 labels. Score = within-U_ms percentile. As a model, B2 is not
  inner-cross-fitted; inside B7 it is, through `ksp.crossfit.crossfit` (A1 eval mode).
- **B7** = L2 logistic regression (C = 1, one configuration). Inputs:
  - standardised features: logit(best-arm Cantley percentile), logit(B2 percentile), PPI direct edge in any screen,
    localisation combined cosine;
  - missing indicators for Cantley, B2 and localisation (NaN filled with neutral values);
  - one-hot kinase intercepts (L2-shrunk, i.e. partial pooling toward the global intercept).

  Sampled negatives carry inverse-sampling-probability weights (normalised to mean 1), so the intercepts reflect each
  kinase's prevalence in the training universe, as M2 needs. Canary runs add the planted CANARY_POS PPI edges
  (`ksp.canary.canary_ppi_partners`) as ordinary PPI edges.
- **Canary runs** train with `load_dev(with_canaries=True)` (B2 fit with canary kinase rows, family/group priors of
  their own) and score only canary candidates. Main runs are trained without canaries, so real-kinase numbers are
  unaffected.
- Runtime = wall-clock on <=8 cores while the node was shared (load ~20-30), so it is indicative only. "fit+predict"
  covers all 5 folds over ~71M candidates. "total" adds `ksp.eval` evaluation, bootstrap and I/O.
- **Canary v2 (installed 15:34) and re-evaluation.** The harness also changed after the first runs (`src/ksp/eval.py` modified 15:31). All main-run metrics in this file were *re-evaluated* from the
  saved OOF predictions with `prepare('cv')` under canary v2 (runs `simple_<m>_reeval__*`, `--from-run`). The
  first-pass numbers, from the earlier harness, were M1 about 0.002 higher for B0/B1/B2. The B7 run already used the
  15:31 harness and did not change on re-evaluation, so the shift comes from the eval.py change, not from canary v2.
  Main models were not refitted, because their training never used canaries.
- **Canary findings (B2) - two pipeline hazards found with canary v2:**
  1. `ksp.canary` gives every canary Manning group `CANARY`. The learned-PSSM hierarchical prior therefore pools the
     12 negative canaries' windows (v2: real labelled sites) into every canary's prior. With alpha = 128-256, that
     generic "labelled-site" prior dominates the positive canary's ~88 training sites, and its E-hard AUROC falls to
     0.83-0.85 (fail). E-hard negatives are labelled sites, so they look labelled-site-like. Fix used for the official
     canary runs: canary rows get an empty group/family (no pooling; prior = background), via
     `ksp.baselines.simple.pssm_kinase_table`. Positive canary is then 0.953.
  2. Shared hyperparameters: `fit_pssm` tunes alpha/beta by macro log-likelihood over all kinases, including canaries.
     Official canary runs tune on real kinases only (`_fit_pssm_canary_safe`) and then fit the canaries with those
     values fixed. Here it changed alpha/beta on 1 of 5 folds (fold 1: 256/256 with canaries vs 128/1024 without), and the effect on the positive canary was +0.017 E-hard; the rule
     matters.

  **Recommendation (motif agent / R2 / R3):** never let canary rows share group/family priors or enter shared
  hyperparameter tuning. Ideally `ksp.canary` should give canaries unique `manning_group` values. The diagnostic
  variants are listed in the canary table.
- Negative canaries: E-all AUROC is 0.58-0.67. This is expected under v2 (generic labelled-site-likeness) and is not
  gated. E-hard is 0.53-0.55, which passes the 0.57 gate but sits above 0.5.
