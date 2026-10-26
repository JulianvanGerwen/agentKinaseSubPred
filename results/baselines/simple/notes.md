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
- **Official B1 = `b1_published`, but it is a near-tie.** It beats (ii) on M1 by 9e-7. That gap is a tie-handling
  artefact of the 2-dp published percentile; the paired CI on M1 is [-5e-6, 8e-6]. The pre-declared rule (max M1)
  is applied as written. Arm (ii), the universe percentile, is significantly better on M2 (+0.0135, CI
  [0.007, 0.020]). **The PI should decide** whether the M2 superiority test (A19) uses arm (i) as declared, or arm (ii)
  as the stronger and more conservative comparator. B7 uses arm (i).
- **Tyrosine kinases are weak for every motif baseline.** Over 65 Y-class kinases, E-hard AUROC is about 0.53-0.54
  (Cantley, B2, B7). Their E-hard negatives are mostly other TKs' substrates with near-identical motifs. Priming arms
  help the priming kinases (0.71 -> 0.79) but lose overall.

## TEST-batch deliverables (G4; `scripts/baselines/simple_test.py`; entries in `freeze_entries/ALL.json`)
- **Common to all files.** All models trained on all DEV (L2; no canaries; sealed/ never read). Every file scores the
  exhaustive `ksp.eval_test.test_candidate_pairs()` set. Each passed `--check-predictions` with `ok: true`, and each is
  pinned with `eval_test.pin_files`.
- **B1 arms (i) and (ii).** Own atlas matrices only: 396 kinases scored, 107 unscored and not imputed (the imputed arm
  was not the CV-selected arm). 10 STY kinases have a matrix for one residue type only: MAP2K1/2/4/5/6/7, TTK, TESK1,
  NEK10, PKMYT1. Their other residue type gets percentile -1, so the kinase is fully scored (398,166 pairs). This
  reproduces CV's `missing='worst'` ranking; the alternative is to leave these 10 kinases unscored.
- **B2** = `fit_pssm` (family prior) on all DEV L2, scored as the U_ms percentile.
- **B7** = the same configuration as on CV, all-DEV fit.
  - The B2 feature for training rows is cross-fitted (inner 5-fold). TEST rows come from a second `crossfit` call
    with no training rows (A1 inner-fit mode). Its 5 inner refits are counted in `inference_wall_s`.
  - `score` = within-kinase linear score (no kinase intercept).
  - `score_site` = A17 softmax over all class-compatible kinases at the site, of (score + kinase intercept)/T.
    T = 0.770 was fitted on labelled DEV sites by OOF log-likelihood (CV run `simple_b7_lr__20261026-153330__s0`).
- **B0-prior.** `score` = DEV L2 count. `score_site` = n_k / sum of n_j over class-compatible kinases. B0-pop and
  B0-miss have no `score_site`.
- **T3 (B2, B7).** For each LKO fold, the labels of that fold's kinases are removed from training, and those kinases'
  TEST candidates are scored by that model. For B7, the held-out kinases have no training rows, so they get the global
  intercept. B2 held-out kinases fall back to group or background priors. Exception: 2 kinases in Manning
  "*-Unique" pseudo-families that span LKO folds (RNASEL/STK31, MAP3K8/MAP3K14) keep a "family" prior.
- **CV caveat.** The CV numbers above are pre-A24 (pooled raw scores). A24 within-fold rank-normalised values were not
  recomputed here.
