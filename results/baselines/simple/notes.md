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
