# Review r1: statistics and ML evaluation (DESIGN.md v0)

Reviewer role: biostatistician / ML-evaluation expert. Scope: PU framing, negative sets, endpoints, minimum positives and variance, cluster bootstrap, multiplicity and model selection, calibration and precision for coverage, baseline comparison, and the test-split stratification algorithm.

## 0. Quick quantitative checks (read-only, approximate)

Scripts are in the reviewer scratchpad (`/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/stat/{prep,split,sim,sim2,ana,pop,var}.py`). They are not part of the pipeline, and the PI should port the logic into `scripts/` once the real split exists.

**Approximations:**
- **U_obs** = PSP human sites ∪ compendium sites, Swiss-Prot canonical, with the residue checked against UniProt 2026_03. This gives 303,218 sites on 18,497 proteins; CPTAC and Ochoa are not included.
- **Labels** = L2 (HMR, SGID mapping plus sequence fallback), canonical only, with the §5 residue class applied: 18,297 pairs, of which 222 out-of-class pairs are dropped.
- **Test positives** = human kinase × human substrate.
- **Substrate grouping:** `mmseqs easy-cluster` (40% identity, 50% coverage, cov-mode 0; 26 s on 4 threads) plus union-find on identical centred 15-mers. This gives 11,766 groups. The largest group holds only 1.1% of label pairs, so there is no giant-component problem.
- **Kinase paralogue proxy:** mmseqs at 50% identity on full-length kinase sequences, not on domains.

| Check | Result |
|---|---|
| Kinases with L2 labels / with ≥5 L2 labels | 431 / 292 |
| **20% cluster split, balanced optimiser (5 seeds): kinases with ≥1 / ≥3 / ≥5 / ≥10 / ≥20 human-native test positives** | **283–301 / 159–166 / 102–108 / 63–66 / 27–29** |
| Same, random 20% cluster split (10 seeds, mean [min–max]) | 248 / 147 [132–164] / 107 / 63 / 32 |
| Ideal case (exactly 20% of each kinase's native labels) | ≥3: 143, ≥5: 99, ≥10: 62 |
| ≥3 test positives, by L2-label bin 1–4 / 5–19 / 20–99 / ≥100 | 4/139, **23/152**, 93/100, 40/40 |
| Same pipeline, grouped-CV out-of-fold (OOF) on DEV: kinases with ≥3 / ≥5 / ≥10 OOF native positives | 297 / 247 / 161 |
| Test universe | 60,632 sites (S 37.6k, T 14.8k, Y 8.2k); **E-hard universe 1,951 sites (Y only 265)** |
| Per-kinase prevalence (160 evaluable kinases), median [10–90%] | E-all 2.2e-4 [6e-5, 1.2e-3]; E-hard 6.8e-3 [1.8e-3, 3.8e-2] |
| Positives per evaluable kinase: median 7; distinct substrate clusters: median 4 | **28% of kinases with ≥3 positives have them on ≤2 clusters** |
| Popularity-only (PSP LT_LIT) macro AUROC, E-all / E-hard | **0.966 / 0.594** (MS count: 0.622 / 0.538) |
| Popularity-only macro AUPRC on E-all | 0.047 (median lift 50×; macro recall in top 1% = 0.42) |
| E-hard negatives labelled to a same-paralogue-cluster kinase | mean 5%, but **55–58% for BLK, FGR, YES1, HCK, LCK, LYN**; 8–12% for CDK2/3/5 |
| Rodent-only L2 pairs on test sites (would be false "negatives" if not masked) | **568 pairs across 84/160 kinases** (PRKACA 73, MAPK1 49, CDK5 34) |
| L2 pairs whose site is not in U_obs | 1,234 / 18,045 (6.8%) |
| Cantley atlas coverage of the 160 evaluable kinases | 158 (the 2 missing, CSNK2B and PKM, are excluded from K anyway) |

**Variance (binormal simulation; per-kinase AUROC 0.80):**
- **Per-kinase SD:**
  - E-hard AUROC has SD 0.21 / 0.12 / 0.10 / 0.07 / 0.05 at n_pos = 1 / 3 / 5 / 10 / 20.
  - E-all AUPRC (AUROC 0.90, 50k negatives) has SD ≈ 2.8× its mean at n = 3 (0.054 vs 0.019) and ≈ 1.2× at n = 10.
- **Macro level:**
  - Macro E-hard AUROC over the evaluable kinases has SD ≈ 0.006 from sampling alone.
  - A paired Δ = 0.02 between two correlated models (ρ = 0.8) has SD ≈ 0.004. The minimum detectable difference (MDD, ≈ 2.8 SD) is therefore ≈ 0.012, **almost the same for thresholds ≥3 (K = 159), ≥5 (K = 105) and ≥10 (K = 63)**.
  - The simulation assumes independent sites. Clustering (0.56 clusters per positive) inflates this by roughly √1.8, giving an MDD of about 0.016.

Conclusion: the test set has enough power for macro comparisons at Δ ≳ 0.015–0.02. However, it **cannot evaluate understudied kinases**: 68% of labelled kinases have fewer than 20 labels, and these are almost all absent from the ≥3-positive set.

---

## 1. Proposed changes (prioritised)

### STAT-1 (critical): Understudied kinases are structurally absent from the test endpoints. Add pre-registered endpoints that include them. §6 (Kinase-centric, Primary, Stratification); §4.3
**Change:**
- (a) Add **secondary endpoint S1**: macro E-hard AUROC over **all kinases with ≥1 test positive** (≈280–300 kinases), reported overall and by training-label bin.
  - Macro AUROC with n_k = 1 is unbiased (it is the kinase-equal mean of per-positive normalised ranks). It is only noisy per kinase, and the macro over about 290 kinases stays precise.
  - Never interpret per-kinase values below n ≈ 10.
- (b) Pre-declare the **DEV grouped-CV OOF** macro metrics for the frozen final configuration as the main evidence for the 5–19 label bin: 247 kinases have ≥5 OOF positives, compared with 1 kinase on TEST.
  - Label these as optimistic, because selection used the same CV.
  - To remove that optimism, run the *whole* selection procedure inside each outer fold (nested, selection-aware CV) for the final comparison. This is cheap for LightGBM.
- (c) Optional and clearly secondary: after freezing the pipeline, run a "frozen-pipeline full-data grouped CV" over DEV ∪ TEST, with no further selection. Report the TEST-cluster and DEV-cluster parts separately.

**Why:** with a 20% cluster split, only 23 of 152 kinases with 5–19 labels, and 4 of 139 kinases with 1–4 labels, reach ≥3 test positives. As written, the primary endpoint is effectively "kinases with ≳20 labels". That contradicts the PROMPT goal of harmonious performance including lowly studied kinases.

### STAT-2 (critical): Remove paralogue positives from E-hard negatives. §5 (Evaluation negatives)
**Change:**
- Redefine **E-hard(k)** as test sites with ≥1 annotation, of k's residue class, that are not labelled to k **and not labelled to any kinase in k's paralogue cluster**. Use the same cluster definition as the §4.5 audit (≥50% kinase-domain identity, or the same Manning family for the SRC/CDK/AKT/PKC groups). Make this the primary negative set.
- Keep the current definition as **E-hard-incl-paralog**, a sensitivity analysis.
- Report E-hard separately for the S/T universe (≈1.7k sites) and the Y universe (≈265 sites).

**Why:**
- For SRC-family kinases, 55–58% of E-hard negatives are sites annotated to another SRC-family member.
- With near-identical specificity, many of these are true but untested positives (worst-case PU contamination). E-hard then measures "which paralogue the paper happened to test".
- That contamination depresses AUROC for exactly the kinases where paralogue separation is unidentifiable. It also favours models that memorise annotation habits.

### STAT-3 (major): Mask, rather than count as negatives, every (k,s) pair with any non-test evidence. §3 (Test labels); §5
**Change:** for kinase k on test sites, the following pairs are **excluded from both positives and negatives**:
- rodent-transferred (L2) and other-species (L3) pairs;
- external curated pairs (UniProt, SIGNOR, OmniPath);
- out-of-class labels;
- paralogue labels under STAT-2.

State this in §5 and test it in `tests/`, for example with a toy kinase whose rodent label sits on a test site: the assert fails if that site counts as a negative.

Also state how U_obs treats labelled sites. 6.8% of L2 pairs fall on sites not present in the PSP-site ∪ compendium union. Recommended default: include all PSP KS-dataset sites in U_obs, since they are low-throughput-observed phosphosites, and run a sensitivity analysis that excludes them. Their count, and the kinases they push across the ≥3 threshold, must be reported.

**Why:** 568 rodent-only pairs on test sites affect 84 of the 160 evaluable kinases. Counting them as negatives systematically penalises models that generalise across species, which is the very pooling the design adopts.

### STAT-4 (major): Popularity is not uninformative on E-hard. Make the popularity-only model a reported floor on both endpoints, and add an ascertainment-matched E-all. §5 (Ascertainment control)
**Change:**
- Delete "where popularity is uninformative by construction".
- Pre-declare the popularity-only model (LT_LIT plus MS count, logistic, fit on DEV) as baseline B0 on every endpoint. Report Δ(model − B0).
- Add **E-all-matched** as a secondary analysis: per kinase, inverse-probability-weight E-all negatives so their LT_LIT-stratum × residue distribution matches the positives (strata: 0, 1, 2–4, ≥5).

**Why:** on the actual split, LT_LIT alone gives macro AUROC 0.966 on E-all and 0.594 on E-hard, and macro E-all AUPRC 0.047 (about 50× the random level). E-hard reduces ascertainment but does not cancel it, because heavily studied sites are positives for many kinases. Features such as PPI degree, abundance and compendium coverage can proxy popularity even though LT_LIT itself is excluded.

### STAT-5 (major): Replace the macro E-all AUPRC co-primary with a top-rank, prevalence-invariant metric. Keep AUPRC as secondary. §6 (Co-primary; AUPRC lift)
**Change:**
- Make the co-primary the **macro per-kinase recall in the top 1% of E-all candidates** (TPR at FPR = 1%). This matches the §6 deployment threshold.
- Report the **McClish-standardised partial AUROC (FPR ≤ 1%)** as its smooth companion.
- Keep AUPRC as secondary, always with the stated prevalence.
- Do not macro-average raw AUPRC lift, which is unbounded at about 1/prevalence and dominated by the rarest kinases. Report its median or geometric mean instead.
- Report S/T and Y universes separately, because candidate counts differ about 6×.

**Why:** per-kinase E-all prevalence is about 2e-4 with a median of 7 positives. At n = 3–5, per-kinase AUPRC has SD about 2.5–2.8× its mean, so the macro is driven by a handful of high-AP kinases. AUPRC values also depend on the arbitrary size of the test universe, which changes with the 20% fraction and the U_obs definition. Recall at a fixed top fraction is directly interpretable ("of PSP substrates, x% fall in the kinase's top 1%"), is far less noisy, and is invariant to universe size.

### STAT-6 (major): Base kinase eligibility on effective sample size; restrict per-kinase claims. §6 (Kinase-centric)
**Change:**
- Eligibility for the primary macro: **≥3 test positives on ≥2 distinct substrate clusters**.
- Report a "≥5 positives on ≥3 clusters" subset as a sensitivity analysis.
- Per-kinase plots and tables show empirical-Bayes-shrunk AUROCs (beta-binomial or normal-normal shrinkage toward the group mean) next to the raw values. Macro metrics always use raw values.
- No per-kinase "wins" are claimed for n < 10.
- Evidence-class stratified macros are compared only on the kinases evaluable in every class. Kinases with ≥3 positives: both 83, in vitro only 76, in vivo only 61.

**Why:**
- 28% of kinases with ≥3 positives have them on only one or two clusters, so they contribute one or two independent observations, and the cluster bootstrap degenerates for them.
- Simulation shows the threshold barely affects macro power (MDD about 0.012 for ≥3, ≥5 and ≥10), so ≥3 is kept for coverage.
- Per-kinase SD is 0.12 at n = 3, so individual-kinase numbers are noise.

### STAT-7 (major): Fully specify the bootstrap and CV estimators. §6 (Uncertainty); §4.3; §4.4
**Change:**
- **(a) T1 bootstrap.** Resample substrate clusters (the split unit), stratified by labelled vs unlabelled cluster, with **B = 2,000** for p-values (1,000 is fine for CIs).
  - Fix the evaluated kinase set from the observed test set.
  - In each replicate, drop kinases with 0 resampled positives from the macro and report the mean number dropped.
  - Duplicated clusters enter as weights.
  - Use percentile CIs on paired differences, with the same replicate indices for all models.
- **(b) Training stochasticity.** The evaluated final predictions are the **mean over the 5 seed models** (pre-declared ensemble). The bootstrap CI is conditional on that ensemble, and the seed SD is reported separately.
- **(c) Kinase-level generalisation.** Report a secondary crossed (kinase-family × substrate-cluster) bootstrap (Owen 2007) for generalisation to "kinases like these".
  - For the LKO tiers T2/T3, the resampling unit is the **kinase family**, since the population of inference is unseen kinases.
  - With 5 family folds, also report fold-level spread.
- **(d) CV per-kinase AUROC.** Compute it **fold-stratified**: only positive–negative pairs within the same fold are compared. This avoids mixing score scales across fold models.
- **(e) Close CV decisions.** Repeat grouped CV with 3 fold-assignment seeds; training seeds alone do not capture split variance.

**Why:** the design states "cluster bootstrap, 1,000" without handling for varying kinase eligibility, seeds, LKO units or cross-fold pooling. Each of these changes CI width or bias materially.

### STAT-8 (major): Pre-register a small confirmatory family with multiplicity control and selection discipline. §6; §8 (Final model); §12
**Change:**
- **Confirmatory family:**
  - Comparisons: final model vs B1 (Cantley, priming-aware), B3 (Phosformer, retrained) and B4 (SELPHI2.0, retrained). As-published B3/B5/B6 are excluded because their training overlaps the test set.
  - Testing order within each baseline: primary first, then co-primary (fixed-sequence gatekeeping). Holm correction across the 3 baselines; one-sided α = 0.05; paired cluster bootstrap.
- **Everything else** (strata, ablations, label variants, LKO) is descriptive, with CIs and no significance claims. If p-values are shown, apply BH-FDR within each table.
- **On CV:**
  - Log every configuration evaluated, and report the count.
  - Choose with the one-SE rule on the primary endpoint.
  - Treat paired differences below about 0.01 (the CV MDD) as ties and prefer the simpler model.
  - Report CV-to-test optimism (best CV − test) for the final model and each baseline.

**Why:** the plan compares at least 6 baselines, 4 regimes, an α grid, negative mixes and 5 label variants on the same CV. Without a declared confirmatory set and selection rule, the headline is a maximum over many noisy comparisons (winner's curse).

### STAT-9 (major): Fair comparison when baselines do not cover every kinase or site. §6 (new "Baseline comparison rules"); §8 (R0)
**Change:**
- (1) The evaluated kinase set is fixed independently of any method.
- (2) Confirmatory comparisons are computed on the **intersection** of kinases and sites both methods can score.
- (3) Also report a **coverage-penalised full-set** metric. Unscored pairs get a tied lowest score, so an unscored kinase contributes AUROC 0.5 and zero recall.
- (4) Report each method's coverage: the share of evaluable kinases and of test positives it can score.
- (5) Apply identical residue-class rules, test masks (STAT-3) and the site universe to every baseline.
- (6) Give retrained baselines (B2–B4) the same CV tuning budget as R2. Cap it at N configurations each, and report N.

**Why:** coverage differences can make either side look better by choice of denominator. On the evaluable set Cantley covers 158 of 160 kinases, so the issue bites mainly for NetworKIN, KolossuS, window-truncated sites and the S1 tail endpoint.

### STAT-10 (major): Specify the test-split stratification algorithm and its acceptance criteria. §4.2
**Change:** replace "greedy … where cluster structure permits" with a deterministic optimiser:
- **Objective:** minimise Σ_k (t_k − 0.2·n_k)² / n_k over kinases with ≥5 human-native labels, where n_k is native positives and t_k is test positives.
- **Hard constraints:**
  - U_obs sites 20 ± 0.5%;
  - total native pairs 20 ± 0.5%;
  - Y sites 20 ± 1%;
  - each kinase group's pairs 20 ± 3%.
- **Search:** single-cluster flip and swap local search, seed 0.
- **Committed output:** the achieved per-kinase distribution, plus a DEV-vs-TEST balance table (LT_LIT, kinases per site, evidence class, has-primed-neighbour, autophosphorylation share).
- **Stratify on human-native positives,** since those are the test labels, not on L2.

**Why:** a naive greedy over clusters (largest first, site budget only) **put 25.5% of native pairs into TEST** and kept only 57.5% of eligible kinases within 10–30%. Random cluster splits keep 47%. The constrained optimiser above, tried here, gives exactly 20.0% of pairs, 77–80% of eligible kinases within [0.10, 0.30] (65–69% within [0.15, 0.25]), and 159–166 kinases with ≥3 positives, compared with 147 for random. The greedy's apparent gain (198 kinases with ≥3) is bought by moving training data into TEST.

### STAT-11 (major): Precision for coverage is not identifiable under PU for dark sites. Redefine and validate it orthogonally. §6 (Coverage); §8 (Calibration)
**Change:**
- (a) Rename the quantity "estimated precision" to **label-precision**: the precision against PSP labels on E-all. It is a lower bound only if labels are correct and labelling does not depend on score.
- (b) Report Priority-2 coverage through a decomposition that is estimable:
  - **Assignment probability** P(k | s is a target of some kinase in K, s) comes from a site-centric softmax calibrated on CV, with reliability evaluated on annotated test sites (Hit@1 calibration).
  - Coverage = the number of dark sites whose top kinase has assignment probability ≥0.5 and ≥0.8, stated as **conditional on the site being a K-kinase target**.
- (c) Validate dark-site predictions orthogonally:
  - non-PSP curated pairs on test proteins at sites without PSP KS labels;
  - inhibitor-response concordance in held-out perturbation studies excluded from features.
- (d) Replace per-kinase isotonic with **hierarchical Platt** (kinase intercept and slope shrunk toward the family), fit on OOF, for kinases with <50 OOF positives.
  - Check calibration on TEST with calibration intercept and slope (with bootstrap CIs) and reliability diagrams.
  - Note that the calibrator is fit on 80%-data fold models but applied to the full-DEV model.

**Why:** "dark" sites have label frequency 0 by definition. No prevalence correction from PSP can estimate their precision, so a "precision ≥0.8" claim for them would be unverifiable. Isotonic regression on fewer than about 50 positives overfits badly.

### STAT-12 (minor): Make the site-centric metric fair and comparable. §6 (Site-centric)
**Change:**
- Rank the true kinase among kinases in K of the matching class **after per-kinase percentile normalisation** (as the atlas does). Report the raw-score variant too.
- Compute the metric over a fixed candidate kinase set, the same for every method (intersection rule, STAT-9). Hit@k depends on the candidate set's size.
- Report both the kinase-frequency prior and the popularity floor.
- Exclude paralogue-cluster kinases of the true kinase from the competitor list in a sensitivity analysis ("Hit@k up to paralogue").

**Why:** raw scores are not comparable across kinases, and paralogue ties are unidentifiable from PSP.

---

## 2. Risks I could not resolve
- **Label noise.** In vivo-only PSP labels include indirect relations (69% of S/T-on-TK labels are in vivo-only). This caps attainable AUROC, and its kinase-specific magnitude is unknown. The in vitro-supported stratum is the cleanest check.
- **Dark-site precision.** True precision for dark sites remains unidentified even after STAT-11. Any absolute coverage number rests on a stated assumption.
- **Approximations in my checks:**
  - full-length rather than domain identity for paralogues;
  - U_obs without CPTAC or Ochoa;
  - set-cover mmseqs clustering plus 15-mer union-find.

  The kinase counts should be stable to about ±10%; re-run them on the real split.
- **Tail evidence.** Evidence for kinases with <5 labels will remain weak on every endpoint. LKO tiers and S1 help, but conclusions there should be framed as exploratory.
- **No temporal validation.** The KS file has no PMIDs or dates. A future PSP release would be the cleanest prospective test; record the current release hash so new pairs can be scored later without leakage.
