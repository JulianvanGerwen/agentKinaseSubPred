# Review r2: statistics (checking DESIGN.md v1)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Reviewer role: biostatistics / ML evaluation. Scope as briefed: §4–§6 and §8.5. I checked endpoints M1–M5, the superiority rule, the matched-universe rule, prior correction, tiered calibration and the coverage estimate. I read r1_statistics.md and compared v1 against it.

## 0. Quick checks on the actual labels (read-only, approximate)

Scripts are in `<scratchpad>/r2stat/` (`dom.py`, `ehard.py`, `more.py`, `keval.py`). Per-kinase output is in `ehard_per_kinase.tsv` in the same folder.

**Inputs:**
- `cache/labels/ks_pairs.parquet` (31,605 rows).
- `cache/index/sites.parquet` (350,736 sites, with `in_U_ms` and `window15`).
- `results/labels/kinase_list.tsv` (503 included kinases).

**Approximations:**
- **Split:** the r1 constrained-optimiser TEST split (r1 mmseqs groups plus identical-15-mer union). The 484 sites not covered by it were assigned randomly at 20%. TEST = 20.1% of sites and 20.8% of native pairs; U_ms TEST sites = 69,736.
- **Par(k)** follows §5 exactly:
  - kinase-domain identity ≥60% (mmseqs, KD coordinates from `kinase_list.tsv`, alignment ≥50% of the shorter domain), OR
  - same Manning subfamily, or same family where no subfamily exists.
- **Masks:**
  - KP(k) = every ks_pairs row for k (all sources and species), plus identical 15-mers;
  - paralog-label mask from all sources;
  - Hamming ≤2 on 15-mers.
- **E-hard membership** uses L2 in-class, non-autophosphorylation labels.
- **Test positives** are human-native, in-class and not autophosphorylation.

### Evaluable sets
| Set | Count | DESIGN v1 says |
|---|---|---|
| K_any (≥1 test positive) | **261** | "~290" (§6.3) |
| ≥3 test positives | 155 | – |
| **K_eval (≥3 positives on ≥2 clusters)** | **135** | "about 160" (§13) |
| ≥5 positives on ≥3 clusters (sensitivity) | 94 | – |
| K_cv (≥5 DEV L2 positives on ≥2 clusters) | 239 | – |
| Kinases with ≥30 DEV positives (per-kinase Platt tier) | 86 | – |

### Does paralog masking make E-hard tiny? No.

**E-hard size per K_eval kinase, after all v1 masks:**

| Class | n kinases | E-hard v0 (median) | E-hard v1 (median) | v1 min | Retained (median / min) | v1 clusters (median) |
|---|---|---|---|---|---|---|
| ST | 104 | 2,047 | 2,028 | 1,774 (CDK2/3) | 99.6% / 86% | 520 |
| STY | 2 | 2,378 | 2,377 | 2,373 | 99.9% | 581 |
| Y | 29 | 319 | 301 | **116** (SRC family) | 95% / **36%** | 138 (SRC family: 66) |

**Key points:**
- **SRC family (SRC, YES1, FYN, LCK, HCK, LYN, BLK, PTK6; 10 paralogs each).** Masking keeps 36–42% of E-hard, which leaves **116–117 negatives on 66–67 clusters**. All eight share almost the same negative set. Masking also removes positives of near-identical specificity, which is the intended effect.
- **Other Y kinases:** 269–324 negatives. ABL2 is the next lowest at 84% retained.
- **Smallest S/T cases:** CDK2/3 (86–90% retained) and PKCs (92%).
- **The Hamming ≤2 mask is negligible** on E-all: median 0 sites, maximum 20 (PRKACA). Low-complexity windows did not cause mass masking.
- **No K_eval kinase has fewer than 100 E-hard negatives.** With n_pos median about 7, AUROC variance is dominated by positives (r1: SD 0.12 at n_pos = 3). The 116 negatives for the SRC family add little, so masking does not threaten M1.

### Other findings that drive the changes below
- **20.5% of test positive pairs fall on sites outside U_ms.** Within K_eval kinases, 2,008 of 2,537 are in U_ms; the share is similar across evidence classes (80–82%). E-all negatives are U_ms only. **13 of the 135 K_eval kinases would have fewer than 3 U_ms positives:** CDC7, PLK4, ERBB2, EIF2AK2, FLT3, ABL2, BLK, MAP2K1, IKBKE, MAP3K5, PINK1, PIM2 and AKT3.
- **Labelled prevalence in U_ms** per kinase (L2 sites / class-matched U_ms sites):
  - 10% / 50% / 90% quantiles: 6.6e-6 / **4.6e-5** / 6.2e-4.
  - A label-scale probability ≥0.5 therefore needs a likelihood ratio ≥ (1−π)/π, which is **≈21,600 for the median kinase** (≥1,600 even at the 10th percentile).
  - Under a binormal model, the fraction of positives above that threshold is about 3e-15 at AUROC 0.80, 2e-6 at AUROC 0.90, and 9e-4 at AUROC 0.95.
- **Coverage by random scores:**
  - Candidate kinases are 415 for S/T sites and 101 for Y sites.
  - With per-kinase top-q selection, the expected share of sites with ≥1 kinase in its top q is 1−(1−q)^K:
    - **q = 1%: 98.5% of S/T sites and 64% of Y sites**;
    - q = 0.5%: 88% / 40%;
    - q = 0.1%: 34% / 10%.
- **Hard-negative pool versus the §5 grid.** The DEV hard pool is roughly 8,000 S/T sites, about 1,200 Y sites and about 470 for the SRC family.
  - **SRC** (564 DEV positives) needs 2,820 hard negatives at the *smallest* grid point (1:20 × 25%).
  - At 1:50 × 50%, 56 kinases need more than 1,200 hard negatives and 10 need more than 8,000.
- **Paralog-set sizes in K_eval:** 21 of 135 have ≥5 paralogs. Some come from the family fallback at low identity: **NEK2 and NEK6 each get 10 paralogs at 33–42% domain identity**; CK1 family (no subfamily) 6; MAP3K5 (STE11) 7.

---

## 1. Proposed changes (prioritised; 8)

### STAT2-1 (critical): The §6.6 coverage measures are uninformative as defined. Replace them with an estimable quantity and a null-referenced curve. §6.6, §1 deliverables
**Problem:**
- **"≥1 kinase in that kinase's top-1%"** is met by **98.5% of S/T sites and 64% of Y sites under random scores**, so it measures the size of the kinome, not the model.
- **"Corrected probability ≥0.5 / ≥0.8"** is on the PSP-label scale (median prevalence 4.6e-5). It needs a likelihood ratio of about 2×10⁴, which an honest model essentially never reaches (see §0). Any non-trivial count would come mainly from over-confident extrapolation (miscalibration), not from accuracy.
- r1 STAT-11(b) is listed as accepted in §15, but the assignment-probability decomposition is not in v1.

**Change:**
- (a) **Primary coverage quantity: site-conditional assignment probability** q(k | s) = P(k ∈ kinases(s) | s is an observed phosphosite).
  - Fit it on CV out-of-fold predictions over *labelled* sites, using the site-grouped model (§8.2) or a per-(k,s) logistic regression on the site-centric score margin. Calibrate it globally plus by residue class.
  - Coverage = the number of U_ms sites (and of sites with no PSP kinase) with max_k q ≥ 0.5 and ≥ 0.8.
  - State the assumption: labelled sites are representative of observed sites conditional on features.
  - Check reliability on TEST labelled sites (reliability diagram; calibration intercept and slope with cluster-bootstrap CIs), stratified by LT_LIT bin {0, 1, 2–4, ≥5} to expose ascertainment drift.
- (b) **Coverage–accuracy curve.**
  - For q ∈ {0.1%, 0.5%, 1%}, report the number of sites covered by ≥1 per-kinase top-q call, next to its **random-score null** 1−(1−q)^K_class and the TEST macro recall at that q (M3 at q).
  - Report only coverage in excess of the null as model-driven.
- (c) Keep the marginal label-scale probability as a deliverable, named "P(PSP-type annotation)" (a lower-bound scale). Report its threshold counts, but do not headline them.

### STAT2-2 (major): The analytic prior offset is wrong under the stratified hard/easy negative sampling, and the negative grid is infeasible for large kinases. §5 training negatives, §6.6, §8.5
**Problem:**
- Training negatives sample the hard stratum (about 3% of U_ms) at 25–50% of all negatives, which oversamples it 10–30×.
- A single case–control offset log ρ_k assumes uniform negative sampling, so the "corrected" probability is biased in a direction that depends on whether a site is labelled to another kinase.
- The grid also cannot be realised: SRC needs 2,820 hard negatives from a pool of about 470 at 1:20 × 25%, and at 1:50 × 50%, 56 kinases need more than the Y-sized pool of about 1,200.

**Change:**
- (i) **Drop the analytic offset.** Fit all calibrators on **exhaustive OOF predictions over every masked U_ms validation site of each fold** (natural labelled prevalence). The calibration intercept then absorbs ρ_k and the hard mix exactly. If an offset is still wanted, use stratum-specific offsets (hard vs. easy), since stratum membership is known at deployment.
- (ii) Define the hard share as **min(requested, available pool)**, sampled without replacement. Log the realised ratio and hard fraction per kinase in the run manifest.
- (iii) **Tiered calibration:**
  - For the 86 kinases with ≥30 DEV positives, fit per-kinase Platt on OOF pooled over the 5 folds. Each fold has only about 6 positives.
  - Below that tier, replace "family-level Platt" with **hierarchical Platt: a kinase-specific intercept with offset log(n_k / N_class) plus a family-shared slope**, shrunk toward a global slope. A pooled family *intercept* would mis-calibrate by up to about 100×, because labelled prevalence varies from 6.6e-6 to 6.2e-4.
- (iv) **Pre-declare an acceptance check:** the TEST calibration slope's 95% CI overlaps [0.8, 1.2]. If it does not, report the miscalibration and do not recalibrate on TEST. This also addresses the mismatch from fitting the calibrator on fold models (80% data) and applying it to the full-DEV model.

### STAT2-3 (major): Fix the M3 universe asymmetry, define the top-1% threshold and ties, and fix the B0-prior floor. §6.2 M3, §6.5
**Problem:**
- 20.5% of test positive pairs are not in U_ms, while E-all negatives are U_ms only.
  - Non-U_ms positives have missing phospho features, so models that use MS features are penalised relative to motif-only baselines (Cantley), which biases the head-to-head.
- "Top 1% of E-all" is ambiguous: positives in or out of the denominator, ties.
- B0-prior is constant within each kinase, so it is 0.5 / 0.01 by construction on M1, M3, M4 and M5.

**Change:**
- (a) M3 positives = test positives **∩ U_ms**.
  - M3 eligibility = K_eval kinases with ≥3 U_ms positives on ≥2 clusters (at most 122 here, since 13 kinases have fewer than 3 U_ms positives; they are listed in §0).
  - Report all-positive M3 as a sensitivity analysis.
- (b) Define **threshold = the per-kinase 99th percentile of E-all negative scores**, i.e. TPR at FPR = 0.01, as r1 STAT-5 intended. A positive tied with the threshold gets fractional credit, so a constant-score model scores exactly 0.01.
- (c) State that B0-prior is a floor **for M2 only**. For the kinase-centric metrics the floors are B0-pop and B0-miss.

### STAT2-4 (major): Make the superiority rule executable. §6.8
**Problem:**
- "Best arm" of B1 chosen after seeing TEST is a selection on TEST.
- Holm needs p-values, not CIs.
- At one-sided α = 0.025 across 9 tests, the smallest Holm threshold is 0.0028. With B = 1,000, a test passes only if at most 1 replicate has Δ ≤ 0, which is Monte Carlo noise.
- Moving the per-test α from 0.025 to 0.0028 raises the MDD from about 2.8 to 3.6 SD, i.e. from about 0.016 to **about 0.021 macro AUROC** (r1 clustering-inflated scale). This should be stated.

**Change:**
- (a) Choose the B1 arm on DEV CV (M1-CV) and record it in FREEZE.json.
- (b) p_j = (1 + #{Δ*_b ≤ 0}) / (B + 1), one-sided. Apply Holm at α = 0.025 across the 9 (comparison × endpoint) tests. Report Holm-adjusted percentile CIs.
- (c) Use **B = 10,000** for the 9 confirmatory tests and 1,000 elsewhere.
  - Cost: about 8.5M weighted-rank operations per replicate per model pair (M1 about 135 × 2k, M3 about 135 × 61k, M2 about 2k sites). That is 10–30 ms in numpy, so a few minutes in total.
- (d) Carry over r1 STAT-7a:
  - K_eval fixed;
  - kinases with 0 resampled positives dropped per replicate, with the mean number dropped reported;
  - duplicated clusters as weights;
  - M3 threshold re-estimated per replicate with the weights.
- (e) Evaluated final predictions = the **mean of the seed models [0..4]**, pre-declared, so that seed noise is not part of the confirmatory Δ.

### STAT2-5 (major): Make M2 consistent with the paralog rule. §6.2 M2
**Problem:**
- M1 treats paralog labels as unidentifiable and masks them, but M2 ranks the true kinase against its paralogs.
- 21 of 135 K_eval kinases have ≥5 paralogs; the SRC family has 10 each.
- A model that correctly ranks the whole SRC family on top has an expected filtered Hit@5 of only about 5/11 ≈ 0.45 on SRC-family sites. M2 therefore partly measures which paralogue the paper tested, which r1 STAT-2 removed from M1.

**Change:**
- **Primary M2:** the filtered rank also removes members of Par(k_true) that are not labelled to the site ("paralog-collapsed Hit@5").
- The current definition becomes secondary ("Hit@5 incl. paralogs"), together with family-level Hit@5.
- Report M2 by residue class: the candidate sets are 415 for S/T and 101 for Y, so random Hit@5 is 1.2% vs 5.0%.

### STAT2-6 (minor): Bound the family fallback in Par(k) by identity. §5
**Problem:**
- Families without a subfamily are masked as a whole even at low identity: Par(NEK2) and Par(NEK6) each hold 10 NEKs at 33–42% domain identity, with distinct specificities.
- The effect on E-hard is small (≤1% for S/T), but under STAT2-5 the fallback also drives M2, and it removes informative training negatives.

**Change:** Par(k) = ≥60% domain identity, OR same subfamily, OR (same family AND ≥50% domain identity).
- This keeps all SRC-family pairs (0.53–0.89) and CK1α/δ/ε/γ (0.56–0.92), and drops the NEK pairs.
- Commit a per-kinase |Par(k)| table at G1.

### STAT2-7 (minor): Pin down E-hard membership and add a size guard and a family-weighted sensitivity. §5, §4.6
**Change:**
- (a) **E-hard membership** is defined by **L2 in-class, non-autophosphorylation labels for every model and label variant**. L1, L2_vitro and L2+ext ablations must not change the evaluation negatives (PROMPT: never change tests).
- (b) **K_eval also requires ≥50 E-hard negatives on ≥20 clusters.** All 135 kinases pass (minimum 116 on 66), so this is a guard against future label or mask changes, not a cut. Commit per-kinase n_neg and n_neg_clusters.
- (c) **Add a sensitivity macro with equal family weight.** The 8 SRC-family kinases share nearly identical negatives and take 6% of M1's weight. The sensitivity macro gives each Manning family weight 1, split equally among its K_eval members.

### STAT2-8 (major): Carry over the missing r1 CV estimator and selection rules, and correct the stale counts. §6.3 M5, §6.8, §13
**Problem:** v1 omits r1 STAT-7(d)/(e) and the STAT-8 selection discipline. Selecting on M1-CV alone can pick a configuration that loses on M2 or M3; the site-grouped objective (§8.2) exists precisely for M2.

**Change:**
- (a) Compute **M1-CV, M5 and every CV per-kinase AUROC fold-stratified**: only positive–negative pairs within a fold are compared, and kinases are weighted by OOF positives across folds. Pooled OOF mixes score scales across fold models.
- (b) **Selection rule:**
  - Maximise M1-CV using the one-SE rule.
  - A candidate must also be within 0.01 of the best on M2-CV and within 0.02 on M3-CV.
  - Ties go to the simpler model.
  - For close calls (Δ < 0.01), repeat CV with 3 fold-assignment seeds.
- (c) **Update the counts** (approximate split, to be re-run at G1): K_eval ≈ 135 (not about 160), K_any ≈ 261 (not about 290), K_cv ≈ 239.

---

## 2. Risks I could not resolve
- **The assignment-probability calibration (STAT2-1a) still assumes** labelled sites are representative of observed sites given features. PSP ascertainment (popular and in-vitro-tested sites) violates this to an unknown degree. The LT_LIT-stratified reliability check only detects drift; it does not fix it.
- **The SRC family's 116 shared E-hard negatives** make its 8 AUROCs highly correlated. The cluster bootstrap captures this only partly, because it resamples substrate clusters, not kinases.
- **The numbers here use the r1 approximate split and my domain-identity proxy** (mmseqs on UniProt KD coordinates, not KinCore). Counts should be stable to about ±10%. Re-run `ehard.py`-style checks on the frozen G1 split, since FAM20A/C, FASTK, PAN3, PEAK3, PLK5 and TSSK6 lack KD coordinates and default to family-only Par(k).
- **In vivo-only labels include indirect pairs.** This caps attainable AUROC by a kinase-specific, unknown amount. The in-vitro stratum is the only check.
