> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Independent audit of the single TEST batch (721dc1977366)

**Auditor:** an independent results reviewer. I produced none of the audited work.
**When:** 2026-10-27, 01:00–02:00 CET, after the COMPLETED ledger block (00:57:29).
**Scope:** `results/test/LEDGER.md`, `results/test/721dc1977366/`, `FREEZE.json`, DESIGN §6 and §16–17,
`docs/G4_RUNBOOK.md`, `docs/audit/audit_cv_r2.md`, `docs/model_r2.md` §6, `docs/strata_cv.md`.

**What I did:**
- I recomputed the batch from the pinned prediction files and `sealed/test_labels.parquet` (read-only), using the
  frozen `ksp.eval` / `ksp.eval_test` functions.
- I changed no committed file. I ran no git commands. Memory stayed under 12 GB; at most 8 cores were used.
- This is a verification pass. Every TEST number below is either a batch number or a recomputation that checks
  one. Rows marked **(audit diagnostic)** are re-cuts of batch quantities (subsets of the batch's per-kinase
  tables, a matched-kinase restriction, or a corrected truth set). They are not new endpoints and must not be
  quoted as TEST results.
- The CV comparisons use the committed CV per-kinase tables (`results/runs/<run>/metrics_a24.json`, and
  `per_kinase.tsv` for R3).

**Scripts:** `docs/audit/audit_test_scripts/`. They write to a session-scratchpad `OUT` path; change it before re-running.

| Script | What it does |
|---|---|
| `recompute.py` | Endpoints, strata, ext tiers and reliability per config |
| `compare.py` | The 9 primary matched comparisons; B = 10,000 bootstrap for B1_i and B4 |
| `m2cv.py` | CV vs TEST on common kinase sets |
| `cands.py` | M2 competitor-set sizes, CV vs TEST |
| `matched_strata.py` | Strata restricted to B1-scored kinases; reliability with full site truth |

## 0. Verdicts

| # | Headline claim | Verdict | Basis / required wording |
|---|---|---|---|
| 1 | Final model TEST **M1 0.794** [0.774, 0.811], **M2 0.363** [0.324, 0.397], **M3 0.363** [0.309, 0.401], M4 0.799 | **TRUSTED** | Reproduced to 0e+00 (§1). CV M1 0.786 lies inside the TEST CI. Kinase-resampling CIs: M1 [0.767, 0.821], M2 [0.310, 0.414]. |
| 2 | Final **beats Cantley** on M1 (+0.075 [0.046, 0.099]), M3 (+0.205) and M2 (+0.179 vs arm ii); Holm p 0.0018 | **TRUSTED** | Point estimates, CIs and p reproduced exactly at B = 10,000 (§3). 97 of 116 matched kinases favour the final model on M1. Also holds under paired kinase resampling (audit diagnostic: M1 [0.057, 0.094]). Comparator = arms (i) and (ii), *not* the priming-aware arm. |
| 3 | Final beats **Phosformer-retrained (B3)** on M1, M2 and M3 | **TRUSTED WITH CAVEATS** | Reproduced. Caveat: B3 is our faithful retrain on our split (fidelity check passed). On CV it already sat far below its published numbers. Say "beats our retrained Phosformer". |
| 4 | Final beats **SELPHI2.0-retrained (B4)** on M1 and M2 | **TRUSTED WITH CAVEATS** | Reproduced. Caveat: 66.6% of B4's TEST pairs use imputed (NA-default) site features. The fold-0 CV check puts the imputation cost at about 0.013 M1, so the gap is not explained by it, but it must be stated. |
| 5 | Final beats B4 on **M3** (+0.059 [0.0096, 0.120], p 0.0178) | **TRUSTED WITH CAVEATS: soften** | It is the only primary test off the p floor. The CI lower bound is 0.01; only 52 vs 27 kinases favour the final model (40 ties). It survives kinase resampling (audit diagnostic: [0.021, 0.097]). Report it as "modest", not "large". |
| 6 | "p = 0.0002" for 8 of 9 tests | **TRUSTED, reword** | 0.0002 = 2/(B+1) is the **resolution floor**: not one of 10,000 resamples had Δ ≤ 0. Write "p < 0.0002 (bootstrap floor); Holm p ≤ 0.0018", not an exact p. |
| 7 | **M2 rose from CV to TEST** (final 0.263 → 0.363) | **TRUSTED; explained** | K_eval composition (§2). On the same 117 kinases, CV M2 = 0.348 vs TEST 0.363 (paired Δ +0.014 ± 0.021). No leak signal. REPORT must not present 0.363 as an improvement over CV, and must quote the K_cv-comparable CV value (0.348) next to it. |
| 8 | CV→TEST M1 stable (0.786 → 0.794) | **TRUSTED** | On common kinases: CV 0.792 vs TEST 0.793 (paired Δ +0.001 ± 0.012). Matches the audited optimism ≲ 0.005. Caveat in §4: ST kinases improved (+0.014) while Y kinases dropped (−0.072 ± 0.033). |
| 9 | Final beats Cantley **in every stratum** / on understudied kinases | **TRUSTED WITH CAVEATS** | True for M1/M4 in every matched stratum and training bin (§5). The exceptions and near-ties are listed in §5. No stratum has a paired CI in the batch; all strata statements are descriptive. |
| 10 | Tyr kinases: final 0.595 vs Cantley 0.496 | **NOT TRUSTED as a skill claim** | Final is barely above the **site-popularity floor B0-pop (0.585)** and beats it on only 10 of 23 Y kinases. T3 on Y is 0.531. Report it as "no demonstrated Tyr-kinase specificity on TEST; Cantley is at chance". |
| 11 | Primed sites: final 0.897 vs Cantley 0.704 | **TRUSTED WITH CAVEATS** | Only vs arm (i), which is not priming-aware. On CV the priming-aware arm (B1prim) beat R2 on primed positives (0.945 vs 0.870–0.893), and it is not in FREEZE, so TEST cannot say we beat priming-aware Cantley. Small: 6 kinases, 65 pairs. Unprimed positives are a near-tie (0.740 vs 0.725). |
| 12 | T3 (unseen family + unseen protein): final M1 **0.752** | **TRUSTED WITH CAVEATS** | Consistent with CV LKO-T2 (0.749). It exceeds label-free Cantley's 0.714 (audit diagnostic on B1-scored kinases: 0.755 vs 0.719, 79 of 116 kinases), but **no paired T3 test exists**, so the comparison is descriptive. On Y kinases T3 = 0.531. |
| 13 | A17 assignment probability honest on labelled TEST sites (ECE 0.026; coverage@0.5 = 5,773, coverage@0.8 = 0) | **TRUSTED WITH CAVEATS** | ST at ≥ 0.5: 0.615 claimed vs 0.603 observed (n = 345, Wilson [0.55, 0.65]). No ≥ 0.8 tier exists, and no Y site can exceed 0.445. Y sites at [0.2, 0.5) are mildly overconfident (0.41 vs 0.34, n = 195). The temporal/external ECEs (0.12 / 0.22) are mostly a truth-set artefact (§6). |
| 14 | External-curated and temporal strata: final ahead of Cantley | **TRUSTED WITH CAVEATS** | Tier A: 0.788 vs 0.684 (41 kinases, 134 pairs). Tier B: 0.804 vs 0.661, but tier B is popularity-heavy (B0-pop 0.725). New-since-2021: matched 0.771 vs 0.713. The A16 calibration on these strata (A29) is **not in the batch**. |
| 15 | Final > R2 and final > R3 (fusion helps) | **TRUSTED** | Secondary family; all 8 tests have Holm p ≤ 0.013. M2 vs R3 (+0.041, p 0.013) is the weakest. |

## 1. Reproduction (claim 1)

**Pins and integrity.**
- All 17 prediction/T3 files match their FREEZE.json SHA-256 pins.
- FREEZE.json's own SHA (721dc197…) equals the ledger value.
- The sealed `test_labels.parquet` SHA matches SEAL.json.
- The current `src/ksp/{eval,eval_test,strata,labels,splits,index,kinases}.py` match the ledger's `code_sha256`.
- TEST and DEV share 0 proteins and 0 substrate clusters.

**Recomputed endpoints.** I used `E.prepare('test', labels_df=sealed, with_canaries=False)` and then
`E.evaluate(..., missing='worst')`, exactly as `run_batch` does. Every difference from the batch is 0.

| config | M1 | M2 | M3 | M4 | batch diff |
|---|---|---|---|---|---|
| final | 0.794046 | 0.362706 (score_site) | 0.363249 | 0.798675 | 0 |
| B1_i | 0.713684 | 0.162761 | 0.160749 | 0.698867 | 0 |
| B1_ii | 0.713659 | 0.188829 | 0.160749 | 0.698855 | 0 |
| B3 | 0.707901 | 0.163372 | 0.104037 | 0.687780 | 0 |
| B4 | 0.704693 | 0.232656 | 0.304402 | 0.700819 | 0 |
| R2 / R3 / B0-prior | 0.7783 / 0.7748 / 0.5000 | 0.2637 / 0.3223 / 0.1088 | 0.3118 / 0.2920 / 0.0100 | | 0 |

**Other checks reproduced exactly (max |diff| = 0):**
- the strata tables of final, B1_i and B1_ii;
- `final_reliability.tsv`;
- coverage@0.5 = 5,773 and @0.8 = 0;
- the label-free isotonic caps: the maximum top-kinase probability over all TEST sites is 0.756 for S/T and 0.445 for Y.

**Final score columns.** I reconstructed the final file from the pinned R2 and R3 TEST files:
- `score` = 0.6·z(R2) + 0.4·z(R3) and `score_site` = 0.5·z(R2) + 0.5·z(R3) (global z over TEST rows), with max |diff| = 0;
- `assignment_prob` sums to 1 per site.

So no label-dependent quantity enters any final column.

## 2. The M2 jump (claim 2)

**Answer: it is the K_eval composition, not a leak or a harness difference.**
- K_eval is the 119 kinases with ≥ 3 TEST positives on ≥ 2 clusters.
- By training bin (5–19 / 20–99 / ≥ 100), K_eval has 13 / 73 / 31 kinases, plus 2 in bin 0. K_cv has 83 / 89 / 31.
- Since M2 is a macro over kinases and Hit@5 rises steeply with label count, dropping 70 low-label kinases raises it.
- Selection also acts *within* bins: the K_eval members of a bin are its better-studied ones.

| | final | R2 | R3 | B1 (ii) | B1 (i) | B0-prior |
|---|---|---|---|---|---|---|
| CV M2 over K_cv (203) | 0.263 | 0.200 | 0.223 | 0.171 | 0.157 | 0.064 |
| **CV M2 over K_eval ∩ K_cv (117)** | **0.348** | **0.246** | **0.320** | 0.180 | 0.160 | **0.111** |
| CV per-bin M2 reweighted to K_eval's bin mix | 0.336 | 0.242 | 0.302 | 0.176 | 0.160 | 0.109 |
| TEST M2 over K_eval ∩ K_cv (batch per-kinase) | 0.363 | 0.260 | 0.328 | 0.184 | 0.157 | 0.111 |
| TEST M2 over K_eval (batch) | 0.363 | 0.264 | 0.322 | 0.189 | 0.163 | 0.109 |
| paired per-kinase TEST − CV (117 kinases), mean ± SE | +0.014 ± 0.021 | +0.014 ± 0.020 | +0.008 ± 0.020 | +0.004 ± 0.016 | −0.003 ± 0.016 | −0.000 ± 0.004 |

**How to read the table:**
- Of the final model's +0.100 CV→TEST change, +0.085 is kinase selection. The remaining +0.014 is well inside one SE.
- **B0-prior is the clean control.** It ranks kinases by DEV count only, so it contains no site information. It goes from 0.064 to 0.109 on TEST, and its CV value on the same kinases is 0.111.
- Cantley is nearly prior-free and flat across bins (CV 0.16–0.19), so it barely moves. That is why it "moved less".

**Other mechanisms checked:**

| mechanism | finding |
|---|---|
| Competitor-set size (`cands.py`) | The M2 candidate set is the same real-kinase set on CV and TEST; canaries are excluded on CV by `_site_centric`. Mean competitors per positive pair (kinase-macro): 345.7 on TEST vs 345.7 on CV for the same kinases. TEST sites carry slightly more labelled kinases (5.1 vs 4.2 per site), so filtering removes about 1.3 more competitors per pair (6.7 vs 5.4). That is a second-order effect, bounded by the +0.014 residual. |
| Label-source mix | Same in TEST and DEV label tables (PSP 36% vs 35%, OmniPath, mask rows), so the "any source" filter behaves the same. |
| score_site vs score | Final uses `score_site` (w = 0.5) on both CV and TEST; R2/R3 use `score` on both. No switch. |
| Leak | TEST/DEV are disjoint at protein and cluster level. score_site is an exact label-free function of the pinned R2/R3 files. The prior-only control reproduces its CV value. **I found no leak signal.** |

**By-bin TEST M2 (batch strata) vs CV on all K_cv:**

| bin | final TEST | final CV | Cantley (ii) TEST | Cantley (ii) CV |
|---|---|---|---|---|
| 5–19 | 0.208 (13 kinases) | 0.142 (83) | 0.165 | 0.162 |
| 20–99 | 0.303 | 0.250 | 0.190 | 0.171 |
| ≥ 100 | 0.568 | 0.621 | 0.177 | 0.193 |

These by-bin TEST numbers are themselves selected kinase subsets: the 13 TEST 5–19 kinases are the best-labelled of that bin. The CV finding that the model is weak on M2 for understudied kinases (it loses to Cantley at 5–19 on CV) is **not refuted** by TEST. TEST simply has too few such kinases.

## 3. The primary tests (claims 2–6)

**Matched universe.**
- B1 scores 386 of 503 kinases: 117 are unscored, including the 10 single-residue STY kinases of A30.
- The common pairs number 19,215,814 of 25,922,452, so frac A kept = 0.7413 ✓.
- Three K_eval kinases are unscored by B1 (TTK, MAP2K1, LIMK2), so the matched M1/M3 are over 116 kinases ✓.
- M2 is matched-kinase (A27): the competitors are restricted to the 386 kinases, which raises the final model's M2 to 0.372.
- B3 and B4 are exhaustive (503 kinases, frac 1.0) ✓.

**Bootstrap.**
- B = 10,000, seed 0, paired cluster bootstrap with a two-sided p = 2·min(·)/(B+1).
- I re-ran final vs B1_i (M1, M3) and final vs B4 (M1–M3) at B = 10,000. **CI bounds and p equal the batch to the last digit.**
- Point estimates for all 9 tests reproduce exactly.

**Holm.**
- The Holm step-down over the 9 batch p-values gives 0.0018 × 8 and 0.0178 ✓.
- p_boot 0.0002 = 2/10,001 is the floor. The true p is smaller, and "p = 0.0002" must not be quoted as exact.
- All 9 meet "beats" (Holm p < 0.05 and CI lower bound > 0).

**Sensitivity (audit diagnostic, from batch per-kinase tables).** I ran a paired kinase-resampling bootstrap, which treats the kinase panel as random (the cluster bootstrap holds it fixed). All 9 still exclude 0. The weakest is B4 M3: [0.021, 0.097].

| test | Δ | cluster CI (batch) | kinase-resampling CI (audit) | kinases A > B / A < B |
|---|---|---|---|---|
| vs B1 (i) M1 | 0.075 | [0.046, 0.099] | [0.057, 0.094] | 97 / 19 |
| vs B1 (i) M3 | 0.205 | [0.165, 0.235] | [0.167, 0.243] | 78 / 4 |
| vs B3 M1 | 0.086 | [0.059, 0.105] | [0.065, 0.110] | 97 / 22 |
| vs B4 M1 | 0.089 | [0.063, 0.117] | [0.069, 0.111] | 98 / 20 |
| vs B4 M2 | 0.130 | [0.091, 0.156] | [0.089, 0.172] | 70 / 16 |
| **vs B4 M3** | **0.059** | **[0.010, 0.120]** | **[0.021, 0.097]** | **52 / 27** |

**Wording for the baseline comparisons:**
- **Cantley gap.** The TEST gap vs Cantley M1 (0.075) is larger than on CV (0.063 [0.048, 0.074]). On the 115 kinases scored by B1 in both, Cantley fell from CV 0.741 to TEST 0.717 (−0.025 ± 0.013), while the final model was flat (−0.003). Almost all of that fall is on Y kinases (0.607 → 0.496). Part of the TEST margin is therefore a weak-baseline draw; quote the CV gap alongside it.
- **B3/B4.** Say "our faithful retrained reimplementations under this benchmark's split, labels and negatives". For B4, add that 66.6% of its TEST pairs use NA-default site features.

## 4. CV→TEST gap for M1 (claim 3)

| | CV (K_cv) | CV on K_eval ∩ K_cv | TEST on the same kinases | paired Δ ± SE |
|---|---|---|---|---|
| final | 0.786 | 0.792 | 0.793 | +0.001 ± 0.012 |
| R2 | 0.776 | 0.786 | 0.778 | −0.007 ± 0.013 |
| R3 | 0.754 | 0.763 | 0.773 | +0.010 ± 0.012 |
| B1 (i) | 0.717 | 0.737 | 0.713 | −0.024 ± 0.013 |
| B3 | 0.677 (A24) | 0.715 | 0.709 | −0.007 ± 0.012 |

**Reading:**
- All movements are composition effects plus noise. B3's +0.029 "gain" is entirely composition: CV on the same kinases is 0.715.
- The final model shows no optimism, which is consistent with the audit's ≲ 0.005.
- M3 is equally stable (final: CV 0.358 on K_eval ∩ K_cv vs TEST 0.361).

**One red flag, in the residue split** (final, same kinases):

| | CV | TEST | paired Δ ± SE |
|---|---|---|---|
| ST (92 kinases) | 0.827 | 0.841 | +0.014 ± 0.012 |
| Y (23 kinases) | 0.667 | 0.595 | −0.072 ± 0.033 |

- On TEST, Y kinases are at the B0-pop floor: 0.595 vs 0.585.
- SYK (0.85 → 0.56), IGF1R (0.82 → 0.38) and LYN (0.68 → 0.49) collapse. Most Y kinases have only 3–4 TEST positives, so per-kinase values are noisy.
- The ST-dominated macro hides this. REPORT must state the Y result separately.

## 5. Strata (claim 5)

**How the table was built.**
- The batch strata are standalone (missing = 'worst'). There, B1's unscored kinases count as AUROC 0.5. This makes B1's STY row exactly 0.500 (TTK and MAP2K1, per A30) and deflates its bin-0/1–4 M4 (12 of 35 such K_any kinases are unscored).
- The table therefore restricts every row to **kinases B1 scores** (audit diagnostic; M1/M3/M4 are then exact matched values).
- M2 is the batch's standalone value (final competitors = 503 kinases), vs arm (ii).
- The B0-pop column is the batch's standalone value (all K_eval / K_any kinases in the stratum); it serves as a floor reference.
- **No stratum has a CI in the batch.** Differences under about 0.03, or with fewer than 5 kinases or 30 pairs, are not interpretable.

| stratum (TEST) | kin / pairs | final M1 | Cantley (i) M1 | B0-pop M1 (batch, standalone) | final M2 | Cantley (ii) M2 | final M3 | Cantley M3 | final ahead? |
|---|---|---|---|---|---|---|---|---|---|
| **all (matched)** | 116 | **0.794** | 0.719 | 0.578 | 0.372* | 0.194* | 0.370 | 0.165 | yes |
| ST kinases | 93 | 0.843 | 0.775 | 0.576 | 0.367 | 0.207 | 0.418 | 0.190 | yes |
| **Y kinases** | 23 | **0.595** | 0.496 | **0.585** | 0.362 | 0.140 | 0.174 | 0.060 | vs Cantley yes; **vs popularity floor no** |
| bin 5–19 (M1) | 13 | 0.759 | 0.646 | 0.639 | 0.208 | 0.165 | 0.262 | 0.095 | yes (small n; CV M2 said no) |
| bin 20–99 | 71 | 0.797 | 0.716 | 0.577 | 0.302 | 0.195 | 0.403 | 0.182 | yes |
| bin ≥ 100 | 31 | 0.795 | 0.748 | 0.542 | 0.568 | 0.177 | 0.318 | 0.127 | yes |
| M4, bin 0 (K_any) | 8 | 0.793 | 0.704 | | | | | | yes (8 kinases) |
| M4, bin 1–4 (K_any) | 15 | 0.756 | 0.634 | | | | | | yes (mostly 1–2 positives each) |
| M4, bin 5–19 (K_any) | 95 | 0.803 | 0.694 | | | | | | yes |
| priming kinases | 8 | 0.786 | 0.726 | 0.449 | 0.415 | 0.120 | 0.211 | 0.051 | yes |
| primed positives | 6 / 65 | 0.897 | 0.704 | 0.517 | 0.637 | 0.117 | 0.345 | 0.000 | vs arm (i) only; priming-aware arm not tested |
| unprimed positives (priming kinases) | 8 / 194 | 0.740 | 0.725 | 0.445 | 0.278 | 0.105 | 0.137 | 0.062 | **near-tie on M1** |
| PIKK (S/T-Q PIKK) | 4 (3) | 0.868 (0.903) | 0.796 (0.834) | 0.517 | 0.605 | 0.436 | 0.407 | 0.230 | yes (small) |
| hub sites (A26) | 84 / 270 | 0.862 | 0.714 | **0.788** | 0.368 | 0.205 | 0.458 | 0.174 | yes; mostly popularity |
| non-hub | 114 / 1,898 | 0.783 | 0.727 | 0.544 | 0.368 | 0.182 | 0.354 | 0.165 | yes |
| in vivo only | 85 / 523 | 0.791 | 0.721 | 0.589 | 0.401 | 0.171 | 0.379 | 0.157 | yes |
| in vitro only | 94 / 794 | 0.827 | 0.735 | 0.615 | 0.362 | 0.215 | 0.382 | 0.152 | yes |
| both | 105 / 851 | 0.777 | 0.722 | 0.539 | 0.362 | 0.179 | 0.359 | 0.178 | yes |
| new since 2021 | 60 / 249 | 0.771 | 0.713 | 0.556 | 0.348 | 0.155 | 0.357 | 0.181 | yes (smaller margin, 0.058) |
| pre-2021 | 115 / 1,919 | 0.799 | 0.725 | 0.580 | 0.379 | 0.202 | 0.375 | 0.162 | yes |
| external tier A (batch, unmatched) | 41 / 134 | 0.788 | 0.684 | 0.531 | | | | | yes |
| external tier B (batch, unmatched) | 19 / 49 | 0.804 | 0.661 | **0.725** | | | | | yes; popularity-heavy |
| Manning STE (M2) | 7 | 0.910 | 0.844 | | 0.244 | 0.224 | 0.600 | 0.270 | M2 near-tie |
| Manning TKL | 2 | 0.884 | 0.707 | | 0.155 | 0.292 (arm i) / 0.000 (arm ii) | 0.393 | 0.333 | M2 inconclusive (2 kinases) |
| **T3 (unseen family), all** | 116 | 0.755 | 0.719 | | | | | | yes, descriptive (no paired test) |
| T3, Y kinases | 23 | 0.531 | 0.496 | 0.585 | | | | | no skill |

\* Matched M2 from `superiority.tsv` (competitors restricted to the 386 B1-scored kinases).

**Summary.**
- **Where the final model beats Cantley clearly:**
  - every S/T stratum, evidence class, temporal set, external tier and training bin, on M1, M3 and M4;
  - M2 almost everywhere (≥ 100 labels: 0.57 vs 0.18).
- **Where it does not, or where the comparison is uninformative:**
  - Tyr kinases: both are at the popularity floor (§4).
  - Unprimed positives of priming kinases: a near-tie.
  - STE M2: a near-tie.
  - TKL M2, and the 0-label bin of K_eval (2 kinases; 1 matched): too small to read.
  - Primed sites against the *priming-aware* Cantley arm: not tested on TEST. CV says Cantley-primed is ahead.
  - **Hub-site gains** are largely popularity (B0-pop 0.788 on hubs).
- **Not delivered by the batch:**
  - The REPORT §7.5 row "Y, hub sites excluded" (only the hub/non-hub split of positives exists).
  - A16 calibration on the external/temporal strata (A29). Neither can be filled without another sealed access.

## 6. A17 assignment reliability on TEST (claim 6)

Batch `final_reliability.tsv` on labelled TEST sites (n = 1,458), ECE 0.026:

| tier | n | claimed | observed |
|---|---|---|---|
| [0, 0.2) | 261 | 0.12 | 0.14 |
| [0.2, 0.5) | 852 | 0.34 | 0.32 |
| [0.5, 0.6) | 236 | 0.563 | 0.564 |
| [0.6, 0.7) | 37 | 0.679 | 0.676 |
| [0.7, 0.8) | 72 | 0.755 | 0.694 |
| ≥ 0.8 | 0 | — | — |

The first two tiers merge the batch's 0.1-wide bins.

**By residue (audit diagnostic, same sites):**

| residue | ECE | ≥ 0.5: claimed / observed (n) | [0.2, 0.5): claimed / observed (n) |
|---|---|---|---|
| ST | 0.021 | 0.615 / 0.603 (345; Wilson [0.55, 0.65]) | — |
| Y | 0.069 | no site ≥ 0.5 (the cap is 0.445) | 0.41 / 0.34 (195) |

**Verdict.** The A31(b) isotonic map is **honest on TEST for S/T** and replicates the CV result (0.62 / 0.61).
- There is no ≥ 0.8 tier, and coverage@0.8 = 0 is the honest statement.
- Coverage@0.5 = 5,773 TEST sites (5,428 unlabelled), all S/T, and means about 60% top-1 precision **on labelled sites**. Under PU selection, that is not a precision on unlabelled sites.
- Y probabilities are mildly overconfident.

**Temporal and external subsets.** The batch reports ECE 0.122 (new since 2021) and 0.218 (external), with ≥ 0.5-tier accuracy of 0.36 and 0.21. Those numbers **must not be quoted as miscalibration.**
- `_site_truths` scores those subsets against only the new (or external) kinases of each site. A top call that is the site's older PSP kinase counts as wrong.
- With the site's full truth set (audit diagnostic), the ≥ 0.5 tier is 0.55 observed vs 0.63 claimed (n = 53, new since 2021) and 0.54 vs 0.61 (n = 28, external). ECE becomes 0.062 and 0.068.
- So there is mild overconfidence at most, and it is within the Wilson intervals.

**Missing from the batch:** A17 requires the random-score and Cantley-only coverage baselines next to coverage. Neither is in the batch, because only the final model supplied `assignment_prob`. The CV versions must be cited instead.

## 7. Surprises and claims to soften (claim 7), ranked

1. **M2 0.363 is not an improvement over CV.** It is K_eval selection (§2). Quote CV-on-K_eval 0.348 next to it. Keep the CV finding that understudied kinases are weak on M2.
2. **Tyr kinases show no specificity on TEST** beyond site popularity (0.595 vs B0-pop 0.585; −0.07 vs CV; T3 0.531). Cantley is at chance (0.496). "Beats Cantley on Y" is true but empty.
3. **The B4 M3 win is modest** (CI lower bound 0.010; p 0.018). Say "small".
4. **The margin over Cantley is partly a weak-baseline draw.** Cantley fell 0.025 from CV to TEST on the same kinases, almost entirely on Y kinases. Quote the CV margin (+0.063) alongside the TEST margin (+0.075).
5. **p-values.** "p = 0.0002" is the floor. Write "p < 0.0002".
6. **Priming.** The primed-site advantage is only against non-priming-aware Cantley. On unprimed positives it is a tie.
7. **Hub-site** and **external tier B** gains are mostly popularity.
8. **Temporal/external A17 ECEs** are a truth-definition artefact (§6). Do not report them as miscalibration.
9. **REPORT count drift.** REPORT §7.6 quotes "263 temporal positives" and "tier A 158". The batch evaluates 262 temporal pairs (K_any) and 134 tier-A pairs. Use the batch counts.
10. **Gaps the batch cannot fill:** hub-excluded Y, A16 calibration on the temporal/external strata, A17 coverage baselines on TEST, and a paired T3 test. They must be listed as not done, not back-filled.
