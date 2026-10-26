# CV strata report (pooled OOF, frozen ksp.eval with A24)

Owner: strata analyst. Generated 2026-10-26 22:10 CET.

Code:
- `src/ksp/strata.py`: canonical definitions.
- `scripts/strata_definitions.py`: writes the definition tables.
- `scripts/strata_cv_report.py`: CV report.
- `scripts/strata_priming_diag.py`: priming split.
- `tests/test_strata.py`: 10 toy tests, all passing.

Results are in `results/strata/`:

| file | contents |
|---|---|
| `cv_strata.tsv` | every stratum × run: n, M1/M2/M3/M5 with CIs |
| `cv_strata_tables.md` | every table, every run |
| `cv_strata_site_per_kinase.tsv` | per-kinase values for the site strata |
| `cv_priming_diag.tsv` | priming split per kinase |
| `cv_strata_sanity.tsv` | sanity checks |
| `kinase_strata.tsv`, `stratum_definitions.tsv`, `cantley_phospho_preference.tsv`, `dev_hub_sites.tsv`, `dev_site_strata_counts.tsv` | definition tables |

## 1. Definitions (fixed before any model metric was looked at)

All definitions live in `ksp.strata`, and every input is DEV-safe:
- the kinase index;
- the site index (window15, U_ms);
- kinase-library 1.8.0 matrices;
- DEV L2 counts from `splits.load_eval_kinase_sets().n_dev_L2`.

Label-defined site strata (hub, evidence) take the label table as an argument, so `eval_test` can pass the TEST labels it already holds. The module never reads labels itself.

- **Kinase group, residue class, training bin (0/1–4/5–19/20–99/≥100).** These are identical to the `ksp.eval` attributes, and a test asserts the equality.
  - On CV, bin 0 is empty, because CV positives are DEV labels.
  - The 1–4 bin has no K_cv kinases (K_cv requires ≥5), so it is reported only as M5 over K_any.
- **Priming, PRIMARY (curated).** These are DESIGN §9, identical to `ksp.eval` 'priming'.
  - Kinases and their priming offsets:

    | kinases | priming offset |
    |---|---|
    | GSK3A/B | +4 (Fiol 1987; Frame & Cohen 2001) |
    | Manning family CK1: CSNK1A1/A1L/D/E/G1/G2/G3 | −3 (Flotow 1990) |
    | CK2α: CSNK2A1/A2/A3 | +1 or +3, where pS replaces the acidic determinants (Meggio & Pinna 2003) |
    | PLK1 | none: PBD docking on a distal S-pT-P (Elia 2003) |

  - That is 13 index kinases, 9 of them in K_cv.
- **Priming, SECONDARY (data-driven, `priming_cantley`).**
  - Rule: in the S/T Cantley log2 matrix, at some position −5..+4, pS or pT is the single top-scoring residue **and** its log2 is ≥ 2 (≥4-fold).
  - The offsets are the positions that qualify.
  - This gives 25 index kinases, 13 of them in K_cv. Besides CK1 (all 7), GSK3A/B and CSNK2A1/A2, these include GRK2/3/4/5, CAMK2A, UHMK1, CLK2, DYRK4, BMPR1B, ACVR2A/B, STK17A and STK32B/C (YANK2/3). PLK1 and CSNK2A3 do not qualify.
  - Tyr matrices are excluded: many TK matrices favour pY at −1/+1/+2, which reflects adjacent multisite phosphorylation, not priming.
  - Curated is primary because DESIGN §9 pre-registered it.
- **PIKK.**
  - All six PIKKs are in the index. TRRAP is a pseudokinase with 0 labels, and SMG1 has 4 DEV positives, so SMG1 is in K_any but not K_cv.
  - K_cv contains ATM, ATR, PRKDC and MTOR.
  - `pikk_sq` excludes MTOR, which has no S/T-Q motif.
- **Secretory kinases:** FAM20C, FAM20A (pseudokinase, 0 labels) and PKDCC (0 labels). In practice this stratum is **FAM20C alone, with 6 positives.**
- **Hub sites (A26).**
  - A hub site has evaluable-positive labels to ≥5 mutually non-paralogous kinases, counted as connected components of the Par(k) graph.
  - DEV has 88 hub sites, 84 in the universe, carrying 719 positive pairs.
  - Two views are reported. "hub/non_hub" splits the positives. "hub_excluded" re-evaluates with hub positives removed **and** hub sites removed from E-hard.
- **Primed sites.**
  - A positive of a priming kinase is primed when its offset residue is an S/T that is itself a U_ms site.
  - DEV counts: 234 primed vs 824 unprimed positives. The unprimed count includes all 213 PLK1 positives.
  - A control applies the same rule (−3 or +4) to non-priming S/T kinases.
- **Size flag (pre-declared).** A stratum is "small" if it has fewer than 5 kinases with defined M1, or fewer than 30 positives.

## 2. Method

**Runs.** These are the leaderboard runs, re-evaluated:

| label | run |
|---|---|
| R2 | `r2_final_knn_s0_cv5__20261026-181711__s0` |
| R3 | `r3_final_distill_cv_pooled__20261026-173414__s0` |
| B1i | `simple_b1_published_reeval_…155352` |
| B1ii | `simple_b1_universe_reeval_…155703` |
| B1prim | `simple_b1_primed_hard_reeval_…155750`, Cantley with hard priming, added because it is the priming-aware arm |
| B2 | `simple_b2_pssm_reeval` |
| B4 | `b4_selphi2_noleak` |
| B7 | `simple_b7_lr_reeval` |
| B0 | `simple_b0_{prior,pop,miss}_reeval` |

Why this R2 run:
- It is the latest R2 pooled 5-fold A24 run listed on LEADERBOARD (generated 18:54).
- `r2_final_ens5` (19:12) is not on the board.
- `r2_alarm_compknow` scores higher but is an alarm/control run.

**Evaluation path.**
1. Load the run's predictions.
2. Apply `E.map_predictions` and then `E.fold_rank_normalise` (A24), blocked by kinase, which gives identical groups.
3. Run `E.evaluate(..., fold_normalise=False, missing='worst')`.

No frozen file was edited.

**Kinase strata.** These are per-kinase metrics, macro-averaged:
- M1, M2 and M3 over K_cv ∩ stratum;
- M5 over K_any ∩ stratum, following the leaderboard "M5bin" convention.

**Site strata.** These subset the positives against unchanged E-hard negatives.
- This is identical to re-evaluating with subset labels, because removed positives stay masked through KP(k).
- Against `evaluate()` on subset labels, the maximum difference is 0 for M1 and Hit@5.
- M3 uses the full ranking's top-1% cut-off. A subset re-evaluation differs by at most 0.002 macro.

**CIs.** 95% CIs come from `E._boot_matrix`: B = 1000 substrate-cluster resamples, seed 0. The 'all' row equals `E.bootstrap` to 4e-16.

**A27 matched-kinase view.** R2 is restricted to the kinases B1 scores (199/203 of K_cv).

**Runtime and memory.**
- Runtime was 9,000 worker-seconds over 13 evaluations, about 2 h wall-clock (20:04–22:00) on at most 3 workers. Per run it was 6.5–19 min, of which evaluate took 2–3.5 min.
- Peak memory was about 13 GB per worker.
- The node's job cgroup (128 GB, shared) reached 135 GB at 20:50, with one OOM kill of another process. I then dropped to 2 workers.

## 3. Sanity checks (`cv_strata_sanity.tsv`)

| check | result |
|---|---|
| full-run M1/M2/M3 vs leaderboard `metrics_a24.json`, 11 runs | difference **0** for every run |
| partition recombination (338 checks): kinase partitions, n-weighted M1/M5 | all ≤ 2e-16 |
| partition recombination: site partitions per kinase (evidence, hub, primed, primed_cantley), M1/M2/M3 | all ≤ 2e-16 |
| B0-prior | M1 = 0.500 exactly, overall and in every stratum |
| B0-pop / B0-miss overall | 0.552 / 0.539 |
| B0-pop / B0-miss, non-small strata | range 0.41–0.73 |

The B0-pop range is driven by label-free popularity, not a harness error:
- On hub sites, B0-pop reaches 0.73, so hub-site AUROCs are partly a popularity effect for every method.
- On primed sites it is 0.41.

## 4. Key results (M1 = macro E-hard AUROC, K_cv; R2 with 95% CI)

| stratum | n kin | n pos | R2 | R3 | B1i | B1prim | B7 | B0pop |
|---|---|---|---|---|---|---|---|---|
| all | 203 | 8416 | **0.774** [0.756, 0.786] | 0.754 | 0.717 | 0.714 | 0.726 | 0.552 |
| ST kinases | 161 | 7228 | 0.808 [0.793, 0.822] | 0.783 | 0.751 | 0.747 | 0.763 | 0.538 |
| Y kinases | 38 | 1142 | **0.629** [0.580, 0.659] | 0.624 | 0.594 | 0.594 | 0.569 | 0.592 |
| Y, hub-excluded | 38 | 1114 | 0.627 [0.577, 0.662] | 0.623 | 0.586 | 0.589 | 0.560 | 0.583 |
| 1–4 labels (M5, K_any) | 79 | 147 | 0.741 [0.690, 0.787] | 0.647 | 0.624 | 0.628 | 0.651 | 0.565 |
| 5–19 labels | 83 | 669 | 0.758 [0.722, 0.782] | 0.748 | 0.689 | 0.683 | 0.705 | 0.571 |
| ≥100 labels | 31 | 5400 | 0.815 [0.803, 0.828] | 0.786 | 0.766 | 0.768 | 0.775 | 0.543 |
| priming kinases | 9 | 1049 | 0.837 [0.812, 0.857] | 0.795 | 0.749 | 0.792 | 0.812 | 0.474 |
| primed positives | 8 | 231 | 0.870 [0.826, 0.910] | 0.827 | 0.662 | **0.945** [0.924, 0.962] | 0.832 | 0.412 |
| unprimed positives (priming kinases) | 9 | 818 | 0.817 [0.787, 0.839] | 0.781 | 0.775 | 0.729 | 0.802 | 0.490 |
| PIKK (small) | 4 | 449 | 0.880 [0.853, 0.907] | 0.846 | 0.833 | 0.822 | 0.876 | 0.511 |
| FAM20C (1 kinase) | 1 | 6 | 0.980 | 0.931 | 0.974 | 0.967 | 0.882 | 0.327 |
| in vivo only | 172 | 2080 | 0.766 [0.743, 0.782] | 0.768 | 0.707 | 0.704 | 0.720 | 0.583 |
| in vitro only | 181 | 3292 | 0.781 [0.762, 0.802] | 0.750 | 0.729 | 0.720 | 0.734 | 0.529 |
| hub positives | 151 | 664 | 0.801 [0.774, 0.827] | 0.804 | 0.769 | 0.757 | 0.787 | 0.732 |

M2 (Hit@5) by training bin, K_cv:

| bin | R2 | R3 | B1ii | B7 | B0prior |
|---|---|---|---|---|---|
| 5–19 | **0.099** [0.069, 0.122] | 0.092 | 0.162 | 0.121 | 0.000 |
| 20–99 | 0.164 | 0.216 | 0.171 | 0.158 | — |
| ≥100 | 0.410 | 0.592 | 0.193 | 0.507 | 0.367 |

Under the A27 matched-kinase view, R2 scores 0.119 and B1ii 0.170 on the 5–19 bin.

B2 and B4 are in `cv_strata_tables.md`:
- B4 has 33% pair coverage, ranked worst, so its numbers mostly measure coverage.
- B1i and B1ii are identical on M1/M3 (A23).

## 5. Answers

**Priming kinases and primed sites.**

R2 handles priming kinases well:
- M1 0.837 against 0.771 for non-priming kinases.
- It is ahead of published Cantley (0.749).
- Within each kinase, primed positives score higher than unprimed ones: +0.078 on average, with 6 of 7 kinases positive.
- The non-priming control shows −0.005, so this is not a generic effect of phospho-dense regions.

There is still a clear, actionable gap:
- On primed positives, the **priming-aware Cantley arm B1prim reaches 0.945 [0.924, 0.962] against R2's 0.870 [0.826, 0.910]**. The CIs are disjoint.
- B1prim is ahead for every large priming kinase (`cv_priming_diag.tsv`):

  | kinase | B1prim | R2 |
  |---|---|---|
  | GSK3B (102 primed) | 0.939 | 0.883 |
  | GSK3A | 0.948 | 0.835 |
  | CSNK1A1 | 0.951 | 0.899 |
  | CSNK1D | 0.895 | 0.761 |
  | CSNK1E | 0.967 | 0.946 |
  | CSNK2A1 | 0.947 | 0.887 |

- On unprimed positives B1prim is worse (0.729 against 0.817).
- So R2 carries the primed Cantley scores (F1 primed_hard/soft) but does not use them selectively. A global pair model has to learn that the primed score matters only for priming kinases at their own offset.

Unprimed GSK3B positives are hard for every method: R2 scores 0.684 on them. This holds both when the +4 residue is not S/T (106 positives, 0.679) and when it is S/T but never observed (22 positives, 0.708).

**PIKK and FAM20C.**
- PIKKs are among the best-handled kinases. R2 M1 is 0.880, and 0.926 for the S/T-Q PIKKs, ahead of B1 (0.833). M2 is 0.580 and M3 0.512.
- MTOR is the weak member: M1 0.742, Hit@5 0.265.
- FAM20C (6 positives, 1 kinase) is anecdotal: M1 0.980.
- Both strata are flagged small.

**Understudied kinases.**
- On M1 and M5, R2 holds up well in the low-label bins:

  | bin | R2 | R3 | Cantley |
  |---|---|---|---|
  | 1–4 (M5) | 0.741 | 0.647 | 0.624 |
  | 5–19 (M1) | 0.758 | — | 0.689 |

- **On M2, R2 is weak for understudied kinases.** In the 5–19 bin its Hit@5 is 0.099 against 0.162 for Cantley (CIs disjoint), while in the ≥100 bin it is 0.410.
- This means the cross-kinase ranking at a site is dominated by well-studied kinases. The B0-prior floor alone reaches 0.367 in the ≥100 bin. This is the "dominated by CDK1/PKA" failure mode.

**Tyr vs S/T.**
- Y kinases are much worse than S/T kinases: 0.629 against 0.808.
- Y kinases are barely above the B0-pop floor (0.592) for every method.
- Excluding hub sites changes nothing (0.627), so for the models the Y gap is not a hub artefact.
- M2 for Y is similar to S/T (0.166 against 0.177).

**Evidence class.** In-vivo-only positives are slightly harder for R2 (0.766 against 0.781). The CIs overlap.

## 6. Suggested fixes (for the R2 owner; nothing was implemented)

1. **Priming.** Add the kinase-conditioned priming-offset features that DESIGN §9 already lists. I found no implementation in `src/ksp/models/r2`:
   - for kinase k, whether the residue at k's own offset (GSK3 +4; CK1 −3; CK2 +1/+3; or the `priming_cantley` offsets) is S/T and in U_ms;
   - the Cantley primed_hard − unprimed delta for k's own matrix;
   - alternatively, gate F1 to primed_hard for `priming`/`priming_cantley` kinases.

   The expected headroom on primed positives is about +0.07 AUROC, which is the B1prim − R2 gap.
2. **M2 for understudied kinases.** Supply an A20 `score_site` whose per-kinase prior is removed, for example a per-kinase calibrated assignment probability with the log training-count offset removed. Alternatively, blend toward the Cantley percentile for kinases with fewer than 20 labels. Cantley's flat M2 across bins (0.16–0.19) shows what a prior-free score achieves.
3. **Tyr.** Hubs are not the cause. The features_tyr diagnosis (atlas anti-correlation for FGFR/Eph/VEGFR/TAM) stands, so nearest-atlas imputation remains the candidate.

## 7. Using `ksp.strata` on TEST (for `eval_test`)

- **Kinase strata.** Use `kinase_table()` / `kinase_sets()` and subset `res.per_kinase`.
- **Site strata.**
  - `primed_mask(pos[kinase_acc, site_id], priming_curated())`;
  - `hub_sites(test_eval_positive_rows, labels.paralog_sets())`;
  - `evidence_class(pos)`.
- Then subset `res.items` / `res.site_level`, or pass subset positives to `evaluate()`. `subset_result()` in `scripts/strata_cv_report.py` is a worked example.

## 8. Caveats

- CV training bins use the full DEV L2 count. Each fold actually trains on about 80% of it.
- M3 for site strata uses the full-ranking cut-off. This is documented above and is exact for recombination.
- M2 and M3 bootstraps keep the cut-off fixed, as `ksp.eval` does.
- B1 leaves 4 K_cv kinases unscored, which receive 0.5 under the 'worst' rule. See the matched view.
- Unrelated to this work: `tests/test_no_test_access.py` currently fails on `src/ksp/models/regime/variants.py` (pattern 'omnipath'), which belongs to another owner.

## 9. A32 follow-ups (run 22:18–23:05; rules in DESIGN §17 A32, applied literally)

**(b) Prior-debiased score_site for M2** (`scripts/models/m2_prior_debias.py`, `results/models/priming/m2_debias/`).
- Base: OOF fusion `r2_fusion_r3__20261026-191901__s0` (score_site = w 0.5 z-blend of the R2 5-seed ensemble and R3).
- Prior: log1p of the DEV L2 training count of k in each fold.
- Scale check, fixed before any M2 was computed (`scale_check.json`):
  - within-site SD of score_site, median over positive sites, is 0.83;
  - SD of the prior is 1.50, a ratio of 1.8;
  - the pre-declared window was 0.25–4, so log1p was used as is, with no rescaling.
- M2 was computed with the frozen `E._site_centric` + `_endpoints` logic. At β = 0 it matches the run's metrics.json exactly.

| β | 0 | 0.1 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| pooled M2 | **0.263** | 0.236 | 0.165 | 0.069 | 0.029 | 0.009 |
| M2, 5–19 bin | 0.142 | 0.146 | 0.115 | 0.079 | 0.051 | 0.023 |
| M2, ≥100 bin | 0.621 | 0.526 | 0.321 | 0.075 | 0.014 | 0.001 |

- β* = 0 by pooled M2.
- Nested leave-one-fold-out selection picks β = 0 in all 5 folds; held-out M2 is 0.269.
- At β = 0.1, pooled M2 falls by 0.027, more than the 0.003 tolerance.
- **Decision: not adopted (β = 0).** The global log-count prior is too blunt. Even β = 0.1 costs well-studied kinases far more than understudied ones gain.

**(a) Priming features.**
- Files: `scripts/models/priming_r2.py`, `priming_compare.py`, `results/models/priming/priming_compare.json`.
- Run: `priming_r2_s0_cv5__20261026-224631__s0`, logged in REGISTRY.
- Pre-check: F1 **already contains** the priming-aware Cantley percentiles (`f1_primed_hard_pct_universe` and `f1_primed_soft_pct_universe`), so the experiment reduces to three added columns:
  - `fp_offset_obs`: the residue at k's best pS/pT offset is an S/T U_ms site;
  - `fp_offset_score`: preference-weighted observed flanks;
  - `fp_strength`: the kinase's maximum pS/pT log2.
- Comparison against `r2_final_s0_cv5`, seed 0, with paired cluster bootstrap B = 1000:

| criterion | result | rule met? |
|---|---|---|
| overall M1 Δ | +0.0037 [−0.001, 0.009] | yes (≥ −0.002) |
| M2 Δ | +0.010 [−0.000, 0.023] | yes (≥ −0.005) |
| M4 Δ (not a criterion) | +0.006 [0.001, 0.010] | — |
| canary | 0.566 | passes, narrowly (≤ 0.57) |
| primed-positives M1 Δ | **−0.006 [−0.022, +0.008]**: 0.893 → 0.888, 8 kinases, 231 positives | **no** (rule is ≥ +0.02 with CI > 0) |
| primed-positives M1 Δ, Cantley-defined offsets | +0.009 [−0.006, 0.018] | not a criterion |

- The new columns rank 86–167 of 175 by gain. The model already exploits priming through F1 primed_hard/soft, which rank 11–53.
- **Decision: not adopted.**
- Note: the baseline r2_final_s0 already scores 0.893 on primed positives, against 0.870 for the knn run in §4. The gap to B1prim (0.945) is therefore smaller than §4 suggests.
