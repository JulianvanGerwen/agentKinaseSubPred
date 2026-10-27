## 7. TEST results (single batch, FREEZE 721dc197; attempt 1; COMPLETED 2026-10-27 00:57:29)

**About the batch.**
- One sealed access. The batch took 2,838 s for 12 configurations at B = 10,000, with harness 974d520 (`results/test/LEDGER.md`).
- Outputs are in `results/test/721dc1977366/`.
- **Independent review:** `docs/audit/audit_test.md`.
  - Every endpoint, stratum and reliability table reproduces exactly (max |diff| = 0).
  - All 17 prediction pins match FREEZE.
  - TEST and DEV share 0 proteins and 0 clusters.
  - The final score columns are an exact label-free function of the pinned R2 and R3 files.
- **Wording rules.** The review's verdicts and caveats bind the wording below.
- **Audit-diagnostic values.** Values marked "(audit diagnostic)" are the review's re-cuts of batch quantities: kinase-resampling, matched subsets, corrected truth sets. They are not TEST endpoints.

### 7.1 T1: seen kinase, unseen protein (K_eval = 119, K_any = 256)

Standalone metrics use missing = 'worst' and are not head-to-heads; see §7.2 for those. CIs are cluster bootstrap, B = 1,000.

| Config | M1 [95% CI] | M2 Hit@5 [CI] | M3 [CI] | M4 | Kinases scored | Train s | Infer s |
|---|---|---|---|---|---|---|---|
| **final (R2+R3)** | **0.794** [0.774, 0.811] | **0.363** [0.324, 0.397] (score_site) | **0.363** [0.309, 0.401] | **0.799** | 503 | 1,136.6 | 316.9 |
| R2 alone | 0.778 [0.758, 0.795] | 0.264 [0.224, 0.293] | 0.312 [0.261, 0.348] | 0.781 | 503 | 417.9 | 298.9 |
| R3 alone | 0.775 [0.752, 0.793] | 0.322 [0.290, 0.362] | 0.292 [0.242, 0.330] | 0.769 | 503 | 718.7 | 3.1 |
| B7 simple-strong | 0.725 [0.700, 0.751] | 0.243 | 0.200 | 0.717 | 503 | 89.1 | 287.0 |
| B1 (i) Cantley as published | 0.714 [0.685, 0.743] | 0.163 | 0.161 [0.129, 0.194] | 0.699 | 386 | 0 | 4.4 |
| B1 (ii) universe percentile | 0.714 | 0.189 [0.155, 0.214] | 0.161 | 0.699 | 386 | 0 | 1.0 |
| B3 Phosformer (our retrain) | 0.708 [0.686, 0.733] | 0.163 [0.145, 0.190] | 0.104 [0.077, 0.129] | 0.688 | 503 | 3,152.4 | 6,478.8 |
| B4 SELPHI2.0 (our retrain; 66.6% of pairs imputed) | 0.705 [0.673, 0.733] | 0.233 [0.204, 0.266] | 0.304 [0.233, 0.351] | 0.701 | 503 | 737.0 | 2,043.4 |
| B2 learned PSSM | 0.685 [0.647, 0.717] | 0.155 | 0.138 | 0.664 | 503 | 3.3 | 39.7 |
| B0-pop (site popularity) | 0.578 [0.543, 0.608] | 0.020 | 0.127 | 0.588 | 503 | 0 | 0.4 |
| B0-miss | 0.565 | 0.037 | 0.066 | 0.579 | 503 | 19.0 | 30.4 |
| B0-prior | 0.500 | 0.109 | 0.010 | 0.500 | 503 | 0 | 0.8 |

Sources: `results/test/LEDGER.md` and `results/test/721dc1977366/<config>_ci.tsv`.

**Coverage and runtime.**
- B1 leaves 1,199,434 of 13,230,299 required pairs unscored (117 kinases).
- Runtimes are the FREEZE manifest values; B3's training time excludes its about 4.3 h of MLM pretraining.

**Kinase-resampling CIs for final** (`final_ci_kinases.tsv`): M1 [0.767, 0.821], M2 [0.310, 0.414], M3 [0.318, 0.407], M4 [0.775, 0.820].

### 7.2 The 9 primary tests (paired cluster bootstrap B = 10,000, matched universe, Holm over 9)

| # | Comparison | Endpoint | Δ [95% CI] | p (bootstrap) | Holm p | Beats | Kinase-resampling CI (audit diagnostic) | Kinases favouring A / B (audit) |
|---|---|---|---|---|---|---|---|---|
| 1 | final − B1 (i) | M1 | +0.075 [0.046, 0.099] | < 0.0002 (bootstrap floor) | 0.0018 | yes | [0.057, 0.094] | 97 / 19 |
| 2 | final − B1 (i) | M3 | +0.205 [0.165, 0.235] | < 0.0002 (floor) | 0.0018 | yes | [0.167, 0.243] | 78 / 4 |
| 3 | final − B1 (ii) | M2 | +0.179 [0.147, 0.213] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 4 | final − B3 | M1 | +0.086 [0.059, 0.105] | < 0.0002 (floor) | 0.0018 | yes | [0.065, 0.110] | 97 / 22 |
| 5 | final − B3 | M2 | +0.199 [0.151, 0.229] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 6 | final − B3 | M3 | +0.259 [0.213, 0.291] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 7 | final − B4 | M1 | +0.089 [0.063, 0.117] | < 0.0002 (floor) | 0.0018 | yes | [0.069, 0.111] | 98 / 20 |
| 8 | final − B4 | M2 | +0.130 [0.091, 0.156] | < 0.0002 (floor) | 0.0018 | yes | [0.089, 0.172] | 70 / 16 |
| 9 | final − B4 | M3 | **+0.059 [0.010, 0.120], modest** | 0.0178 | 0.0178 | yes | [0.021, 0.097] | 52 / 27 (40 ties) |

Sources: `results/test/721dc1977366/superiority.tsv`; `docs/audit/audit_test.md` §3. The audit re-ran tests 1, 2 and 7–9 at B = 10,000 and reproduced the batch to the last digit.

**Matched universe.**
- Against B1: 386 of 503 kinases, 74.1% of pairs, and 116 K_eval kinases for M1/M3. TTK, MAP2K1 and LIMK2 are not scored by B1.
- M2 against B1 is matched-kinase (A27); there the final model's M2 is 0.372.
- B3 and B4 are exhaustive (503 kinases).

**Required wording.**
- B1 is Cantley arms (i) and (ii). The priming-aware arm was not a TEST comparator.
- B3 and B4 are "our faithful retrained reimplementations under this benchmark's split, labels and negatives".
  - B3's pipeline reproduces the published AUROC on its own split (0.888 vs 0.885; `results/baselines/phosformer/fidelity_check.json`).
  - 66.6% of B4's TEST pairs use NA-default site features. The fold-0 CV check puts the imputation cost at about 0.013 M1, which does not explain the gap.
- The B4 M3 win (#9) is **modest**: it is the only test off the p floor, and its CI lower bound is 0.010.

### 7.3 Secondary comparisons (Holm over 17, separately; `results/test/721dc1977366/secondary.tsv`)

| Comparison | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| final − R2 | +0.016 [0.009, 0.022] | +0.099 [0.073, 0.128] | +0.051 [0.034, 0.071] | +0.017 [0.007, 0.025] |
| final − R3 | +0.019 [0.008, 0.030] | +0.041 [0.007, 0.060] (Holm p 0.013; weakest) | +0.071 [0.041, 0.097] | +0.030 [0.017, 0.040] |
| final − B7 | +0.069 [0.044, 0.088] | +0.120 [0.083, 0.146] | +0.163 [0.123, 0.198] | +0.082 [0.051, 0.098] |
| final − B2 | +0.109 [0.076, 0.140] | +0.208 [0.168, 0.239] | +0.225 [0.171, 0.263] | +0.135 [0.098, 0.157] |
| final − B1 (i) | — | — | — | +0.087 [0.057, 0.103] |

All 17 meet the "beats" rule; Holm p ≤ 0.0034 except final − R3 M2 (0.013).

**Fusion helps on TEST.** Both components are beaten on every endpoint, which replicates the CV choice of late fusion.

### 7.4 T3: unseen kinase family and unseen protein

Each kinase is scored by the model trained without its LKO family (`results/test/LEDGER.md`; `<config>_t3.json`).

| Config | T3 M1 | T3 M3 | T3 M4 | Missing kinases |
|---|---|---|---|---|
| **final** | **0.752** | 0.285 | 0.753 | 0 |
| R2 | 0.744 | 0.247 | 0.746 | 0 |
| B7 | 0.708 | 0.171 | 0.699 | 0 |
| R3 | 0.703 | 0.209 | 0.703 | 0 |
| B2 | 0.616 | 0.072 | 0.609 | 0 |

T3 M2 is not reported, because it mixes LKO models.

**Reading.**
- **Consistent with CV.** CV LKO-T2 for the fusion was 0.749 (`docs/model_r2.md` §6.3).
- **Relative to Cantley.** The final T3 M1 sits above label-free Cantley's T1 M1 (0.714). On B1-scored kinases (audit diagnostic) it is 0.755 vs 0.719, with 79 of 116 kinases ahead. There is **no paired T3 test**, so this is descriptive.
- **Y kinases.** T3 on Y kinases is 0.531 (23 kinases; `final_t3_per_kinase.tsv`), below the B0-pop floor of 0.585. There is no skill there.

### 7.5 TEST strata (M1 macro E-hard over K_eval ∩ stratum; standalone; no CIs in the batch, so descriptive)

| Stratum | Kinases / pairs | final | B1 (i) | B3 | B4 | B0-pop (floor) |
|---|---|---|---|---|---|---|
| all | 119 | 0.794 | 0.714 | 0.708 | 0.705 | 0.578 |
| ST kinases | 94 | 0.841 | 0.772 | 0.751 | 0.736 | 0.576 |
| **Y kinases** | 23 | **0.595** | 0.496 | 0.534 | 0.575 | **0.585** |
| STY | 2 | 0.865 | 0.500† | 0.695 | 0.733 | 0.604 |
| bin 5–19 | 13 | 0.759 | 0.646 | 0.610 | 0.694 | 0.639 |
| bin 20–99 | 73 | 0.799 | 0.710 | 0.706 | 0.701 | 0.577 |
| bin ≥100 | 31 | 0.795 | 0.748 | 0.757 | 0.703 | 0.542 |
| M4, bin 0 (K_any) | 16 | 0.807 | 0.602† | 0.508 | 0.676 | 0.677 |
| M4, bin 1–4 (K_any) | 19 | 0.757 | 0.606† | 0.712 | 0.636 | 0.548 |
| priming kinases | 8 | 0.786 | 0.726 | 0.740 | 0.635 | 0.449 |
| primed positives | 6 / 65 | 0.897 | 0.704 | 0.788 | 0.665 | 0.517 |
| unprimed positives (priming kinases) | 8 / 194 | 0.740 | 0.725 (near-tie) | 0.729 | 0.644 | 0.445 |
| PIKK (S/T-Q PIKK) | 4 (3) | 0.868 (0.903) | 0.796 (0.834) | 0.859 | 0.739 | 0.517 |
| FAM20C (M4, K_any) | 1 | 0.999 | 0.969 | 0.855 | 0.766 | 0.371 |
| hub sites (A26) | 85 / 272 | 0.863 | 0.712 | 0.731 | 0.802 | **0.788** |
| non-hub | 117 / 1,932 | 0.784 | 0.721 | 0.707 | 0.682 | 0.544 |
| in vivo only | 88 / 535 | 0.788 | 0.713 | 0.713 | 0.706 | 0.589 |
| in vitro only | 96 / 807 | 0.827 | 0.730 | 0.715 | 0.751 | 0.615 |
| both | 107 / 862 | 0.779 | 0.718 | 0.701 | 0.676 | 0.539 |

Source: `results/test/721dc1977366/{final,B1_i,B3,B4,B0_pop}_strata.tsv`.

† Standalone B1 counts its unscored kinases as AUROC 0.5. This deflates the STY row (TTK and MAP2K1, per A30) and the bin-0 and bin-1–4 rows. Restricted to B1-scored kinases (audit diagnostic, `audit_test.md` §5), the values are M4 bin 0: 0.793 vs 0.704, and M4 bin 1–4: 0.756 vs 0.634.

**M1 by Manning group.**

| Group | Kinases | final | B1 (i) |
|---|---|---|---|
| AGC | 24 | 0.816 | 0.741 |
| CAMK | 15 | 0.829 | 0.774 |
| CMGC | 23 | 0.880 | 0.823 |
| STE | 8 | 0.913 | 0.801 |
| Other | 15 | 0.804 | 0.721 |
| CK1 | 4 | 0.817 | 0.725 |
| Atypical | 4 | 0.868 | 0.796 |
| TKL | 3 | 0.819 | 0.638 |
| TK | 23 | 0.595 | 0.496 |

**M2 (Hit@5) by training bin, final vs B1 (ii).**

| Bin | final | B1 (ii) |
|---|---|---|
| 5–19 | 0.208 | 0.165 |
| 20–99 | 0.303 | 0.190 |
| ≥100 | 0.568 | 0.177 |

These are selected subsets: the 13 TEST kinases with 5–19 labels are the best-labelled of that bin. The CV M2 weakness for understudied kinases is not refuted (audit §2).

**Where the final model does not clearly beat Cantley, or the comparison is uninformative:**
- Tyr kinases: both sit at or near the popularity floor.
- Unprimed positives of priming kinases: a near-tie.
- STE M2: a near-tie (0.255 vs 0.196, arm ii).
- TKL M2: 3 kinases, uninformative.
- Primed sites against the *priming-aware* Cantley arm: not tested.

Hub-site gains are mostly popularity (B0-pop 0.788 on hub positives).

### 7.6 Temporal and external curated sets

| Set | K_eval kinases / pairs | final | R2 | R3 | B1 (i) | B3 | B4 | B0-pop |
|---|---|---|---|---|---|---|---|---|
| new since 2021 | 61 / 262 | **0.766** | — | — | 0.710 | 0.660 | 0.696 | 0.556 |
| pre-2021 | 118 / 1,942 | 0.800 | — | — | 0.719 | 0.721 | 0.711 | 0.580 |
| external tier A | 41 / 134 | **0.788** | 0.771 | 0.766 | 0.684 | 0.721 | 0.679 | 0.531 |
| external tier B (popularity-heavy) | 19 / 49 | 0.804 | 0.752 | 0.842 | 0.661 | 0.778 | 0.754 | **0.724** |
| external, all (M1 / M4) | 49 / 130 (K_any 75 / 183) | 0.800 / 0.800 | | | 0.699 / 0.691 | 0.734 / 0.710 | 0.696 / 0.726 | 0.577 / 0.589 |

Sources: `<config>_strata.tsv` (temporal, external all) and `<config>_ext_tiers.tsv`.

- **Discoveries since 2021:** the margin over Cantley is smaller (+0.056) but still positive. M3 is 0.351 against 0.178.
- **Count drift:** the batch evaluates 262 temporal pairs and 134 tier-A pairs. The earlier draft quoted 263 and 158, from the split report's pre-filter counts (audit §7 #9).

### 7.7 A17 assignment reliability on labelled TEST sites (final model only)

Overall: n = 1,458 sites, top-1 accuracy 0.353, ECE 0.026 (`final_reliability.tsv`).

| p_top tier | n | claimed | observed |
|---|---|---|---|
| [0, 0.2) | 261 | 0.12 | 0.14 |
| [0.2, 0.5) | 852 | 0.34 | 0.32 |
| [0.5, 0.6) | 236 | 0.563 | 0.564 |
| [0.6, 0.7) | 37 | 0.679 | 0.676 |
| [0.7, 0.8) | 72 | 0.755 | 0.694 |
| ≥ 0.8 | 0 | — | — |

By residue (audit diagnostic, `audit_test.md` §6):
- **ST:** ECE 0.021; the ≥0.5 tier is 0.615 claimed vs 0.603 observed (n = 345; Wilson [0.55, 0.65]).
- **Y:** ECE 0.069; no site can exceed 0.445; the [0.2, 0.5) tier is 0.41 claimed vs 0.34 observed (n = 195), mildly overconfident.

**Coverage over TEST candidate sites.** 5,773 at 0.5 (5,428 unlabelled), all ST; 0 at 0.8.

The A31(b) isotonic map is honest for S/T and replicates CV (0.62 / 0.61). There is no ≥0.8 tier. The ≥0.5 tier means about 60% top-1 precision on *labelled* sites; under PU selection, that is not a precision on unlabelled sites.

**Temporal and external ECE.** The batch reports ECE 0.122 (new since 2021) and 0.218 (external), with ≥0.5-tier accuracy of 0.36 and 0.21. **These are not miscalibration.** They are mostly a truth-set artefact:
- The batch scores each subset only against that subset's new or external kinases at the site, so a top call that is the site's older PSP kinase counts as wrong.
- With each site's full truth set (audit diagnostic), the ≥0.5 tier is 0.55 observed vs 0.63 claimed (n = 53, new since 2021) and 0.54 vs 0.61 (n = 28, external). ECE becomes 0.062 and 0.068.
- That is mild overconfidence at most, within the Wilson intervals (`audit_test.md` §6).

### 7.8 CV→TEST gap (leakage alarm)

**M1 on common kinases** (K_eval ∩ K_cv; `docs/audit/audit_test.md` §4):

| Config | CV (K_cv) | CV on common kinases | TEST on the same kinases | paired Δ ± SE |
|---|---|---|---|---|
| **final** | 0.786 | 0.792 | 0.793 | **+0.001 ± 0.012** |
| R2 | 0.776 | 0.786 | 0.778 | −0.007 ± 0.013 |
| R3 | 0.754 | 0.763 | 0.773 | +0.010 ± 0.012 |
| B1 (i) | 0.717 | 0.737 | 0.713 | −0.024 ± 0.013 |
| B3 | 0.677 | 0.715 | 0.709 | −0.007 ± 0.012 |

- **No optimism in the final model.** This matches the audited bound of ≲0.005 over about 90 CV configurations. M3 is equally stable: CV 0.358 on common kinases vs TEST 0.361.
- **Residue split, final model on the same kinases:**

  | | CV | TEST | paired Δ ± SE |
  |---|---|---|---|
  | ST (92 kinases) | 0.827 | 0.841 | +0.014 ± 0.012 |
  | **Y (23 kinases)** | 0.667 | 0.595 | **−0.072 ± 0.033** |

  SYK (0.85 → 0.56), IGF1R (0.82 → 0.38) and LYN (0.68 → 0.49) collapse; most Y kinases have only 3–4 TEST positives.

**M2 composition.** "On the same 117 kinases, CV M2 = 0.348 vs TEST 0.363 (paired Δ +0.014 ± 0.021). No leak signal."

| | final | R2 | R3 | B1 (ii) | B0-prior |
|---|---|---|---|---|---|
| CV M2 over K_cv (203) | 0.263 | 0.200 | 0.223 | 0.171 | 0.064 |
| CV M2 over K_eval ∩ K_cv (117) | **0.348** | 0.246 | 0.320 | 0.180 | 0.111 |
| TEST M2 over K_eval | **0.363** | 0.264 | 0.322 | 0.189 | 0.109 |

Source: `audit_test.md` §2.

- Of the +0.100 CV→TEST change, about +0.085 is kinase selection. K_eval holds 13 / 73 / 31 kinases in the 5–19 / 20–99 / ≥100 bins, against 83 / 89 / 31 for K_cv, and Hit@5 rises steeply with label count.
- The label-free B0-prior control moves the same way (0.064 → 0.109; CV on the same kinases 0.111).
- **TEST M2 0.363 must not be presented as an improvement over CV.**

**Margin over Cantley, CV vs TEST.**

| | Δ M1, final − B1 |
|---|---|
| CV, matched (`results/models/r2/compare__r2_fusion_r3__vs__simple_b1_universe.tsv`) | +0.063 [0.048, 0.074] |
| TEST, matched (`superiority.tsv`) | +0.075 [0.046, 0.099] |

- Part of the larger TEST margin is a weak-baseline draw. On the 115 kinases B1 scores in both, Cantley fell from 0.741 to 0.717 (−0.025 ± 0.013), almost entirely on Y kinases (0.607 → 0.496). The final model was flat (−0.003) (`audit_test.md` §3).

### 7.9 Not computed in the batch (cannot be back-filled without a second sealed access, A29)

Source: `docs/audit/audit_test.md` §5, §6, §7 #10. The review has no §8; this list collects its "not delivered" items.

1. **Y kinases with hub sites excluded.** The planned §7.5 row; only the hub/non-hub split of positives exists.
2. **A16 calibration on the external and temporal strata** (A29).
3. **A17 coverage baselines on TEST** (random scores and Cantley-only). Only the final model supplied `assignment_prob`. The CV versions are in `docs/model_r2.md` §6.1; they predate A31, so they are not like-for-like.
4. **A paired T3 test** against any baseline. The T3 comparisons above are descriptive.
5. **CIs for strata.** Every §7.5–7.6 stratum comparison is descriptive.
6. **The priming-aware Cantley arm** as a TEST comparator. It was not in FREEZE, so no TEST claim can be made against it.
