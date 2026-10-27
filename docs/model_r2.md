# R2 — global kinase–site pair model (LightGBM)

Owner: R2 agent. Status (21:25): CV work done; final DEV/TEST runs are prepared (§6.4) and waiting for the G4 go.

**Headline.** Pooled 5-fold DEV CV, ksp.eval A24 (harness 974d520), K_cv = 203 kinases.

| model | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| final candidate = R2 (v3, 5-seed ensemble) late-fused with R3 | **0.786** | **0.263** | **0.361** | 0.774 |
| R2 alone, 5-seed ensemble | 0.776 | 0.200 | 0.312 | 0.767 |
| B1 Cantley (universe percentile) | 0.719 | 0.171 | 0.133 | |

- Per-fold M1 of the final candidate: 0.767 / 0.788 / 0.782 / 0.794 / 0.773.
- Paired vs B1-universe on the matched universe: ΔM1 +0.063 [0.048, 0.074], ΔM2 +0.105 [0.082, 0.125],
  ΔM3 +0.227 [0.189, 0.251].
- Pooled negative canary (v2) 0.539: passes. Positive canary 1.00.
- LKO-T2: R2 M1 0.751; the fusion 0.749.
Code: `src/ksp/models/r2/{xfit,assemble,sample,train}.py`, scripts `scripts/models/r2_*.py`. Caches: `cache/r2/`.
Results: `results/models/r2/`, runs `results/runs/r2_*`, registry rows `results/cv/REGISTRY.tsv` (model = R2).

## 1. Pipeline
1. **Cross-fitted label-derived features** (`xfit.py`; DESIGN §7, A1). F2a (learned PSSM), F3-label (PPI substrate
   neighbourhood), F6a (co-regulation; plus an F6-clean variant). All are built through `ksp.crossfit.crossfit()`,
   with one query row per U_ms site:
   - training rows = all U_ms sites on outer-training proteins;
   - eval rows = all U_ms sites of the eval fold.

   The driver assigns inner folds by substrate cluster and raises `LeakageError` if a fitted label shares a protein or
   cluster with a query site. Eval rows get a randomly assigned inner-fold fit (A1 default).

   A thin wrapper (`DenseExpand`) expands each query site to all 516 kinases (503 real + 13 canaries) and writes the
   result into a dense K × S cache. Expanding over kinases at fixed sites adds no label information, because the
   disjointness guarantee concerns sites, proteins and clusters.

   **Toy leak test:** `crossfit.LeakyLabelFlag` passed through the same wrapper gives 0 of 9,721 fold-0 training
   positives flagged. Without cross-fitting, all of them would be flagged. Run it with
   `scripts/models/r2_xfit.py --toy`.

   The families are subclassed only to (a) include canary kinases and (b) run F6a on CPU.
2. **Pair-feature assembler** (`assemble.py`). Features are gathered per kinase: row *k* of every kinase × site and
   kinase × protein store, site vectors, and kinase descriptors. Nothing is ever materialised for all pairs. Feature
   groups F0–F8 are listed in the module docstring. Forbidden inputs (A6/A11) are not loaded:
   - centre-site observation flags;
   - observation counts;
   - HPA reliability;
   - `loc_n_sources`;
   - COMPARTMENTS knowledge (ablation only).
3. **Training set** (`sample.py`). Positives are the fold-training DEV L2 labels (U_ms, residue class matched),
   including canaries. Negatives are drawn per kinase from training-fold U_ms sites, after:
   - the NegativeMasks KP / Par / Hamming ≤ 2 masks, built from **training-fold label rows only**;
   - a hard/E-all mix as in A15.

   Each row's sample weight is n_k^−α, normalised. W_k and ρ_k are logged per kinase in `train_stats_fold*.tsv`.
4. **Training** (`train.py`). LightGBM with 16 threads. Early stopping uses an inner validation split of training
   clusters (`splits.train_valid_split_inner`, 20%), scored by macro per-kinase AUROC. The model is then refitted on
   all training rows with the early-stopped number of rounds.
5. **OOF prediction.** Exhaustive over eval-fold U_ms sites × every kinase of matching class (~20 M rows per fold).
   Evaluation uses `ksp.eval.evaluate`: `prepare('cv', fold)` for fold-0 search and pooled `prepare('cv')` for 5-fold
   CV.

### Scripts
| script | purpose |
|---|---|
| `scripts/models/r2_xfit.py --family f2a f3nb f6a f6a_clean --folds 0 1 2 3 4` (`--toy` = leak test) | dense cross-fitted features |
| `scripts/models/r2_run.py --name N --folds F --set k=v ...` | one configuration on the given folds (train, OOF predict, eval, REGISTRY) |
| `scripts/models/r2_cv5.sh N threads stage "k=v ..." [fold0_run]` | 5 per-fold runs, then `r2_pool.py` (pooled OOF evaluation + bootstrap + strata, registered) |
| `scripts/models/r2_queue.sh queue.q threads logprefix` | sequential queue of `r2_run.py` configurations |
| `scripts/models/r2_compare.py --a RUN --b RUN [--fold F]` | matched-universe paired cluster bootstrap (`ksp.eval.paired_bootstrap`) |
| `scripts/models/r2_canary_diag.py --runs ...` | density-stratified negative-canary check |
| `scripts/models/r2_combine_m2.py`, `r2_calib.py`, `r2_partial_pool.py` | site-grouped M2 combination, A16/A17 calibration and assignment, §8.3 λ sweep |

## 2. Assumptions / deviations (stated, not silent)
- **F8 (Sugiyama 2019):** no feature store existed when R2 was built, so the group is empty. This is counted as a
  manifest fallback (`F8_unavailable`).
- **F1 for kinases without an own atlas matrix** uses the nearest-atlas imputed matrix in the same column, with
  `f1_impute_dist` set to 0 for an own matrix.
- **Canary kinases** have no rows in the kinase × site / kinase × protein stores, so those features are NaN. The
  exception is the planted PPI flag of CANARY_POS, exposed as `f3_edge_any`, as `ksp.canary` specifies. F7 slots use
  `canary_descriptors()`.
- **Ranking objectives** use no per-row weights.

## 3. Runtime (measured, fold 0; 8 LightGBM threads, shared machine)
| step | time |
|---|---|
| xfit F2a / F3-label / F6a / F6a-clean per outer fold (5 inner fits + eval rows, all 516 kinases × 278k sites) | ~6.5 / ~8.5 / ~3 / ~3 min (CPU; under contention) |
| FoldData (labels, masks) | 7 s |
| training-row sampling (1:50, 496k rows) | 4 s |
| training-row feature gather (172 features) | 58 s (single thread; now threaded) |
| LightGBM early-stopped fit (63 leaves, ~210 rounds) + refit on all rows | 23 s + 14 s |
| exhaustive OOF prediction, 21.3 M rows (gather + predict) | 28 s |
| ksp.eval fold-0 evaluation | ~1–2 min |

## 4. Staged search (A14) — fold 0, then 5 folds
Selection metric: M1 on CV fold 0 (K_cv kinases with fold-0 positives, ~180), as pre-declared. M2, M3 and the
canaries are reported alongside. Every configuration is a REGISTRY row.

**Stage (a): tree size and learning rate.** All feature groups, α 0.5, 1:50, 50% hard, binary objective, inner
early stopping. The canary columns are **canary v1** (sampler with the protein-density bias; see §5.1).

| config | rounds (ES) | M1 | M2 | M3 | neg canary hard/all (v1) |
|---|---|---|---|---|---|
| leaves 31, lr .05 | 212 | 0.7503 | 0.154 | 0.261 | 0.545 / 0.561 |
| leaves 63, lr .05 | 213 | 0.7531 | 0.186 | 0.257 | 0.549 / 0.560 |
| leaves 127, lr .05 | 208 | 0.7545 | 0.190 | 0.254 | 0.566 / 0.562 |
| leaves 127, lr .05, ES on hard negatives | 199 | 0.7558 | 0.189 | 0.256 | 0.566 / 0.562 |
| leaves 255, lr .03, min_leaf 50 | 560 | 0.7547 | 0.228 | 0.278 | 0.611 / 0.597 |
| **leaves 127, lr .03** (chosen) | 364 | **0.7591** | 0.221 | 0.269 | 0.577 / 0.572 |
| leaves 255, lr .03, min_leaf 100 | 797 | 0.7578 | 0.219 | 0.279 | 0.578 / 0.565 |

The M1 differences are within fold-0 noise (about ±0.01). The chosen configuration has the best M1, and it is
about 4× cheaper to predict than 255 leaves (83 s vs 353 s for 21 M rows).

For reference, fold 0 vs B1-universe (Cantley) on the matched universe (180 kinases), leaves 63:
- ΔM1 = +0.021 [−0.002, 0.046]
- ΔM2 = +0.031 [0.005, 0.083]
- ΔM3 = +0.080 [0.040, 0.125]
- ΔM4 = +0.031 [0.005, 0.056]

**Canary v2 (15:31).** Core-infra replaced the negative canaries. They are now drawn from real DEV-labelled sites,
and the gate is on E-hard only. Real labels are unchanged (fold 0: 9,433 real training pairs before and after).

Because the canary rows of the cross-fitted features had been fitted on v1 labels, all xfit stores were rebuilt with
v2. The v1 stores are archived in `cache/r2/xfit/cv_canaryv1/`.

**Correction (audit, 17:30).** The rebuild started at 15:33. Core-infra rewrote `canary_dev.parquet` once more at
15:34, so six stores were fitted on an intermediate canary table: fold 0 F2a/F3-label/F6a/F6a-clean (331 canary
training pairs, against 348 in final v2) and fold 3 F2a/F3-label (329, against 350). Each store's json records
`n_canary_labels`, which is how this was found.

Those stores were rebuilt at 17:30 with `--overwrite`. Real-kinase rows of F3-label/F6a are per-kinase computations,
so they are unaffected. Every final-configuration run (`r2_final_s*`) was started after the rebuild.

Runs completed before 17:30 carry fold-0/fold-3 canary rows from the intermediate table: G2, full, ablations and
alarms. Their **canary** numbers are therefore approximate (the real-kinase numbers are not affected beyond tree
randomness). The canary gate is reported from the final runs.

Real-kinase features are unaffected, except for a possible tiny shift in F2a's internally tuned α/β, which are tuned
over all kinases including canaries. Stage-(a) selection therefore stands. Every run after stage (a) uses v2.

The chosen configuration re-run with v2 gives fold-0 M1 0.7537 (0.7591 with v1). A small perturbation of the F2a
features alone moves fold-0 M1 by 0.005, which is a useful noise floor for the search.

**Stages (b)–(d)**, fold 0, canary v2. Base: leaves 127, lr 0.03, α 0.5, 1:50, 50% hard, binary.

| stage | config | M1 | M2 | M3 | neg canary E-hard |
|---|---|---|---|---|---|
| b | **α 0.5 (base)** | **0.7537** | 0.223 | 0.262 | 0.635 |
| b | α 0 | 0.7498 | 0.159 | 0.266 | 0.573 |
| b | α 1 | 0.7487 | 0.225 | 0.263 | 0.593 |
| c | 1:20, 25% hard | 0.7468 | 0.195 | 0.285 | 0.592 |
| c | 1:20, 50% hard | 0.7496 | 0.168 | 0.265 | 0.590 |
| c | 1:50, 25% hard | 0.7425 | 0.202 | 0.254 | 0.574 |
| c | 1:50, 75% hard (extra arm) | 0.7507 | 0.193 | 0.263 | 0.604 |
| c | 1:200, 25% hard | 0.7415 | 0.203 | 0.274 | 0.588 |
| c | 1:200, 50% hard | 0.7460 | 0.201 | 0.274 | 0.631 |
| d | lambdarank grouped by kinase (≤1000-row queries) | 0.7535 | 0.086 | 0.261 | 0.598 |
| d | site-grouped lambdarank (M2 objective) | 0.7297 | 0.189 | 0.184 | 0.400 |
| d | site-grouped, 600 fixed rounds | identical: inner NDCG@5 peaks before round 50 | | | |

Decision: keep the base configuration (binary, α 0.5, 1:50, 50% hard). Nothing beat it on M1 by more than the noise
floor.
- The site-grouped objective does not beat the binary model on M2 (0.189 vs 0.223), so M2 will use the calibrated
  cross-kinase assignment score of the binary model (A20, §7).
- lambdarank by kinase matches M1, but its scores are not comparable across kinases (M2 collapses).

Configurations tried on fold 0: 7 (a) + 1 (a, v2 rerun) + 2 (b) + 6 (c) + 3 (d) = 19, plus one label-free
control.

## 5. Ablation curve and leakage alarms
Setup: seed 0, all 5 folds, chosen tree configuration, pooled DEV OOF. Metrics come from `ksp.eval` **A24**
(within-(fold, kinase) rank normalisation; harness frozen at commit 974d520) and are re-evaluated with
`scripts/models/r2_reeval.py`; the full table is `results/models/r2/cv_table.tsv`. CIs are 95% cluster bootstrap,
B = 1000. These runs predate the final v3 settings (F2a fitted jointly with canaries, descriptor v1), and their canary
numbers are approximate (see §4).

| step | M1 [95% CI] | per-fold mean M1 | M2 | M3 | M4 | neg/pos canary E-hard |
|---|---|---|---|---|---|---|
| F1 (Cantley) | 0.723 [0.707, 0.737] | 0.725 | 0.136 | 0.125 | 0.703 | 0.50 / 0.50 |
| + F0 (site) | 0.733 [0.714, 0.747] | 0.729 | 0.169 | 0.236 | 0.717 | 0.55 / 0.45 |
| + F2a (learned PSSM) | 0.738 [0.718, 0.751] | 0.733 | 0.094 | 0.239 | 0.720 | 0.53 / 0.79 |
| + F3 (PPI + F3-label) | 0.748 [0.730, 0.761] | 0.745 | 0.143 | 0.260 | 0.728 | 0.55 / 0.99 |
| + F4 (localisation) | 0.753 [0.734, 0.767] | 0.750 | 0.141 | 0.269 | 0.742 | 0.54 / 1.00 |
| + F6a (co-regulation) | 0.755 [0.737, 0.768] | 0.751 | 0.149 | 0.274 | 0.746 | 0.53 / 0.99 |
| + F6b (inhibitor) | 0.755 [0.738, 0.767] | 0.752 | 0.154 | 0.281 | 0.742 | 0.52 / 1.00 |
| + F6c (CPTAC) = G2 run | 0.754 [0.736, 0.767] | 0.751 | 0.146 | 0.272 | 0.735 | 0.54 / 1.00 |
| + F5 (co-expression) | 0.759 [0.741, 0.773] | 0.757 | 0.160 | 0.283 | 0.746 | 0.53 / 1.00 |
| + F8 (Sugiyama) | not available (no store) | | | | | |
| + F7 (kinase descriptors) = full | **0.773** [0.755, 0.785] | 0.768 | **0.208** | **0.312** | 0.770 | 0.54 / 1.00 |
| + knn_sim (final v3 settings; §6) | 0.774 | 0.769 | 0.175 | 0.311 | 0.761 | 0.55 / 1.00 |

Reading:
- The largest single gains come from F0 (M3 0.125 → 0.236), F3 (+0.010 M1) and F7 (+0.014 M1, +0.05 M2).
- F2a alone lowers M2: the per-kinase learned-PSSM scales are not comparable across kinases. F7 repairs that.
- F6 adds about +0.002 M1 in total and +0.01 M3.

**Leakage alarms** (all relative to the full model, 0.773):

| alarm | M1 | M2 | M3 | verdict |
|---|---|---|---|---|
| (i) + COMPARTMENTS knowledge channel | 0.776 (+0.003) | 0.204 | 0.320 | small gain; the curated channel adds little beyond HPA/OpenCell/HCM/Hein; stays out |
| (ii) F6-clean (18 label-source PMIDs removed from F6a/F6b) | 0.773 (−0.000) | 0.194 | 0.302 | no dependence on label-source studies |
| (iii) B0-miss (LightGBM on F6 missingness only, same pipeline) | 0.522 | 0.018 | 0.054 | missingness carries ~no within-kinase E-hard signal; the F6 gain is not a missingness artefact |
| label-free control F0+F7 (fold 0) | 0.663 (fold 0) | 0.152 | 0.183 | negative canary 0.54 |

### 5.1 Negative-canary caveat (found 15:15; for PI / core-infra)
On fold 0, R2's negative canaries score 0.545–0.566 (E-hard) and 0.56 (E-all). That passes the ≤ 0.57 gate, but
it is consistently above 0.5. I believe the cause is the canary **design**, not leakage.

`ksp.canary` samples random proteins first and then 1–4 sites per protein. Canary positives are therefore enriched
on proteins that carry few U_ms sites. A label-free, protein-level score shows this directly. On pooled DEV CV:
- −(number of U_ms sites on the protein) gives a **mean negative-canary AUROC of 0.74 (E-all) / 0.76 (E-hard)**
  (range 0.54–0.88);
- −(protein length) gives 0.70 / 0.63;
- random scores give 0.49.

So any model with protein-level features can exceed 0.5 on the negative canaries without seeing a label. R2's
0.55–0.57 is far below what this proxy allows, which is reassuring, but the gate is not a clean null.

Suggested fix (core-infra): sample canary sites uniformly over U_ms sites, or match the protein site-count
distribution of real labels.

**Density-stratified check** (`scripts/models/r2_canary_diag.py`, `results/models/r2/canary_density_diagnostic.tsv`).
This is the fold-0 negative-canary AUROC on E-all, recomputed *within* deciles of protein site count and averaged
with weights by positives:

| run | raw E-all | density-stratified | proxy −n_sites (raw) |
|---|---|---|---|
| r2_a_l63 | 0.560 | **0.482** | 0.686 |
| r2_a_l127 | 0.562 | **0.488** | 0.686 |
| r2_a_l255 (lr 0.03, min_leaf 50) | 0.597 (fails the 0.57 gate) | **0.515** | 0.686 |

Within density strata the negative canaries sit at chance. The excess AUROC is therefore the sampler bias, which
larger trees learn from the canaries' training labels. It is not label leakage. The positive canary is 0.999
(stratified 0.998).

Selection policy: the official A5 gate stays binding. A configuration that fails it is reported, but it is not
selected unless the PI accepts this diagnostic.

### 5.2 Canary v2 results and diagnosis (fold 0; for PI / core-infra)
| run | neg canary E-hard (gated ≤ 0.57) | neg E-all | pos canary |
|---|---|---|---|
| full model, chosen config (α 0.5) | **0.635** (fails) | 0.802 | 0.999 |
| same, α 0 | 0.573 (fails marginally) | 0.756 | 0.999 |
| same, α 1 | 0.593 (fails) | 0.760 | 0.999 |
| **label-free control**: F0 (site) + F7 (kinase descriptors) only | **0.539** (passes) | 0.734 | 0.810 |

How the diagnosis was done:
- Stratifying the E-hard AUROC by the number of labelled (E-hard-eligible) sites per protein does not remove the
  excess (0.627).
- Univariate AUROCs of the canary's own feature columns, canary positives vs E-hard negatives (11 negative
  canaries):
  - cross-fitted F3-label (proximity of the site's protein to the canary's training substrate proteins) scores
    0.39, i.e. 0.61 inverted: canary positives sit on proteins *less* central than E-hard sites;
  - site structure and ESM features reach 0.59 (RSA, SASA, several ESM PCs).

So canary-v2 positives and E-hard negatives are different populations even without labels. The v2 sampler picks
proteins first, then 1–4 labelled sites per protein. E-hard is site-weighted and so dominated by heavily studied
hub proteins.

A model with kinase-specific capacity learns "canary k prefers less-central, more exposed sites" from the canary's
training labels and carries that over to the eval fold. The cross-fitted F3-label is the main channel for this, and
its eval rows are verified disjoint by protein and cluster in `ksp.crossfit`. This looks like legitimate
generalisation of a sampler artefact, not a label leak.

I have **not** proven that, though. A decisive gate needs negative-canary sites drawn **site-uniformly from the
E-hard-eligible pool**, so that positives and E-hard negatives share one distribution. I recommend this to
core-infra.

Until then the gate failure is reported on every run, and selection uses real-kinase M1 only.

## 6. Results (pooled 5-fold DEV OOF; ksp.eval A24 = harness commit 974d520; K_cv = 203 kinases)
All numbers in this section use the A24 eval. Pre-A24 numbers in earlier progress notes are superseded; in
particular full-model M4 0.782 is now 0.770.

**G2 gate: met (16:27).** R2 with F1+F0+F2a+F3+F4+F6 (`r2_g2_main_cv5`):
- M1 0.754 [0.736, 0.767]; per-fold 0.737 / 0.756 / 0.765 / 0.760 / 0.737;
- paired vs B1-universe on the matched universe (pre-A24 eval): ΔM1 +0.023 [0.009, 0.034], ΔM3 +0.130 [0.098, 0.153],
  ΔM2 −0.009 (not significant).

**Final R2 configuration v3.** Binary objective, 127 leaves, lr 0.03, inner early stopping on real kinases only,
α 0.5, 1:50 with 50% hard. Groups F1+F0+F2a+F3+F4+F6abc+F5+F7, with F2a fitted with canaries separate (`f2a_sep`)
and descriptor v2 (A25).

Seeds 0–4 (`r2_final_s{0..4}_cv5`):

| | M1 | per-fold mean M1 | M2 | M3 | M4 |
|---|---|---|---|---|---|
| mean ± sd over seeds | 0.7711 ± 0.0028 | 0.7669 ± 0.0025 | 0.187 ± 0.008 | 0.300 ± 0.002 | 0.761 ± 0.001 |
| seed ensemble (mean logit; `r2_final_ens5`) | **0.776** | 0.771 (0.756 / 0.782 / 0.772 / 0.787 / 0.761) | 0.200 | 0.312 | 0.767 |

- Negative canary (pooled E-hard) per seed: 0.547 / 0.547 / 0.541 / 0.550 / 0.525, all passing ≤ 0.57.
  Ensemble: 0.548.
- Positive canary: 0.999.

**+knn_sim** (`Fknn`, 16 cross-fitted columns built through `ksp.models.r2.xfit`). They agree with feat-knn's own
stores to float16 precision (r > 0.99999). Seeds 0–4: M1 0.7736 / 0.7738 / 0.7748 / 0.7732 / 0.7748, mean 0.7740 ±
0.0008.
- vs v3 without knn: **+0.003 mean M1**. The seed-0 pair gave +0.0054 [0.000, 0.012], but seed 1 only +0.0004.
- ΔM3 +0.012 [0.002, 0.024].
- **ΔM2 −0.018 [−0.030, −0.006]**: the per-kinase kde/log_ntrain columns break cross-kinase comparability.

Decision under the PI rule (M1 gain ≥ 0.005 with the canary gate passing): **not adopted**. The seed-averaged gain
is 0.003, and M2 drops significantly.

**F9 (Tyr)** ablation (`r2_final_f9_s0_cv5`, seed 0). Columns: sh2_dock_{max,local,spec}, itam_gated,
kinase_desc{is_tm,is_receptor,has_sh2_matrix} and itam_order.
- M1 0.7699 vs 0.7682 (+0.0017); M2 0.195 vs 0.192; M3 0.310 vs 0.298; canary 0.537 (passes).
- The overall M1 gain is below the 0.005 adoption threshold, so F9 is **not adopted**. The Y-stratum test was not
  needed.

**Late fusion with R3 (A13)**, `r2_fusion_r3__20261026-191901__s0`. The blend is w·z(R2 seed ensemble) +
(1−w)·z(R3 final OOF), with z = standardised logits per fold. The weight is chosen on pooled-CV M1; the score_site
weight (M2) is chosen on M2 (A20).

| w_R2 | 0 (R3 alone) | 0.3 | 0.4 | 0.5 | **0.6** | 0.7 | 1 (R2 alone, ensemble) |
|---|---|---|---|---|---|---|---|
| M1 | 0.754 | 0.780 | 0.784 | 0.786 | **0.786** | 0.785 | 0.776 |
| M2 | 0.223 | 0.258 | 0.261 | **0.263** | 0.257 | 0.246 | 0.200 |
| M3 | 0.294 | 0.353 | 0.363 | 0.360 | 0.361 | 0.351 | 0.312 |

The R3-alone column is from the seed-0 sweep. The ensemble sweep covered w ∈ {0.3, …, 0.8}, and both sweeps agree
to within 0.002.

**Candidate final model = fusion (score w = 0.6, score_site w = 0.5):**
- M1 **0.786**; per-fold 0.767 / 0.788 / 0.782 / 0.794 / 0.773 (mean 0.781);
- M2 **0.263**, M3 **0.361**, M4 0.774;
- negative canary E-hard 0.539 (passes); positive canary 1.000.

Paired vs B1-universe on the matched universe: ΔM1 **+0.063 [0.048, 0.074]**, ΔM2 **+0.105 [0.082, 0.125]**,
ΔM3 **+0.227 [0.189, 0.251]**, ΔM4 +0.071 [0.053, 0.080] (p = 0.002 for all).

### 6.1 Calibration (A16) and site-level assignment (A17) — final fusion model
`scripts/models/r2_calib.py`; outputs in `results/models/r2/calib/r2_fusion_r3__20261026-191901__s0__{calibrated,score_site}/`.

Everything here is cross-fitted over CV folds: the calibrator and the temperature are fitted on 4 folds' OOF rows and
applied to the 5th. Calibration is exhaustive over DEV U_ms, with E-all negatives subsampled to 20k per kinase with
inverse-probability weights.

Calibrator tiers:
- per-kinase Platt for kinases with ≥ 30 OOF positives (61 kinases);
- otherwise a family-level slope with a per-kinase intercept (261), or a global slope (17);
- kinases with no positives (164) get an intercept set so that their expected number of positives is 0.5. Borrowing
  the family intercept instead over-predicted the 0-label bin by about 2,000 positives in a first version.

Calibration-in-the-large, out of sample (expected / observed curated positives) by training bin:

| bin | expected / observed |
|---|---|
| 1–4 | 173 / 147 |
| 5–19 | 863 / 848 |
| 20–99 | 2,558 / 2,361 |
| ≥ 100 | 5,739 / 5,400 |
| 0 | 84 / 0 |

ECE < 1e-4 in every bin, which is trivial at this prevalence; see `reliability_by_bin.tsv`. Under PU labels these
probabilities estimate P(curated human-native PSP positive | score), a **lower bound** on the true probability.

**Assignment probability** p(k | s) is a softmax over the **calibrated** kinase logits of the real kinases of matching
class. The temperature T ≈ 0.93 is fitted on labelled DEV sites by OOF log-likelihood.

Reliability on labelled DEV sites (out of sample):

| p_top bin | sites | mean p_top | top kinase is a true kinase |
|---|---|---|---|
| [0, 0.2) | 1,835 | 0.13 | 0.18 |
| [0.2, 0.5) | 3,021 | 0.32 | 0.40 |
| [0.5, 0.8) | 794 | 0.62 | 0.61 |
| [0.8, 1] | 167 | 0.87 | **0.54** (overconfident) |

Overall top-1 accuracy is 0.365.

Cantley-only, same procedure, for comparison:
- top-1 accuracy 0.282;
- on labelled sites only 152 reach p ≥ 0.5 (observed 0.66), and none reach 0.8.

Its large unlabelled coverage@0.5 comes from a smaller candidate set and from confident motif-family calls on
unlabelled sites. That cannot be verified, so **coverage counts are not comparable across methods unless they share a
candidate set**. Reliability on labelled sites is the comparable quantity: R2 puts 961 labelled sites at ≥ 0.5 with
0.59 observed precision, Cantley 152.

**Coverage over DEV U_ms OOF (277,924 sites):**

| | coverage@0.5 | coverage@0.8 |
|---|---|---|
| this model | **14,235** (13,274 unlabelled) | **1,674** (1,507 unlabelled) |
| random scores | 5 | 0 |
| Cantley-only (B1-imputed; only the 339 kinases with ≥ 1 DEV positive, so fewer candidates per site) | 23,867 (23,715 unlabelled) | 477 |

Coverage is conditional on the site being a target of one of the modelled kinases.

The ≥ 0.8 tier is overconfident: 0.87 claimed vs 0.54 observed. So coverage@0.8 should be read as "≥ 0.5 precision
on labelled sites", not 0.8.

An uncalibrated softmax (raw fused logits, T ≈ 0.47) is much worse: [0.8, 1] observed 0.13, and coverage@0.5 of
10,722 against 9,757 for random scores. The reason is that pseudokinases and 0-label kinases (PRAG1, TTK, ALPK2,
BUB1B) get extreme raw logits.

**Correction (audit round 2, R2-M1; PI 22:12).** The tiers above are **not trusted**:
- Y sites at p ≥ 0.8 are 0/41 correct, and the ST top tier is 0.85 claimed vs 0.72 observed.
- The adopted fix is the auditor's isotonic recalibration of p_top, fitted separately for ST and Y on labelled DEV
  OOF sites and cross-fitted by fold on CV (`docs/audit/audit_cv_r2_scripts/calib_check3.py`).

Result:
- ST: [0.5, 0.8) 0.62 claimed / 0.61 observed (n = 1,272), and only 2 sites at ≥ 0.8.
- Y: no site reaches 0.5.
- **Coverage@0.5 = 21,115 DEV U_ms sites (19,841 unlabelled), all ST; coverage@0.8 = 9.**

So the model has no reliable ≥ 0.8 tier and no reliable Y-site assignment at ≥ 0.5. The coverage numbers in the
table above are superseded.

For TEST, `assignment_prob` uses the same chain with maps fitted on all DEV OOF (`scripts/models/r2_assign_test.py`):
1. the tiered calibrators;
2. a softmax with T = 0.926 (mean of the fold temperatures);
3. the ST/Y isotonic map applied to the top kinase;
4. the other kinases rescaled to fill the remaining probability mass.

`score_site` (M2) stays the raw fused score.

**M2 component (A20).** The calibrated logit gives M2 0.205. The raw fused score_site (w = 0.5) gives **0.263**.
Each component is selected on its own endpoint, so M2 uses the raw fused score_site and coverage uses the calibrated
assignment probability.

### 6.2 benchmarKIN-clean variant (DESIGN §6.7)
Phospho features without the 22 `clean_bench` PMIDs and without any CPTAC-derived feature:
- F6b `*_clean_bench` stores;
- F6a refitted on compendium blocks minus those PMIDs, with the CPTAC columns NaN (xfit family `f6a_bench`);
- F6c dropped.

Otherwise it is the v3 configuration, seed 0 (`r2_final_cleanbench_s0_cv5`):
- **CV M1 0.771**, M2 0.189, M3 0.300, M4 0.768;
- canary 0.551 / 0.999.

v3 seed 0 gives 0.768, so the benchmark studies and CPTAC contribute nothing measurable to M1.

The DEV-trained U_obs library is `r2_final_cleanbench_uobs__20261026-204513__s0`:
`cache/runs/<id>/uobs_predictions.parquet`, 131 M rows = all U_obs × the 503 real kinases of matching class.
- Training took 274 s and inference 729 s.
- The full 131 M-row parquet OOM-killed the downstream harness. It was cut to each kinase's top 2·|PSP_all(k)| + 200
  rows (`uobs_predictions_topk.parquet`, made by `scripts/models/r2_topk_filter.py`), which is all the arms need.
- `scripts/downstream/run_predictions.py --tag r2_clean` was then run; outputs are in
  `results/downstream/r2_clean_*`.

benchmarKIN results (z-score, native, U_obs; mean AUROC, with the number of evaluable kinases in brackets):

| benchmark | r2_clean top50 | top50 −PSP | topPSP | PSP_all | PSP_dev | Cantley top50 | shuffled PSP |
|---|---|---|---|---|---|---|---|
| perturbation (all) | 0.775 (50) | 0.720 (48) | 0.745 (38) | 0.783 (41) | 0.784 (37) | 0.620 (13) | 0.48–0.55 |
| tumour (CPTAC) | 0.653 (149 kinases) | 0.653 | 0.692 | 0.664 | 0.666 | 0.540 | 0.51–0.56 |

- The clean R2 library matches PSP on the perturbation benchmark with more evaluable kinases (50 vs 41).
- It beats PSP on the tumour benchmark with the topPSP arm (0.692 vs 0.664).
- With all PSP pairs removed (the −PSP arm) it still clearly beats Cantley and the shuffled libraries.

### 6.2b "nophospho" variant (PI 21:45; a proxy for leave-drug-out)
This is the clean configuration with **all** F6 removed (no compendium, no CPTAC):
F1+F0+F2a+F3+F4+F5+F7, v3 settings, seed 0.
- CV (`r2_final_nophospho_s0_cv5`): **M1 0.771** (per-fold mean 0.767), M2 0.187, M3 0.296, canary 0.549/0.999.
  The equivalents are v3 0.768 and clean 0.771, so F6 adds nothing measurable to CV M1.
- DEV-trained U_obs library: `r2_final_nophospho_uobs__20261026-213950__s0` (train 152 s, inference 251 s). Its top-k
  file `uobs_predictions_topk.parquet` passes `pred_topk.py` (no ties at the cut, all arms identical to the superset),
  and `run_final.py --tag r2_nophospho` has been run.

Matched-kinase deltas vs PSP_all (z-score, U_obs, paired bootstrap; ours − PSP_all):

| subset | arm | clean Δ [95% CI] | nophospho Δ [95% CI] |
|---|---|---|---|
| perturbation, all | topPSP | −0.017 [−0.054, 0.018] | −0.007 [−0.047, 0.032] |
| perturbation, all | top50 | −0.013 [−0.043, 0.012] | −0.043 [−0.087, 0.004] |
| perturbation, all | top100 | −0.019 [−0.053, 0.012] | −0.050 [−0.088, −0.013] |
| drug unseen | topPSP | −0.058 [−0.123, 0.003] | −0.047 [−0.116, 0.026] |
| drug unseen | top50 | −0.050 [−0.096, −0.001] | −0.077 [−0.149, −0.011] |
| drug unseen | top100 | −0.041 [−0.099, 0.017] | −0.081 [−0.140, −0.017] |
| tumour (CPTAC7) | topPSP | +0.016 [−0.018, 0.051] | +0.036 [−0.006, 0.077] |
| tumour (CPTAC7) | top50 | +0.011 | +0.042 |

Reading:
- Removing all phospho features does **not** close the gap to PSP on unseen-drug experiments. topPSP is about the
  same; top50/top100 are worse (−0.08).
- So the clean library's parity with PSP on the perturbation benchmark does not come from response data for the same
  compound. If anything, the phospho features (from other drugs) help on unseen drugs.
- The remaining unseen-drug gap of about 0.04–0.06 is shared by both variants.
- On the tumour benchmark the nophospho variant is slightly better vs PSP. CPTAC features were excluded in both
  variants anyway.

### 6.3 Leave-kinase-out T2 (kinase folds `results/splits/kinase_folds.tsv`, pooled; evaluated like R3's LKO pool)
For each kinase fold, the held-out kinases' labels are removed from positives, negative masks and all label-derived
features (xfit tags `lko<j>`). The held-out kinases are then scored on all DEV U_ms sites.

| model | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| R2 v3, seed 0 (`r2_final_lko`) | **0.751** | 0.133 | 0.295 | 0.749 |
| R2 + R3 fusion (w 0.6 / 0.5) | 0.749 | 0.114 | **0.337** | 0.748 |
| R3 final distilled (reference, R3 agent) | 0.686 | | | |

On unseen kinase families with seen proteins, R2 loses only 0.02 M1 relative to CV (0.771 → 0.751). The label-free
Cantley atlas (F1, with nearest-atlas imputation) and the kinase descriptors carry the specificity.

### 6.4 Final DEV retrain and TEST-candidate predictions (DELIVERED 23:02 on PI go 22:12; G4 runbook format)
**Final deliverable: `r2r3_final__20261026-230138__s0`.** Both files pass the `eval_test --check-predictions`
pre-flight (`"ok": true`: 503 kinases full, 0 partial, 0 rows outside the candidate set).

- `test_predictions.parquet`: sha256 `1d95b271672b1db20ae63c2918a32a911e6800dc969b53b6514a3dceff7ce9c4`; 25,922,452 rows.
  Columns: kinase_idx, site_idx, score (fused, w 0.6), score_site (fused, w 0.5) and assignment_prob (A17 with the
  audit's ST/Y isotonic fix; sums to 1 per site).
- `test_t3_predictions.parquet`: sha256 `d62195b8f3a66df64cb0f91f5261469d1c80ed1204922bd016470aaafc8ad1b1`. This is R2
  T3 with **5 seeds** per LKO fold, fused with R3's T3 file.
- Runtime: train_wall_s 1136.6 = R2 418 + R3 719; inference_wall_s 316.9 = R2 299 + R3 3 + fusion/assignment 15.
- Components:
  - `r2_final_test__20261026-221039__s0`: R2 alone, 5 seeds, sha `2aa381c9…`;
  - `r2_final_t3__20261026-221039__s0`: 5 seeds; R2 train 2,520 s, inference 406 s;
  - `r3_final_distill_test__20261026-220424__s0`: pinned by the PI;
  - `r2r3_final_test__20261026-222344__s0` (fusion), `r2r3_final_test_assign__20261026-222401__s0` (assignment_prob),
    `r2r3_final_t3__20261026-230121__s0` (T3 fusion).
- TEST assignment summary: T = 0.926; coverage@0.5 = 5,773 TEST sites, all ST; coverage@0.8 = 0. The isotonic maxima
  are ST 0.756 and Y 0.445.

The original preparation notes follow.
The stores are built:
- `cache/r2/xfit/final/all/{f2a_sep,f3nb,f6a,f6a_bench}`: rows = DEV U_obs, eval rows = TEST-protein U_obs sites,
  A1 inner-fold fits;
- `cache/r2/xfit/t3lko{0..4}/all/*`: being pre-built.

Commands, with V3 = `num_leaves=127 learning_rate=0.03 f2a_family=f2a_sep desc_v=2 es_real_only=true`:

    # 1. R2 final, 5-seed ensemble, TEST candidates (U_ms on TEST proteins x 503 real kinases of matching class)
    PYTHONPATH=src R2_THREADS=16 .venv/bin/python scripts/models/r2_final.py --name r2_final_test --target test --seeds 0 1 2 3 4 --set $V3
    # 2. R2 T3: each LKO fold's kinases scored by the model trained without that fold's labels
    PYTHONPATH=src R2_THREADS=16 .venv/bin/python scripts/models/r2_final.py --name r2_final_t3 --target t3 --seeds 0 1 2 3 4 --set $V3
    # 3. Final = late fusion with R3's TEST predictions (CV weights: score 0.6, score_site 0.5); likewise for T3
    PYTHONPATH=src .venv/bin/python scripts/models/r2_fuse_test.py --name r2r3_final_test --r2 cache/runs/<r2_final_test>/test_predictions.parquet --r3 <R3 test_predictions.parquet>
    PYTHONPATH=src .venv/bin/python scripts/models/r2_fuse_test.py --name r2r3_final_t3 --kind t3 --r2 cache/runs/<r2_final_t3>/test_t3_predictions.parquet --r3 <R3 test_t3_predictions.parquet>
    # 4. pre-flight
    .venv/bin/python -I src/ksp/eval_test.py --check-predictions cache/runs/<run_id>/test_predictions.parquet

Expected runtime:
- (1) about 5 × (≈ 5 min train + ≈ 2 min inference of 25.9 M rows) ≈ 35 min at 16 threads;
- (2) about 5 folds × 5 seeds × ≈ 4 min ≈ 100 min, or ≈ 20 min with `--seeds 0`.

**Plan at the PI's go (~01:00):** run (1) and (2) concurrently at 12 threads each, with 5 seeds for both. Expected
finish is ~02:00 for (1) and ~03:00 for (2). If (2) is not past fold 2 by 02:15, it is restarted with `--seeds 0`
(T3 = 1 seed, stated in the manifest notes). Then (3) fusion with R3's files, and (4) pre-flight.

The pre-built stores (`final/`, `t3lko*/`) use label variant L2. If the regimes agent changes the label variant,
they must be rebuilt first: about 25 min for final and about 60 min for T3.

The manifests record `train_wall_s` and `inference_wall_s`.

(further results appended below)


## 7. G5 deployment (prepared 00:15; waiting for "GO deploy")
The frozen final configuration (v3 R2, 5 seeds, + R3, fusion weights 0.6 / 0.5) is retrained on DEV+TEST labels
through `ksp.deploy.load_deploy_labels`. That loader refuses until the TEST batch is COMPLETED; it has not been called.

Canary pseudo-kinases stay in R2 training, as in the frozen code path (`with_canaries=True`). Only the 503 real
kinases are scored and released.

1. `PYTHONPATH=src .venv/bin/python scripts/models/r2_deploy.py --stage xfit --families <fam>`, once per family:
   f2a_sep, f3nb, f6a, each in parallel at 4 threads, about 6 min each. These are the cross-fitted stores with tag
   `deploy`: rows = all 350,736 U_obs sites; inner folds over substrate clusters of DEV+TEST proteins.
2. `PYTHONPATH=src R2_THREADS=16 .venv/bin/python scripts/models/r2_deploy.py --stage train --seeds 0 1 2 3 4`
   writes `cache/release/<run>/r2_uobs_scores.parquet` (131 M rows, kinase-chunked), with train/inference wall-clock
   in the manifest. It takes about 5 × (4 min train + 5 min inference).
3. `PYTHONPATH=src .venv/bin/python scripts/models/r2_release.py --r2 <run> --r3 <R3 deploy parquet>` does the
   following:
   - fusion: z-score per model over all rows;
   - `prob`, the per-kinase calibrated probability;
   - `pct`, the within-kinase percentile;
   - `assignment_prob`, from DEV-OOF calibrators, T and the ST/Y isotonic map. Limitation: there is no OOF for
     DEV+TEST;
   - outputs in `cache/release/<run>/release_scores.parquet`;
   - coverage in `results/release/coverage.tsv`: Coverage@0.5/0.8 over all U_obs sites and over no-PSP sites, plus
     the number of sites in some kinase's top 1%; for the model, random scores and Cantley-only, the last with its
     own DEV-OOF recipe;
   - release files `results/release/kinase_top200.tsv.gz` and `site_top5.tsv.gz`, plus `RELEASE_MANIFEST.json`
     (SHA-256, licence: PSP-derived, CC BY-NC-SA, **non-commercial**, limitations).

Smoke-tested end to end on DEV labels with 6 kinases (`r2_deploy_smoke__20261027-000833__s0`,
`r2r3_release_smoke__20261027-001416__s0`; outputs under cache/release/, not results/). Expected memory is about
20 GB at peak, under the 50 GB cap.


### 7.1 G5 deployment: DONE (01:29)
Inputs: `load_deploy_labels` (DEV+TEST, L2, FREEZE 721dc197 COMPLETED), giving 15,202 training label pairs including
canaries (code path as frozen). Only the 503 real kinases are released.

| step | run | wall-clock |
|---|---|---|
| cross-fitted stores `deploy` (rows = all 350,736 U_obs sites) | `cache/r2/xfit/deploy/all/{f2a_sep,f3nb,f6a}` | 322 / 208 / 297 s, in parallel (≈ 5.4 min) |
| R2, 5 seeds, U_obs × 503 kinases | `r2_deploy__20261027-010348__s0` | train 310 s, inference 985 s (131,144,410 rows) |
| R3 (R3 agent) | `r3_deploy__20261027-005817__s0` (sha 2f62988c…, verified) | train 447 s, inference 12 s |
| fusion + calibration + assignment + coverage + release files | `r2r3_release__20261027-012616__s0` | 169 s, peak RSS 18.7 GB |

Early stopping chose 53 / 235 / 249 / 543 / 299 rounds for seeds 0–4. Seed 0 stopped at the minimum-rounds floor
(inner validation 0.830 vs 0.85–0.86 for the others). This is the frozen procedure, and the other seeds dominate the
ensemble mean.

**Release files**, licensed CC BY-NC-SA (PSP-derived, non-commercial):

| file | sha256 | size |
|---|---|---|
| `cache/release/r2r3_release__20261027-012616__s0/release_scores.parquet` (kinase_idx, site_idx, score, score_site, prob, pct, assignment_prob; 131 M rows) | `cceb32c73be9c092923e8d123bf8dd74a676afe2bea845e002f33f296de839f5` | 4.1 GB; not for commit, the manifest pins it |
| `results/release/kinase_top200.tsv.gz` (503 × 200) | `98f348ed330784c8dbe52f9d327a14226975a74ab07b1e7612ebb6cf58d64567` | 1.9 MB |
| `results/release/site_top5.tsv.gz` (350,736 × 5) | `61e20f51b4eb043365ef492a0f20c31168f30bc6bacd637d40154fd0a8c3729f` | 34.8 MB |
| `results/release/RELEASE_MANIFEST.json`, `results/release/coverage.tsv` | | |

**Coverage** (`results/release/coverage.tsv`). U_obs = 350,736 sites, of which 338,525 have no PSP kinase.

| method | Cov@0.5 (no-PSP) | Cov@0.8 (no-PSP) | sites in some kinase's top 1% (no-PSP) |
|---|---|---|---|
| R2+R3 final (deployment) | 29,886 (26,762) | **0** (0) | 162,783 (152,377) |
| random scores | 6 (6) | 0 | 329,349 (318,032) |
| Cantley-only (B1, own DEV-OOF recipe) | 30,258 (29,220) | 1,828 (1,780) | 170,065 (162,611) |

Limitations:
- All maps are fitted on DEV OOF; there is no OOF for the DEV+TEST model.
- Probabilities are PU lower bounds.
- assignment_prob is conditional on the site being a target of one of the 503 modelled kinases.
- The model has no ≥ 0.8 tier, and no Y-site assignment reaches 0.5.
- Cantley's ≥ 0.8 tier comes from the tail of its own isotonic map, fitted on few labelled sites: only 152 labelled
  sites reached p ≥ 0.5 on CV. It is not evidence of higher precision.
- Coverage counts compare candidate sets of equal size here (503 kinases for both), but they are not
  precision-verified on unlabelled sites.

Spot checks of site_top5: TP53 S15 → ATM (0.70) ✓; AKT1 S473 → PDPK1 (0.51), AKT1, MTOR (3rd);
EIF4EBP1 T37 → GSK3B, CDK1, MTOR (3rd); EIF2S1 S52 → PKC family (a miss: the true kinases are EIF2AK1–4).
