# Label variants and training regimes (regimes scientist)

Owner: regimes scientist. Status: **complete** (CV exploration, DEV only; `sealed/` was never read).

- **Harness:** frozen `ksp.eval` (commit 974d520). M1, M3, M4 and M5 use A24 within-(fold, kinase) rank
  normalisation; M2 uses the raw cross-kinase score.
- **Scope of numbers:** pooled 5-fold DEV OOF over K_cv = 203 kinases, seed 0, unless marked otherwise.
- **Code** (all new; R2 code is imported read-only):
  - `src/ksp/models/regime/{variants,run,perkinase}.py`
  - `scripts/models/regime_{xfit,run,perkinase,perkinase_eval,partial_pool,partial_pool_eval,resboost,resboost_eval,pergroup,pergroup_eval,report,queue}.py`
- **Results:** `results/models/regime/` holds:
  - `summary.tsv` and `summary_table.md`;
  - `evidence_strata.tsv`;
  - `ppool_lambda_{fold0,cv5}.tsv` and `resboost_lambda_{fold0,cv5}.tsv`;
  - `perkinase_subset_cv5.tsv`;
  - paired bootstraps `compare__*.tsv`.
- **Runs:** `results/runs/regime_*`. REGISTRY rows use model = `REGIME` (37 rows, 24 distinct configurations, plus 2
  extra seeds of L2ext).

## 0. Set-up and sanity check

### Base configuration
The base is the frozen final R2 configuration v3 (`r2_final_s0_cv5__20261026-181435__s0`):
- binary objective, 127 leaves, lr 0.03, α 0.5;
- negatives 1:50 with 50% hard;
- inner early stopping on real kinases;
- feature groups F1+F0+F2a(`f2a_sep`)+F3+F4+F6abc+F5+F7, kinase descriptor v2;
- seed 0.

knn_sim (Fknn) is **not** in the final configuration, because R2 did not adopt it, so it is not refitted here.

### Label-derived features
Three families in the final configuration are label-derived: F2a (learned PSSM, `f2a_sep`), F3-label (substrate
neighbourhood, `f3nb`) and F6a (co-regulation, `f6a`). All other groups are label-free.

Refitting **all three** per variant turned out to be affordable:
- they are rebuilt from the variant's training labels through `ksp.models.r2.xfit.build` → `ksp.crossfit.crossfit`;
- 54 stores, in `cache/regime/xfit/<variant>/`;
- 272 job-minutes in total, about 1.5 cores and 7 GB per job.

As a confound diagnostic, L1 and L2_vitro were also run with features **fixed at L2** (`L2fixed`).

### How variants are wired in
- Variants reach `sample.FoldData` and `crossfit.outer_split` through a process-local wrapper of
  `ksp.labels.load_dev` (`variants.patch()`). No frozen file and no R2 file was edited.
- Evaluation positives and masks are **always** the frozen `ksp.eval.prepare('cv')` sets: human-native PSP pairs,
  autophosphorylation excluded. Only the training labels change.

### Sanity check (passed exactly)
- **xfit stores:** my rebuild of the L2 stores on fold 0 is identical to R2's. F2a and F3: max|Δ| = 0. F6a: 3e-8.
- **Full run:** L2 through my driver gives M1 **0.768153**, which equals R2's `r2_final_s0_cv5` value of 0.768153. M2,
  M3, M4 and the best_iteration of every fold are identical too.

### Training-label variants (DEV, U_ms)
| variant | definition | pairs | kinases | kinases with ≥20 pairs |
|---|---|---|---|---|
| L2 | human + mouse/rat→human PSP (default) | 11,795 | 361 | 109 |
| L1 | human-native PSP | 9,234 | 344 | 91 |
| L2_vitro | L2 pairs with in vitro evidence (incl. "both") | 8,987 | 342 | 92 |
| L2_vivo | L2 pairs with in-vivo-**only** evidence (= L2 − L2_vitro); my reading of "in-vivo-only" | 2,808 | 257 | 31 |
| L2ext | L2 + filtered external curated pairs on DEV proteins (see below) | 12,972 (+1,177) | 375 | 120 |

L2ext adds the `ext_test` pairs:
- DL-5 near-copy filter, A9 PhosphoSite-source drop, cap of 20 pairs per PMID;
- 827 OmniPath and 350 UniProt pairs;
- tier A 750, tier B 427.

Canary labels are identical in every variant (all 428 carry every flag), so the canary gate is comparable across
variants.

## 1. Results (pooled 5-fold CV)
Δ is against the **matched-seed** R2 baseline:
- seed 0: R2 seed-0 run, reproduced bit-for-bit;
- seeds 1 and 2: R2's own `r2_final_s1_cv5` (M1 0.7734) and `r2_final_s2_cv5` (M1 0.7690).

M5 is macro E-hard AUROC over K_any, by DEV training-count bin. Bin 0 is empty for E-hard. M1 ST / M1 Y are K_any
macros by residue class. "adopt?" applies the mechanical pre-declared rule.

| config | M1 | ΔM1 | fold-mean | M2 | ΔM2 | M3 | M4 | M5 1-4 | 5-19 | 20-99 | ≥100 | M1 ST | M1 Y | canary neg/pos | rule |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **L2 (R2 global, baseline)** | 0.7682 | — | 0.7645 | 0.1923 | — | 0.2982 | 0.7602 | 0.7306 | 0.7591 | 0.7719 | 0.8063 | 0.7893 | 0.6255 | 0.547/0.999 PASS | — |
| L1 refit | 0.7711 | +0.0030 | 0.7688 | 0.1630 | −0.0294 | 0.2947 | 0.7569 | 0.7127 | 0.7610 | 0.7717 | 0.8072 | 0.7943 | 0.5873 | 0.527/0.999 PASS | no |
| L1, features fixed at L2 | 0.7676 | −0.0005 | 0.7656 | 0.1738 | −0.0186 | 0.2882 | 0.7576 | 0.7276 | 0.7544 | 0.7723 | 0.8048 | 0.7922 | 0.6016 | 0.541/0.999 PASS | no |
| L2_vitro refit | 0.7647 | −0.0034 | 0.7648 | 0.1453 | −0.0470 | 0.2760 | 0.7532 | 0.7261 | 0.7438 | 0.7726 | 0.8069 | 0.7910 | 0.5856 | 0.545/0.999 PASS | no |
| L2_vitro, features fixed at L2 | 0.7683 | +0.0001 | 0.7669 | 0.1731 | −0.0193 | 0.2820 | 0.7623 | 0.7491 | 0.7502 | 0.7754 | 0.8102 | 0.7954 | 0.6109 | 0.546/0.999 PASS | no |
| L2_vivo refit | 0.7392 | −0.0289 | 0.7394 | 0.1028 | −0.0896 | 0.2246 | 0.7217 | 0.6841 | 0.7155 | 0.7409 | 0.7881 | 0.7544 | 0.5693 | 0.547/0.999 PASS | no |
| L2ext refit, seed 0 | 0.7740 | **+0.0058** | 0.7692 | 0.2180 | +0.0257 | 0.3085 | 0.7584 | 0.7178 | 0.7603 | 0.7741 | 0.8063 | 0.7969 | 0.5853 | 0.549/0.999 PASS | met (see §4) |
| L2ext refit, seed 1 (replicate) | 0.7694 | −0.0040 | 0.7672 | 0.1995 | +0.0178 | 0.3066 | 0.7594 | 0.7323 | 0.7563 | 0.7723 | 0.8035 | 0.7928 | 0.6057 | 0.548/0.999 PASS | no |
| L2ext refit, seed 2 (replicate) | 0.7712 | +0.0022 | 0.7674 | 0.1979 | +0.0017 | 0.3062 | 0.7608 | 0.7376 | 0.7564 | 0.7720 | 0.8059 | 0.7962 | 0.6000 | 0.547/0.999 PASS | no |
| per-kinase LightGBM (N≥20, else global) | 0.7563 | −0.0119 | 0.7493 | 0.1735 | −0.0188 | 0.2592 | 0.7529 | 0.7306 | 0.7591 | 0.7492 | 0.7932 | 0.7805 | 0.6219 | 0.549/0.999 PASS | no |
| per-kinase L2-logistic (N≥20, else global) | 0.7553 | −0.0128 | 0.7481 | 0.1711 | −0.0212 | 0.2606 | 0.7526 | 0.7306 | 0.7591 | 0.7473 | 0.7962 | 0.7817 | 0.6165 | 0.542/1.000 PASS | no |
| partial pooling, ridge residual (λ = 300, selected on fold 0) | 0.7684 | +0.0002 | 0.7647 | 0.1990 | +0.0066 | 0.2984 | 0.7604 | 0.7308 | 0.7592 | 0.7721 | 0.8069 | 0.7895 | 0.6256 | 0.546/0.999 PASS | no |
| partial pooling, per-family residual boosting (λ = 300, selected on fold 0) | 0.7683 | +0.0002 | 0.7646 | 0.1945 | +0.0022 | 0.2981 | 0.7603 | 0.7306 | 0.7591 | 0.7722 | 0.8064 | 0.7894 | 0.6253 | 0.548/0.999 PASS | no |
| per-Manning-group LightGBM | 0.7607 | −0.0074 | 0.7606 | 0.1701 | −0.0223 | 0.2978 | 0.7496 | 0.7276 | 0.7354 | 0.7716 | 0.8028 | 0.7764 | 0.6253 | 0.547/0.999 PASS | no |

For per-kinase and per-group models, M2 is shown with score_site = **qmap**. This is a label-free, rank-preserving
map of each own model's scores onto the global model's score distribution for the same (fold, kinase). M2 on raw
logits:
- per-kinase LightGBM 0.164;
- per-kinase logistic 0.182;
- per-group 0.170.

Paired cluster bootstrap against the matched R2 baseline (`scripts/models/r2_compare.py`, B = 1000, matched
universe; files `results/models/regime/compare__*.tsv`). Values are Δ [95% CI] (p):

| comparison | ΔM1 | ΔM2 | ΔM3 |
|---|---|---|---|
| L2ext s0 vs L2 s0 | +0.0058 [−0.0006, +0.0121] (0.076) | +0.0257 [+0.0124, +0.0351] (0.002) | +0.0103 [−0.0007, +0.0214] (0.064) |
| L2ext s1 vs R2 s1 | −0.0040 [−0.0097, +0.0025] (0.21) | +0.0178 [+0.0050, +0.0278] (0.006) | +0.0038 [−0.0089, +0.0181] (0.56) |
| L1 refit vs L2 | +0.0030 [−0.0019, +0.0097] (0.19) | −0.0294 [−0.0462, −0.0148] (0.002) | −0.0035 [−0.0147, +0.0102] (0.74) |
| per-kinase LightGBM vs global | −0.0119 [−0.0168, −0.0085] (0.002) | −0.0188 [−0.0283, −0.0124] (0.002) | −0.0390 [−0.0495, −0.0326] (0.002) |
| L2_vitro refit vs L2 | −0.0034 [−0.0102, +0.0059] (0.53) | −0.0470 [−0.0627, −0.0283] (0.002) | −0.0221 [−0.0361, −0.0012] (0.038) |
| ridge partial pooling λ = 300 vs global | +0.0002 [−0.0001, +0.0005] (0.19) | +0.0066 [+0.0038, +0.0096] (0.002) | +0.0002 [−0.0027, +0.0022] (0.78) |

## 2. A — label variants

**L1 vs L2.** Human-only training labels give M1 +0.003 with refitted features and −0.0005 with L2-fixed features.
Both are inside seed noise (R2 seed SD 0.0028). M2 drops 0.019–0.029.
- Mouse and rat transferred pairs do not hurt M1 and help cross-kinase M2.
- Keep L2.

**L2_vitro / L2_vivo.**
- In-vitro-only training costs M1 −0.003 (refit) / 0.000 (fixed) and M2 −0.02 to −0.05.
- In-vivo-only training (2,808 pairs, 31 kinases with ≥20) is much worse everywhere: M1 −0.029, M2 −0.09.
- The two variants are size-confounded (8,987 vs 2,808 pairs), so L2_vivo mainly measures label count.

**Does training on in vitro labels help in vitro test positives?** Modestly, yes. Pair strata are macro E-hard AUROC
over K_cv. Test positives are human-native PSP pairs by native evidence class. Stratum sizes:
- both: 3,044 positives, 195 kinases;
- in vitro only: 3,292 positives, 181 kinases;
- in vivo only: 2,080 positives, 172 kinases.

| training labels | both | in vitro only | in vivo only | external curated |
|---|---|---|---|---|
| L2 (baseline) | 0.7784 | 0.7720 | 0.7595 | 0.7673 |
| L1 refit / fixed | 0.7860 / 0.7846 | 0.7820 / 0.7753 | 0.7567 / 0.7580 | 0.7654 / 0.7588 |
| L2_vitro refit / fixed | 0.7741 / 0.7848 | **0.7806 / 0.7804** | 0.7524 / 0.7603 | 0.7552 / 0.7612 |
| L2_vivo refit | 0.7406 | 0.7509 | 0.7363 | 0.7319 |
| L2ext s0 / s1 / s2 | 0.7813 / 0.7784 / 0.7801 | 0.7837 / 0.7804 / 0.7788 | 0.7638 / 0.7611 / 0.7619 | 0.7718 / 0.7727 / 0.7777 |

- L2_vitro raises the in-vitro-only stratum by **+0.008** with either feature mode.
- It leaves in-vivo-only test positives unchanged (fixed, +0.001) or lowers them (refit, −0.007).
- L2_vivo shows **no** matching in-vivo specificity: its loss is −0.023 on in-vivo-only test positives vs −0.021 on
  in-vitro-only.
- So evidence-matched training gives a small, one-directional gain (in vitro → in vitro) at no net M1 benefit. It does
  not justify changing the default.

**L2+ext.**
- On seed 0, adding the 1,177 filtered external DEV pairs gives M1 **+0.0058** (p = 0.076), M2 +0.026 (p = 0.002) and
  M3 +0.010. M1 is better on all 5 folds (+0.007 / +0.000 / +0.008 / +0.001 / +0.008).
- The pre-registered replicates do not reproduce the M1 gain: seed 1 gives −0.0040, seed 2 +0.0022, so the mean ΔM1
  over seeds 0–2 is **+0.0013**.
- R2's seed 0 is its lowest seed (0.768, against 0.769–0.775 for seeds 1–4), which accounts for part of the seed-0
  gap.
- The M2 (+0.026 / +0.018 / +0.002) and M3 (+0.010 / +0.004 / +0.006) gains are positive in all three seeds.
- Y-kinase M1 is lower in every seed (0.585–0.606 vs 0.61–0.63).

**Refit vs fixed features.** The confound is real but small, ±0.004 M1, and its sign is not consistent:
- L1: refit better (+0.0035);
- L2_vitro: fixed better (+0.0036; the PSSM/F3/F6a features fitted on fewer labels are weaker).

Refitting is necessary for an unconfounded label comparison. Features fixed at L2 flatter label-poor variants.

## 3. B — training regimes

### Per-kinase models
N = 20 training positives was chosen a priori:
- it equals the eval's own bin edge (20–99 / ≥100);
- with 1:50 negatives it gives ≥1,000 rows, and ≥4 inner-validation positives for early stopping.

Kinases below N fall back to the global model. The same rule applies to canaries, so the gate tests the regime.
Own models were fitted for 89–95 real kinases per fold; 78 K_cv kinases had ≥ N in all 5 folds.

| subset | global R2 M1 | per-kinase LightGBM | per-kinase logistic |
|---|---|---|---|
| ≥N in all folds (78 kinases) | **0.8046** | 0.7788 (−0.026) | 0.7794 (−0.025) |
| all K_cv (fallback to global for kinases below N) | **0.7682** | 0.7563 | 0.7553 |

| subset | global R2 M3 | per-kinase LightGBM M3 |
|---|---|---|
| ≥N in all folds | 0.319 | 0.232 |

Per-kinase models lose most exactly where they can be fitted. The global model's cross-kinase sharing is worth
about +0.025 M1 even for data-rich kinases.

### Partial pooling (§8.3)
`scripts/models/r2_partial_pool.py` was reviewed and found sound:
- training-row offsets are out of sample, from inner 3-fold cluster refits at the outer round count;
- the per-kinase ridge residual is shrunk to the Manning-family mean.

It was reused through per-fold wrappers. The pre-declared grid was λ ∈ {0, 3, 10, 30, 100, 300, ∞}, with selection on
fold-0 M1. A second form, per-family LightGBM residual boosting with the same λ shrinkage, was added.

| λ (pooled CV) | 0 | 3 | 10 | 30 | 100 | 300 (selected) | ∞ (global) |
|---|---|---|---|---|---|---|---|
| ridge residual M1 | 0.7610 | 0.7649 | 0.7671 | 0.7682 | 0.7684 | 0.7684 | 0.7682 |
| ridge residual M2 | 0.1993 | **0.2322** | **0.2316** | 0.2193 | 0.2062 | 0.1990 | 0.1923 |
| residual boosting M1 | 0.7655 | 0.7668 | 0.7676 | 0.7682 | 0.7683 | 0.7683 | 0.7682 |
| residual boosting M2 | 0.1940 | 0.2005 | 0.2013 | 0.1992 | 0.1978 | 0.1945 | 0.1923 |

Fold 0 selected λ = 300 for both forms, with fold-0 M1 0.7507 / 0.7505 against 0.7505 for the global model. Pooled
gains are +0.0002 M1, which is nothing.

**Where do per-kinase terms start to pay off?** Nowhere, for M1.
- Small λ hurts bins 5–19 and 20–99.
- Ridge λ = 0 helps bin 1–4 (0.748 vs 0.731), but those few kinases have fewer than 5 labels each, so this is noise.

**Unplanned finding (not selected; hypothesis for R2-main's calibration).** The per-kinase ridge **intercept** at
λ = 3–10 raises M2 by **+0.04**, from 0.192 to 0.232, for M1 −0.001 to −0.003. This is a cross-kinase calibration
effect: the residual mostly re-offsets kinase logits. It suggests per-kinase intercept calibration of `score_site`
(A16/A20) could help M2. It should be tested by R2-main on its own CV endpoint.

### Per-Manning-group models
- One LightGBM per Manning group, with the final R2 parameters. Groups with <200 training positives fall back to the
  global model (TKL), and so do canaries.
- Result: M1 0.7607 (−0.0074), M2 −0.02.
- Splitting the data by group loses more than group-specific capacity gains.

## 4. Adoption verdict (pre-declared rule)
**Rule:** pooled-CV ΔM1 ≥ +0.005 over the matched R2 baseline (same seed), and M2 not dropping by more than 0.01, and
the canary gate passing.

- **Label variants L1, L2_vitro, L2_vivo:** not adopted. They are reported as findings.
- **Regimes (per-kinase, both partial-pooling forms, per-group):** not adopted. The global R2 model stays the default.
- **L2ext:**
  - **The seed-0 comparison formally meets the rule:** ΔM1 +0.0058, ΔM2 +0.026, canaries pass.
  - Because the margin is at the threshold and inside seed noise, I ran seeds 1 and 2 against R2's own matched seeds
    before reporting. They give ΔM1 −0.0040 and +0.0022; the mean over three seeds is +0.0013.
  - **My recommendation to the PI is not to adopt L2ext as the default on M1 grounds:** the seed-0 gain does not
    replicate.
  - It is a candidate on M2/M3, which are positive in all three seeds: M2 +0.002 to +0.026, M3 +0.004 to +0.010.
    That is the PI's call (§8.6: M1 first, M2 as tie-break).
  - If adopted, note DESIGN §3: the external TEST pairs stay out of training, because they are on TEST proteins. But
    the training labels would then share curation sources with that secondary test set.

## 5. Runtime and resources
Wall time: 19:57 → 00:03 (about 4 h), on ≤ 24 cores and ≤ 3 concurrent LightGBM jobs.

| component | time |
|---|---|
| xfit refits (54 stores = 4 variants × 5 folds × 3 families, plus the 3 L2 fold-0 sanity stores; 2–4 concurrent, ~1.5 cores and ~7 GB each) | 272 job-min; wall 20:00–21:42 |
| label-variant training, 6 threads per job | 1.5–7 min per fold, 6–29 min per 5-fold run |
| pooled A24 evaluation per run | ~6 min |
| per-kinase models, per fold | 1.3–2.5 min (all 5 folds: 10 min) |
| ridge partial pooling (inner 3-fold R2 refits), per fold | 3.7–5.3 min (23 min in total) |
| residual boosting, per fold | 2.7–4.5 min (18 min in total) |
| per-group, per fold | 2.4–3.5 min (15 min in total) |
| each λ-sweep pooled evaluation (7 λ) | ~25–40 min |
| fold-0 pilot (L2 through my driver, including fold-0 eval) | 5.3 min train+predict, +1.5 min eval |

## 6. Deviations, incidents, caveats

**Memory.** The Slurm job cgroup is shared by all agents and capped at 128 GiB (`memory.max`). Its OOM killer fired
during this work:
- my per-kinase 5-fold evaluation was killed at about 20:40 and re-run;
- two other kills, at about 20:14 and 21:05 (~30 GB and ~21 GB anon), may have hit other agents' processes, because my
  footprint was ~60–67 GB at those times.

After that I capped my concurrency (xfit 2–3, ≤ 2 training jobs, evaluations sequential).
- **PI note:** "≤ 60 GB per agent" is close to the whole job's headroom once other agents run.

**Per-kinase segfault.** Fitting LightGBM models from several Python threads segfaulted in `lib_lightgbm.so`
(folds 1 and 4). Those folds were re-run sequentially. Every model uses `num_threads = 1`, so results do not depend on
the worker count.

**Residual boosting at λ = 0.** For kinases with n_k = 0 the weight n_k/(n_k+λ) is undefined; it is set to w = 0.

**qmap.** The scoring variant for per-kinase/per-group M2 was implemented after the fold-0 pilot. It is label-free.
Its speed fix (`np.interp` on the sorted reference, identical to `np.quantile`, max|Δ| 1e-16) changed no number.

**Evidence strata.** Test positives use **native** evidence (`evidence_native`). L2_vitro and L2_vivo are defined on HMR
evidence, as `ksp.labels` does.

**Selection optimism.** I tried 24 configurations. The best of 4 label variants on one seed (L2ext, +0.0058) is
exactly the kind of result that shrinks on replication, and it did.
