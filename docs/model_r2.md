# R2 — global kinase–site pair model (LightGBM)

Owner: R2 agent. Status: **in progress** (updated incrementally; latest numbers at the top of §6).
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

Real-kinase features are unaffected, except for a possible tiny shift in F2a's internally tuned α/β, which are tuned
over all kinases including canaries. Stage-(a) selection therefore stands. Every run after stage (a) uses v2.

## 5. Ablation curve and leakage alarms
(filled as run)

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

## 6. Results
(filled as run)
