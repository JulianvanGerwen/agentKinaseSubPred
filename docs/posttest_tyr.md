# Post-TEST exploratory: tyrosine-kinase specificity (DESIGN A34, target (c))

Owner: feat-tyr. **DEV CV only.** The released/evaluated model is unchanged (the frozen final
configuration, LESSONS 14); nothing here touches `sealed/`, `FREEZE.json`, the ledger or any existing
run artefact. These results belong in REPORT §11 as post-TEST exploratory, CV evidence only.

| Path | What |
|---|---|
| `scripts/models/posttest_tyr/tyr_fusion_explore.py` | the whole experiment (read-only use of the R2/R3 CV runs) |
| `results/models/posttest_tyr/explore.json` | all numbers reported below |
| `results/models/posttest_tyr/per_kinase_Y.tsv` | per-Y-kinase E-hard AUROC, baseline vs A vs B |
| `results/models/posttest_tyr/selection_trace.tsv` | every inner-fold selection score, for audit |

## 1. Motivation

- TEST: Y-kinase M1 = 0.595 against the B0-pop floor of 0.585, beating it on 10 of 23 kinases;
  T3 Y = 0.531; Cantley Y = 0.496 (`docs/audit/audit_test.md` §4).
- CV: Y M1 ≈ 0.63 against B0-pop 0.592 (`docs/strata_cv.md` §4).
- `docs/features_tyr.md` §1.5 found that the Cantley Tyr matrices of FGFR1/2/3, BTK and the Eph,
  VEGFR and TAM/Axl families are *anti-correlated* with their own labels (row means 0.39–0.49),
  while §1.4 found the Y E-hard pool is dominated by promiscuous pTyr docking hubs.

Two cheap interventions follow directly, both operating on the existing pooled-CV out-of-fold
predictions of the two fusion components, so neither needs a retrain.

## 2. Setup

Baseline is the frozen final fusion reproduced exactly: `score = 0.6 · z(R2) + 0.4 · z(R3)`, where
`z()` standardises each component's logits within CV fold over **all** rows, as in
`scripts/models/r2_fusion.py`. Components:

- R2 = `r2_final_ens5__20261026-191213__s-1` (the R2 side of the frozen fusion),
- R3 = `r3_w7_nce05_cpt_aux_cv_pooled__20261026-162205__s0`.

106,457,574 pooled OOF rows (503 kinases + 13 canaries × class-matched DEV U_ms sites).

**Variant A — Y-specific fusion weight** (brief (ii)+(iii)). On Y-class kinase rows only,
`score = w_Y · z(R2) + (1 − w_Y) · z(R3)`, with `w_Y ∈ {0, 0.2, …, 1.0}`. ST rows keep w = 0.6.

**Variant B — Tyr motif repair** (brief (i)). On Y-class kinase rows only,
`score' = score_fused + γ · z(C_k)`, where `C_k` is kinase *k*'s own Cantley Tyr percentile row if
that row ranks *k*'s positives above its E-hard negatives, and otherwise the row of *k*'s nearest
kinase-domain neighbour that owns an in-class Tyr matrix; if neither row carries signal, *k* is left
alone. `γ ∈ {0, 0.05, 0.1, 0.2, 0.3, 0.5}`, shared across kinases.

### 2.1 Nesting (how labels are kept out)

For each outer CV fold *f*, **every** tunable — `w_Y`, `γ`, and for B the per-kinase own/neighbour/none
decision, i.e. the "is this matrix anti-correlated?" judgement itself — is chosen using only the four
inner folds, then applied to fold *f*. The held-out fold's labels never enter the choice. The
per-fold choices are logged in `explore.json` and `selection_trace.tsv`.

**Inner-fold selection metric:** the mean over the four inner folds of the macro per-kinase
*within-fold* E-hard AUROC over Y kinases. Within a single fold the A24 mid-rank transform is
monotone, so a within-fold AUROC on raw scores equals the official one exactly; only the aggregation
differs from pooled-OOF, and it matches the `M1_fold_mean` convention already used in the fusion
sweep. This was chosen so that the selection criterion cannot drift from the reporting metric by a
rank-normalisation artefact.

### 2.2 Reporting metric

The official path: `ksp.eval.prepare('cv')` → `ksp.eval.evaluate(..., missing='error')` with A24
fold-rank normalisation on, then the Y stratum read off `per_kinase`. CIs and the paired test use
`ksp.eval._boot_matrix` with B = 1000 substrate-cluster resamples, seed 0 — the same route
`docs/strata_cv.md` uses for stratum CIs. The paired bootstrap is restricted to the Y kinases of
K_cv, on identical resamples for baseline and variant.

Variant B modifies **all** of a Y kinase's candidate rows, not only its positives and E-hard
negatives, because A24 ranks within (fold, kinase) over the whole candidate set; modifying a subset
would have made the pooled ranks inconsistent.

Both variants touch Y-class kinase rows only, so ST per-kinase metrics are unchanged by
construction. The script asserts this rather than assuming it (`ST_drift` in `explore.json`).

## 3. Results

Headline Y stratum = the 38 `residue_class == 'Y'` kinases of K_cv, matching `docs/strata_cv.md` §4.
The interventions act on all 42 Y+STY-class kinases, so that wider set is given as a secondary row.
Paired substrate-cluster bootstrap, B = 1000, seed 0, on identical resamples.

| | overall M1 | ST M1 | **Y M1 (38)** | Y Δ vs baseline [95 % CI] | p_boot | Y+STY M1 (42) |
|---|---|---|---|---|---|---|
| baseline (frozen fusion, w = 0.6) | 0.7864 | 0.8208 | **0.6387** [0.5925, 0.6674] | — | — | 0.6544 |
| A. Y-specific fusion weight | 0.7859 | 0.8208 | 0.6347 | **−0.0039** [−0.0097, +0.0010] | 0.128 | 0.6520 |
| B. Tyr motif repair | 0.7873 | 0.8208 | 0.6429 | **+0.0042** [−0.0103, +0.0167] | 0.549 | 0.6585 |

ST M1 is **bit-identical** across all three (`ST_drift` = 0.0), as the construction requires.
Overall M1 moves by −0.0005 (A) and +0.0009 (B), both inside the A34 tolerance of 0.002.
For reference, B0-pop on this stratum is 0.592 (`docs/strata_cv.md` §4): the baseline's lower CI
bound, 0.5925, only just clears the popularity floor.

### 3.1 What the nested selection chose

- **A:** `w_Y` = 0.6 (i.e. no change) in 3 of 5 outer folds, 0.4 and 0.3 in the other two.
- **B:** γ = 0.5 in all 5 folds. Repair source per kinase in fold 0: own matrix 16, nearest
  neighbour 14, nearest cross-family neighbour 5, dropped 7 (similar in the other folds). So the
  inner folds replace the matrix for **roughly 45 % of Y kinases**.

### 3.2 The decisive diagnostic: inner gains do not transfer

`results/models/posttest_tyr/inner_outer_transfer.tsv`, per outer fold, inner-fold gain of the
selected setting vs the gain actually realised on the held-out fold:

| variant | mean inner gain | mean held-out gain | folds improved |
|---|---|---|---|
| A | +0.0006 | −0.0049 | 0 / 5 |
| B | +0.0187 | +0.0056 | 3 / 5 |

B's inner-fold gain is positive and sizeable in **all five** folds (+0.012 to +0.027) and shrinks by
about 70 % out of fold, with two folds going the wrong way (−0.013, −0.015). With 38 kinases and
~1,150 Y positives, this stratum cannot support the selection of even one or two parameters: the
apparent gain is mostly selection noise. Had the own/neighbour choice and γ been fitted on all folds
and reported, variant B would have looked like a ~+0.019 Y gain. It is not one.

### 3.3 Per-kinase and per-family (variant B)

B improves 18 of 38 Y kinases. The gains land where `features_tyr.md` §1.5 predicted — the families
whose own atlas matrix is anti-correlated with their labels:

| family | n | n_pos | baseline | B | Δ |
|---|---|---|---|---|---|
| FGFR | 3 | 28 | 0.490 | 0.550 | **+0.060** |
| Abl | 2 | 139 | 0.590 | 0.628 | +0.038 |
| Jak | 3 | 43 | 0.601 | 0.625 | +0.024 |
| Src | 9 | 601 | 0.667 | 0.662 | −0.005 |
| Tec | 3 | 37 | 0.574 | 0.543 | −0.031 |
| Axl | 2 | 17 | 0.647 | 0.594 | −0.053 |

Largest single gains: TNK2 +0.119 (cross-family neighbour), JAK3 +0.118 (neighbour), FGFR2 +0.109,
FGFR3 +0.082. Largest losses: MERTK −0.098, HCK −0.067, JAK2 −0.059. Every one of these kinases has
5–21 positives, so none of the individual numbers is meaningful on its own.

### 3.4 A bug worth recording

The first run of variant B silently tested nothing: `ksp.labels.kinase_domain_identity()` returns a
**long edge table** (`acc_a, acc_b, identity, method`), not a square matrix, so the neighbour lookup
returned `None` for every kinase and only the "own matrix or drop" arm existed. With that bug the
result was Y Δ = −0.0047; after the fix, with genuine neighbour arms, it is +0.0042. The first
number is in the git-less history only; the artefacts in `results/models/posttest_tyr/` are from the
fixed run. Anyone reusing `kinase_domain_identity()` should note the shape.

Related: the nearest kinase-domain neighbour is usually in the *same* Manning family (FGFR1→FGFR2,
AXL→MERTK, KDR→FLT4), which tends to share the same defect. The cross-family arm (`nbr_xfam`) was
added for that reason and is what carries the FGFR2/FGFR3/TNK2 gains.

## 4. Verdict against the A34 bar

A34: promising = nested-CV gain on the target stratum with a cluster-bootstrap CI above 0, and no
loss on overall M1 beyond 0.002.

- **A (Y-specific fusion weight): fails.** The point estimate is negative (−0.0039) and the inner
  folds mostly re-chose the frozen w = 0.6 anyway. R2 and R3 are close enough on Y (0.629 vs 0.624)
  that there is nothing for a weight to exploit. Also covers brief idea (iii): a TK/TKL-group
  weight is the same instrument, and a Y-only weight is its most favourable case.
- **B (Tyr motif repair): fails the bar, but is the one worth keeping on the list.** +0.0042
  [−0.0103, +0.0167] is the right sign, lands on the predicted families, costs nothing at inference
  and does not touch ST — but the CI straddles 0 and the inner→outer shrinkage shows most of the
  apparent effect is selection noise.

**Neither is adopted.** The released model is unchanged.

### 4.1 What I would do next, and what I would not

- **Not** more score-level Tyr post-processing. Two reasonable instruments both landed inside the
  noise, and §3.2 shows the Y CV stratum is too small to select on. The next honest test of the
  motif-repair idea needs a *pre-registered* rule (e.g. "use the nearest cross-family atlas for any
  kinase whose own matrix has inner-fold AUROC < 0.5") fixed before looking, evaluated on new data
  or a new split — not another CV search.
- The more promising direction remains the one from `features_tyr.md` §1.6, which is about the
  endpoint rather than the model: the Y E-hard pool is 810 sites dominated by promiscuous pTyr
  docking hubs, and a pre-registered hub/A-loop sensitivity analysis would tell us how much of the
  0.639 is even addressable. `ksp.strata` already has `hub_sites`; the Y hub-excluded row in
  `strata_cv.md` (0.627 vs 0.629 for R2) suggests the answer is "not much", which would mean the Y
  ceiling is lower than the ST one for reasons no model can fix.

## 5. Reproduction, cost, limits

```
/workspace/.venv/bin/python -I scripts/models/posttest_tyr/tyr_fusion_explore.py --threads 8
```
693 s wall, peak ~16 GB RSS, 8 threads (within the ≤8 thread / ≤30 GB budget). Nothing was
committed; no frozen artefact, `FREEZE.json`, ledger or `sealed/` file was read or written.

Limits to keep in mind when reading §3:
- The inner-fold selection metric aggregates as a fold mean, while the reported metric is pooled
  OOF; these differ slightly by construction (§2.1).
- B's γ is shared across kinases while the matrix choice is per kinase; a per-kinase γ was not
  tried, and would only add selection freedom to a stratum that already cannot support it.
- The A24 rank transform is applied by `evaluate`; variant B shifts all of a Y kinase's candidate
  rows so the transform stays well defined, but the E-all-based endpoints (M3) were not examined
  here at all — only M1.
