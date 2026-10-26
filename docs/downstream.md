# Downstream: kinase-activity inference on benchmarKIN (DESIGN §6.7, §16 A21)

Code: `src/ksp/downstream/` (`bench_data`, `libraries`, `inference`, `evaluate`, `api`). Scripts: `scripts/downstream/`
(`run_validation.py`, `run_predictions.py`, `selftest.py`). Results: `results/downstream/`. Caches of the mapped benchmark:
`cache/downstream/`. Run time: about 6 min for the full validation (4 methods, 22 libraries, 1000 iterations, 1 core).

Reference: Müller-Dott et al. 2025 *Nat Commun* 16:4771 (doi:10.1038/s41467-025-59779-y). Data: Zenodo 12566560
(`external/benchmarkin`) plus the 18-experiment tyrosine set, which ships only in the benchmarKIN GitHub repo
(`external/benchmarkin_github`, see DATA.md).
- The protocol was taken from the paper's Methods (Europe PMC full text) and from the benchmarKIN R sources and
  decoupler-py 1.6 `get_performances` / `mc_perm`. Those sources were read as text only.
- Nothing downloaded was executed. `.Rds`/`.rda` files are read with the pure-Python `rdata` parser (added to
  requirements.txt). `pyreadr` cannot read R lists.

## 1. Benchmark data → canonical IDs
| set | experiments | sites mapped | rule |
|---|---|---|---|
| Hernandez-Armenta (2017 compilation) | 94 with metadata (103 columns) | 90.7% of rows (55,304 sites) | `GENE\|S103\|ENSG\|ENSP`, no window. Candidate accession from the ENSG (via HGNC) or else the symbol; accepted if exactly one canonical has that residue at that position. This is the same symbol + position rule the paper uses. Of the losses, 9.1% are residue mismatches, i.e. old Ensembl isoforms. |
| Hijazi 2020 | 118 (HL60, MCF7) | 96.1% | symbol + position. Targets are the curated, TTD-based `;` lists. |
| tyrosine | 18 | 87.4% | symbol + position. |
| CPTAC tumour (BRCA, CCRCC, GBM, HNSCC, LSCC, LUAD, UCEC) | 95–122 tumours per cohort | 99.8% | 15-mer window on the ENSG candidates (`CanonicalStore.map_window`). |

- **Totals:** 230 perturbation experiments, the same as the paper's 230.
- **Duplicates:** several source rows mapping to one site_id are nan-averaged.
- **Kinase names:** resolved with `resolve_kinase_names(source='benchmarkin')`.
  - All 80 perturbation target names resolve.
  - In the tumour gold standard (top/bottom 5% kinase protein level), 111 non-protein kinases (AK*, HK*, DGK*, PIK3C2*, …) are dropped. That leaves 13,232 GS pairs on 400 kinases.
- **Mapping logs:** `results/downstream/benchmark_mapping.json`.

## 2. Libraries (`libraries.py`)
- **PSP_all:** human-kinase × human-substrate PSP rows, mapped as in `ksp.labels`. 14,351 pairs, 384 kinases.
  - **It includes TEST-protein pairs.** That is allowed for this benchmark only, because the benchmark evaluates libraries.
  - It is built in memory, never written to disk and never used by a model.
  - Residue class is not enforced, matching the paper.
- **PSP_dev:** the same definition restricted to DEV proteins (`load_dev('all')`, provenance PSP, native_human). 11,433 pairs, 355 kinases. It is a strict subset of PSP_all.
- **Cantley_top{N}:** top-N sites per kinase by the as-published percentile (`motif/cantley_unprimed_pct_published`), ties broken by the raw log2 score.
  - N ∈ {per-kinase PSP_all size, 50, 100}. For the PSP-size arm, kinases without PSP pairs are absent.
  - Universe is `_Uobs` (all 350,736 index sites) or `_measured` (the 92,054 sites measured in any benchmark).
- **KL15site_measured:** the paper's Kinase-Library convention, the top-15 kinases per measured site. **PSP_all+KL15site** is its union with PSP_all.
- **Shuffled:** `shuffled(lib, seed)` applies a random permutation of the library's distinct sites. Every kinase keeps its set size and every site keeps its degree. Seeds are [0..4] on PSP_all.
- **Predictions:**
  - `from_predictions` / `topn_from_long` take the top-N per kinase from a long parquet of `kinase_acc|kinase_idx`, `site_id|site_idx`, `score`.
  - `union(a, b)` builds the "ours ∪ PSP" arm.
  - `topn_from_long(..., exclude=psp_all())` builds the **"ours minus PSP"** arm: PSP pairs are removed *before* ranking, so the top N are all novel. `minus(a, b)` is the plain set difference.
- Library statistics, including the number of kinases with ≥5 measured targets: `results/downstream/library_stats.tsv`.

## 3. Methods (`inference.py`, checked against reference implementations in `selftest.py`)
Kinase k is scored in sample j only if it has ≥5 measured targets there (benchmarKIN `minsize`). Notation: x = the sample's observed values; S = sd(x), ddof 1; T = the kinase's measured targets; n = |T|; N = number of observed sites.

| method | definition |
|---|---|
| `zscore` | RoKAI z-score as in benchmarKIN `run_zscore`: Σ_T x / (S·√n). No centering. |
| `ksea` | (mean_T − mean_all)·√n / S |
| `wmean` | the target mean (decoupler raw wmean with unit weights; no permutation normalisation) |
| `ulm` | t-value of the slope in x ~ 1 + 1[T], computed in closed form |

ULM and KSEA give AUROCs that agree to 3 decimals here: with binary weights and n ≪ N they are near-monotone transforms of each other. For the tumour benchmark, sites on the kinase's own protein are removed from its target set, as in the paper and vignette.

## 4. Evaluation (`evaluate.py`)
**Perturbation AUROC** (benchmarKIN `run_perturbBench` → decoupler `get_performances`):
1. Divide activities per experiment by the RMS: R `scale(center=FALSE)` = sqrt(Σx²/(n−1)).
2. Multiply by the perturbation sign and set 0 to NaN.
3. Flatten over experiments. Positives are (experiment, perturbed kinase); negatives are every other scored cell.
4. `mcauroc`: 1000 draws of n_pos negatives *with replacement*, seed 42. The mean is reported (median also stored).

`rm_bg=True` (the headline setting) keeps only experiments in which a perturbed kinase was scored; this is the paper's "experiments with ≥5 targets of the perturbed kinase". `rm_bg=False`, the package default, is also stored.

**Scaled rank** (`run_rank`, average=T):
1. Per experiment, rank by sign-adjusted activity (1 = best).
2. Divide by the number of scored kinases.
3. Average per perturbed kinase.
4. Report the mean (paper) and the median (vignette) over kinases.

`PHit(k)` is the fraction of kinases whose averaged rank is ≤ k.

**Tumour AUROC** (`benchmarkROC`):
- Per cohort, keep kinases with ≥30 scored tumours and z-score each across tumours.
- Pool GS+/GS− values over cohorts, for kinases with ≥1 GS+ and ≥1 GS− in that cohort.
- 1000 draws of ⌊0.8·n⌋ GS+ and ⌊0.8·n⌋ GS− values *without replacement*; report the mean.

**Modes:**
- `native` = each library on its own evaluable set. This is the paper protocol.
- `matched:<group>` = the cells or kinases scored by every member of the group (DESIGN matched-universe rule). Use this for head-to-head claims.

## 5. Validation (z-score; `results/downstream/validation_summary.tsv`)
| library (native) | pert. AUROC | scaled rank (mean) | PHit 5/10/20 | pert. kinases / exps | tumour AUROC | tumour kinases |
|---|---|---|---|---|---|---|
| PSP_all | **0.782** | **0.252** | 0.29/0.46/0.83 | 41 / 147 | **0.664** | 73 |
| PSP_dev | 0.784 | 0.264 | 0.35/0.54/0.92 | 37 / 140 | 0.666 | 61 |
| PSP_all+KL15site | **0.650** | 0.390 | — | 77 / 217 | **0.583** | 297 |
| KL15site_measured | 0.628 | 0.406 | — | 77 / 217 | 0.581 | 297 |
| Cantley_topPSP U_obs / measured | 0.623 / 0.639 | 0.44 / 0.44 | — | 22 / 45 | 0.552 / 0.558 | 64 / 109 |
| Cantley_top50 U_obs / measured | 0.620 / 0.503 | 0.40 / 0.49 | — | 13 / 55 | 0.540 / 0.546 | 126 / 254 |
| Cantley_top100 U_obs / measured | 0.529 / 0.548 | 0.48 / 0.49 | — | 37 / 66 | 0.544 / 0.555 | 166 / 291 |
| PSP_all shuffled (5 seeds) | 0.518 ± 0.028 | 0.502 | 0.01/0.06/0.40 | 34 / 128 | 0.531 ± 0.019 | 86 |

**Agreement with the paper:**
- **PSP with z-score.** The paper reports AUROC 0.79, mean scaled rank 0.24, PHit 0.27/0.43/0.89 (k = 5/10/20), 37 kinases, and 135 experiments. We get 0.782, 0.252, 0.29/0.46/0.83, 41 kinases and 147 experiments.
  - On Hernandez + Hijazi only we get 0.792, 0.225, **37 kinases**, 132 experiments.
  - The remaining difference is plausibly the PSP release: theirs was accessed 2023-04, ours is 2025.
- **Tumour benchmark.** The paper reports 0.662 for z-score + PSP on unnormalised data; we get 0.664.
- **PSP + Kinase Library (top-15 kinases per site).** The paper reports 0.65 (perturbation) and 0.58 (tumour); we get 0.650 and 0.583.
- **Cantley alone.** The paper did not report it. Our Cantley-only arms are consistently below PSP, and so is the paper's finding that adding KL targets lowers performance.
- **Shuffled control.** It sits near chance. Its residual excess, about 0.52–0.53, comes from library size: |z| grows with √n, and well-studied perturbed kinases have large sets.

Per dataset (PSP_all): Hernandez 0.839, Hijazi 0.756, tyrosine 0.656 (5 kinases, 15 experiments; small).

**Method ranking** (PSP_all, perturbation / tumour):

| method | perturbation | tumour |
|---|---|---|
| zscore | 0.782 | 0.664 |
| ksea | 0.776 | 0.654 |
| ulm | 0.776 | 0.654 |
| wmean | 0.765 | 0.663 |

This matches the paper's finding that z-score is best.

**Matched comparisons:**

| group | cells | PSP_all | other members |
|---|---|---|---|
| psp | 20 kinases / 82 experiments | 0.754 | PSP_dev 0.768; shuffled 0.524 |
| cantley_measured, perturbation | 31 kinases / 80 experiments | 0.730 | KL15 0.660, PSP+KL15 0.689, Cantley_topPSP 0.598 |
| cantley_measured, tumour | 68–70 kinases | 0.662 | KL15 0.595 |

**Evaluable kinases** (z-score, native; full lists in `evaluable_kinases.tsv`; per-kinase ranks in `validation_per_kinase.tsv`):
- **Perturbation, PSP_all (41):** ABL1, AKT1–3, ALK, ATM, ATR, AURKA/B, BRAF, CDK1/2/4/6, CHEK1, CSNK2A1, EGFR, FYN, GSK3A/B, LCK, MAP2K1, MAPK1/3/8/9/12/14, MTOR, PAK1, PLK1, PRKACA, PRKCA/B/G/Z, PRKD1, ROCK1, RPS6KB1, SRC, ZAP70.
- **Perturbation, other libraries:** up to 77 perturbed kinases are evaluable with the KL15 libraries, which add e.g. CSK, DDR1, EPHA2, FLT1/3/4, JAK1–3, KDR, MET.
- **Tumour:** 73 kinases with PSP_all, and up to 297 with dense predicted libraries.

## 6. Running a model library (later)
```
python scripts/downstream/run_predictions.py --pred <clean_preds.parquet> --tag <run_id> [--universe U_obs|U_ms|measured] \
       [--methods zscore,ulm] [--experiments-file <exps.txt>]   # leave-drug-out: restrict to the held-out drugs' experiments
```
Python API: `api.prediction_arms(path, tag)` gives the top{PSP,50,100}, `+PSP` and `-PSP` arms; then `api.run_benchmarks(libs, methods, matched={...})`. Outputs go to `results/downstream/<tag>_{summary,per_kinase,library_stats}.tsv`. Experiment ids look like `hernandez:102_117`, `hijazi:HL60.MK2206`, `tyrosine:QSAV1`; metadata is in `bench_data.load_perturbation().experiments`.

**Leakage rule (DESIGN §6.7):**
- The parquet must come from the **clean** model variant: phospho features without the `clean_bench` PMIDs (`results/features/phospho/exclude_lists.json`: 22 Hernandez + Hijazi compendium PMIDs) and without any CPTAC-derived feature. The tumour benchmark uses 7 of our 10 CPTAC cohorts.
- The 6 tyrosine-set PMIDs are not in the compendium, so the exclusion list needs no update.
- The harness cannot verify provenance. Record the model run id in `--tag`.

**The "ours" library also contains PSP_dev pairs.** A DEV-trained model ranks its training labels highly. Compare against PSP_dev, and report the `-PSP` arm, which removes all PSP pairs before ranking. PSP_all, the `+PSP` arm and the PSP reference arms all contain TEST labels; that is allowed for this downstream benchmark.

**Universe caveat.** Top-50 over U_obs (350k sites) leaves few *measured* targets. Only 13 perturbed kinases are evaluable for Cantley_top50_Uobs. Report U_obs (the pre-registered library) and `measured` side by side. `measured` uses only which sites were quantified, never their values.

## 7. Assumptions and limitations
- Hernandez, Hijazi and tyrosine sites are mapped by symbol + position, with no window; roughly 9% of Hernandez rows are lost to isoform or version shifts. Hijazi FDR filtering is as distributed.
- PSP_all uses the 2025 PSP release; the paper used the 2023 one. The paper's "curated" library also had PTMsigDB and GPS-gold, which we don't build.
- Benchmark experiments may be among the HTP evidence that PSP curates. This affects every PSP-based arm, as it did in the paper.
- The scaled rank uses stable ordering for ties. The mcauroc seed is 42 (decoupler) and the tumour seed is 1 (benchmarKIN). The numpy RNG differs from R's, so subsampled values differ from R in the third decimal.
- "Leave-drug-out" (A21) is supported only by restricting the experiments plus a matching predictions file. The leave-drug-out model itself belongs to the R2 team.
