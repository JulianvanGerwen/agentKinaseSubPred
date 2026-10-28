# Downstream: kinase-activity inference on benchmarKIN (DESIGN §6.7, §16 A21)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Code: `src/ksp/downstream/` (`bench_data`, `libraries`, `inference`, `evaluate`, `api`, `paired`, `ldo`). Scripts: `scripts/downstream/`
(`run_validation.py`, `run_predictions.py`, `pred_topk.py`, `run_final.py`, `selftest.py`). Results: `results/downstream/`. Caches of the mapped benchmark:
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

## 7. Final pre-registered analysis: R2 clean library (2026-10-26 22:10)
**Inputs.**
- **Model:** `r2_final_cleanbench_uobs__20261026-204513__s0` (docs/model_r2.md §6.2). It is DEV-trained; its phospho features exclude the 22 `clean_bench` PMIDs and every CPTAC-derived feature; CV M1 is 0.771.
- **Library source:** its exhaustive U_obs predictions (131,144,410 rows; `cache/runs/<id>/uobs_predictions.parquet`).

**Pipeline.**
1. `scripts/downstream/pred_topk.py --tag r2_clean` streams the parquet with pyarrow `iter_batches` (4 M rows per batch, 30 s, a few GB of RAM).
2. `scripts/downstream/run_final.py --tag r2_clean --methods zscore`, and the same with `ulm`. Each takes about 4 min.
3. Outputs are in `results/downstream/final_r2_clean_{zscore,ulm}_{native,paired,coverage,per_kinase,ldo_strata}.tsv` and `r2_clean_topk_check.json`.

**Truncation check** (`r2_clean_topk_check.json`). R2's truncated file (`uobs_predictions_topk.parquet`, the per-kinase top 2·|PSP|+200 rows) contains every row every U_obs arm needs:
- There are 0 kinases with ties at the cut.
- All six arms (top{PSP,50,100}, with and without the `-PSP` exclusion) are identical whether built from the truncated file or from the streamed full ranking.
- The `-PSP` arms need at most rank 2·|PSP|+100, i.e. K−100.
- The `measured`-universe arms cannot come from the U_obs truncation. They were built from the stream as the top-K within the 92,054 benchmark-measured sites.

**Arms.**
- Pre-registered universe U_obs:
  - ours_top{PSP,50,100};
  - ours_top{N}+PSP (∪ PSP_all);
  - ours_top{N}-PSP (PSP_all pairs removed before ranking);
  - Cantley_top{N};
  - PSP_all and PSP_dev;
  - shuffled controls: PSP_all, ours_topPSP and ours_top50, 5 seeds each.
- The same model arms and Cantley arms are repeated in the `measured` universe as a sensitivity analysis.
- Primary method: z-score. ULM is the sensitivity method and gives the same conclusions to within ±0.02 (`final_r2_clean_ulm_*`).

**Leave-drug-out (A21).** No leave-drug-out *model* exists: it would need F6a/F6b with `exclude_drugs` and an R2 retrain, which belongs to the feature and R2 teams. As an approximation, `ksp.downstream.ldo` classifies every benchmark experiment by whether its drug appears anywhere in the clean model's compendium.
- **How drugs are matched:** the F6b union-find drug identity (Klaeger/Reinecke/Davis synonyms plus `drug_aliases.tsv`), with a substring fallback over Condition/Control text, applied to compendium conditions outside `clean_bench`.
- **drug_unseen** (114 experiments: 24 Hernandez, 88 Hijazi, 2 tyrosine): the clean model had no phosphoproteomic response data for that drug. This is the leave-drug-out stratum.
- **drug_seen:** 57 experiments.
- **no_drug:** 59 experiments (stimuli, cell cycle, DNA damage, genetic).
- **Caveat:** drug_unseen is dominated by Hijazi, so a drug effect cannot be separated from a dataset effect. The per-dataset rows are reported as well.

**Paired comparisons on matched kinase sets** (`ksp.downstream.paired`). A matched set contains only what both libraries score:
- perturbation: (experiment, kinase) cells, keeping experiments in which a perturbed kinase is scored by both;
- tumour: per-cohort kinases with ≥30 scored tumours in both, and GS pairs present in both.

Bootstrap: B = 1000 paired draws. The unit is the experiment for perturbation AUROC, the perturbed kinase for scaled rank, and the kinase for tumour AUROC. CIs are 95% percentile intervals. AUROC is the full AUROC, whose expectation equals the paper's balanced-subsample mcauroc. The CIs are unadjusted: downstream is a secondary, descriptive endpoint, and no superiority claim is pre-registered for it.

**Results** (z-score, U_obs, Δ = ours − reference with 95% CI, followed by the matched AUROC of the reference → ours and the matched n):

| comparison | perturbation, all 230 exps | leave-drug-out (drug_unseen) | tumour (CPTAC, 7 cohorts) |
|---|---|---|---|
| ours_topPSP vs PSP_all | −0.017 [−0.054, 0.018]; 0.772→0.755 (35 kin, 129 exp) | −0.058 [−0.123, 0.003] | +0.016 [−0.018, 0.051]; 0.680→0.696 (57 kin) |
| ours_top50 vs PSP_all | −0.013 [−0.043, 0.012]; 0.769→0.756 (33 kin, 106 exp) | −0.050 [−0.096, −0.001] | +0.011 [−0.027, 0.046] (54 kin) |
| ours_top100 vs PSP_all | −0.019 [−0.053, 0.012] (39 kin, 126 exp) | −0.041 [−0.099, 0.017] | +0.029 [−0.010, 0.071] (70 kin) |
| ours_top{PSP,50,100} vs PSP_dev | −0.023 / −0.021 / −0.014 (all CIs include 0) | −0.057 / −0.063* / −0.021 | +0.024 / +0.021 / +0.039 (all CIs include 0) |
| ours_top{PSP,50,100}+PSP vs PSP_all | +0.005 / +0.002 / +0.001 (all CIs include 0; 41 kin, 147 exp) | −0.015 / −0.006 / −0.011 | +0.015 / +0.015 / **+0.034 [0.005, 0.068]** (73 kin) |
| ours_top{PSP,50,100}−PSP vs PSP_all | **−0.045*** / −0.052 / **−0.068*** | −0.075* / −0.125* / −0.107* | +0.013 / +0.018 / +0.020 (all CIs include 0) |
| ours_topPSP vs Cantley_topPSP | **+0.126 [0.059, 0.190]** (22 kin) | −0.002 [−0.106, 0.102] (15 kin) | **+0.164 [0.095, 0.240]** |
| ours_top100 vs Cantley_top100 | **+0.176 [0.091, 0.263]** (36 kin) | +0.095 [−0.019, 0.208] | **+0.114 [0.078, 0.149]** (146 kin) |
| ours_top{N}−PSP vs Cantley_top{N} (PSP / 100) | **+0.103 / +0.119** | +0.013 / +0.079 | **+0.161 / +0.114** |
| ours_topPSP vs its shuffle (s0) | **+0.260 [0.180, 0.334]** | **+0.184 [0.052, 0.305]** | **+0.200 [0.136, 0.261]** |

Notes on the table:
- `*` = the 95% CI excludes 0. Bold marks the main CI-excludes-0 results.
- ours_top50 vs Cantley_top50 at U_obs is not interpretable: only 9 kinases and 25 experiments are matched, because Cantley's U_obs top-50 is almost never measured.
- In the `measured` universe the picture is the same. ours vs PSP_all on perturbation is +0.013 / −0.027 / −0.024 (n.s.). On tumour it is +0.037 / +0.018 / +0.038, with CIs excluding 0 for topPSP; versus PSP_dev it is +0.048 [0.012, 0.085]. ours vs Cantley is +0.13 to +0.22 on perturbation and +0.10 to +0.13 on tumour, all with CIs excluding 0.

**Per dataset** (ours_top50 vs PSP_all):

| dataset | Δ AUROC [95% CI] |
|---|---|
| Hernandez | +0.014 [−0.031, 0.058] |
| Hijazi | −0.042 [−0.087, 0.001] |
| tyrosine | +0.041 [−0.026, 0.114] |
| drug_seen | +0.016 |
| no_drug | +0.013 |

**Coverage** (`final_r2_clean_zscore_coverage.tsv`; native, each library on its own):

| library | evaluable perturbed kinases / experiments | native perturbation AUROC | evaluable tumour kinases | native tumour AUROC |
|---|---|---|---|---|
| PSP_all | 41 / 147 | 0.782 | 73 | 0.664 |
| PSP_dev | 37 / 140 | 0.784 | 61 | 0.666 |
| ours_topPSP | 38 / 141 | 0.744 | 59 | 0.692 |
| ours_top50 | 50 / 123 | 0.775 | 149 | 0.653 |
| ours_top100 | 65 / 152 | 0.749 | 274 | 0.645 |
| ours_top50+PSP | 60 / 168 | 0.787 | 172 | 0.659 |
| ours_top100+PSP | 66 / 181 | 0.770 | 282 | 0.653 |
| ours_top100−PSP | 62 / 141 | 0.707 | 275 | 0.646 |
| Cantley_top50 (U_obs) | 13 / 31 | 0.620 | 126 | 0.540 |
| Cantley_top100 (U_obs) | 37 / 85 | 0.529 | 166 | 0.544 |
| shuffled PSP_all (5 seeds) | 34 / 128 | 0.518 ± 0.028 | 86 | 0.531 ± 0.019 |

Overlap with PSP: 28% of ours_topPSP pairs, 13% of top50 and 9% of top100 are PSP_all pairs.

**Interpretation (pre-registered arms, z-score, U_obs):**
1. **Ours vs PSP.** On matched kinase sets the clean R2 library is statistically indistinguishable from PSP. Perturbation Δ runs from −0.013 to −0.019; tumour Δ from +0.011 to +0.029; every CI includes 0.
   - The R2 team's native contrasts (top50 0.775 on 50 kinases vs PSP 0.783 on 41; tumour topPSP 0.692 vs 0.664) are not like-for-like.
   - Matched, the tumour contrast is 0.680 → 0.696, +0.016 [−0.018, 0.051]. It is not significant.
2. **Coverage.** The model's libraries make many more kinases evaluable: 50–65 perturbed kinases vs 41, and 149–274 tumour kinases vs 73. At those extra kinases, performance stays well above chance.
3. **Adding predictions to PSP.** This leaves the perturbation benchmark unchanged (Δ ≤ +0.005) and helps the tumour benchmark slightly (+0.015 to +0.034; the top-100 CI excludes 0). The paper found the opposite for Kinase Library and Phosformer targets: adding them lowered performance, to 0.65 on perturbation and 0.57–0.58 on tumour. We reproduce their KL result (§5).
4. **Novel predictions only (−PSP).** These are worse than PSP on perturbation (−0.045 to −0.068) but equal to PSP on tumour. They beat Cantley (+0.10 to +0.12 perturbation, +0.11 to +0.16 tumour) and the shuffles.
5. **Leave-drug-out.** In experiments whose drug the clean model never saw (mostly Hijazi), ours trails PSP by 0.04–0.06; the top50 CI just excludes 0. On drug_seen and no_drug experiments it is level with PSP.
   - So part of the perturbation-benchmark parity depends on drug-response data for the same compounds from other studies.
   - On unseen drugs at U_obs, the model's advantage over Cantley also shrinks to n.s. (+0.00 to +0.10). In the `measured` universe it stays significant (+0.09 to +0.12).
   - A true leave-drug-out retrain would be needed to separate the drug effect from the Hijazi-dataset effect.

## 8. Assumptions and limitations
- Hernandez, Hijazi and tyrosine sites are mapped by symbol + position, with no window; roughly 9% of Hernandez rows are lost to isoform or version shifts. Hijazi FDR filtering is as distributed.
- PSP_all uses the 2025 PSP release; the paper used the 2023 one. The paper's "curated" library also had PTMsigDB and GPS-gold, which we don't build.
- Benchmark experiments may be among the HTP evidence that PSP curates. This affects every PSP-based arm, as it did in the paper.
- The scaled rank uses stable ordering for ties. The mcauroc seed is 42 (decoupler) and the tumour seed is 1 (benchmarKIN). The numpy RNG differs from R's, so subsampled values differ from R in the third decimal.
- Leave-drug-out (A21) is approximated by the drug_unseen stratum (§7), which is confounded with the Hijazi dataset. A true leave-drug-out run needs a model whose F6a/F6b exclude the benchmark drugs (`exclude_drugs`), evaluated with `--experiments-file`; that model belongs to the R2 team.
- The paired CIs are not adjusted for the many comparisons. Read them descriptively; downstream is a secondary endpoint.
