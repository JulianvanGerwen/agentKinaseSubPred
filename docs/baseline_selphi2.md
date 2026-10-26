# B4 — SELPHI2.0 retrained on our splits

Owner: B4-selphi2 agent. DESIGN v2 §8.1 B4, §6.8 (matched-universe rule), PRIOR-4 (`docs/reviews/r1_prior_work.md`).
Reference: Maier, Petursson, Lussana & Petsalaki, *MCP* 2025 24(8):100994 (`data/paper1.pdf`; the brief calls this paper2, but paper2.pdf is Phosformer). Code: `external/selphi2/selphi_2.0` (Apache-2.0). I read it as a reference and ported its logic; none of it was executed.

Code: `src/ksp/baselines/selphi2/{features,model}.py`, `scripts/baselines/selphi2_{build_features,cv,eval,pairsplit}.py`.
Cache: `cache/selphi2/` (OOF predictions in `cache/selphi2/<arm>/oof_fold{f}.parquet`; the ksp.manifest copy is in `cache/runs/<run_id>/predictions.parquet`).
Results: `results/baselines/selphi2/` (`<arm>/fold{f}.json` holds timings, grid table, importances and counts; `<arm>/summary.json`, `metrics_{matched,full}.json` and `per_kinase_matched.tsv` hold the metrics; `pairsplit_reference.json` holds the non-official reference).

## 0. Results: CV, DEV out-of-fold, 5 substrate-cluster folds, pooled; `ksp.eval`, 1,000-resample cluster bootstrap

| arm / universe | M1 AUROC E-hard (K_cv) | M2 Hit@5 | M3 top-1% recall | M4/M5 | AUROC E-all |
|---|---|---|---|---|---|
| **faithful, matched U_B4** (missing='drop'; the claim-bearing view) | **0.708** [0.685, 0.728] | 0.176 [0.158, 0.196] | 0.316 [0.272, 0.356] | 0.714 | 0.887 |
| faithful, full universe (unscored = worst) | 0.587 [0.568, 0.613] | 0.112 [0.100, 0.124] | 0.297 [0.263, 0.336] | 0.581 | 0.701 |
| noleak (drops functional_score, PWM_max_mss, netpho_max_KIN, isKinaseCoreg), matched | 0.724 [0.700, 0.742] | 0.190 [0.171, 0.211] | 0.265 [0.224, 0.302] | 0.725 | 0.875 |
| noleak, full universe | 0.589 | 0.119 | 0.272 | 0.583 | 0.698 |
| *reference:* B1 Cantley universe percentile (unprimed), same U_B4 pairs | 0.728 [0.711, 0.749] | 0.182 | 0.123 | 0.733 | 0.763 |

- **K_cv** = 203 kinases, all evaluated. M1 is defined for 202 of them on U_B4; 16 have fewer than 3 positives inside U_B4.
- **Re-filtered matched set (PRIOR-1).** Requiring ≥ 3 positives on ≥ 2 clusters inside U_B4 leaves 184 kinases: faithful M1 0.709, M2 0.185, M3 0.317; noleak 0.722 / 0.197 / 0.262.
- **Paired comparison against Cantley** on the 26.7M pairs both methods score: SELPHI − Cantley M1 = −0.020 [−0.048, −0.001], M2 +0.001 [−0.026, 0.024], M3 +0.187 [0.137, 0.224]. These are unadjusted bootstrap CIs, for orientation only; they are not part of the pre-registered superiority tests.
- **Reading.** Retrained honestly, SELPHI2.0 does **not** beat the label-free Cantley percentile on kinase-specific ranking (M1, M2). It does win on M3 / E-all, because its site-level Ochoa features (functional score, abundance, conservation, NetPhos) separate "well-studied phosphosite" from "any site". This is a popularity effect, which is why dropping the four leakage-suspect columns raises M1 by 0.016 and lowers M3 by 0.051.
- **Strata (faithful, matched).** M1 by group: AGC 0.78, CAMK 0.77, STE 0.75, CMGC 0.72, TKL 0.72, Other 0.69, CK1 0.65, **TK 0.60** (Tyr is SELPHI's weakest class; Y-class M1 0.61). By training bin: 1–4: 0.70, 5–19: 0.73, 20–99: 0.69, ≥100: 0.73. The full strata tables are in `faithful/summary.json`.
- **Grid winners per fold:** max_depth 20 in every fold; min_samples_split 8/10/12/8/12; n_estimators 1500/1000/1500/1500/1500. Inner AUROCs over the grid are within ±0.005, except depth 10, which is about 0.005 lower.
- **Feature importance (faithful, fold 0):** functional_score 0.15, score (learned PSSM) 0.15, pssm_ser_thr 0.12, rna.tiss / rna.cell / sel.1 about 0.04 each. Each coreg column is about 0.015.

### Coverage (`ksp.eval.coverage`, pooled DEV)
- **Sites:** B4 scores 89,419 of 277,924 DEV U_ms sites (32.2%).
- **Kinases:** all 503 index kinases.
- **Required pairs:** 23.5M of 69.0M scored (positives plus E-all candidates of every kinase with a positive) = **34.1%**.
- **Positive pairs:** 5,135 of 8,756 scored = **58.6%** (labelled sites are enriched in the Ochoa table).
- **Gap reason:** every gap is "site absent from the Ochoa 2020 feature table". There are no kinase gaps.

### Non-official reference: SELPHI's own pair-split protocol, one fold (`pairsplit_reference.json`)
**Clearly NON-OFFICIAL; not an endpoint.**
- **Setup.** All DEV L2 positives on U_B4 (6,971) plus 10× random masked pairs. StratifiedKFold(5) over **pairs**, fold 0. Their metric: pooled AUROC of positives against random negatives.
- **Result.** Pair split: AUROC **0.947**, AP 0.764. On the pair split, 73% of the held-out positive sites also occur in training rows. Our cluster split, fold 0, same metric: AUROC **0.947**, AP 0.760. **The inflation from pair splitting is ≈ 0.000 here.**
- **Caveat.** In this reference the label-derived columns were still cross-fitted, as the project rule requires. The extra optimism SELPHI's protocol gets from computing `score` / `coreg` from all labels (the target pair sits in its own PSSM and regulon) is therefore *not* included. So 0.947 is a lower bound on what their protocol would report; they report 0.97.
- **Where the gap comes from.** The gap from their headline (AUC 0.97) to our M1 (0.71) comes from the **evaluation design**, not the split. Random-pair negatives are mostly unstudied sites, which the site-level features recognise. Per-kinase E-hard negatives (sites labelled to other kinases) remove that shortcut. On our data the same model gives 0.95 on "positives vs random pairs" and 0.71 on M1.

### Runtime (12 cores, no GPU)
- **Faithful arm:** 109 min wall for 5 folds, about 22 min per fold:
  - masks and loading, about 1.1 min;
  - negative sampling, 0.4 min;
  - cross-fitted label features, about 1.8 min;
  - grid search (72 configurations × 3 inner folds, realised as 9 × 3 fits of 1,500 trees), about 5.6 min;
  - 10 bag fits, about 3.5 min;
  - prediction of about 7M pairs per fold with 10 × 1,000–1,500 trees, about 10 min (about 11k pairs/s).
- **noleak arm:** about 6 min per fold (5 bags × 500 trees, grid winner reused).
- **Feature build:** about 1 min (`selphi2_build_features.py`).
- **Evaluation:** about 11 min per arm (single thread; includes the bootstrap and the Cantley reference).
- **Pilot:** fold 0, 1 bag, 300k eval rows, 3 min.
- **Whole B4 effort:** about 2 h 50 min of wall time (14:36–17:25).

### Run manifests (`ksp.manifest`)
- `results/runs/b4_selphi2_faithful__20261026-164243__s0/` (predictions: `cache/runs/<same>/predictions.parquet`, 35.6M rows)
- `results/runs/b4_selphi2_noleak__20261026-171245__s0/`

The PI can recompute every number from these predictions with `ksp.eval`.

### Deviations from the published method (all deliberate; details in §2)
1. Our substrate-cluster 5-fold CV replaces their pair-level StratifiedKFold(10). The inner grid-search CV is GroupKFold(3) over clusters.
2. 10 bags instead of 100; the grid search runs once per outer fold, not once per bag.
3. Negative ratio 10× (the paper says 10/20/50× are equivalent; their code uses 1× when there are > 10k positives), with our masks applied.
4. The label-derived `score` and `coreg*` columns are cross-fitted. `score` is our F2a learned PSSM (substituted for their `pssm.py`). The KSEA permutation test is replaced by the KSEA z-score. Spearman is approximated by ranking each profile over its own observed conditions.
5. `coreg293` uses CPTAC BRCA (BCM v1.2, protein-residualised) in place of Mertins 2016. GTEx v8 was downloaded. HPA is v25 (theirs was v23).
6. Universe U_B4 = U_ms ∩ Ochoa. SELPHI's own scoring of non-Ochoa sites with imputed defaults was not run.
7. Labels are L2 (human + mouse/rat) on U_ms, not PSP human-only.
8. Canaries were not run through B4. The PI's note on canary v2 arrived mid-run. Including canaries would need the B4 feature builders to support the 13 extra kinase rows, which the time box did not allow.
9. The noleak arm uses 5 bags × 500 trees, not 10 × 1,000–1,500. Its bags 0–4 use the same negative draws as the faithful arm.


## 1. Feature table: 44 model columns (residue + 43)

SELPHI's released feature table (`data/psp/pos_set_features_20240910.tsv`) has **44 model columns**. The paper says "48 generated, 45 selected". Its Supplemental Table 1, which would list the selected features, is not in the repo, and the code (`learn.py`) trains on exactly these 44 columns. So we reproduce the 44 and record the paper-versus-table discrepancy of one column. PRIOR-4 counted "43" because it left out `residue`. It also called `score`, `sel.1` and `sel.2` "undocumented". They are resolved here from the code:

- `score` is the PSP-derived PSSM built in `pssm.py`.
- `sel.1` and `sel.2` are the skewness of the kinase's and the substrate's HPA tissue profiles (`coexpr.r`).

Every definition below was traced in `src/feature_table_generation/*.r` / `pssm.py`. Imputation follows their code: `sign-test-set.r` applies the replace_na defaults, then `merge_kin_sites.r` sets every remaining NA to 0. Factor columns get R's `as.integer(as.factor())` codes; `residue` uses learn.py's coding {T: 0, Y: 1, S: 2}.

**Status of the 44 columns: 39 regenerated (same definition and source data, with the minor version notes in the table); 3 regenerated with a method substitution (`coreg.ntera/hl/mcf7`: same Hijazi data, KSEA z-score, cross-fitted); 2 substituted (`score` becomes our learned PSSM, `coreg293` uses CPTAC BRCA instead of Mertins 2016); 0 dropped.** Label-derived columns (5): `score`, `coreg293`, `coreg.ntera`, `coreg.hl`, `coreg.mcf7`. All of them go through `ksp.crossfit`.

| # | column | SELPHI definition (source) | status here |
|---|---|---|---|
| 1 | residue | S/T/Y code | regenerated |
| 2 | score | PSSM built from PSP KS pairs with family fallback / SDR assignment (`pssm.py`). **Label-derived** | **substituted**: our F2a learned PSSM (`motif_learned.fit_pssm`, empirical-Bayes family prior, no phospho channel). Cross-fitted through `ksp.crossfit` (`SelphiPSSM`) |
| 3 | coreg293 | −log10 p of the Spearman correlation between KSEA kinase activity (PSP regulon) and the site, in Mertins 2016 BRCA. **Label-derived** | **substituted**: CPTAC BRCA (BCM v1.2, 122 tumours, our protein-residualised block) in place of Mertins 2016. Regulon = training labels. Cross-fitted (`SelphiKSEACoreg`) |
| 4 | functional_score | Ochoa 2020 MOESM5, median-imputed | regenerated (a leakage-suspect feature; see the noleak arm) |
| 5 | gtex | Spearman correlation of kinase and substrate GTEx tissue median TPM | regenerated. GTEx v8 median TPM was downloaded to `external/gtex_v8` and logged in DATA.md |
| 6 | rna.tiss | Spearman correlation of HPA tissue profiles (v23) | regenerated from our store `context/coexpr/coexpr_hpa_tissue_spearman` (HPA v25 consensus; Spearman, so the log transform does not matter) |
| 7 | rna.cell | Spearman correlation of HPA cell-line profiles | regenerated (HPA `rna_celline`, nTPM, 1,206 lines) |
| 8–9 | sel.1, sel.2 | skewness of the kinase's / substrate's HPA tissue profile | regenerated (HPA tissue consensus nTPM) |
| 10–12 | coreg.ntera, coreg.hl, coreg.mcf7 | as `coreg293`, in Hijazi 2020 (PXD015943) per cell line. **Label-derived** | regenerated from the compendium's Hijazi block (PMID 31959955, 61 × 60-min conditions per cell line, robust-SD scaling undone). Cross-fitted |
| 13–35 | is_disopred, disopred_score, exp3d_ala/acid_ddG_effect, log10_hotspot_pval_min, isHotspot, isInterface, adj_ptms_w21, netpho_max_all, netpho_max_KIN, paxdb_abundance_log10, w0_mya, w3_mya, quant_top1, quant_top5, PWM_max_mss, ACCpro, SSpro, SSpro8, sift_{min,mean,ala,acid}_score | Ochoa 2020 MOESM4 columns with SELPHI's imputation | regenerated (joined through our 15-mer mapping; 111,141 of 116,258 Ochoa rows map to U_ms sites) |
| 36–42 | isProteinDomain, isProteinKinaseDomain, isUniprotRegion, isCytoplasmic, isELMLinearMotif, isEV_ala_prediction_epistatic5, isKinaseCoreg | flags derived from Ochoa columns (`sign-test-set.r`: domain_uniprot, "inase" in the domain, region_uniprot, topology_uniprot == Cytoplasmic, elm_* \| isELMkinaseMotif, EV_ala < −5, !is.na(ptmdb_maxcoreg_kinase)) | regenerated. They come from the Ochoa table, not from UniProt/COMPARTMENTS as PRIOR-4 assumed |
| 43–44 | pssm_ser_thr, pssm_tyr | Cantley (kinase-library) score, min–max normalised per kinase (`norm-pssm.r`) | regenerated from our store `motif/cantley_unprimed_score`, min–max per kinase and residue type over U_B4. A missing matrix becomes 0 (their rule) |

**Check against SELPHI's own table.** For 4,050 positive sites that both tables share, our site columns match theirs:

- Pearson r ≥ 0.95 for 24 of the 31 site columns.
- `is_disopred` and `isKinaseCoreg` are 1/2 codes inverted. A random forest is invariant to that.
- `SSpro` and `SSpro8` have an offset of 1, which is also RF-invariant.
- Their `isUniprotRegion` treats the literal string "NA" as non-empty, an R/openxlsx artefact. We keep the intended meaning (r = −0.68 against theirs).
- They factor-code `log10_hotspot_pval_min` lexicographically, an R artefact. We keep it numeric (r = 0.78).
- `isEV_ala_prediction_epistatic5` is nearly constant (93% equal).

Leakage-suspect columns (PRIOR-4), dropped together in the **noleak** ablation arm: `functional_score` (Ochoa's model was trained on PSP regulatory sites), `PWM_max_mss` and `netpho_max_KIN` (outputs of other predictors' PSSMs), and `isKinaseCoreg` (ptmDB kinase co-regulation, derived from PSP).

## 2. Model, negatives, CV: faithful parts and deviations

| element | SELPHI2.0 | here |
|---|---|---|
| model | sklearn RandomForest, bootstrap=True | same (sklearn 1.x, max_features='sqrt' default) |
| grid | max_depth {10,20,50,70,80,100} × min_samples_split {8,10,12} × n_estimators {150,500,1000,1500}; GridSearchCV(cv=3, roc_auc) in every run | **same full grid, 72 configurations.** Run once per outer fold (on bag 0) and reused by the other 9 bags. Inner 3-fold is **GroupKFold over substrate clusters**, not their pair-level CV. Two shortcuts give exactly the same results: n_estimators is scored from prefixes of a 1,500-tree forest (sklearn seeds trees sequentially), and depths ≥ the deepest realised tree (41–50) are not refitted |
| bagging | 100 runs, each with a new random negative set; final probability = mean | **10 bags** (time box), each with a new negative draw; mean probability |
| negatives | random kinase–site pairs. Paper: 10×/20×/50× ("no differences"; 50× presented). Code: 1× when > 10,000 positives | **10×** uniform random (kinase, site) pairs over training-fold U_B4 sites, residue-class matched. **Adapted to our masks**: minus KP(k), 15-mer Hamming ≤ 2 and Par(k) (`ksp.labels.NegativeMasks`, built from training-fold label rows only), minus positives. Kinases = those with ≥ 1 training positive (their kinase list = PSP kinases) |
| positives | PSP 14,542 human pairs | our L2 training labels (`load_dev('L2')`, U_ms) on U_B4 training sites. About 5.5k per fold (PRIOR-4 deviation 3) |
| outer validation | StratifiedKFold(10) over pairs | **our 5-fold substrate-cluster CV** (`results/splits/cv_folds.tsv`) |
| label-derived features | computed from all PSP labels (the target pair is in its own PSSM / regulon) | `ksp.crossfit`: inner 5-fold over substrate clusters for training rows, all outer-training labels for eval rows. Raises LeakageError on any protein/cluster overlap |
| KSEA | evocellnet `ksea_batchKinases` (permutation, trial = 1000) | Casado z-score per condition, (mean regulon − mean all)·√m / sd (speed) |
| Spearman | `cor.test` pairwise-complete | ranks over each profile's own observed conditions, then Pearson over shared conditions (≥ 5). An approximation |
| prediction universe | 421 kinases × 238k PSP sites (non-Ochoa sites median/zero-imputed) | all 503 index kinases × **U_B4** (U_ms sites with Ochoa features), residue-class matched; about 7M pairs per fold, 35M in total |

**U_B4** is the universe where SELPHI's features are available: U_ms ∩ Ochoa 2020 = 111,141 sites (32% of U_ms; S 81,819 / T 23,256 / Y 6,066). PRIOR-4 declared this restriction in advance. SELPHI itself also scored non-Ochoa sites with imputed defaults (`cache/selphi2/site_features_all.parquet` supports this), but we did not run that arm within the time box.

Canaries (AIFM-3) were **not** run through B4. B4 is a baseline, not one of our models, and the shared crossfit driver's LeakageError check covers the label-derived columns.
