# B4 — SELPHI2.0 retrained on our splits

Owner: B4-selphi2 agent. DESIGN v2 §8.1 B4, §6.8 (matched-universe rule), PRIOR-4 (`docs/reviews/r1_prior_work.md`).
Reference: Maier, Petursson, Lussana & Petsalaki, *MCP* 2025 24(8):100994 (`data/paper1.pdf`; the brief calls this paper2, but paper2.pdf is Phosformer). Code: `external/selphi2/selphi_2.0` (Apache-2.0). I read it as a reference and ported its logic; none of it was executed.

Code: `src/ksp/baselines/selphi2/{features,model}.py`, `scripts/baselines/selphi2_{build_features,cv,eval,pairsplit}.py`.
Cache: `cache/selphi2/` (OOF predictions in `cache/selphi2/<arm>/oof_fold{f}.parquet`; the ksp.manifest copy is in `cache/runs/<run_id>/predictions.parquet`).
Results: `results/baselines/selphi2/` (`<arm>/fold{f}.json` holds timings, grid table, importances and counts; `<arm>/summary.json`, `metrics_{matched,full}.json` and `per_kinase_matched.tsv` hold the metrics; `pairsplit_reference.json` holds the non-official reference).

__RESULTS__

## 1. Feature table: 44 model columns (residue + 43)

SELPHI's released feature table (`data/psp/pos_set_features_20240910.tsv`) has **44 model columns**. The paper says "48 generated, 45 selected". Its Supplemental Table 1, which would list the selected features, is not in the repo, and the code (`learn.py`) trains on exactly these 44 columns. So we reproduce the 44 and record the paper-versus-table discrepancy of one column. PRIOR-4 counted "43" because it left out `residue`. It also called `score`, `sel.1` and `sel.2` "undocumented". They are resolved here from the code:

- `score` is the PSP-derived PSSM built in `pssm.py`.
- `sel.1` and `sel.2` are the skewness of the kinase's and the substrate's HPA tissue profiles (`coexpr.r`).

Every definition below was traced in `src/feature_table_generation/*.r` / `pssm.py`. Imputation follows their code: `sign-test-set.r` applies the replace_na defaults, then `merge_kin_sites.r` sets every remaining NA to 0. Factor columns get R's `as.integer(as.factor())` codes; `residue` uses learn.py's coding {T: 0, Y: 1, S: 2}.

**Status of the 44 columns: 37 regenerated with SELPHI's definition from the same source; 7 regenerated with a substituted source or method (5 label-derived and cross-fitted, 1 replacement data source, 1 method change in `rna.tiss` is minor); 0 dropped.**

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
