# DESIGN — Human kinase–substrate (site-level) predictor

Version: v0 draft (2026-10-26). Revisions are logged in §13. Inputs: `LITERATURE.md`, `docs/lit/*`, `docs/data_profile.md`, `docs/data_acquisition.md`.

## 1. Goal and outputs (fixed by PROMPT.md)
Predict, for every human protein kinase *k* and every human phosphoacceptor site *s* (S/T/Y), a score that *k* directly phosphorylates *s*.
- **Priority 1:** accuracy of KS prediction, harmonious across many kinases, including understudied ones.
- **Priority 2:** coverage of as many human phosphosites as possible.
- **Deliverables:**
  - per-(k,s) calibrated score and a within-kinase percentile for the prediction universe;
  - per-kinase substrate lists for kinase activity inference;
  - code and reports.
- **Runtime** is reported for every method: train time, and inference time for the full universe.

## 2. Entities and universe
- **Kinases (K):** human protein kinases from UniProt (release 2026_03, `src/ksp/uniprot.py`).
  - Atypical protein kinases are included: PIKKs (ATM, ATR, MTOR, PRKDC, SMG1), FAM20C, PDHKs, BCKDK, alpha-kinases, RIOKs, TRPM6/7 (BRD4 and other proteins with disputed kinase activity are excluded).
  - **Excluded:** lipid kinases (PI3K, PI4K, PIP5K…), non-kinase regulatory subunits (CSNK2B, PRKAB/G, CCNs), fusions (BCR-ABL, which is mapped to nothing), metabolic enzymes with disputed protein-kinase activity (PKM, PGK1), and STK19/WHR1.
  - The list and rationale are written to `results/labels/kinase_list.tsv`.
  - Non-human kinases in labels are mapped to their human ortholog by gene symbol (99.3% map).
- **Sites (U):** human S/T/Y residues on Swiss-Prot canonical sequences. Two coverage tiers:
  - **U_obs** (main universe): observed phosphosites, defined as PSP `Phosphorylation_site_dataset` human (240,698), the compendium (189,962), CPTAC and Ochoa 2020, mapped to canonical by 15-mer.
  - **U_all:** all 1.87M S/T/Y residues, scored for coverage with sequence/context features only and flagged as "not observed".
- **Site identity** = (UniProt canonical accession, position) after 15-mer-based mapping. Isoform-only sites with no canonical match are dropped and their count is reported.

## 3. Labels
- **Gold source:** PSP Kinase_Substrate_Dataset (2025-05-19).
  - Rows are kept if they have an S/T/Y residue, a kinase in K, and a substrate mappable to human canonical.
  - Mapping for non-human substrates: SITE_GRP_ID to human site, falling back to ortholog sequence alignment.
  - Metadata kept per pair: evidence (in vivo / in vitro / both), species of origin (human-native vs rodent-transferred), and autophosphorylation flag.
- **Label variants** (all built; one default chosen a priori, the others are ablations evaluated on the *same* test definition):
  - L1 = human-native only.
  - **L2 (default)** = human-native + mouse/rat transferred.
  - L3 = L2 + all other species.
  - L2-vitro = L2 restricted to pairs with in vitro evidence (in vitro only or both).
  - L2-noauto = L2 with autophosphorylation removed.
- **Test labels** are always human-native PSP pairs on test proteins (all evidence types). Results are reported stratified by evidence class (in vivo-only / in vitro-only / both). Rodent-transferred pairs are never used as test labels, so the test definition is independent of the pooling choice.
- **Independent external test** (never used for training in the main model): curated non-PSP pairs on test proteins from UniProt ECO:0000269 "by <kinase>" and the OmniPath curated non-PSP rows (SIGNOR, HPRD, NCI-PID, Reactome, Phospho.ELM).
  - Filtered to single-kinase gene names, S/T/Y, and canonical-mappable.
  - These same sources on *dev* proteins may augment training (ablation L2+ext).
- **Autophosphorylation** pairs are kept in training labels (default L2). They are excluded from per-kinase test positives in a sensitivity analysis, because protein-level features (substrate = kinase) trivially predict them.

## 4. Data splits (leakage control)
1. **Substrate grouping.**
   - mmseqs2 clusters all human canonical proteins that carry a U_obs site, at 40% identity and 50% coverage (cov-mode 0).
   - Union-find additionally merges proteins that share any identical 15-mer window centred on a site in U_obs, which captures paralog sites and duplicates.
   - The resulting groups ("substrate clusters") are the unit of splitting.
2. **Outer test split (locked):**
   - 20% of substrate clusters go to TEST, 80% to DEV.
   - Assignment is greedy: it is stratified so that every kinase with ≥5 L2 labels has ≈20% of its human-native positives in TEST where cluster structure permits.
   - Uses seed 0, is computed once, and the test list is written to `results/splits/test_proteins.tsv` and committed.
   - Test labels are only read by `src/ksp/eval_test.py --final`; every invocation is logged to RESEARCH_TREE.md. Development uses DEV only.
3. **Development CV:** 5-fold grouped CV over DEV substrate clusters, using the same stratification. All model selection and ablations use the CV out-of-fold metrics.
4. **Leave-kinase-out (LKO) tiers.**
   - Kinases are clustered by kinase-domain sequence (mmseqs2 at 50% identity on KinCore/UniProt domain boundaries) and grouped with the Manning family.
   - 5 kinase folds are assigned by family, so no family is split.
   - Tier T2 (unseen kinase family, seen substrate proteins): the model is trained on DEV labels of kinases outside the fold and scored on DEV sites for the held-out kinases. This uses out-of-kinase-fold predictions.
   - Tier T3 (unseen kinase family and unseen protein): the model is trained on DEV × kinases outside the fold and evaluated on TEST sites for held-out kinases (final only).
   - Tier T1 = standard (seen kinase, unseen protein).
   - An additional "near-zero-shot" LKO uses random kinase folds, where paralogs can be in train; it is reported for contrast only.
5. **Paralog-kinase redundancy audit:** for kinase pairs with ≥80% domain identity (e.g. MAPK1/3, CDK1/2, AKT1/2), report the overlap of their label sets. T1 results are not counted as "generalisation" for kinases whose substrates are shared with a paralog.

## 5. Negatives and the PU framing
- There are no true negatives. For kinase *k*, **"unlabelled" sites** are those in U_obs, of a residue class *k* acts on, that are not labelled to *k*.
- **Kinase residue class** is fixed a priori:
  - TK group → Y;
  - others → S/T;
  - dual-specificity kinases (MAP2K1–7, DYRK1A/1B/2/3/4, TTK, CLK1–4, TESK1/2, PKMYT1; listed in `kinase_list.tsv`) → S/T/Y; non-TK tyrosine kinases (WEE1, WEE2) → Y.
  - Labels outside a kinase's class (e.g. S/T sites for TKs, which are 69% in vivo-only and likely indirect) are kept in training as positives only if the kinase is dual-specificity. Otherwise they are dropped, and the count is reported.
- **Training negatives:** per kinase, sampled from unlabelled DEV-fold sites of the matching residue class.
  - Mixed: 50% "hard" (sites labelled to other kinases) and 50% "random observed".
  - Default ratio is 1:20 positive:negative per kinase (min 50 negatives). The ratio and mix are tuned on CV.
- **Evaluation negatives** (both pre-declared):
  - **E-all:** all unlabelled test sites of matching residue class within the evaluated split. This is realistic, but confounded by ascertainment.
  - **E-hard:** only test sites with ≥1 kinase annotation that are not labelled to *k*. This cancels the "studied site" bias and measures kinase specificity.
- **Ascertainment control:** observation-count / literature-count / number-of-conditions-measured are NOT used as features. A diagnostic "popularity-only" model reports how much of E-all performance a popularity score alone achieves. Our model must beat it on E-hard, where popularity is uninformative by construction.

## 6. Evaluation metrics (pre-registered; not changed after seeing results)
- **Kinase-centric** (per kinase *k* with ≥3 test positives):
  - AUROC and AUPRC under E-all and E-hard.
  - AUPRC is reported raw and as **AUPRC lift** = AUPRC / prevalence, to compare across kinases.
- **Primary endpoint:** macro-average over kinases of per-kinase **AUROC under E-hard** (kinase-equal weighting).
  - Justification: it is prevalence-invariant, insensitive to ascertainment bias, and weights understudied kinases equally.
- **Co-primary:** macro per-kinase **AUPRC under E-all.**
  - Justification: it is the realistic discovery setting, at the precision that matters to users.
- **Stratification:** all macro metrics are also reported by training-label bins (0 / 1–4 / 5–19 / 20–99 / ≥100), by kinase group (AGC, CAMK, CMGC, CK1, STE, TK, TKL, Atypical, Other), by evidence class, and for "priming-dependent" kinases (Cantley `pp` flag + GSK3A/B, CSNK1*).
- **Site-centric:** for test sites with ≥1 label, the rank of each true kinase among all kinases in K of matching class, summarised as Hit@1, Hit@5, Hit@10 and MRR. This is compared against a kinase-frequency prior baseline.
- **Coverage (Priority 2):** the number of U_obs sites (and of sites without any PSP annotation) with ≥1 kinase prediction at an estimated precision ≥ 0.5 and ≥ 0.8. Precision is estimated from per-kinase calibration on CV (E-all, prevalence-corrected to the full U_obs), and also at fixed per-kinase percentile thresholds (top 1%).
- **Downstream (secondary):** kinase activity inference with benchmarKIN (perturbation + CPTAC tumour benchmarks), using top-N predicted substrates per kinase vs PSP-only and Cantley-only libraries.
  - Because the compendium contains perturbation studies, this evaluation uses a model variant whose phosphoproteomics features exclude the benchmark's studies (matched by PMID/PXD).
- **Uncertainty:** 95% CIs by cluster bootstrap over substrate clusters (1,000 resamples). Model comparisons are paired, on the same bootstrap samples. Training stochasticity uses seeds [0,1,2,3,4] (mean ± sd).
- **Metric code** is verified on toy data with known answers before use (§11).

## 7. Features (all computed inside each fold from training labels only, if label-derived)
Each feature family is a module in `src/ksp/features/` with a provenance flag (label-derived? source? leakage risk), recorded in `docs/feature_provenance.md`.

**F0 Site-intrinsic** (kinase-independent; these enter only as interactions in pair models):
- residue type; one-hot ±7 window;
- AlphaFold pLDDT, RSA, disorder, secondary structure (structuremap-like, computed ourselves from AFDB PDB/CIF);
- ESM-2 650M per-residue embedding at the site (later 3B or ESM-C 600M).

**F1 Motif (published):** Cantley PSSM score and percentile for (k,s), for kinases in the atlases (S/T 303, Y 93).
- **Priming-aware:** flanking residues known to be phosphorylated (in U_obs) are scored with the pS/pT/pY matrix rows. Both "unprimed" and "primed" scores are kept.
- For kinases not in the atlases, the score of the nearest atlas kinase (by domain identity) is used, with a distance feature.

**F2 Motif (learned from labels, fold-aware):**
- (a) per-kinase PSSM / log-odds from training positives, with a Dirichlet prior from the kinase's family;
- (b) a neural sequence model (§8 R3) score.

**F3 PPI (leak-filtered):** only from systematic screens: BioPlex 3.0, HuRI (no Lit-BM), OpenCell, Buljan 2020 kinase AP-MS (to acquire), Human Cell Map.
- direct-interaction flag;
- random-walk-with-restart proximity of kinase to substrate protein;
- **fold-aware:** proximity of substrate protein to the kinase's *training* substrates.
- **Excluded:** STRING (all channels), BioGRID biochemical-activity/PTM, IntAct enzymatic, curated pathway DBs. This is justified by the leakage audit (`results/leakage/`).

**F4 Localisation:** cosine similarity of kinase vs substrate compartment vectors from HPA IF, OpenCell, COMPARTMENTS (knowledge+experiments only) and Human Cell Map, plus missingness flags.

**F5 Co-expression / co-dependency:** kinase–substrate correlation across:
- GTEx/HPA RNA tissues;
- CPTAC protein abundance (per cohort, averaged);
- DepMap Chronos co-essentiality.

**F6 Phosphoproteomics, kinase-specific** (all study-aware):
- (a) **Co-regulation with training substrates (fold-aware, leave-site-out):** for site *s* and kinase *k*, the correlation of *s*'s profile with the profile of each training substrate of *k*, aggregated (mean of top-m and median Fisher-z). It is computed per data block (compendium per-study blocks; each CPTAC cohort after regressing out protein abundance), then combined with study-size-aware weighting.
- (b) **Inhibitor response:** mean / min log2FC of *s* under conditions where *k* is inhibited (Hijazi matrix; decryptM drugs mapped to targets via Kinobeads/DTC if acquired), versus conditions where *k* is not inhibited (a target-specificity contrast).
- (c) **Correlation with kinase state in CPTAC:** correlation of *s* with *k*'s own protein abundance and with *k*'s activation-loop phosphosites (not label-derived).
- (d) **Kinase-activity profile correlation:** per condition, estimate *k* activity from training substrates (mean z-score, leave-site-out), then correlate with *s*.
- **Missing data:** pairwise-complete with minimum n (≥8 conditions, ≥2 studies), otherwise NA with a flag. Missingness itself is not given as a numeric "count" feature.

**F7 Kinase descriptors (for sharing / zero-shot):** Manning group/family, ESM-2 kinase-domain embedding (PCA-reduced), Cantley PSSM vector (PCA), and log number of training labels (used only for calibration, not as a ranking feature; ablated).

**F8 Auxiliary high-throughput evidence:** Sugiyama 2019 in vitro lysate hits (kinase-specific); PhosphoNetworks protein-level array hits. These are features only, never labels. They are ablated, since they are in vitro data that may overlap PSP in vitro labels.

**Excluded as features:**
- Ochoa 2020 functional score, PSP Regulatory_sites, LT_LIT/MS counts (trained on or proxies of PSP annotation);
- CoPheeMap/CoPheeKSA networks (trained on PSP);
- any third-party predictor output (NetworKIN, GPS, KolossuS, Phosformer), which are used only as baselines.

## 8. Models / training regimes (compared on CV, E-hard + E-all)
- **R0 Baselines:**
  - B1 Cantley percentile (as published; priming-naive and priming-aware);
  - B2 per-kinase PSSM learned from training labels;
  - B3 Phosformer re-implemented (RoBERTa-like kinase-domain + peptide encoder; we re-run their recipe and pretraining on kinase domains, or initialise from ESM-2, documented) and trained on our folds; plus the released Phosformer weights as-published (inference only, with the caveat that its training data overlap our test);
  - B4 SELPHI2.0 retrained on our folds with its feature set regenerated fold-aware (Apache-2.0 code from GitLab EBI, adapted to our features where its raw inputs are unavailable; deviations documented);
  - B5 NetworKIN predictions (benchmarKIN Zenodo) as-published, with coverage reported;
  - B6 KolossuS as-published if inference fits (ESM-2 15B across 2 GPUs), else skipped.
- **R1 Per-kinase models:** LightGBM/logistic regression per kinase with ≥20 training labels; fallback to the family-level model.
- **R2 Global pair model:** one LightGBM over (k,s) pairs with F1–F8 + kinase descriptors F7; binary loss vs lambdarank grouped by kinase; sample weights equalising kinases (∝ 1/n_k^α, α ∈ {0, 0.5, 1} tuned).
- **R3 Neural sequence model:** a site encoder (ESM-2 embedding of the ±15 window, or a fine-tuned small ESM) and a kinase encoder (ESM-2 kinase-domain embedding + Cantley vector), combined by a bilinear / MLP head and trained multi-task across all kinases with hard negatives. Its score feeds R2 (stacked; out-of-fold within DEV to avoid leakage).
- **R4 Hierarchical:** global R2 + per-kinase residual models (R1) blended with weight depending on n_k.
- **Calibration:** per-kinase isotonic / Platt on out-of-fold CV predictions, with prevalence correction; percentile ranks within kinase over U_obs.
- **Final model** = the best CV configuration retrained on all DEV labels → evaluated once on TEST. A deployment model is then retrained on DEV+TEST for the released predictions, and is clearly separated from the evaluated one.

## 9. Unorthodox kinases
- **Priming-dependent** (GSK3A/B +4 pS; CK1 −3 pS; CK2 acidic/primed; PLK1 polo-box docking on pS/pT-P): priming-aware F1 scoring, plus explicit features for a phosphorylatable residue at the priming offset and whether it is an observed phosphosite / predicted target of a common priming kinase.
- **Atypical** (PIKK S/T-Q, FAM20C S-x-E secreted pathway, alpha-kinases, PDHK mitochondrial): handled by the global model with localisation and descriptor features; reported as a separate group.
- **Docking-motif kinases** (MAPKs D-motif; CDKs RxL/Cks1; AGC PIF): substrate-protein-level docking motif features (ELM class regexes, distance to site).
- Each of these groups is reported separately in all stratified results.

## 10. Compute plan (to be measured early; LESSONS records actual times)
- Features: CPU (≤32 cores); co-regulation is the heaviest, vectorised on GPU if needed.
- ESM-2 650M embeddings of all U_obs windows: <1 GPU-h. ESM-2 3B: optional.
- LightGBM global: about 1–5M training rows, minutes per fold on CPU.
- R3 neural: ≤1 GPU-h per fold on one 4090.
- Inference over U_obs × K ≈ 300k × 430 ≈ 130M pairs: chunked, feasible in about 1 h CPU.
- **GPUs:** 0–1 for neural models, 2 for embeddings, 3 for baselines (Phosformer/KolossuS).

## 11. Verification and controls (AI-agent failure-mode mitigations)
- **Toy-data tests** (`tests/`):
  - metric functions on hand-computed examples;
  - split code (assert no protein/cluster/identical-window overlap between TEST and DEV, or between CV folds);
  - fold-aware features (planted-signal toy where leakage would give a perfect score: the assert fails if the feature uses held-out labels).
- **Negative controls:**
  - label-permutation within kinase → macro AUROC ≈ 0.5;
  - random-feature model;
  - a popularity-only model (diagnostic).
- **Positive controls:** a planted-motif synthetic kinase is recovered by R2/R3; known canonical motifs are recovered by the learned PSSMs (e.g. PKA R-R-x-S, CDK S/T-P, CK2 S-x-x-E/D).
- **Leakage alarms:**
  - the curated-vs-screen PPI ablation;
  - comparing the per-kinase CV metric against the test metric (a large CV>test drop is investigated);
  - a feature-importance audit for implausible dominance.
- **Independent review:** every key result is reviewed by an agent that did not produce it, before it enters RESEARCH_TREE as [x]/[!].
- **Locked test:** only one final evaluation run per final model; any re-run is logged with a reason.
- **Reproducibility:** fixed seeds [0,1,2,3,4]; every result table is produced by a committed script from cached inputs; DATA.md has checksums; `make`-style pipeline `scripts/run_all.sh`.

## 12. Experiment plan (priority order)
1. Labels L1–L3, kinase list, site universe, splits + tests.
2. Baselines B1, B2 on CV, plus the evaluation harness (fast; establishes the bar).
3. F0–F1 and F2a features + R2 → first real model.
4. Add F3–F6 one family at a time (ablation curve).
5. R3 neural → stacking into R2; B3 Phosformer retrain in parallel on a separate GPU.
6. B4 SELPHI2.0 retrain.
7. Training-regime comparison: R1 vs R2 vs R4, kinase weighting α, negative mix.
8. Label variants (L1/L2/L3/vitro/+ext) and the unorthodox-kinase analysis.
9. LKO tiers T2 (CV), plus T3 at final.
10. Freeze → final TEST evaluation → downstream benchmarKIN → coverage predictions → REPORT.md.

## 13. Revision log
- v0 (2026-10-26): initial draft by PI.
