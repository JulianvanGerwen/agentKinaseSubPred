# 01 — Prior methods and benchmarks for kinase–substrate (site-level) prediction

Author: literature agent (methods/benchmarks). Date: 2026-10-26.
Scope: (1) deep read of the two provided papers; (2) Cantley kinome atlases as a baseline; (3) survey of other methods, emphasis 2020–2026; (4) synthesis + concrete benchmark/baseline recommendations.
All claims cite DOIs/URLs. Numbers quoted from papers are *as reported by the authors* — most are not comparable across papers because the test sets differ. Criticism is explicit and separated from reporting.

---

## 0. TL;DR for our design

- No published benchmark in this field is leakage-controlled in the way our brief requires (group by substrate protein / homologous proteins, split on sites, leave-kinase-out). Every method below splits on kinase–substrate **pairs** at random, or on peptides, or evaluates on a self-generated label set. We therefore cannot reuse any published number; we must retrain.
- The two provided "benchmarks" are of very different kinds: Phosformer (sequence-only PLM, retrainable in hours on our hardware) and SELPHI2.0 (a 45-feature random forest over co-regulation/co-expression/PSSM features, retrainable in principle but feature regeneration is the expensive part and some inputs are themselves PSP-derived → leakage risk).
- The Cantley atlases (Johnson 2023 Ser/Thr, Yaron-Barir 2024 Tyr) are the only *non-PSP-trained* baseline and are the honest reference point. They need no retraining; the only PSP contact is in their *evaluation* (and, for the Ser/Thr percentile reference distribution, the Ochoa 2020 phosphoproteome — not PSP K–S labels).

---

## 1. Paper 1 (provided): SELPHI2.0 — Maier, Petursson, Lussana, Petsalaki

Maier BD*, Petursson B*, Lussana A*, Petsalaki E. "Data-Driven Extraction of Human Kinase-Substrate Relationships From Omics Datasets." *Mol Cell Proteomics* 24(8):100994, 2025. doi:[10.1016/j.mcpro.2025.100994](https://doi.org/10.1016/j.mcpro.2025.100994). Local PDF: `/workspace/data/paper1.pdf`.

### 1.1 Task and scope
Site-level probability for every (kinase, phosphosite) pair: **421 kinases × 238,374 phosphosites** (199,262 S/T, 39,112 Y) on 17,469 proteins = 73,207,860 predictions. Sites = PhosphoSitePlus sites that also have an Ochoa 2020 functional score. Dual-specificity handling: 33 kinases annotated for both S/T and Y in PSP (or dual in NCBI gene summaries) are scored on all 238,374 sites; others are restricted by residue type.

### 1.2 Inputs / features (45 of 48 kept)
Direct evidence from the code (`data/psp/pos_set_features_20240910.tsv` header at <https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0>) gives the actual feature list:

`prot1 (kinase), prot2 (substrate), psites, residue, score, coreg293, functional_score, gtex, rna.tiss, rna.cell, sel.1, sel.2, coreg.ntera, coreg.hl, coreg.mcf7, is_disopred, disopred_score, exp3d_ala_ddG_effect, exp3d_acid_ddG_effect, log10_hotspot_pval_min, isHotspot, isInterface, adj_ptms_w21, netpho_max_all, netpho_max_KIN, paxdb_abundance_log10, w0_mya, w3_mya, quant_top1, quant_top5, PWM_max_mss, ACCpro, SSpro, SSpro8, sift_min_score, sift_mean_score, sift_ala_score, sift_acid_score, isProteinDomain, isProteinKinaseDomain, isUniprotRegion, isCytoplasmic, isELMLinearMotif, isEV_ala_prediction_epistatic5, isKinaseCoreg, pssm_ser_thr, pssm_tyr`

Grouped:
- **Kinase-specificity (pair-level, only 2–3 features)**: `pssm_ser_thr`, `pssm_tyr` — min–max-normalised-per-kinase PSSM scores from the Cantley Ser/Thr and Tyr atlases, computed by the companion Nextflow pipeline <https://github.com/alussana/selphi2.0-pssms> (`bin/compute_pssm_scores.py` implements Johnson Supplementary Note 2 scoring: multiply normalised PSSM values across positions −5..+4, scale by random-peptide probability). Earlier Invergo-derived PSSMs also feature in the README lineage. PSSMs are downloaded from the Springer supplementary XLSX directly (URLs in `nextflow.config`, see §2.4).
- **Co-regulation / co-phosphorylation (pair-level)**: `coreg293` (Mertins 2016 breast-cancer proteogenomics), `coreg.ntera/hl/mcf7` — Spearman correlation between *estimated kinase activity* (KSEA over PSP substrate sets, `ksea` R package) and target-site phosphorylation across 63 inhibitor conditions in NTERA2/HL60/MCF7 from Hijazi 2020; `isKinaseCoreg`.
- **Co-expression (pair-level)**: `gtex`, `rna.tiss`, `rna.cell` (GTEx v8, HPA v23.0 tissue + single-cell).
- **Site-intrinsic (substrate-only, from Ochoa 2020 functional landscape)**: functional score, disorder (DISOPRED3), structural ΔΔG, hotspot, interface, neighbouring PTM density, NetPhos scores, PaxDb abundance, conservation (`w0_mya`, `w3_mya`), MS quantification frequency, secondary structure/accessibility (SSpro/ACCpro), SIFT, domain/ELM/localisation flags.

**Critique**: only ~5 of 45 features are kinase-specific. The rest are substrate/site properties shared by all kinases for that site, so the model largely learns "is this site phosphorylatable/functional/well-measured" plus a weak kinase term. That inflates pooled AUC enormously and is exactly the failure our benchmark must expose (per-site ranking across kinases, not pooled pair classification).

### 1.3 Training data, species, in vivo/in vitro, negatives
- Positives: **14,542 kinase–phosphosite pairs from PhosphoSitePlus** (downloaded 2024-09-03). Human-only in the activity scripts (`KIN_ORGANISM == "human" & SUB_ORGANISM == "human"`). No in vivo/in vitro stratification; no quality filtering beyond PSP itself.
- Negatives: **random (kinase, site) pairs** at 10×, 20×, 50× the positive set ("no differences between results"; the 50× version is reported). No hard negatives, no decoy matching by residue/motif/abundance.
- Feature selection: 100 train/test datasets at 10× negatives, RF + RFE-CV (10-fold), majority vote >50% → 45 features.
- Final model: **scikit-learn RandomForestClassifier**, grid search over `max_depth ∈ {10,20,50,70,80,100}`, `min_samples_split ∈ {8,10,12}`, `n_estimators ∈ {150,500,1000,1500}`, bootstrap=True; trained **100×** with resampled negatives; final probability = mean over the 100 models.

### 1.4 Splits and leakage — serious problems
- Splits are **stratified 10-fold CV over kinase–substrate pairs** (`StratifiedKFold`, `src/model_training/cross-val.py`). No grouping by substrate protein, no sequence-identity control, no leave-kinase-out.
- **Leak 1 (substrate protein)**: many features are substrate-protein-level (co-expression, abundance, disorder, localisation, functional score). A protein with several annotated sites appears in both folds → the model memorises "substrates of well-studied proteins".
- **Leak 2 (kinase activity features)**: `coreg.*` features are correlations with KSEA kinase activities computed *from PSP substrate sets*. The scripts (`calc-cuttilas-act.r`, `cutillas_act-cor.r`) build the regulon from the same PSP table used for the positive labels and, as far as the published code shows, do **not** remove the target site from its kinase's regulon before correlating. For a kinase with few substrates this is close to correlating a site with itself. (Contrast CoPheeKSA, which explicitly removes the site — §3.8.) This alone can explain much of the 0.97 CV AUC.
- **Leak 3 (negatives)**: random pairs drawn from the same 73M grid include true-but-unannotated positives; with 50:1 imbalance and pooled AUC this mostly inflates the metric.
- The authors are aware of the pooled-metric problem only indirectly: they note residual study bias (more predictions for well-cited kinases) and that only 5 of their 33 "dual-specificity" kinases look truly dual, blaming PSP annotation errors (e.g. serine sites attributed to EGFR from indirect evidence) — a useful warning about PSP label quality for us.

### 1.5 Evaluation and reported numbers
- Internal CV: **AUC-ROC 0.97** (mean of 100 runs), PR curve also shown. 721,707 edges with p>0.5 (1.7% of them in PSP; >90% of known PSP interactions get p>0.5).
- External: **Sugiyama 2019** (in vitro LC-MS/MS kinome-wide substrate screen, *Sci Rep* 9:10503, doi:[10.1038/s41598-019-46385-4](https://doi.org/10.1038/s41598-019-46385-4)) and **Hijazi 2020** (inhibitor phosphoproteomics-derived KS edges, *Nat Biotechnol* 38:493, doi:[10.1038/s41587-019-0391-9](https://doi.org/10.1038/s41587-019-0391-9)), after removing pairs already in PSP. Pairwise comparisons were restricted to the intersection of predictions each method can make, with 100 validation sets at 10:1 negatives:
  | Comparator | SELPHI2.0 AUC | comparator AUC |
  |---|---|---|
  | GPS 6.0 | 0.83 | 0.81 (349 kinases) |
  | KinomeXplorer | 0.79/0.81 | 0.61 |
  | NetPhos 3.1 | 0.75 | 0.69 (17 kinases only) |
  | GPS 5.0 | 0.77 | 0.56 |
  | PhosphoPICK | 0.76 | 0.47 |
  | LinkPhinder | 0.77 | 0.49 |
  Per-family on the external set: Atypical 0.81, CMGC 0.79, CAMK 0.76, AGC 0.72, Other 0.72, TK 0.69, STE 0.68, CK1 0.65, TKL 0.65.
- Downstream: PCSF fitting of the network to 436 Ochoa 2016 phosphoproteomics conditions raises F1 for recovering known edges from 0.031 → 0.19 (precision 0.016 → 0.22).

**Critique of the external evaluation**: the drop from 0.97 (internal) to 0.75–0.83 (external) is the honest signal and confirms the internal number is inflated. But (a) the "negatives" in the external set are all non-reported pairs, so AUC is again dominated by easy negatives; (b) Sugiyama is an **in vitro** screen, so a model scoring high there is being judged on intrinsic specificity, not cellular targeting; (c) comparators were run via web servers/standalone with "no threshold", and several (NetPhos 3.1 with 17 kinases) are so restricted that the comparison is near-meaningless; (d) Hijazi edges were also used to *build features* (`coreg.hl/mcf7/ntera` come from the Hijazi inhibitor data) — so the Hijazi-based validation is partially circular even though the labels differ from the features.

### 1.6 Code, weights, data, compute
- Code: <https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0> (**Apache-2.0**, 1 commit, created 2024-09-20; R + Python; LSF `bsub` shell wrappers under `modules/shell_scripts/`, container via `env/Dockerfile`). PSSM feature pipeline: <https://github.com/alussana/selphi2.0-pssms> (Nextflow + Docker/Apptainer, **no license file** — treat as all-rights-reserved for redistribution, fine to read/reimplement). Radial trees: <https://github.com/alussana/radialtree>.
- Data: Zenodo doi:[10.5281/zenodo.13821320](https://doi.org/10.5281/zenodo.13821320), **CC BY 4.0**, 27.1 GB total: `feature_table.tsv` (10.7 GB), `prediction_matrix.csv` (8.3 GB), `k_p_ser_thr_pssm_scores.tsv` (4.1 GB), `k_p_tyr_pssm_scores.tsv` (3.9 GB), `prediction_matrix_high_conf_curated.csv` (87 MB), `positive-set-20240910.tsv` (224 KB), citation matrices. **Do not download the big ones** — the high-confidence matrix + positive set + the repo's `data/psp/pos_set_features_20240910.tsv` are enough to inspect the schema; full feature table only if we decide to reuse their features verbatim.
- Web server: <https://selphi2.com> (R Shiny; Docker image in the GitLab repo).
- Compute: RF on ~14.5k positives × 50 negatives × 45 features, ×100 runs with grid search — CPU-only, embarrassingly parallel; on our 64 cores this is hours, not days. The expensive part is **feature regeneration** (PSSM scoring of 421×238k pairs; KSEA activity correlations across Hijazi/CPTAC).

### 1.7 What a fair retrain requires
1. Regenerate the 45 features ourselves, on **our** site universe. Reuse `alussana/selphi2.0-pssms` logic for the PSSM features (it is clean and documented) and Ochoa 2020 supplementary tables for the site features.
2. **Fix the activity-feature leak**: when computing `coreg.*` for pair (K, s), exclude s (and ideally all sites of the same substrate protein) from K's KSEA regulon, and compute the regulon only from *training-fold* labels. This is the single most important change.
3. Replace stratified-pair CV with our site/protein-grouped splits and leave-kinase-out folds; keep their negative ratio as one arm but add motif/abundance-matched hard negatives.
4. Keep their model class (RF, 100 bagged runs) as given — this is a feature-engineering baseline, not an architecture baseline.
5. Report per-kinase and per-site-ranking metrics, not pooled AUC.
Feasible verdict: **yes, retrainable**, but it is a feature-pipeline project (~1–2 agent-days) rather than a model-training one.

---

## 2. Paper 2 (provided): Phosformer — Zhou, Yeung et al. (Kannan lab)

Zhou Z*, Yeung W*, Gravel N, Salcedo M, Soleymani S, Li S, Kannan N. "Phosformer: an explainable transformer model for protein kinase-specific phosphorylation predictions." *Bioinformatics* 39(2):btad046, 2023. doi:[10.1093/bioinformatics/btad046](https://doi.org/10.1093/bioinformatics/btad046). Local PDF: `/workspace/data/paper2.pdf`.

### 2.1 Inputs / model
- Inputs: **unaligned kinase domain sequence** + **11-mer substrate peptide** (centre S/T/Y; `<pad>` if within 5 residues of a terminus). No other features at all — no PPI, no family labels, no structure.
- Architecture: RoBERTa-style encoder, **6 attention layers, 12 heads, 768 dims**; tokenised `<s>peptide</s><s>kinase</s>`; the embedding at the central peptide residue feeds a feed-forward + binary head (softmax). Framed as context-based question answering (peptide = context, kinase = question).
- Pre-training: MLM (15% masking) on 295,320 kinase domain sequences from NCBI nr (18,832 organisms; 86.9 M aa) + 5,338 full-length substrate proteins from UniProt (29 organisms; 3.85 M aa). 50 epochs, **8× NVIDIA A5000**, batch 48/device, AdamW lr 5e-5.
- Fine-tuning: focal loss (γ=2) to handle 1:16 positive:negative imbalance, AdamW lr 2e-5, early stopping patience 5.

### 2.2 Training data: sources, species, negatives
- Positives: Phospho.ELM + PhosphoNetworks + PhosphoSitePlus, cross-checked against UniProt → **24,534 unique kinase–substrate pairs, 800 kinases, 13 organisms**. The HF `curated/phosphosites_11mer_kinase_specific.tsv` I inspected has 25,729 rows / 800 kinase accessions with per-source flags (`source|PhosphoSitePlus`, `|PhosphoELM`, `|PhosphoNetworks`). No in vivo/in vitro split.
- Negatives, two tiers (their main conceptual contribution):
  - **Hard negatives**: peptides with evidence of phosphorylation, paired with a kinase not reported to phosphorylate them.
  - **Easy negatives**: 11-mers with central S/T/Y and no phosphorylation evidence, from the same substrate proteins.
  - Training uses n easy + m hard per positive (best: total 1:16); val/test use exactly **1 hard negative per positive**; a separate test set of **8,504 easy negatives** is used only for FPR.
- Augmentation: upsample under-represented kinases (kinases >50 positives → 1,118; <50 positives → 50), random ±5-residue shifts of kinase domain boundaries, 5% kinase token masking, multilevel negative sampling.

### 2.3 Split and leakage — measured, not assumed
Stated split: positives 70:15:15 **at random over pairs**, with kinases having <50 positives confined to training; **106 kinases** in val/test.

I downloaded the released splits (<https://huggingface.co/datasets/waylandy/phosformer_curated>, files `splits/test_pos_neg.csv`, `splits/val_pos_neg.csv`, `curated/phosphosites_11mer_kinase_specific.tsv`) and measured the overlap directly (script kept at `/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/phosformer_leak.py`):

- test set: 5,318 rows, 2,659 positives, 106 kinase accessions; test positives are 2,368 human, 176 mouse, 85 rat, 30 bovine → **cross-species orthologue pairs sit in the same split**.
- **51.1%** of test-positive peptides also occur as a *positive* 11-mer in the training portion (paired with some other kinase).
- **29.9%** of test positives have the exact (kinase **family**, peptide) combination present as a training positive.
- **2.7%** of test positives have the identical (kinase domain sequence, peptide) pair in training — i.e. literal duplicates via orthologous kinases with identical domain sequence.
- **15.7%** of test hard-negative peptides are training positives (that is by construction, and is fine).

So: no substrate-protein grouping, no sequence-identity clustering of peptides or kinases, no species separation. AUROC 0.885 is therefore an upper bound on in-distribution performance and says nothing about novel substrates or novel kinases.

### 2.4 Metrics and reported numbers
Metrics: AUC-ROC and AUC-PRC on the pos + hard-negative test set (balanced 1:1); FPR on the easy-negative-only test set. Headline: **AUC-ROC 0.885** pooled; **mean 0.860** across the 106 test kinases; "majority of kinases above 0.800". Group-level selections from their Table 1 (Phosformer column): TK 0.931/0.910/FPR 0.049; AGC 0.898/0.900/0.018; CMGC 0.898/0.882/0.006; CAMK 0.830/0.813/0.009; CK1 0.808/0.822/0.000; TKL 0.605/0.612/0.011. Family level: MAPK 0.924, Src 0.927, CDK 0.907, PKC 0.884, CK2 0.881, "other families (38)" 0.857. Comparators (MusiteDeep, DeepPhos, PhosIDN, PhosIDNSeq, EMBER) were run with their own extra required inputs; several cannot predict for most groups.

**Critique**: a 1:1 positive:hard-negative test set is not the deployment regime (one site vs ~300 candidate kinases, prevalence ≪1%); AUC-PRC on a balanced set is almost meaningless. The FPR-on-easy-negatives test is a nice idea but is a *general* phosphosite task, not kinase assignment. Per-kinase averaging is good practice and worth copying. TKL at 0.605 shows the 11-mer-only model fails where specificity is docking/proximity-driven.

### 2.5 Code, weights, compute
- Code: <https://github.com/esbgkannan/phosformer> (fork of `waylandy/phosformer`). **LICENSE = CC BY-NC-ND 4.0** → non-commercial, no derivatives. For an academic benchmark this is usable, but *retraining* is arguably a derivative work; we should note this and, if needed, reimplement the architecture from the paper (it is a 6-layer RoBERTa with a trivial head — a half-day reimplementation) rather than redistribute modified weights.
- Repo contents: `Phosformer/{modules.py, roberta_modified.py}`, `environment.yaml`, `example.ipynb`, `data/reference_human_kinases.csv`, `data/example_input_peptide.csv`. **Inference only — no training script.**
- Weights: HF `waylandy/phosformer` (`RobertaForSequenceClassification`, loads via `from_pretrained`).
- Data: HF `waylandy/phosformer_curated` (~7.9 GB total; `splits/train.csv` alone is 715 MB — the big file is just the expanded negative pairs; the small split files and the curated TSVs are a few MB each and are what we need).
- Compute: reported **220 predictions/s on a single RTX 2080** for inference. Pre-training was 8×A5000 × 50 epochs, but we do **not** need to redo pre-training — the released MLM checkpoint (or plain ESM-2) can be the starting point. Fine-tuning ~25k positives × 16 negatives ≈ 400k examples on a 6-layer/768 model: **well under an hour per epoch on one RTX 4090**; a full retrain with 5 seeds is a few GPU-hours.

### 2.6 What a fair retrain requires
1. Rebuild the pair table from **our** PSP-derived labels (our species policy, our in vivo/in vitro policy), not their 13-organism union.
2. Map each kinase to its **kinase domain sequence** (Pfam PK/PK_Tyr boundaries, or reuse their `reference_human_kinases.csv` domain definitions for human).
3. Regenerate easy/hard negatives *within training folds only*, using our site universe; keep their n:m knob as a hyperparameter (they chose 1:16 total).
4. Keep the augmentations (resample, ±5 boundary shift, 5% mask) — they are cheap and reported to help.
5. Train on our site/protein-grouped splits + leave-kinase-out; evaluate with per-kinase and per-site-ranking metrics at realistic prevalence.
6. Also run the **released checkpoint zero-shot** on our test set as "Phosformer-as-published", and report both numbers (published checkpoint has seen much of any PSP-derived test set, so this number is an optimistic ceiling and must be labelled as contaminated).
Feasible verdict: **yes, cleanly retrainable**, cheapest of all baselines. Highest priority.

### 2.7 Successor: Phosformer-ST (important, and its benchmark is circular)
Zhou Z, Yeung W, Soleymani S, Gravel N, Salcedo M, Li S, Kannan N. "Using explainable machine learning to uncover the kinase–substrate interaction landscape." *Bioinformatics* 40(2):btae033, 2024. doi:[10.1093/bioinformatics/btae033](https://doi.org/10.1093/bioinformatics/btae033).
- **Labels are not literature annotations**: positives/negatives come from the **Johnson 2023 atlas percentile scores** mapped onto 89,784 S/T phosphosites: percentile >90 = positive, <50 = hard negative; "easy negatives" are 884,203 S/T sites with no phosphorylation evidence, randomly paired. 300 kinases (BCKDK, PDK1, PDK4 dropped), 86,043 peptides.
- **Leakage control is the best in the field**: CD-HIT at 70% identity over kinases (187 clusters from 300); greedy clustering of the central 11 residues of peptides at Hamming distance ≤5 (39,077 Ser clusters, 34,470 Thr clusters); the representative dataset keeps one pair per (kinase cluster, peptide cluster) combination → 1,010,423 positives, 3,080,095 hard negatives, 8,681,574 easy negatives; 60:20:20 split with 1:1:1 pos:hard:easy in val/test.
- Model: ESM-2 fine-tuned, 15-mer peptide + unaligned kinase domain, cross-attention block from kinase→peptide residues, multitask (classification + kinase MLM on the 295,320-sequence corpus). 3 epochs on **8× NVIDIA A6000**, effective batch 64/device.
- Reported: **AUROC 0.992 / AUPRC 0.992 / FPR 0.019**, Hit@1 group 0.683, family 0.822. Ablation with **~10% of kinases entirely held out**: random-init 0.833, ESM-pretrained-no-multitask 0.863, full model **0.886** AUROC on unseen kinases.
- **Critical reading**: the 0.992 is "can the model reproduce the Kinase Library PSSM percentile?" — it is a *distillation* metric, not a biology metric. It is not comparable to Phosformer's 0.885 on PSP labels (their Table 1 lists Phosformer at 0.581 on the Phosformer-ST test set, which simply reflects a different label definition). Predictions are explicitly in vitro. The **unseen-kinase ablation (0.886)** is the genuinely informative number and the closest published analogue to our leave-kinase-out benchmark.
- Code: <https://github.com/gravelCompBio/Phosformer-ST> (**CC BY-NC-ND 4.0**; inference notebook + architecture files; weights via HF/Zenodo link file). Web: <https://esbg.bmb.uga.edu/phosformer/>.
- For us: Phosformer-ST is best used as (a) an architecture template (cross-attention + multitask + cluster-aware splitting), and (b) a *pretraining* source — distilling the atlas is a reasonable way to inject specificity priors before fine-tuning on literature labels. Its headline number should never be quoted as a KS-prediction benchmark.

---

## 3. Cantley lab kinome atlases (Kinase Library) — the no-retraining baseline

### 3.1 Johnson et al. 2023 (Ser/Thr)
Johnson JL, Yaron TM, Huntsman EM, et al. "An atlas of substrate specificities for the human serine/threonine kinome." *Nature* 613:759–766. doi:[10.1038/s41586-022-05575-3](https://doi.org/10.1038/s41586-022-05575-3). PMC9876800.
- **Data**: positional scanning peptide arrays (PSPA), combinatorial library substituting each of **22 residues (20 aa + pThr + pTyr)** at **9 positions** (−5..−1, +1..+4) around a central 1:1 Ser/Thr acceptor; X positions randomised over all aa except Ser/Thr/Cys. **303 Ser/Thr kinases** profiled with purified recombinant kinase. Purely **in vitro, peptide-level, intrinsic specificity** — no cellular context.
- **Matrix normalisation**: densitometry column-normalised at each position by the sum of the 17 randomised aa (excluding S, T, C; 16 for PDHK1/PDHK4 which also exclude Tyr); Cys row scaled so its median is 1/17; S and T values at each peripheral position set to the median of that position. Central position: the Ser-vs-Thr acceptor favourability S₀/T₀ is computed from the summed Ser and Thr densitometry rows corrected for the different central (1:1) vs peripheral composition, then normalised to the larger of the two (Supplementary Note 1). `alussana/selphi2.0-pssms/bin/get_pssm_ser_thr_kinome_2023_suppl_table_2.py` implements exactly this: `S_ctrl = 0.75*S_S − 0.25*S_T`, `T_ctrl = 0.75*S_T − 0.25*S_S`, then divide both by `max(S_ctrl, T_ctrl)`.
- **Scoring**: multiply the PSSM values of the observed residues at their positions, scaled by the probability of a random peptide (Supplementary Note 2). The Kinase Library website/package instead reports a **log2-summed score** (positive = preferred over random, negative = selected against).
- **Percentile rank**: score the substrate against a kinase's *a priori* score distribution over a reference phosphoproteome of **82,735 Ser/Thr sites** (Ochoa et al. 2020, doi:[10.1038/s41587-019-0344-3](https://doi.org/10.1038/s41587-019-0344-3)); the percentile is the substrate's rank in that kinase-specific distribution. Site annotation in the paper covers **89,752 sites** (82,735 HT + 7,017 LT-only). Their operational "biochemically favoured kinase" rule = **kinase in the top 15** for that site (and/or percentile >90).
- **Phospho-priming is handled explicitly**: pThr and pTyr are rows in the matrix (pSer assumed equivalent to pThr). **129 of 303** kinases select a phosphorylated residue as their single most preferred residue at ≥1 position (GSK3, CK1, CK2 families plus, unexpectedly, GRK and YANK families). This is directly relevant to the brief's "unorthodox kinases / priming" requirement: the atlas can score a primed site if we pass lowercase s/t/y in the motif window.
- **Caveat the authors state**: first-order peptide arrays do not measure cooperative/second-order effects, docking, scaffolding, or localisation; ~1/3 of sites match many kinases ("promiscuous" sites), and negative selection (being disfavoured by non-cognate kinases) carries much of the signal.

### 3.2 Yaron-Barir et al. 2024 (Tyr)
Yaron-Barir TM, Joughin BA, Huntsman EM, et al. "The intrinsic substrate specificity of the human tyrosine kinome." *Nature* 629:1174–1181. doi:[10.1038/s41586-024-07407-y](https://doi.org/10.1038/s41586-024-07407-y). PMC11136658.
- **78 catalytically active conventional Tyr kinases** (plus non-canonical kinases and 12 *C. elegans* orthologues; 93 kinases in the clustering, 109 total profiled). Positions −5..+5; X randomised over all aa except Tyr and Cys.
- Normalisation: column-normalised by the sum of the **18** randomised aa (16 for the customised library used for some kinases); Cys median fixed. Scoring: PSSM values scaled by 18 (or 16) so the mean randomised aa = 1, rounded to 1e-4, then multiplied.
- **Percentile reference distribution: 5,431 Tyr sites with localisation probability >0.75** (Ochoa 2020). Annotation covers 7,315 known Tyr sites.
- **Priming is a dominant feature**: **47 of 78** conventional Tyr kinases select a phosphorylated residue as their single most preferred residue overall, most often pTyr at −1, +1 or +2; a phosphoresidue at +3 acts as a "phospho-obstruction". They validate priming-aware scoring on FGFR1 autophosphorylation by scoring earlier-phosphorylated Tyr residues as pTyr.
- ~30% of mapped human Tyr sites are poorly matched by any canonical Tyr-kinase motif → ceiling on motif-only prediction for Tyr.

### 3.3 Are they trained on PSP? (crucial)
**No.** The PSSMs come only from peptide arrays. PSP enters in two places, both evaluation-side:
- Johnson 2023: "Experimentally validated kinase–substrate relationships were obtained from PhosphoSitePlus (July 2021)"; used for the AUCDF percentile-distribution analyses and the "more reports → better percentile" trend (n = 9,073 / 3,945 / 544 / 224 / 201 pairs with 1/2/3/4/5+ reports).
- Yaron-Barir 2024: PSP (April 2022), Tyr sites on human proteins, **filtered out sites with another phosphorylated residue within 5 aa** and kinases not profiled.
- The **percentile calibration** uses the Ochoa 2020 MS phosphoproteome, not PSP K–S labels. So percentiles are calibrated against the *observed phosphoproteome's* composition — a mild but real dependence on which sites get detected by MS (abundance/study bias), and sites in our test set that are also in the Ochoa reference contribute to their own kinase's background distribution. Harmless in practice, but worth one sentence in our methods.
- The trend "pairs with more literature reports get better percentile scores" is itself evidence that **PSP annotation quality correlates with motif-match**, i.e. single-report PSP pairs are enriched in indirect/wrong assignments. Directly relevant to our label-quality policy.

### 3.4 Getting the PSSMs
- **Python package `kinase-library`** (PyPI, v1.8.0; source <https://github.com/TheKinaseLibrary/kinase-library>; **CC BY-NC-SA 3.0**). Ships the matrices in-tree: `src/kinase_library/databases/matrices/{ser_thr,tyrosine}/{raw,densitometry,norm,log2}/*.tsv` — **311 Ser/Thr** and **93 Tyr** matrices in the current release (v1.2.0 data: 311 S/T, 78 Tyr + 15 non-canonical Tyr; v0.1.0 is the "legacy" 303/78/15 exactly as published). Rows are positions (S/T: −5..+4 plus `P`; Tyr: −5..+5 plus `P`), columns are the 20 aa plus **`s`, `t`, `y`** (phospho-residues) — so priming is directly supported by passing lowercase letters in the window.
- Reference distributions for percentiles ship as `databases/substrates/ochoa/scored_phosprots/{ser_thr,tyrosine}_phosphoproteome_scored.parquet` (63 MB / 2.3 MB).
- API surface (from `modules/scoring.py`): `score_protein(seq, kinases=None, phosphoacceptor=['S','T','Y'], pp=False, st_fav=True, non_canonical=False, score_promiscuity_threshold=1, percentile_promiscuity_threshold=90)`; objects `Substrate`, `PhosphoProteomics(...).predict(kin_type=..., st_fav=..., non_canonical=...)`; also binary enrichment, differential phosphorylation and MEA/GSEA enrichment modules. Key flags: **`pp=True`** enables phospho-residues (lowercase s/t/y) in the input window; **`st_fav`** applies the Ser-vs-Thr acceptor favourability; **`non_canonical`** adds the 15 non-canonical Tyr kinases. Docs site is a stub ("Coming Soon") — read the source and `src/notebooks/`.
- Raw supplementary tables (what SELPHI2.0 uses), no package needed:
  - Ser/Thr Supplementary Table 2: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM4_ESM.xlsx`
  - Tyr Supplementary Table 2: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM4_ESM.xlsx`
  - Ser/Thr Supplementary Table 3 = 89,752 scored sites with percentiles/ranks/promiscuity index (useful to sanity-check our reimplementation).
- R alternative: Bioconductor `JohnsonKinaseData` (PWMs for the 303 Ser/Thr kinases + background scores for 85,603 sites).
- Also `psite-annotation` (PyPI, kusterlab) wraps kinase-library motif scoring for phosphosite tables.

### 3.5 How to use as a baseline — and the pitfalls
- Two baseline variants, both zero-training: **(a) raw log2 score**, **(b) percentile rank** (kinase-wise calibrated). Percentile is the fairer cross-kinase comparator; raw score is the fairer within-site comparator. Report both.
- Also report the **"rank of the true kinase among all atlas kinases for this site"** metric — that is the deployment question and the one the atlas papers themselves use (top-15 rule).
- Pitfalls to control:
  - **Name mapping**: atlas kinase names are Cantley-style (AKT1, PKACA, CK2A1, ERK1, P70S6K, …), not HGNC. `src/kinase_library/databases/kinase_data/kinome_information.tsv` carries the mapping; verify manually for ambiguous ones (CK1 family, PDHK/BCKDK, GRKs, MAP2Ks).
  - **Coverage**: only ~311 S/T + ~93 Y kinases; PSP has 421 human kinases with substrates. Any comparison must be restricted to the intersection *and* we must report what fraction of our test kinases the atlas simply cannot score (the main axis on which we can beat it).
  - **Calibration on PSP?** No recalibration is needed or appropriate. If we ever fit a per-kinase threshold or a logistic calibration of percentiles against PSP labels, that **is** training and must be done inside training folds only. Flag it if we do.
  - **Priming**: to score primed sites we must decide which neighbouring residues to mark as phospho. Options: (i) always uppercase (ignore priming, the conservative default); (ii) mark a neighbour lowercase if it is a known PSP/Ochoa phosphosite; (iii) both, as two baseline arms. (ii) is a mild information leak from the same site-annotation universe; declare it.
  - Multi-phosphorylated peptides were excluded from the atlas analyses; we should at minimum stratify results by "has another phosphosite within ±5".

---

## 4. Survey of other methods (emphasis 2020–2026)

Ordered roughly by relevance to us. "Leakage" column = what they did about it.

### 4.1 Motif/PSSM + network era (pre-2020, still the default comparators)
| Method | Modalities | Model | Leakage handling | Availability |
|---|---|---|---|---|
| **Scansite** (2003, doi:[10.1093/nar/gkg584](https://doi.org/10.1093/nar/gkg584)) | peptide-library PSSMs, ~60 motifs | PSSM score | n/a (no training) | web |
| **NetPhorest** (2008) / **NetworKIN** (2007) / **KinomeXplorer** (Horn 2014, doi:[10.1038/nmeth.2968](https://doi.org/10.1038/nmeth.2968)) | motif NN/PSSM **+ STRING network proximity** (NetworKIN = NetPhorest × network score) | ensemble NNs + network diffusion | random splits, family-level | **websites dead**: networkin.science, netphorest.science unreachable; kinomexplorer.info returns 503 (checked 2026-10-26). Only usable via mirrors/KinomeXplorer-DB dumps or the KSTAR-bundled NetworKIN graph. |
| **NetPhos 3.1** (Blom 2004) | sequence only | NN | n/a | DTU web service; only ~17 kinases → near-useless as a comparator |
| **PhosphoPICK** (Patrick 2015, doi:[10.1093/bioinformatics/btu663](https://doi.org/10.1093/bioinformatics/btu663)) | PPI, co-expression, cell-cycle abundance, localisation + motif | Bayesian network | pair-level CV | server bioinf.scmb.uq.edu.au unreachable (checked 2026-10-26) |
| **GPS 5.0** (Wang 2020) / **GPS 6.0** (Chen 2023, doi:[10.1093/nar/gkad383](https://doi.org/10.1093/nar/gkad383)) | sequence only (PSP(30,30) windows, 10 feature encodings) | LightGBM + DNN + penalised logistic regression; general model pretrained on 490,762 sites/71,407 proteins, then **transfer learning** to 577 kinase-specific predictors (group/family/single-kinase) from 30,043 ssKSRs | **CD-HIT 40% clustering** of the general p-site data + timestamp-based test split; kinase-specific models: 10-fold CV for kinases with ≥30 sites; negatives = non-phosphorylated S/T/Y on the same proteins | <https://gps.biocuckoo.cn/download.php> (local PyInstaller packages, Win/Mac/Linux). Reachable. Best-performing comparator in the SELPHI2.0 external benchmark (AUC 0.81). |
| **KinasePhos 3.0** (Ma 2023, doi:[10.1016/j.gpb.2022.06.004](https://doi.org/10.1016/j.gpb.2022.06.004)) | sequence + physicochemical features, SHAP explanations | 771 SVM/ML models (group/family/individual), ≥15 positives each; 41,421 sites from PSP+UniProt+GPS5.0+Phospho.ELM, 1,380 kinases | redundancy removal within clusters; kinases without UniProt IDs placed by phylogenetic tree | <https://github.com/tom-209/KinasePhos-3.0-executable-file> (MIT). Reported "average accuracy 87.2%" across 771 models — accuracy on balanced sets, not comparable. |
| **LinkPhinder** (Nováček 2020, doi:[10.1371/journal.pcbi.1007578](https://doi.org/10.1371/journal.pcbi.1007578)) | knowledge graph over PSP (kinase–site–substrate triples) + MEME-derived motifs | KG link prediction (TransE-family) | 100 random 90:10 splits of PSP triples; negatives = random kinase perturbations | predictions download link dead (checked); splits at figshare doi:[10.6084/m9.figshare.12179925.v1](https://doi.org/10.6084/m9.figshare.12179925.v1) |

**LinkPhinder is the single most instructive cautionary tale in this literature.** On PSP random splits it reports AUPR 0.973 / AUROC 0.968, beating GPS, NetworKIN, NetPhorest, Scansite, NetPhos, PhosphoPredict. On an *independent* dataset (site-specific interactions of 103 kinases in cancer cells) **every method collapses to chance**: AUROC 0.532 (LinkPhinder), 0.509 (GPS), 0.511 (NetworKIN), 0.510 (NetPhorest), 0.502 (Scansite/PhosphoPredict), 0.505 (NetPhos). A 0.97 → 0.53 drop is what random pair splits over a KG buy you: the model learns node degree and study bias, not biochemistry. Our benchmark must be designed so this cannot happen silently.

### 4.2 Deep learning on sequence (2017–2022)
- **MusiteDeep** (Wang 2017, doi:[10.1093/bioinformatics/btx496](https://doi.org/10.1093/bioinformatics/btx496)): CNN + attention, family-level only.
- **DeepPhos** (Luo 2019, doi:[10.1093/bioinformatics/bty1051](https://doi.org/10.1093/bioinformatics/bty1051)): densely connected CNN blocks, 8 group/family models.
- **PhosIDN / PhosIDNSeq** (Yang 2021, doi:[10.1093/bioinformatics/btab551](https://doi.org/10.1093/bioinformatics/btab551)): sequence encoder **+ PPI embeddings** fused; 10 group/family/subfamily/kinase models. CD-HIT 40% on substrate proteins; ~20% random holdout for kinase-specific sites. Code <https://github.com/ustchangyuanyang/PhosIDN> (no license, unmaintained since 2020). Notable as the clearest demonstration that **adding PPI helps over sequence alone** within a family.
- **EMBER** (Kirchoff & Gomez 2022, doi:[10.1093/bioinformatics/btac083](https://doi.org/10.1093/bioinformatics/btac083)): multi-label prediction over 5 kinase families, siamese/metric-learning flavour; one model, family-level. <https://github.com/gomezlab/EMBER> (no license, 2020). Weak in both Phosformer benchmarks (AUPRC 0.25 in Phosformer-ST Table 1).
- **DeepKinZero** (Deznabi 2020, doi:[10.1093/bioinformatics/btaa013](https://doi.org/10.1093/bioinformatics/btaa013)): the first explicit **zero-shot** formulation — BiRNN phosphosite encoder + kinase representation from sequence embeddings and the kinase functional taxonomy, joined by a bilinear compatibility function, so kinases with zero known sites get predictions. Code <https://github.com/Tastanlab/DeepKinZero> (no license, last push 2021). Conceptually the ancestor of Phosformer/KolossuS; its evaluation protocol (hold out *kinases*, measure accuracy on sites of unseen kinases) is the right one and we should copy it.

### 4.3 Co-phosphorylation / phosphoproteomics-driven
- **CoPhosK** (Ayati 2019, doi:[10.1371/journal.pcbi.1006678](https://doi.org/10.1371/journal.pcbi.1006678)): naive-Bayes over co-phosphorylation across MS datasets, seeded with NetworKIN/PSP.
- **Ayati 2023 PSB** (doi:[10.1142/9789811270611_0007](https://doi.org/10.1142/9789811270611_0007)): kinase/phosphosite "functional landscape" features.
- **CoPheeMap / CoPheeKSA** — Jiang W, Jaehnig EJ, Liao Y, Shi Z, Yaron-Barir TM, Johnson JL, Cantley LC, Zhang B. "Deciphering the dark cancer phosphoproteome using machine-learned co-regulation of phosphosites." *Nat Commun* 16, 2025. doi:[10.1038/s41467-025-57993-2](https://doi.org/10.1038/s41467-025-57993-2). PMC11926083. **The best-engineered recent integrative method; read it before designing our features.**
  - CoPheeMap: XGBoost classifier for "are these two phosphosites co-regulated", features = 15-mer sequence similarity, STRING kinase-interaction-profile similarity of host proteins, pairwise site correlations in each of 11 CPTAC cancer types (1,195 tumours). Ground truth = site pairs sharing a kinase (192,926 positives) vs pairs regulated by kinases from different groups (1,956,116 negatives). Independent-test AUROC **0.78**; applied to >3 billion pairs at FPR 0.2% → **764,049 co-regulated pairs over 26,280 sites**.
  - CoPheeKSA: XGBoost over **55 features** = 16 CoPheeMap Node2Vec embedding dims + 16 "KMap" embedding dims + 22 dynamic features (correlations of site abundance with kinase protein abundance and with inferred kinase activity, per cohort) + **1 static feature: the Kinase Library motif score**. Positives: 2,353 KSAs (GPS 5.0 curated set ∩ CoPheeMap). Negatives: substrates of kinase A assigned to kinase B when A and B are in **different families**, excluding sites with any overlapping upstream kinase group → 114,530 negatives; train/val/test at 1:10 pos:neg, 10× Monte-Carlo CV (235 positives per test group). Median AUROC **0.95**; at threshold 0.76 (LLR 5.5) → **24,015 KSAs over 9,399 sites and 104 kinases incl. 26 understudied**.
  - **They explicitly prevent the leak SELPHI2.0 has**: "to avoid information leakage, we computed kinase activity scores for each positive KSA separately, removing the corresponding phosphosite in our computation", and they match missing-value rates between positives and negatives so the model cannot exploit missingness. They also intentionally held out CSNK1A1/CSNK2A1/CSNK2A2 from the ground truth as a sanity check.
  - Code: <https://github.com/bzhanglab/CoPheeMap> (**MIT**). Data: CPTAC pan-cancer portal <https://pdc.cancer.gov/pdc/cptac-pancancer>; HCC matrices at <https://zenodo.org/records/14553766>.
  - **Residual criticism**: still a pair-random split (no substrate-protein grouping), only 104 S/T kinases, negatives assume cross-family pairs are true negatives — which Johnson 2023 explicitly contradicts (kinases from different families often share motifs). Their own analysis shows ground-truth positives sometimes score *worse* on the Kinase Library than their predictions (CAMK2G, CDK7, DYRK1B, PAK4) — read as evidence of PSP/GPS label noise.

### 4.4 Knowledge-graph / embedding (2023–2025)
- **KSFinder** (Anandakrishnan 2023, doi:[10.7717/peerj.16164](https://doi.org/10.7717/peerj.16164)): KG embedding + MLP; predicts **substrate proteins** (not sites) for 432 kinases incl. 68 understudied. Negatives combine experimentally validated non-interacting pairs, different-subcellular-localisation pairs, and random sampling, explicitly to remove entity-representation bias (their critique of PredKinKG, where one substrate appeared in 216 negative triples vs 2 positive, is a good illustration of how negative construction biases supervised KG classifiers). <https://github.com/manju-anandakrishnan/ksfinder>.
- **KSMoFinder** (2025, doi:[10.1093/bioadv/vbaf289](https://doi.org/10.1093/bioadv/vbaf289), PMC12664573): site/motif-level successor. KG of protein context + kinase-specific motifs; positives from iPTMnet + PSP (19,574), **biologically motivated synthetic negatives** (kinase × non-interacting substrate protein; kinase × motifs it is experimentally *unfavourable* toward, taken from the Cantley atlas) (69,934 negatives). ROC-AUC **0.851**, PR-AUC 0.839 on a balanced test set; covers **430 human kinases**. They report their KG embeddings beat ProtT5/ESM2/ESM3 embeddings on this task. <https://github.com/manju-anandakrishnan/KSMoFinder>. Their survey of negative-sampling strategies across tools is the most useful part for us.

### 4.5 PLM-based (2023–2026) — the current frontier
| Method | PLM / encoding | Training labels | Split | Reported | Code |
|---|---|---|---|---|---|
| **Phosformer** (2023) | custom 6-layer RoBERTa, kinase domain + 11-mer | PSP/ELM/PhosphoNetworks, 800 kinases | random pairs (leaky, §2.3) | AUROC 0.885 | HF + GitHub, CC BY-NC-ND |
| **Phosformer-ST** (2024) | ESM-2 + cross-attention, multitask MLM, kinase domain + 15-mer | **Johnson 2023 percentiles** (in vitro) | CD-HIT 70% kinases + Hamming-5 peptide clusters | 0.992 (distillation); **0.886 unseen kinases** | GitHub, CC BY-NC-ND |
| **KolossuS** (Jha, Shonai, …, Singh R, Soderling SH; bioRxiv 2025.04.27.650849, doi:[10.1101/2025.04.27.650849](https://doi.org/10.1101/2025.04.27.650849)) | **ESM-2 15B** (also 650M; ESM-C 6B compared), **full kinase sequence** + 15-mer substrate window, mean-pooled, MLP projection into a shared space, score = **cosine similarity** | pretrain on PSPA atlases (393 kinases), fine-tune on **Sugiyama in vitro kinome** (354 kinases), evaluate on PSP (human/mouse/rat) | **phylogeny-aware** ~80:10:10 (CK1 family held out) | PSP human AUROC **0.755** (S/T), **0.838** (Tyr); atlas-only 0.754/AUPR 0.667; +fine-tune 0.762/0.679; raw ESM-2 zero-shot 0.503 | `pip install kolossus`, <https://github.com/rohitsinghlab/KolossuS> (**no license file**), weights for 650M + 15B in-repo; web <https://singhlab.net/kinase> |
| **Esmaili et al. 2025** (doi:[10.1016/j.csbj.2025.03.003](https://doi.org/10.1016/j.csbj.2025.03.003), PMC11968300) | ESM2-650M encoder + autoregressive decoder, **protein-level** (kinase × substrate protein, not site) | 11,391 pairs after **CD-HIT 0.7** redundancy reduction, one pair per cluster | random 90:10 (+ frequency-stratified and "Novel Kinase Pair" zero-shot sets) | top-k and zero-shot analyses; **full kinase sequence beat kinase-domain-only** | <https://github.com/farz1995/substrate_kinase_prediction> |
| **Prot2Token-style** (bioRxiv 2025.03.12.642740) | ESM-2 650M encoder + autoregressive decoder, multi-task (kinase group, KS interaction, p-site) | PSP-derived | tasks trained separately **because of leakage concerns**; zero-shot set of **29 dark kinases** (8,026 samples) | SOTA claimed | preprint |
| **SAGEPhos** (ICLR 2025, arXiv 2502.07384) | ESM-2 sequence + AlphaFold2 structure fusion, kinase-conditioned | new "refined" dataset | unclear | "significantly outperforms baselines" | <https://github.com/ZhangJJ26/SAGEPhos> — **v2 WITHDRAWN by an author citing author-team disagreements. Do not cite as a result; do not use as a baseline.** |
| **SPINET-KSP** (2026, doi:[10.3390/kinasesphosphatases4010003](https://doi.org/10.3390/kinasesphosphatases4010003)) | GNN + ESM-3 embeddings + AlphaFold3 structures, kinase–substrate–**phosphatase** triads | — | — | claims to beat single-modality | niche journal; treat with caution |

**KolossuS is the most important recent result for us**, and its numbers deserve emphasis: a 15B-parameter PLM, trained on the *in vitro* atlases and the Sugiyama screen, evaluated on PSP literature labels with a phylogeny-aware split, reaches **AUROC 0.755 (S/T) / 0.838 (Tyr)**. That is the realistic ceiling for specificity-only models on literature labels — an order of magnitude away from the 0.95–0.99 numbers reported by methods that split randomly on their own labels. It also shows raw ESM-2 embeddings are useless zero-shot (0.503), i.e. all the signal comes from supervised specificity, not from the PLM prior.

### 4.6 Related but different tasks (do not confuse with KS prediction)
- **KSTAR** (Crowl 2022, doi:[10.1038/s41467-022-32017-5](https://doi.org/10.1038/s41467-022-32017-5)): patient-level *kinase activity* from phosphoproteomics by pruning NetworKIN's dense weighted graph into many sparse binary graphs + hypergeometric enrichment + Mann-Whitney against 150 random experiments. <https://github.com/NaegleLab/KSTAR>. Their key finding for us: **thresholded NetworKIN graphs are "highly problematic"**, and homologous kinases share up to 80% of evidence for sites in an experiment (EGFR/ERBB2), which their pruning reduces to 18%. Any prediction network we ship for activity inference must be audited for this kind of kinase-kinase redundancy.
- **KinPred** (Xue, Jordan, Rizvi, Naegle 2021, doi:[10.1371/journal.pcbi.1008681](https://doi.org/10.1371/journal.pcbi.1008681)): framework to re-run and unify proteome-level KS predictions on a current reference proteome, with analyses of study bias and network properties across three resources. <https://github.com/NaegleLab/KinPred> (MIT). **Directly relevant**: this is the established methodology for making heterogeneous predictors comparable (common reference proteome, common ID mapping, restriction to the human phosphoproteome). We should follow it when we run baselines.
- **benchmarKIN** — Müller-Dott S, et al. "Comprehensive evaluation of phosphoproteomic-based kinase activity inference." *Nat Commun* 16, 2025. doi:[10.1038/s41467-025-59779-y](https://doi.org/10.1038/s41467-025-59779-y). PMC12098709. Benchmarks *activity inference* methods × prior-knowledge libraries (PSP, PTMsigDB, GPS gold, OmniPath, iKiP-DB, NetworKIN). Findings we should internalise: methods matter little, **libraries matter a lot**; the curated combination (PSP+PTMsigDB+GPS gold) gives the best perturbation-benchmark AUROC **0.73**, OmniPath 0.71, iKiP-DB 0.59; NetworKIN gives the best activating-site-based AUROC 0.68 and **adding predicted targets boosts tumour-based performance**; on a common 17–31-kinase subset NetworKIN (0.77) and OmniPath (0.76) actually beat the curated combination (0.75); performance drops for "poor" (less-studied) kinases, particularly for NetworKIN and iKiP-DB. Kinase coverage: OmniPath 467, curated 414, PSP 390, GPS gold 352; interactions: iKiP-DB 26,786, NetworKIN 22,788. **This is the benchmark our final network should be evaluated in** (an orthogonal, downstream-utility evaluation that does not depend on PSP labels as ground truth).
- **iKiP-DB** (Kuznetsova/Bantscheff et al., doi:[10.1021/acs.jproteome.2c00198](https://doi.org/10.1021/acs.jproteome.2c00198)): in vitro kinase-to-phosphosite DB for >300 kinases — a candidate independent evaluation set (and a candidate extra training source, with the in vitro caveat).
- Structure-based: AlphaFold2 models of the active form of all **437** catalytically competent human kinase domains (doi:[10.1101/2023.07.21.550125](https://doi.org/10.1101/2023.07.21.550125); 2026 update doi:[10.1042/bcj20260137](https://doi.org/10.1042/bcj20260137) analysing substrate P−5..P+5 contacts with the activation loop). AlphaFold-Multimer for kinase–peptide complexes is **not** reliable at scale: the authors report it often cannot produce an active activation-loop conformation without templates, and forcing the active state disrupts substrate placement. Global confidence scores (ipTM/pDockQ) miss transient, IDR-mediated kinase docking; local-interaction-score variants (LIS/iLIS, e.g. doi:[10.1101/2025.10.10.681672](https://doi.org/10.1101/2025.10.10.681672)) are proposed as fixes. Structure-guided *peptide design* works (Subtimizer, doi:[10.1101/2025.07.04.663216](https://doi.org/10.1101/2025.07.04.663216)) but that is design, not proteome-scale prediction. **Recommendation: do not invest in per-pair folding** (437 kinases × 10⁵ sites is infeasible anyway); use structure only as cheap substrate-side features (pLDDT/disorder, accessibility, domain context) and possibly kinase-side docking-groove descriptors.

---

## 5. Synthesis

### 5.1 Standard evaluation protocols and their flaws
What the field does, and why each is broken:

1. **Random splits over kinase–substrate pairs** (SELPHI2.0, Phosformer, LinkPhinder, CoPheeKSA, Esmaili, KSFinder/KSMoFinder). A substrate protein's other sites, and the same peptide paired with other kinases, are in the training set. Measured on Phosformer's own released splits: **51% of test-positive peptides are training positives** under another kinase and **2.7% are literal (kinase-domain, peptide) duplicates**. Fix: group by substrate protein and by homologous protein cluster (mmseqs2 at e.g. 30–50% identity), split sites, never pairs.
2. **No homology control on the kinase side**. Only Phosformer-ST (CD-HIT 70%), GPS 6.0 (CD-HIT 40% on substrates), PhosIDN/Esmaili (CD-HIT 40%/70%) do anything. Paralogues (CDK1/CDK2, AKT1/2/3, PRKACA/B, ERK1/2, EGFR/ERBB2) make "unseen kinase" evaluation meaningless unless whole paralogue groups move together. Fix: leave-kinase-**group**-out, not leave-one-kinase-out; cluster kinase domains and hold out clusters.
3. **Cross-species duplication**. Phosformer's test set mixes human/mouse/rat/bovine orthologue pairs of identical peptides. If we pool species for training (as the brief permits), orthologous sites must be assigned to the same group.
4. **Study bias dominates pooled metrics.** ~80% of annotated sites belong to ~20% of kinases; our provided PSP file has CDK1 at 694 rows while 146 kinases have <5. Pooled AUC is then essentially "performance on CDK1/PKA/CK2/SRC". Johnson 2023's own analysis shows percentile scores *improve* with the number of literature reports — meaning well-studied pairs are both over-represented and easier. Fix: macro-average over kinases, report per-kinase distributions, stratify by annotation count (rich/medium/poor, as benchmarKIN does), and never report a single pooled number as the headline.
5. **Negative sampling choices drive the result more than the model.**
   - Random pairs (SELPHI2.0, LinkPhinder): trivially separable; inflates AUC; contains unknown true positives.
   - Hard negatives = phosphosites of other kinases (Phosformer, Phosformer-ST): much better, but conflates "not reported" with "not phosphorylated".
   - Cross-family negatives (CoPheeKSA, LinkPhinder): contradicted by Johnson 2023 (distant kinases share motifs) — KSMoFinder flags this explicitly.
   - Experimentally unfavourable motifs (Phosformer-ST, KSMoFinder): the most defensible, but ties your labels to the atlas.
   - Embedding-distance hard negatives (Esmaili): elegant, but risks selecting actual positives among close paralogues.
   Fix: evaluate on a **fixed, declared negative definition** (we propose: all other detected phosphosites of the same residue type, within the same site universe = "rank the true kinase among all kinases for this site" and "rank the true sites among all sites for this kinase"), plus a sensitivity analysis over negative definitions. Never tune the negative definition to improve a result.
6. **Balanced test sets at 1:1** (Phosformer, Phosformer-ST, KSMoFinder) misstate precision by orders of magnitude relative to the real ~1:300 (per site) or ~1:10⁴ (per kinase) prevalence. Fix: report AUPRC and precision@k at the true prevalence of the evaluation universe, and state the prevalence.
7. **Self-referential labels.** Phosformer-ST's 0.992 is distillation of the Cantley PSSMs; SELPHI2.0's co-regulation features are derived from PSP-defined kinase activities (and partly from the Hijazi data used for validation). Fix: a written feature-provenance audit — for every feature, state whether its computation touched the test-fold labels.
8. **Comparators run out-of-the-box with different coverage.** SELPHI2.0 compared against NetPhos with 17 kinases; Phosformer had to invent inputs (family labels, PPIs) for competitors. Fix: follow KinPred — common reference proteome, common ID mapping, and always report coverage (fraction of test pairs each method can score) alongside accuracy; score uncovered pairs as "no prediction" and report metrics both with and without them.

### 5.2 Which modalities actually buy accuracy
Ranked by the evidence available:

1. **Peptide motif / specificity (PSSM or learned)** — the backbone. Necessary but far from sufficient: motif-only ceiling on literature labels is ~0.75 AUROC (S/T) per KolossuS; ~1/3 of S/T sites and >30% of Tyr sites are motif-ambiguous or motif-suboptimal (Johnson 2023, Yaron-Barir 2024). Also cannot separate paralogues with near-identical motifs.
2. **Kinase identity encoded by sequence (PLM)** — the key enabler of understudied-kinase coverage and leave-kinase-out generalisation (DeepKinZero → Phosformer → Phosformer-ST 0.886 on unseen kinases → KolossuS). ESM-2 embeddings alone are useless (0.503); the gain comes from supervised specificity learning *plus* kinase-sequence conditioning that transfers across the kinome.
3. **Co-regulation / co-phosphorylation from large phosphoproteomics** — the largest *complementary* gain in the integrative papers. CoPheeKSA: network-embedding + dynamic correlation features beat PSSM-only by a wide margin (median AUROC 0.95 vs PSSM-only, their Fig. 3b). SELPHI2.0 likewise leans on it. **But this is also where the leakage lives** — it must be computed leave-site-out and fold-aware.
4. **Context: PPI/STRING, co-expression, localisation, abundance** — consistent moderate gains (PhosphoPICK, NetworKIN, PhosIDN, KinomeXplorer, CoPheeKSA's STRING feature). They also encode study bias (well-studied proteins have more interactions), so they must be evaluated under leave-protein-out.
5. **Perturbation/inhibitor phosphoproteomics** (Hijazi 2020, decryptM) — strong evidence for *directionality and causality*, under-exploited; our provided compendium + `kinhib_inhib_comb__1.csv` is precisely this modality.
6. **Site-intrinsic features** (disorder, accessibility, conservation, functional score) — improve *phosphorylatability* but are kinase-agnostic; they inflate pooled metrics without improving kinase assignment. Include them, but verify via per-site kinase-ranking metrics that they are not carrying the result.
7. **Structure (AF2/AF3 complexes)** — no demonstrated proteome-scale gain yet; AlphaFold-Multimer unreliable for kinase–peptide. Low priority.
8. **Priming / phospho-residues in the motif** — demonstrably important (129/303 S/T kinases, 47/78 Tyr kinases prefer a phospho-residue at some position) and essentially ignored by every ML method except via the atlas matrices. **This is a concrete differentiation opportunity for us**: encode neighbouring known phosphosites as lowercase s/t/y (or as an extra channel) and evaluate specifically on primed sites.

### 5.3 Recommendations for our benchmark design
1. **Split unit**: substrate-protein clusters (mmseqs2, `/home/tools/mmseqs/bin/mmseqs`, e.g. 50% identity / 80% coverage; also merge cross-species orthologues) → assign whole clusters to folds. Sites, never pairs. Seeds [0,1,2,3,4].
2. **Three difficulty tiers**, reported separately:
   - T1 **unseen sites** on seen proteins is *excluded* (leaky) — do not report it.
   - T1' **unseen substrate proteins** (protein-cluster holdout), seen kinases.
   - T2 **unseen kinases** — hold out whole kinase-domain clusters (CD-HIT/mmseqs ~70% on kinase domains, as Phosformer-ST did), so paralogues move together.
   - T3 **unseen kinase × unseen protein** (both held out) — the honest "dark" scenario.
3. **Metrics**: primary = per-kinase AUPRC and per-site "rank of true kinase" (MRR / Hit@k / top-15 rate, matching the atlas convention), macro-averaged over kinases with bootstrap CIs; secondary = pooled AUROC/AUPRC (for comparability with the literature only, clearly flagged). Always stratified by kinase annotation count (rich ≥50, medium 10–49, poor <10) and by residue type (S/T vs Y) and by whether the site has a neighbouring phosphosite (priming stratum).
4. **Negatives**: fix one evaluation universe (all S/T/Y sites in our site set, correct residue type for the kinase) and report the real prevalence. Report a sensitivity analysis over 2–3 negative definitions (random, residue+abundance-matched, atlas-unfavourable) but pre-register the primary one.
5. **Label quality arms**: PSP in vivo vs in vitro vs both; ≥2 reports vs 1 report (Johnson's report-count analysis predicts single-report pairs are noisier). Decide once, document, do not revisit to improve numbers.
6. **External/orthogonal validations** that do not use PSP labels: (a) Sugiyama 2019 in vitro kinome screen, (b) Hijazi 2020 inhibitor-derived edges (only if not used as features in that arm), (c) iKiP-DB, (d) downstream **kinase-activity inference** via benchmarKIN (doi:10.1038/s41467-025-59779-y) — the last is the strongest evidence of practical utility and is independent of our label set.
7. **Feature-provenance audit table** in the design doc: every feature × "does its computation see test-fold labels?" Hard rule: any feature derived from kinase activity or kinase substrate sets must be recomputed per fold with the target site (and ideally the target protein) excluded.
8. **Kinase-redundancy audit** (KSTAR's point): report the overlap of predicted substrate sets between paralogue pairs; a model that cannot separate EGFR/ERBB2 or CDK1/CDK2 is not useful for activity inference even at high AUROC.

### 5.4 Which baselines to retrain, and concretely how
Priority order:

**B1. Cantley Kinase Library (no retraining).**
- `pip install kinase-library==1.8.0` into `/workspace/.venv` (CC BY-NC-SA 3.0 — academic use fine; cite Johnson 2023 + Yaron-Barir 2024). Or read matrices directly from the repo tree `src/kinase_library/databases/matrices/{ser_thr,tyrosine}/{log2,norm,densitometry}/`.
- Arms: (i) log2 score; (ii) percentile vs the shipped Ochoa reference distribution; (iii) priming-aware (lowercase s/t/y for neighbouring known phosphosites) vs priming-naive; (iv) `st_fav` on/off; (v) `non_canonical` Tyr on/off.
- Sanity check: reproduce Johnson Supplementary Table 3 percentiles for a sample of sites; and reproduce the "PHKG1/PHKG2 top-2 for PYGL S15" and "ATM top for TP53 S15" cases.
- Record coverage: which of our test kinases have no atlas matrix (expect ~100+ of 421).

**B2. Phosformer (retrain; cheapest, highest value).**
- `git clone https://github.com/esbgkannan/phosformer` for the architecture (`Phosformer/roberta_modified.py`, `modules.py`); note **CC BY-NC-ND** — prefer reimplementing the ~100-line model wrapper in our `src/` and loading only the published MLM-pretrained encoder (HF `waylandy/phosformer`) as initialisation, or substituting ESM-2 35M/150M.
- Changes needed: (a) replace their dataset builder with ours (our PSP label policy, our site universe, our kinase-domain definitions); (b) negatives (easy+hard) generated **within training folds**; (c) our grouped splits + leave-kinase-cluster-out; (d) keep focal loss γ=2, AdamW 2e-5, resample/shift/mask augmentations, 1:16 pos:neg; (e) log per-kinase metrics.
- Also run the **published checkpoint zero-shot** on our test set, labelled "contaminated upper bound".
- Compute: ≲1 GPU-hour/epoch on one RTX 4090; 5 seeds × ~10 epochs ≈ a few GPU-hours. Use **one** 4090.

**B3. Phosformer-ST / atlas distillation (optional, as a pretraining stage not a baseline).**
- `git clone https://github.com/gravelCompBio/Phosformer-ST`; weights via the link file in the repo (HF/Zenodo). CC BY-NC-ND again → reimplement the ESM-2 + cross-attention head.
- Use as: pretrain on atlas percentiles (their exact label rule: >90 positive, <50 hard negative, with their cluster-dedup), then fine-tune on our literature labels. Report whether this pretraining helps on our splits — that is a genuinely novel, decision-relevant experiment.
- Do **not** quote their 0.992.

**B4. SELPHI2.0 (retrain; the integrative-features baseline).**
- `git clone https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0` (Apache-2.0) and `https://github.com/alussana/selphi2.0-pssms` (Nextflow; no license — reimplement rather than redistribute).
- Reuse: `bin/compute_pssm_scores.py` + `bin/get_pssm_{ser_thr,tyr}_*.py` (clean PSSM implementation), the feature list from `data/psp/pos_set_features_20240910.tsv`, `src/model_training/learn.py` (RF + grid search).
- Must change: (1) **kinase-activity features recomputed per fold with the target site removed from the kinase's regulon** (follow CoPheeKSA's wording); (2) grouped splits; (3) positives/negatives from our label policy; (4) drop or re-derive `netpho_max_KIN` (a NetPhos kinase score — another model's output, and NetPhos covers 17 kinases).
- Feature inputs to fetch (another agent's job; URLs only): Ochoa 2020 functional landscape supplementary tables; GTEx v8; HPA v23.0; Mertins 2016; Hijazi 2020 (we already have the compendium + inhibitor matrix); Cantley supplementary XLSX (URLs in §3.4). Avoid the 10.7 GB Zenodo feature table unless we decide to reuse their exact values.
- Compute: CPU-only; the RF itself is minutes-to-hours on 64 cores; the PSSM scoring of 421 × ~240k pairs is the bulk (hours, parallelisable).

**B5. GPS 6.0 (run as published; do not retrain).** Downloadable local package at <https://gps.biocuckoo.cn/download.php>; it was the strongest comparator in the SELPHI2.0 external benchmark. It is trained on PSP/GPS-gold, so our test sites are contaminated for it — run it, report it, and label it "trained on overlapping labels; upper bound". Retraining it is not practical (577 transfer-learned models, no training pipeline released).

**B6. KolossuS (run as published).** `pip install kolossus`; CLI `kolossus-cli -p pairs.txt -s seqs.fasta -o out` (pairs = `kinase_id substrate_id site_offset`), with `--model-small` for the ESM-2-650M variant; `kolossus-extract` to precompute embeddings. The 15B model needs substantial GPU memory — use the 650M variant first, and only one 4090. No license file → use for benchmarking, do not redistribute. Trained on atlases + Sugiyama (in vitro) and evaluated on PSP, so it is the **cleanest published comparator for our setting** and the number to beat is AUROC ~0.755 (S/T) / ~0.838 (Tyr) on PSP-style labels.

**B7. Not worth retraining**: NetworKIN/NetPhorest/KinomeXplorer (servers dead, code unmaintained), PhosphoPICK (server dead), NetPhos (17 kinases), LinkPhinder (download dead; and its own independent evaluation is at chance), EMBER/MusiteDeep/DeepPhos/PhosIDN (family-level only, superseded). Cite them, run GPS 6.0 as the representative of the classical family, and state why the others were excluded.

**B8. Worth reading before our feature design, not as baselines**: CoPheeKSA (MIT code, leave-site-out activity features, CPTAC-based co-regulation), KSMoFinder (negative-sampling taxonomy), KinPred (MIT, prediction unification), KSTAR (paralogue-redundancy pruning), benchmarKIN (downstream evaluation).

---

## 6. Open questions / requests for the PI
- Licensing: Phosformer and Phosformer-ST are **CC BY-NC-ND 4.0**. Retraining them arguably produces a derivative. Decision needed: reimplement the architecture from the paper (safe, ~half a day) vs train from their released code (fast, but ND). My recommendation: reimplement, and use only the published weights for the "as-published" zero-shot arm.
- `alussana/selphi2.0-pssms` and `rohitsinghlab/KolossuS` have **no license files**. Fine for running locally; we must not redistribute their code/weights in our repo.
- Decide early whether the priming-aware motif encoding uses known phosphosites from the same annotation universe as our labels (mild leak) or only from an independent source.

## 7. Key URLs (no large downloads performed)
- Cantley Ser/Thr PSSMs: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM4_ESM.xlsx`
- Cantley Tyr PSSMs: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM4_ESM.xlsx`
- Kinase Library package/data: <https://github.com/TheKinaseLibrary/kinase-library> · <https://kinase-library.phosphosite.org>
- Phosformer code/weights/data: <https://github.com/esbgkannan/phosformer> · HF `waylandy/phosformer` · HF `waylandy/phosformer_curated`
- Phosformer-ST: <https://github.com/gravelCompBio/Phosformer-ST> · <https://esbg.bmb.uga.edu/phosformer/>
- SELPHI2.0: <https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0> · <https://github.com/alussana/selphi2.0-pssms> · Zenodo doi:10.5281/zenodo.13821320 · <https://selphi2.com>
- KolossuS: <https://github.com/rohitsinghlab/KolossuS> · <https://singhlab.net/kinase>
- CoPheeMap/KSA: <https://github.com/bzhanglab/CoPheeMap> · CPTAC <https://pdc.cancer.gov/pdc/cptac-pancancer> · <https://zenodo.org/records/14553766>
- GPS 6.0: <https://gps.biocuckoo.cn/download.php>
- KinPred: <https://github.com/NaegleLab/KinPred> · KSTAR: <https://github.com/NaegleLab/KSTAR>
- KSFinder/KSMoFinder: <https://github.com/manju-anandakrishnan/ksfinder> · <https://github.com/manju-anandakrishnan/KSMoFinder>
- DeepKinZero: <https://github.com/Tastanlab/DeepKinZero> · EMBER: <https://github.com/gomezlab/EMBER> · PhosIDN: <https://github.com/ustchangyuanyang/PhosIDN>
