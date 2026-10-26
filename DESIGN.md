# DESIGN — Human kinase–substrate (site-level) predictor

**Version: v2 (2026-10-26 ~14:30) — FROZEN for implementation.** v2 = v1 + the amendments in §16, which take precedence over earlier sections where they conflict. This revision follows review round 1 (`docs/reviews/r1_*.md`). §15 maps each accepted or rejected reviewer change to the place it lands.

Supporting material:
- `LITERATURE.md` and `docs/lit/*`
- `docs/data_profile.md`
- `docs/labels_universe.md`
- `docs/data_acquisition.md`

## 0. Time gates (hard; T0 = 2026-10-26 13:20, deadline T0+24h)

| Gate | Deadline | Must exist (otherwise cut scope as listed) |
|---|---|---|
| G1 | T0+3.5h (16:50) | Splits frozen and committed. Evaluation harness, test seal and canaries pass their tests. B0/B1 baselines run on CV. |
| G2 | T0+7.5h (20:50) | R2 LightGBM (F1+F2+F3+F4+F6) beats B1 Cantley on CV. If not: drop R3 and the deep baselines, and debug features. |
| G3 | T0+13h (02:20) | R3 neural model plus stacking; B3 Phosformer and B4 SELPHI2.0 retrained on CV; ablations done. |
| G4 | T0+16h (05:20) | FREEZE.json committed, then the single batch evaluation on TEST. |
| G5 | T0+20h (09:20) | Downstream benchmarKIN, deployment model, coverage predictions, independent result review. |
| G6 | T0+23h | REPORT.md final and committed. |

Deferred, done only if ahead of schedule:
- U_all (unobserved S/T/Y)
- L3 label variant
- docking-motif features
- KolossuS 15B
- AF-Multimer pilot
- 5 seeds for every ablation; ablations use seed 0, and final configurations use seeds [0..4]

## 1. Goal and outputs (fixed by PROMPT.md)
For each human protein kinase *k* and human phosphosite *s* (S/T/Y), output a score that *k* directly phosphorylates *s*.

- **Priority 1:** accuracy, harmonious across kinases, including understudied ones.
- **Priority 2:** coverage of human phosphosites.

Deliverables:
- a per-(k,s) score;
- a within-kinase percentile over the universe;
- a prior-corrected probability with an honest uncertainty statement (§6.6);
- per-kinase substrate libraries for activity inference;
- code and reports.

Runtime (train and inference) is reported for every method.

## 2. Entities and universe (implemented: `src/ksp/{ids,kinases,sites}.py`; `docs/labels_universe.md`)
**Kinases.** `results/labels/kinase_list.tsv` lists 503 included kinases. They are keyed by **UniProt accession** everywhere, and symbols are display-only.
- Every external name is joined through `results/labels/kinase_synonyms.tsv`, and unknown names raise an error.
- Cantley matrices are joined by UniProt accession. Example: the matrix named "PDK1" is PDPK1.
- Pseudokinases are flagged; they are kept in the list but excluded from evaluated sets unless they have ≥3 labels.
- Regulatory-subunit rows (CK2β, AMPKβ/γ, cyclins) are **dropped**, not remapped.

**Residue class (revised).** The class is the Manning-group default (TK → Y; else → ST), with these curated overrides (`config/kinase_curation.tsv`):
- **ST:** LMTK1–3, STYK1. Also DYRK1A/1B/2/3/4, CLK1–4, HIPK1–4 and GSK3A/B: their Y labels are cis-autophosphorylation, so they are ST for trans substrates.
- **STY:** MAP2K1–7, TTK, TESK1/2, PKMYT1, PKDCC, NEK10.
- **Y:** WEE1/2.

Labels outside a kinase's class are dropped, with counts reported. ERN1 and MAP3K21 are reviewed in the same pass.

**Sites.**
- **U_obs:** 350,736 sites, the union of PSP sites, the compendium, CPTAC (BCM v1.2) and Ochoa 2020.
- **U_ms** ⊂ U_obs: sites with ≥1 MS observation (PSP MS_LIT+MS_CST ≥1, or present in the compendium, CPTAC or Ochoa). E-all evaluation uses U_ms (§5).
- **Site ID:** `ACC_ResPos` on the Swiss-Prot canonical sequence (UniProt 2026_03), mapped by 15-mer.

## 3. Labels (implemented: `src/ksp/labels.py`, `cache/labels/ks_pairs.parquet`)
- **Default training labels L2:** human-native plus mouse/rat PSP pairs transferred to human (18,638 pairs; 401 kinases).
- **Ablations:** L1 (human only); L2_vitro (in vitro evidence); L2_noauto (no autophosphorylation); L2+ext (adds the filtered external curated pairs on DEV proteins).
- **Autophosphorylation** is identified by substrate accession = kinase accession, not by symbol. It is kept in training labels by default, and **excluded from all primary evaluation endpoints** (DL-3). It is reported as a separate stratum.
- **Test positives** are human-native PSP pairs only, within class and not autophosphorylation, whatever the training label variant. They are reported stratified by evidence class (in vivo only / in vitro only / both).
- **"Known positive" set KP(k)**, used for negative exclusion: every pair for *k* from any source and any species, plus any site whose 15-mer is identical to a labelled site of *k*. Sources are PSP all species, external curated, and UniProt. **No KP site is ever a negative** in training or evaluation (DL-2, STAT-3).
- **External curated test (secondary):** pairs that are not in PSP, from OmniPath curated sources and UniProt ECO:0000269, restricted to TEST proteins. Filters (DL-5):
  - drop near-copies of PSP, meaning the same protein with the PSP site renumbered within ±5 residues, or the same site assigned to a paralog of a PSP kinase;
  - single-kinase names only.
  - When the main model trains on L2, these pairs are never in training.
- **Temporal test (secondary, PI-4):** the PSP Kinase_Substrate_Dataset bundled in kinase-library 1.8.0 (dated 2021) versus the 2025 release. Human-native pairs **new since 2021** on TEST proteins form a "recent discoveries" stratum. A temporal retrain on 2021-only DEV labels is optional (deferred).

## 4. Splits (to implement: `src/ksp/splits.py`; frozen at G1)
1. **Substrate clusters.**
   - mmseqs2 all-vs-all search over all proteins carrying a U_obs site, keeping edges with ≥40% identity and ≥50% coverage.
   - **Connected components** of that graph (DL-9; not mmseqs greedy clusters), unioned with identical-15-mer edges.
   - If any component holds >3% of labels, its internal edges are re-thresholded at 50% identity, and this is reported.
2. **TEST (locked).** About 20% of substrate clusters, chosen by a **constrained optimiser** (STAT-10: integer-program or simulated annealing, seed 0). The objectives are:
   - each kinase with ≥5 human-native evaluable labels gets ≈20% (15–25%) of its positives in TEST;
   - TEST holds 19–21% of all labels;
   - residue-class balance.
3. **Sealing (AIFM-2).** TEST labels are written to `sealed/test_labels.parquet` with its SHA-256 recorded. Rules:
   - The DEV label table physically excludes TEST proteins.
   - The canonical label loader `ksp.labels.load_dev()` is the only allowed reader.
   - A test (`tests/test_no_test_access.py`) greps the code base for direct reads of `sealed/`.
   - Only `src/ksp/eval_test.py` reads `sealed/`, and only with a committed `FREEZE.json` listing every configuration to be scored. All are scored in **one batch**, and the results are appended to `results/test/LEDGER.md`.
4. **DEV CV.** 5-fold grouped CV over DEV substrate clusters with the same constraints (seed 0). All selection, tuning and ablation happens here.
5. **Leave-kinase-out.**
   - Kinase clusters are connected components at ≥50% kinase-domain identity (KinCore domains), merged with Manning family.
   - 5 kinase folds, with no family split across folds.
   - **T1:** seen kinase, unseen protein (the standard setting).
   - **T2:** unseen kinase family, seen proteins; scored on CV.
   - **T3:** unseen family and unseen protein; scored at the final TEST run.
   - Random-kinase LKO is reported for contrast only.
6. **Evaluable kinase sets** are frozen at G1 with the split (AIFM-6, STAT-6):
   - **K_eval:** kinases with ≥3 test positives on ≥2 distinct substrate clusters.
   - **K_any:** kinases with ≥1 test positive.
   - **K_cv:** kinases with ≥5 DEV positives on ≥2 clusters.

## 5. Negatives (PU framing)
For kinase *k*, candidate negatives are sites of a residue class *k* acts on, minus KP(k), minus the **paralog mask**.

The paralog mask (KB-1, ALG-2, STAT-2, DL-2) removes:
- sites labelled to any kinase in Par(k), defined as ≥60% kinase-domain identity OR the same Manning subfamily (or the same family where no subfamily exists; e.g. SRC family = SRC, YES1, FYN, FGR, LCK, HCK, BLK, LYN);
- sites within Hamming distance ≤2 (in the 15-mer) of a KP(k) site.

The negative sets are:
- **E-hard (primary negative set):** sites in the evaluated split with ≥1 non-autophosphorylation label to a kinase outside {k} ∪ Par(k).
- **E-all:** all U_ms sites in the evaluated split, of matching class, after the masks. Evaluation is **exhaustive** over this set, with no subsampling.
- **Training negatives** use the same masks. Per kinase, they are a mix of E-hard-type and E-all-type sites from training clusters. The ratio and mix are tuned on CV over a grid {1:20, 1:50, 1:200} × {25%, 50% hard}.

## 6. Endpoints (pre-registered; frozen in this file at G1)
### 6.1 Primary
**M1: macro mean over K_eval of per-kinase AUROC on E-hard.** Kinases are weighted equally, and autophosphorylation is excluded.

### 6.2 Co-primary
- **M2, site-centric (user question: "which kinase for this site?").** For every test site with ≥1 label, take the *filtered* rank of each true kinase among all kinases of matching class. Filtered means the site's other true kinases are removed from the ranking. Report **Hit@5** averaged within kinase, then macro over K_eval. Hit@1, Hit@10, MRR and family-level Hit@5 are secondary.
- **M3, discovery (user question: "which substrates for this kinase?").** Macro over K_eval of **recall of test positives within the kinase's top 1% of E-all**, using U_ms sites in the TEST split. This replaces E-all AUPRC, which is unstable at prevalence ~2e-4 (STAT-5). Partial AUROC (FPR ≤ 0.01) and AUPRC lift are secondary.

### 6.3 Understudied-kinase endpoints (STAT-1)
- **M4:** macro AUROC (E-hard) over K_any, the ~290 kinases with ≥1 test positive. For a kinase with one positive, this is the positive's percentile.
- **M5:** CV-based macro AUROC (E-hard) over K_cv on DEV out-of-fold predictions. This includes kinases with 5–19 labels.
- Both are stratified by training-label bins (0 / 1–4 / 5–19 / 20–99 / ≥100).

### 6.4 Strata (all endpoints)
- kinase group;
- residue class;
- evidence class;
- priming-dependent kinases;
- PIKK;
- new since 2021 (temporal);
- external-curated test.

### 6.5 Floors and controls reported on every endpoint (STAT-4, PI-9, DL-4)
- **B0-prior:** the kinase frequency prior, i.e. the kinase's training count.
- **B0-pop:** a site popularity score (LT_LIT, MS counts, compendium observation count). This is diagnostic only and never a feature.
- **B0-miss:** a model trained only on missingness patterns of the phospho features.

### 6.6 Coverage (Priority 2; STAT-11, ALG-1)
Probabilities are **prior-corrected** (case–control offset log ρ_k for the training sampling ratio vs. U_ms prevalence). Calibration is tiered (§8.5).

Reported coverage measures:
- the number of U_obs sites, and of sites with no PSP kinase, having ≥1 kinase with corrected probability ≥ 0.5 and ≥ 0.8;
- the number of sites having ≥1 kinase in that kinase's top-1% percentile.

These are explicitly labelled as **model-based estimates**. Their calibration is checked on the external-curated and temporal sets. Precision on unannotated sites cannot be measured directly from PSP.

### 6.7 Downstream (secondary; KB-2, DL-6)
Kinase activity inference with benchmarKIN on the perturbation (Hernandez-Armenta) and CPTAC tumour benchmarks.
- Libraries use the top-N predicted substrates per kinase, with N ∈ {PSP size, 50, 100} fixed a priori. They are compared against PSP-only, Cantley-only and a shuffled-library control.
- Evaluation uses a **clean model variant**: its phospho features exclude every compendium study that is in the benchmark (matched by PMID or PXD) and all CPTAC-derived features (for the tumour benchmark).

### 6.8 Statistics
- 95% CIs come from a cluster bootstrap over substrate clusters (1,000 resamples); kinases are resampled for macro endpoints in a sensitivity analysis.
- Model comparisons are paired, on identical bootstrap draws.
- **Superiority claims** (PRIOR-8) are pre-declared for 3 comparisons on M1, M2 and M3: final vs B1 Cantley (best arm), vs B3 Phosformer-retrained, and vs B4 SELPHI2.0-retrained. A claim is "beats" if the paired 95% CI of the difference is >0 after Holm correction across the 3 comparisons × 3 endpoints.
- **Matched-universe rule (PRIOR-1):** every head-to-head is computed on the (k,s) pairs that both methods score. Coverage is reported separately.
- **CV model selection** uses M1-CV. To bound optimism, the number of configurations tried is logged in `results/cv/REGISTRY.tsv`.

## 7. Features
Each feature family has a provenance entry in `docs/feature_provenance.md`: source, whether it is label-derived, leakage risk, coverage, and expected value range. Features are stored as kinase × site (or kinase × protein) float16 matrices, indexed by the canonical orders in `cache/index/{kinases,sites,proteins}.parquet` (written once by the core module).

**Cross-fitting rule for label-derived features (critical; DL-1, PI-2, AIFM-4, ALG-6).**
- In CV fold *f*, label-derived features for validation and test rows use all fold-*f* training labels.
- For **training rows**, they are computed with **inner 5-fold cross-fitting over training substrate clusters**: the features of a site in inner fold *i* use labels from the other inner folds only.
- A site's own protein and substrate cluster is therefore never used to compute its own feature.
- The same rule applies in the final model (DEV → TEST).
- Canary kinases (§11) test this rule end-to-end.

### 7.1 Feature families
**F0. Site-intrinsic (not label-derived):**
- residue;
- AlphaFold pLDDT, smoothed pLDDT, exposure and IDR flag;
- ESM-2 650M per-residue embedding (full-protein context; PCA-64 for GBDT, full for R3).

Site-only features enter pair models only through interaction with kinase features. A "site-only" ablation measures how much they act as a popularity proxy.

**F1. Published motif (not label-derived):** Cantley score and **label-free percentile**. The percentile is recomputed against all U_ms sites, or against the Ochoa reference as shipped; the PSP-based reference is not used (DL-10).
- **Priming-aware variant:** flanking residues that are U_ms phosphosites are scored as pS/pT/pY, kept as a *soft* weight that uses the site's observation evidence (KB-10).
- For kinases without an atlas matrix, the nearest atlas kinase is used, with a distance feature; this is not part of the "as-published" baseline.

**F2. Learned motif (label-derived, cross-fitted):**
- (a) Per-kinase position-specific log-odds with an empirical-Bayes prior from family and group matrices.
- (b) R3 neural score, nested out-of-fold (§8).

**F3. PPI (systematic screens only):** BioPlex 3.0, HuRI (excluding Lit-BM), OpenCell, Buljan 2020 kinase AP-MS, Human Cell Map.
- Direct-edge flags per source.
- Random-walk-with-restart (RWR) proximity of kinase to substrate protein.
- **Label-derived, cross-fitted:** the maximum and mean RWR proximity of the substrate protein to the kinase's training substrate proteins, excluding the same protein and cluster.
- **Excluded:** STRING (all channels), BioGRID biochemical activity and PTM, IntAct enzymatic, curated pathway databases, ELM instances, GO BP/MF, COMPARTMENTS text mining.

**F4. Localisation:** cosine similarity of kinase and substrate compartment vectors from HPA IF v25, OpenCell, Human Cell Map, COMPARTMENTS (experiments channel; knowledge channel only in an ablation, because it imports GO/UniProt curation) and Hein 2025.
- Missingness flags.
- Membrane topology compatibility (KB-8): a substrate site that is extracellular or luminal pairs with a FAM20C-like secretory kinase. UniProt topology, if available.

**F5. Co-expression / co-dependency (low priority):** GTEx/HPA RNA, CPTAC protein, DepMap co-essentiality. Only if time permits by G2.

**F6. Phosphoproteomics (study-aware):**
- **(a) Co-regulation with training substrates (label-derived, cross-fitted).**
  - Correlate the site with *k*'s training substrates, per data block (compendium studies; CPTAC cohorts after regressing out protein abundance, purity, plex and other confounders per KB-9).
  - Exclude substrates on the same protein, same cluster, or same peptide.
  - Aggregate (top-m mean Fisher-z, median), with a minimum of 8 shared conditions over ≥2 studies; otherwise NA.
- **(b) Inhibitor response (not label-derived).**
  - Hijazi conditions (60 min) plus decryptM conditions, mapped to drug targets using the provided Hijazi matrix, Reinecke 2024, CurveCurator kinobeads and Davis 2011.
  - Per site, a **target-deconvolution regression**: the site's log2FC across inhibitor conditions regressed on the drug × kinase inhibition matrix (ridge, non-negative for decreases). Increases and decreases are kept separate.
  - Kinases whose inhibitor profile is not separable from another kinase's are flagged, and the flag is passed as a feature (KB-5).
- **(c) CPTAC kinase-state correlation (not label-derived):** correlation of the site with *k*'s protein abundance and *k*'s activation-loop phosphosites (from UniProt / known A-loop positions), per cohort after confounder regression.
- **(d) Kinase-activity-profile correlation (label-derived, cross-fitted):** low priority, because it overlaps with (a).
- **Missingness:** no count features; only NA, which LightGBM handles natively. B0-miss quantifies the residual signal from missingness.

**F7. Kinase descriptors:**
- Manning group and family;
- ESM-2 kinase-domain embedding (PCA);
- Cantley PSSM vector (PCA; S/T and Y separately);
- KinCore-aligned specificity-determining residues (optional).

**F8. Auxiliary high-throughput:**
- Sugiyama 2019 / iKiP in vitro lysate hits (kinase-specific flag plus score);
- Buljan 2020 S6 predicted sites are **excluded**, because they are a prediction.

This family is ablated.

**Excluded as features:**
- Ochoa functional score; PSP Regulatory_sites; LT/MS counts;
- CoPheeMap/CoPheeKSA;
- KinaseCanvas (contains curated edges);
- any third-party predictor (NetworKIN, GPS, KolossuS, Phosformer, KinoPlex percentiles).

## 8. Models
### 8.1 Baselines
**B0 floors:** prior, popularity and missingness (§6.5).

**B1 Cantley**, in three arms (PRIOR-2):
- (i) as published: atlas percentile, atlas kinases only;
- (ii) universe-recalibrated percentile;
- (iii) priming-aware.

**B2:** per-kinase learned log-odds PSSM (F2a alone).

**B3 Phosformer**, faithful reimplementation of paper2.pdf (PRIOR-3, CMP-2):
- 6 layers, 12 heads, d=768;
- MLM pretraining on kinase domains as in the paper (about 6 GPU-h; GPU 0);
- fine-tuning with focal loss (γ=2) at lr 2e-5, using the paper's input format (kinase domain plus 11-mer, ±5), on our folds;
- negative ratio re-tuned on our CV within {1:16 as in the paper, 1:50};
- plus the released weights as-published, inference only, with a caveat that its training data overlaps our test set.

**B4 SELPHI2.0** retrained (PRIOR-4):
- 37/43 features regenerated from our data; coreg/GTEx substituted by our F6a/F5 equivalents, cross-fitted;
- its random forest and hyperparameter grid, on our folds;
- every deviation documented.

**B5 NetworKIN** (benchmarKIN file), as published, with coverage reported. It is demoted: matched-universe comparison only.

**B6 KolossuS-small** (ESM-2 650M head), as published, if it runs within 1 GPU-h. Its weights may have been fine-tuned on human PSP (PRIOR-5), so it is reported as a contaminated comparator.

**B7:** Cantley + learned PSSM + minimal context (logistic regression), a "simple strong" baseline (PRIOR-6).

### 8.2 R2 global pair model
LightGBM over (k,s) rows with F0–F8 and F7, using binary loss. Variants:
- a lambdarank objective grouped by kinase;
- a **site-grouped** objective grouped by site, with candidate kinases = all kinases of matching class (PI-5), which targets M2;
- kinase sample weights ∝ n_k^(−α), α ∈ {0, 0.5, 1}.

Tree-size search is done on one fold first. Jobs run at ≤16 threads, 4 in parallel (CMP).

### 8.3 Partial pooling (ALG-3; replaces R1/R4 as separate regimes)
The final score is the global R2 logit plus a per-kinase residual model (a logistic regression on a small feature set), shrunk with strength λ:
- λ = ∞ gives the global model only;
- λ = 0 gives per-kinase models.

λ is chosen on CV per training-count bin. The result "where per-kinase models start to pay off" is reported.

### 8.4 R3 neural model (ALG-4/5)
- **Site encoder:** frozen ESM-2 650M per-residue embeddings (full-protein context) for the ±7 window, plus a phospho-channel bit marking U_ms phosphosites in the window, for priming. The encoder is a small transformer or CNN.
- **Kinase encoder:** ESM-2 kinase-domain embedding, Cantley vector and family one-hot, as blocks with **descriptor dropout p=0.2**.
- **Head:** bilinear plus MLP.
- **Loss:** BCE with masked hard negatives, plus a **listwise InfoNCE over kinases per positive site**, filtered.
- **Training:** GPUs 0–1, ≤1 GPU-h per fold.
- **Stacking into R2:** nested out-of-fold (inner folds for training rows), the same as the cross-fitting rule.

### 8.5 Calibration (ALG-7)
Prior correction (§6.6), then calibration in tiers:
- kinases with ≥30 DEV positives: per-kinase Platt scaling on CV out-of-fold predictions;
- otherwise: family-level, then global Platt.

### 8.6 Final and deployment models
- **Final:** the best CV configuration retrained on all of DEV, then the single TEST batch (G4).
- **Deployment:** retrained on DEV+TEST, with the same configuration, for the released predictions. It is clearly separated from the evaluated model.

## 9. Unorthodox kinases
**Priming-dependent** kinases: GSK3A/B (+4), CK1 family (−3), CK2 (acidic or primed), and PLK1 (PBD, S-pT-P). They get:
- F1 priming-aware scores;
- explicit priming-offset features: the residue is S/T, its U_ms status, and its predicted CK2/CDK/MAPK score;
- an R3 phospho channel.

Atypical kinases (PIKKs S/T-Q, FAM20C S-x-E with secretory localisation, PDHK/BCKDK mitochondrial, alpha-kinases) rely on localisation and topology features. Strata are reported separately; in practice the "atypical" strata are PIKK and FAM20C (KB-11).

Docking motifs are deferred (KB-4: the naive D-motif regex matches 76% of proteins). If done, they are restricted to disordered and conserved regions, with distance-to-site encoded.

## 10. Compute (measured in r1: `docs/reviews/r1_algorithms_compute.md`)
| Step | Measured cost |
|---|---|
| ESM-2 650M, full proteome | 8–15 GPU-min |
| Co-regulation, full matrix | ~4 s on GPU (masked matmul) |
| LightGBM, 400k×160, 1000 trees | ~15 min at 64 threads |
| LightGBM inference | ~1M rows/s; 129M pairs ≈ 2 min |

The bottleneck is materialising pair features: 129M×150 ≈ 78 GB in float32. Pair features are therefore gathered on the fly per kinase chunk from the kinase × site matrices.

GPU allocation:
- GPU 0: Phosformer B3;
- GPU 1: R3;
- GPUs 2–3: embeddings, then free for R3 folds.

## 11. Verification, controls and process (AIFM-1…11)
1. **Shared modules, mandatory:**
   - `ksp.ids`, `ksp.kinases`, `ksp.labels` (`load_dev`), `ksp.splits`;
   - `ksp.eval`: the only metric code. It takes a predictions parquet (kinase_acc, site_id, score) and returns all endpoints.
   - Unknown IDs raise errors. No `except: pass` and no silent fallbacks: every fallback is counted and logged in the run manifest.
2. **Run manifests:** every training or evaluation run writes `results/runs/<run_id>/manifest.json` with the git hash, config, seed, input checksums, runtime, n_rows and fallbacks, plus `predictions.parquet` (cache) and `metrics.json` (committed). The PI recomputes headline metrics from predictions with `ksp.eval`; agent-reported numbers are not trusted on their own.
3. **Canaries (AIFM-3):**
   - 12 **negative canary kinases**: pseudo-kinases whose labels are random U_ms sites (size-matched to real kinases at 10–100 labels; random descriptors). They run through the full pipeline, and the mean CV AUROC must be ≤0.57.
   - 1 **positive canary**: a planted motif plus planted PPI. Its CV AUROC must be ≥0.9.
   - Canaries are excluded from all reported metrics.
4. **Toy tests:**
   - metrics on hand-computed examples;
   - splits (no cluster, protein or identical-window overlap; constraints met);
   - cross-fitting (planted-signal toy where a leak gives AUROC 1);
   - negative masks.
5. **Biological controls (KB-7):** a frozen panel of about 40 canonical pairs (`config/bio_controls.tsv`, written at G1, not used in tuning). Also motif-fraction checks of top-ranked predictions: CDK P+1, ATM/ATR Q+1, CK2 acidic +3, PKA R−3, MAPK P+1, GSK3 +4 S/T.
6. **Leakage alarms:**
   - screens-only vs +curated-PPI ablation;
   - a large CV-vs-TEST gap;
   - feature-importance audit;
   - comparison of the B0-pop / B0-miss floors.
7. **Independent review:** every key result is reviewed by an agent that did not produce it, before it is marked [x]/[!] in RESEARCH_TREE.
8. **Process:**
   - Only the PI commits; agents never run `git`, because NFS plus concurrency breaks it.
   - Agents get explicit GPU/core budgets, and work and write only in their assigned paths.
   - `docs/DELEGATION.md` registers who owns which files.
   - Time pilots run before any job over 30 min.
   - `scripts/run_all.sh` documents the end-to-end pipeline.

## 12. Experiment plan (ordered by gates)
**G1:**
- splits, seal, harness, canaries and toy tests;
- feature stores F1, F0, F3, F4, F6b/c, which are not label-derived;
- B0, B1.

**G2:**
- cross-fitted F2a, F3-label and F6a;
- R2 v1;
- B2, B7;
- ablation curve by feature family.

**G3:**
- R3 plus stacking;
- B3 Phosformer (pretraining starts at G1 on GPU 0);
- B4 SELPHI2.0;
- regime comparison (α, λ, objectives, negative mix);
- label variants;
- LKO T2.

**G4:** FREEZE, then the TEST batch (T1, T3, external, temporal).

**G5:** downstream benchmarKIN with the clean variant; deployment model; coverage; independent review.

**G6:** REPORT.md.

## 13. Risks acknowledged
- **Test power:** about 160 kinases are in K_eval. Understudied kinases are assessed mainly through M4/M5 (CV), which are noisier.
- **PU labels:** all precision estimates are model-based.
- **Shared Unix user:** the test seal is enforced by convention plus grep and canaries, not by permissions.
- **Possible overlap** between PSP in vivo labels and compendium studies (PSP has no PMIDs per KS row). Mitigated by an in-vitro-label stratum and the study-removal sensitivity analysis (DL-6).

## 14. Revision log
- v0 (14:00): initial draft.
- v1 (14:30): round-1 revisions (see §15).
- v2 (~14:30): round-2 revisions (§16). No round 3: round 2 found no design-level critical issues beyond those fixed here (decision: PI; reason: avoid rework while 6 implementers run, per LESSONS #6/#7).

## 15. Round-1 change map
**Accepted, with location:**

| Reviewer | Changes → section |
|---|---|
| PI | PI-1 → §0; PI-2 → §7; PI-3 → §6.2/§4.6; PI-4 → §3; PI-5 → §8.2; PI-6/7 → §7 F6b, F8, Sugiyama acquired; PI-8 → §8.1; PI-9 → §6.5 |
| STAT | STAT-1 → §6.3; STAT-2/3 → §5/§3; STAT-4 → §6.5; STAT-5 → §6.2 M3; STAT-6 → §4.6; STAT-7–9 → §6.8; STAT-10 → §4.2; STAT-11 → §6.6; STAT-12 → §6.2 M2 |
| DL | DL-1 → §7; DL-2 → §3/§5; DL-3 → §3; DL-4 → §2/§5/§6.5; DL-5 → §3; DL-6 → §6.7/§13; DL-7 → §2; DL-8 → §2; DL-9 → §4.1; DL-10 → §7 F1; DL-11 → mapping already done (98.5%); DL-12 → §7 exclusions |
| ALG | ALG-1 → §6.6; ALG-2 → §5; ALG-3 → §8.3; ALG-4/5 → §8.4; ALG-6 → §7/§8.4; ALG-7 → §8.5; CMP-1–4 → §10/§8.1 |
| AIFM | AIFM-1…11 → §4.3/§11 |
| KB | KB-1 → §5; KB-2 → §6.7; KB-3 → §2; KB-4 → §9 (deferred); KB-5 → §7 F6b; KB-6 → covered by E-hard plus M2; KB-7 → §11.5; KB-8 → §7 F4; KB-9 → §7 F6a/c; KB-10 → §7 F1; KB-11 → §9; KB-12 → acquisition wave 2 done |
| PRIOR | PRIOR-1 → §6.8; PRIOR-2 → §8.1 B1; PRIOR-3 → B3; PRIOR-4 → B4; PRIOR-5 → B6; PRIOR-6 → B0/B7; PRIOR-7 → §4.6/§6.3; PRIOR-8 → §6.8 |

**Modified or rejected:**
- **Paralog identity threshold.** Reviewers proposed 60–80%. I chose ≥60% domain identity OR same subfamily, to be conservative for negatives.
- **KB-6 site-within-protein metric:** not a separate endpoint. E-hard plus M2 already isolate site-level specificity; it may be added as a secondary diagnostic if cheap.
- **External curated pairs as training augmentation** are kept only as the L2+ext ablation, so the external test stays independent.

## 16. v2 amendments (review round 2: docs/reviews/r2_*.md)
**Interfaces and process**
- **A1. Cross-fitting driver** (r2 PI-1, critical). There is exactly one driver, `src/ksp/crossfit.py`, owned by core-infra.
  - Every label-derived feature family is a `LabelDerivedFeature` subclass with `fit(train_pairs, ctx)` / `transform(query, ctx)`, registered with `@register`.
  - The driver raises `LeakageError` if any fitted label lies on a query protein or cluster.
  - **Eval rows** are transformed by a randomly assigned inner-fold fit, not by the all-labels fit (DL2-5). This makes their label coverage match the training rows (≈80%), so count-type features are not shifted between training and eval.
  - A grep test ensures labels enter feature code only through `fit()`.
- **A2. Feature store** (AIFM2-1, critical): every matrix is saved and loaded through `src/ksp/store.py`, which carries the index hash and raises on mismatch. Cached arrays are never indexed by position.
- **A3. Label sealing** (AIFM2-2, critical): at split freeze, every label table that contains TEST moves to `sealed/`. `ksp.labels.load_dev()` is the only reader for modelling, and a grep guard over a pattern list enforces it.
- **A4. Name resolution** (AIFM2-4): one strict resolver, `ksp.kinases.resolve_kinase_names`, with a per-source homonym table (e.g. "PDK1" = PDPK1 in HPRD/SIGNOR/Cantley).
- **A5. Canaries** (r2 PI-2): negative canaries are sampled clustered by protein (1–4 sites per protein) and injected before the split is frozen. The acceptance threshold is checked on both E-hard and E-all.

**Universe and labels**
- **A6. Candidate universe** (AIFM2-3, STAT2-3, r2 PI-4): positives and negatives for training and evaluation are both restricted to **U_ms**.
  - Non-U_ms positives and the 956 positives with no index row are dropped, with counts reported.
  - Centre-site observation flags are forbidden as features (DL2-7).
  - The deployment model still scores all of U_obs.
- **A7. Mask-only positives** (DL2-2): UniProt family-level notes ("by PKA/PKC"), complex-subunit rows and regulatory-subunit rows are added to KP as **mask-only** rows: never positives, never negatives. Text-mining and KEA rows are masked only in a sensitivity analysis.
- **A8. Paralog fallback** (STAT2-6): the "same family" fallback in Par(k) requires ≥50% domain identity.
- **A9. External test tiers** (DL2-3):
  - drop rows that list PhosphoSite among OmniPath sources;
  - cap 20 pairs per PMID;
  - tier A = SIGNOR/HPRD/UniProt-ECO with a single kinase and a named residue.
  - The expected size is small, about 24–32 evaluable kinases, so this set is secondary.
- **A10. Temporal test** (DL2-4): a pair is "new since 2021" if it is absent from the 2021 PSP file in all species. These pairs are reported with their in-vivo share, which is 40% against 23% for older pairs.

**Features**
- **A11. Forbidden sources** (DL2-1, DL2-8):
  - the iKiP-DB combined-with-PSP GMT, the iKiP activity and C-score tables;
  - COMPARTMENTS integrated;
  - the Hein literature tables;
  - HPA reliability as a feature;
  - Buljan publication counts.

  F8 = Sugiyama 2019 Table S2 only. Any feature source overlapping PSP KS pairs at more than 3× chance fails `tests/test_feature_overlap.py`.
- **A12. Label-source study circularity** (DL2-6): build the compendium-study ↔ PSP-label list before G2 and run the "F6-clean" sensitivity analysis (studies cited by PSP removed) at G2.

**Models**
- **A13. R3 integration** (r2 PI-3, DL2-5): nested stacking is replaced by **late fusion**. The R3 out-of-fold score is blended with R2 at score level, with the weight chosen on CV. Six trainings in total (5 CV folds + final), not 36.
- **A14. Cuts** (r2 PI-3): KolossuS (B6) is cut unless time allows, and so is the L2_noauto variant. The regime search is staged: tree size, then α, then the negative mix, then λ.
- **A15. Negative grid** (STAT2-2): where the hard-negative pool is smaller than requested, all available hard negatives are taken and the rest are filled with E-all-type negatives.

**Calibration and coverage**
- **A16. Calibration** (STAT2-2): calibrate on **exhaustive OOF predictions** over U_ms (no sampling, so no single prior-correction offset is needed). The fallback tier uses a family-level slope plus a per-kinase intercept.
- **A17. Coverage, redefined** (STAT2-1, critical).
  - **Site-level assignment probability:** P(k | s is a substrate of a kinase in K), a softmax over calibrated kinase scores at site *s*, with the temperature fitted on labelled DEV sites by OOF log-likelihood.
  - **Coverage@τ** = number of U_obs sites whose top kinase has assignment probability ≥ τ (τ ∈ {0.5, 0.8}).
  - Its reliability is checked on labelled TEST sites, external and temporal sets.
  - It is always reported next to a **random-score baseline** and a Cantley-only baseline.
  - It is conditional on the site being a target of one of the modelled kinases, and is stated as such.

**Endpoints and statistics**
- **A18. M2 paralog filtering** (STAT2-5): for true kinase *k* at site *s*, the ranking excludes Par(k) members that are not labelled to *s*. Family-level Hit@5 stays as a secondary metric.
- **A19. Superiority tests** (STAT2-4):
  - Cantley's "best arm" is chosen on CV, not on TEST.
  - The final superiority tests use B = 10,000 paired cluster-bootstrap resamples with Holm correction over the 9 tests.
  - K_eval counts are as computed at G1; about 135 is expected.
- **A20. Which score each endpoint uses** (r2 PI-6): M1, M3 and M4 use the per-kinase calibrated score. M2 uses the cross-kinase assignment probability (A17). Each component is selected on its own CV endpoint.
- **A21. Downstream additions** (r2 PI-7): leave-drug-out and "ours minus PSP" library arms.
- **A22. Phosformer** (AIFM2-8, r2 PI-8):
  - TEST proteins and clusters are removed from the pretraining corpus;
  - augmentations follow the paper;
  - negatives come from the shared `ksp.labels` generator;
  - if pretraining overruns, steps are reduced, with the change documented.

## 17. Post-freeze clarifications (decisions logged; they apply equally to all methods)
- **A23. Cantley comparator per endpoint** (2026-10-26 16:55; simple-baselines report). For superiority tests, B1 is the Cantley arm that is best on CV *for that endpoint*: arm (i) for M1 and M3, arm (ii) universe-recalibrated for M2.
  - This deviates from A19's single "best arm", in the direction that makes the baseline stronger.
  - Arms (i) and (ii) tie on M1 (difference 9e-7, from rounding), and arm (ii) is significantly better on M2.
- **A24. Pooled-OOF scale invariance.** Scores whose offsets vary by fold (per-kinase intercepts, count features) distort pooled-OOF per-kinase AUROC.
  - Example: B0-prior scored with fold-training counts gave M1 = 0.40 instead of 0.5.
  - **Rule:** CV metrics M1, M3, M4 and M5 are computed on scores rank-normalised **within fold and within kinase** (each kinase's percentile among that fold's candidates). M2 uses the cross-kinase score unchanged.
  - Per-fold mean M1 is reported alongside.
  - TEST has a single model, so the issue doesn't arise there.
- **A25. Canary descriptors.** Each canary gets its own Manning group, so group- or family-level priors cannot pool canary labels.
- **A26. Tyrosine-kinase analyses** (2026-10-26 17:15; `docs/features_tyr.md` §1.6).
  - (a) A **hub-site stratum** is ADDED as a secondary, reporting-only analysis: sites labelled to ≥5 non-paralog kinases (e.g. GLO1 Y136 with 15 TKs). It does not change any endpoint.
  - (b) Treating cross-kinase activation-loop labels as mask-only is **REJECTED** for all endpoints: it would change labels after seeing results. It may be reported as a clearly-labelled sensitivity analysis only.
- **A27. Audit follow-ups** (`docs/audit/audit_cv_r1.md`, 2026-10-26 17:25).
  - A24 is confirmed by the PI as a correctness fix applied to all methods equally. Documents show pre-A24 and A24 values where both exist. TEST is unaffected.
  - **Matched-kinase M2:** for head-to-head M2 comparisons (e.g. vs Cantley), the competitor kinase set is restricted to kinases scored by *both* methods. Unscored kinases are otherwise "worst", which inflates the M2 of low-coverage methods by about 0.02.
  - The canary gate remains on pooled CV; per-fold values are reported.
- **A28. Harness freeze** (2026-10-26 17:35). `src/ksp/eval.py`, `labels.py`, `splits.py`, `crossfit.py` and `canary.py` are frozen at the current commit; the harness version is that commit's hash.
  - Any further change needs a PI decision logged here, and must be applied to all methods.
  - Exploration on CV closes at **01:00**. After that, only the pre-declared final selection rule (§8.6: best M1-CV; ties → M2-CV) and the seed runs happen before FREEZE.json (G4).
- **A29. TEST-batch integrity and contents** (2026-10-26 20:45; independent review of `eval_test.py`).
  - FREEZE.json pins every predictions file by SHA-256. A STARTED ledger block is written before `sealed/` is read, and a COMPLETED block after. Any COMPLETED block blocks further TEST batches, except through a logged override. A crashed run may be resumed only with identical FREEZE, input and code hashes. The run requires `--real`, B = 10,000 and harness commit 974d520.
  - **Contingency:** if B3 or B4 is not ready at FREEZE, the primary family stays m = 9. Absent tests take p = 1 and are reported "NOT TESTED".
  - **In the single batch** (there is no second sealed access): per-kinase and site-level tables; kinase-resampling sensitivity CIs (§6.8); A17 assignment-reliability on labelled TEST sites; calibration on the external and temporal strata. Secondary analyses use B = 1,000.
  - Fold-0 rehearsal note: final vs B1 arm (i) M1 was Holm p 0.018 at B = 1,000 and 0.050 at B = 10,000. Small-B p-values are not trusted.
- **A30. B1 on TEST: single-residue-matrix kinases** (2026-10-26 21:50). Ten dual-specificity kinases have an atlas matrix for only one residue type: MAP2K1/2/4/5/6/7, TTK, TESK1, NEK10 and PKMYT1. On TEST they are left **unscored** in both B1 arms, so they drop out of matched head-to-heads instead of receiving a placeholder worst score.
  - On CV, the uncovered residue type was ranked last.
  - The change can only help the baseline, so it is conservative for superiority claims.
- **A31. Decisions from audit round 2** (2026-10-26 22:12; `docs/audit/audit_cv_r2.md`).
  - (a) **M2 score (clarifies A20).** M2 uses the raw fused score: CV 0.263 against 0.205 for the calibrated assignment probability, chosen on cross-fitted CV only (nested leave-one-fold-out: 0.259). The A17 assignment probability is used only for reliability analysis and coverage.
  - (b) **Assignment probability.** The calibrated softmax is overconfident in the ≥0.8 tier (0.87 claimed, 0.54 observed). Y sites at ≥0.8 are 0/41 correct. Fix (CV-only): isotonic recalibration of the top-kinase probability, separately for ST and Y sites, fitted on DEV OOF. Expected consequence: Coverage@0.8 is close to zero, and Y sites rarely reach 0.5. This is reported as is.
  - (c) **Canary gate.** It stays pooled, per A27 (pooled 0.539, pass). Per-fold values are reported: 2 of 5 folds exceed 0.57 (0.575 and 0.676; fold 4 is driven by R3).
  - (d) **LKO T2.** Reported stratified by whether the held-out site has other-kinase labels: M1 is 0.703 on sites without them (55%).
  - (e) **Clean variant disclosure.** Site existence (U_ms membership) and primed-Cantley evidence weights still draw on the compendium and CPTAC.
- **A32. Pre-declared rules for two late CV experiments** (2026-10-26 22:16; prompted by `docs/strata_cv.md`). These rules were declared before either experiment was run.
  - **(a) Priming features.** The finding: on primed positives the priming-aware Cantley arm scores 0.945, against 0.870 for R2.
    - The experiment adds label-free, kinase-specific priming features to the R2 final config. They are the priming-aware Cantley percentile, plus an indicator that the kinase's preferred phospho-offset residue (from its Cantley phospho-preference) is an observed phosphosite.
    - **Adopt only if all of the following hold** (pooled 5-fold CV, seed 0, against r2_final_s0):
      - overall M1 Δ ≥ −0.002;
      - M2 Δ ≥ −0.005;
      - primed-positives-stratum M1 Δ ≥ +0.02, with a cluster-bootstrap CI above 0;
      - canary pass.
  - **(b) Prior-debiased site score for M2.** The finding: Hit@5 is 0.099 for kinases with 5–19 labels, against 0.162 for Cantley.
    - The experiment computes score_site = fused score − β·log(DEV training count of k, per fold), with β ∈ {0, 0.1, 0.25, 0.5, 0.75, 1}. It is applied post hoc to the OOF fusion and needs no retraining.
    - **Selection:** β is chosen by pooled CV M2 (macro, A24-compatible), with nested leave-one-fold-out selection reported.
    - **Adopt β > 0 only if:**
      - M2 does not fall by more than 0.003 against β = 0;
      - the 5–19-label-bin M2 improves.
      This deliberately favours harmonious performance across kinases (PROMPT), at a bounded cost to the macro endpoint.
- **A33. Label variant and training regime: final decision** (2026-10-27 00:03; `docs/regimes_labels.md`). The default stays: **L2 labels with the global R2 model.**
  - **L2+ext** meets the A32-style rule on seed 0 only: ΔM1 +0.0058. Seeds 1 and 2 give −0.004 and +0.002, a mean of +0.001. Its M2 is positive on all three seeds (+0.026 / +0.018 / +0.002).
  - **Why L2+ext is not adopted:**
    1. DESIGN §16 (round-2 change) pre-registered external curated pairs as **ablation only**, so that the external-curated TEST stratum stays independent. That rule takes precedence.
    2. The M1 gain does not replicate across seeds.
    3. A rebuild of TEST and T3 (about 2.5 h) is not justified by a gain within noise.
  - **Other variants:** L1, L2_vitro and L2_vivo are all worse on M2 (pooling mouse and rat labels helps M2 by +0.029). Per-kinase models (−0.012), per-group models (−0.007) and partial pooling (+0.0002, with λ selected at 300, i.e. essentially global) are not adopted.
  - **Side finding, recorded as a follow-up (not explored, because CV closed at 01:00):** per-kinase ridge intercepts lift R2 M2 to 0.232.
