# LITERATURE — kinase–substrate prediction: prior work, data sources, modalities

This file is the synthesis; full reviews with citations (DOIs/URLs) are in `docs/lit/`:
- `docs/lit/01_methods_benchmarks.md` — prior methods, benchmarks, evaluation flaws, how to retrain baselines
- `docs/lit/02_ks_annotation_sources.md` — sources of experimentally verified KS annotations and their quality
- `docs/lit/03_phosphoproteomics_sources.md` — phosphoproteomics resources (provided compendium, CPTAC, inhibitors, Kinobeads)
- `docs/lit/04_other_modalities.md` — PPI, localisation, co-expression, structure, PLMs, kinase descriptors
Supporting audit scripts: `scripts/phosformer_split_leakage_check.py`, `scripts/lit/ks_source_overlap_stats.py`, `scripts/audit_ppi_ks_leakage.py` (results in `results/leakage/`).

## 1. Prior methods and what their numbers really mean
| Method | Signal | Model | Evaluation & caveat |
|---|---|---|---|
| Cantley Ser/Thr atlas (Johnson 2023 Nature, 10.1038/s41586-022-05575-3) & Tyr atlas (Yaron-Barir 2024 Nature, 10.1038/s41586-024-07407-y) | In vitro positional-scanning peptide arrays; 303 S/T + 93 Y PSSMs | Product-of-PSSM score; percentile vs Ochoa-2020 reference phosphoproteome | Not trained on PSP → usable as-is (no retraining needed). Handles priming natively (pS/pT/pY matrix rows; 129/303 S/T kinases favour a phospho-residue). Motif only: no context. |
| Phosformer (10.1093/bioinformatics/btad046) — paper2.pdf | Kinase-domain sequence + 11-mer | 6-layer RoBERTa, MLM-pretrained on kinase domains, fine-tuned | **Random split over pairs**; we measured 51% of test-positive peptides also occur as training positives under another kinase, 30% share (family, peptide). Reported AUROC 0.885 is an in-distribution upper bound. CC BY-NC-ND, no training script → reimplement to retrain. |
| Phosformer-ST (btae033) | as above, S/T | trained on Cantley percentiles | Distillation of the atlas; informative number is 0.886 AUROC on held-out kinases. Best published split design (CD-HIT-70% kinases + Hamming-5 peptides). |
| SELPHI2.0 (10.1016/j.mcpro.2025.100994, Petsalaki lab) — paper1.pdf | 45 features (only ~5 kinase-specific: 2 PSSM scores, co-expression, co-regulation) | Random forest | Stratified CV **over pairs**; co-regulation features correlate sites against KSEA activities from the same PSP substrate sets (leak). CV AUC 0.97 → 0.75–0.83 on external in vitro/inhibitor data. Apache-2.0; retrainable on CPU with fold-aware feature regeneration. |
| KolossuS (bioRxiv 2025.04.27.650849) | ESM-2 15B kinase + 15-mer | cosine head; trained on atlases + Sugiyama | Phylogeny-aware split, PSP evaluation: AUROC 0.755 S/T / 0.838 Y — the realistic level for sequence-only specificity on literature labels. |
| NetworKIN/KinomeXplorer, GPS 5/6, LinkPhinder, KinasePhos3, DeepKinZero, EMBER, PhosIDN, KSFinder, CoPheeKSA, KSTAR | motif ± network ± co-phosphorylation | various | Mostly random-pair splits on PSP-derived data; LinkPhinder 0.97 (random) → 0.53 (independent). GPS 6.0 ≈0.81 on SELPHI's external set. |

**Lessons for evaluation** (01 §4): random pair splits leak through shared peptides, paralogous substrates and orthologues; protein-level features (abundance, "studied-ness") let models learn annotation popularity instead of kinase specificity; pooled metrics are dominated by ~20 kinases (top 20 kinases = ~50% of PSP pairs). Recommended: substrate-protein clusters (mmseqs2) with orthologues/identical windows merged; three difficulty tiers (unseen protein; unseen kinase cluster; both); per-kinase macro metrics stratified by annotation count; per-site "rank of true kinase"; declared negative sets; feature-provenance audit.

## 2. Kinase–substrate annotations (labels)
- **PhosphoSitePlus** (provided, release 2025-05-19): 23,904 rows; 14,724 human→human; 421 human kinases. Evidence: ~40% in vitro only, 25% in vivo only, 35% both. In vivo-only is enriched for indirect relations (S/T sites "assigned" to tyrosine kinases are 69% in vivo-only vs 25% baseline). 893 rows use isoform accessions → map sites by 15-mer to canonical sequence, not by residue number. Mouse/rat pooling adds ~4.7k (kinase, site-group) pairs; mapping to human needs site groups (now possible: user provided `Phosphorylation_site_dataset`).
- **Additional curated sources (critically filtered):** SIGNOR 4 (direct, residue-specified, single-kinase relations; host unreachable → OmniPath copy), UniProt ECO:0000269 "by KINASE" (~960 pairs not in PSP; filter family names/autophosphorylation/orthology-propagated), OmniPath curated non-PSP rows (~5.8k; HPRD/Phospho.ELM/NCI-PID/Reactome; require agreement). Most other resources (EPSD, dbPTM, qPhos, GPS gold, iPTMnet) re-aggregate PSP/UniProt.
- **Not labels:** Sugiyama 2019 in vitro lysate screen (175k pairs, recovers ~4% PSP; Müller-Dott 2025 found adding it hurts activity inference) → feature only; PhosphoNetworks site assignments are motif predictions; text-mined (RLIMS-P/eFIP/REACH) excluded; predictor outputs (GPS/NetworKIN/PhosphoNET) excluded.
- **Negatives:** no true negatives → positive–unlabelled setting. Residue-matched negatives from observed phosphosites; easy (random observed sites) and hard (sites annotated to other kinases) sets.

## 3. Phosphoproteomics
- **Provided compendium:** 1,937 contrasts / 178 PMIDs / 189,962 sites / 19,233 proteins; only 4.6% dense. Contains decryptM (Zecha 2023 Science; 790 dose contrasts = 41% of columns) and Hijazi 2020 (183 contrasts); 461 contrasts from qPTM. Heavily study-structured (dose series share controls) → study-aware aggregation needed; only singly-phosphorylated peptides (primed sites under-represented).
- **CPTAC:** 10 cohorts, public (LinkedOmics S3 freeze v1.2; PDC pan-cancer harmonised). Co-phosphorylation across tumours is dominated by protein abundance → correlate after protein-abundance adjustment, per cohort.
- **Inhibitor → kinase maps:** provided Hijazi matrix; Klaeger 2017 Kinobeads (blocked, 403), Reinecke 2024, Davis 2011/DTC. Bayer…Kuster bioRxiv 2025 (17M dose–response curves, 133 inhibitors) — not in PRIDE (request to user).
- **Leakage:** CoPheeMap network trained on PSP; SELPHI2.0 coreg; Ochoa 2020 functional score trained on PSP functional sites; PSP in vivo labels may derive from compendium studies. Rule: any feature derived from known substrates is recomputed per fold from training substrates only (leave-site-out).
- **Downstream benchmark:** benchmarKIN (Müller-Dott 2025 Nat Commun, 10.1038/s41467-025-59779-y; Zenodo 12566560): kinase-activity inference on perturbation (Hernandez-Armenta) and tumour (CPTAC) benchmarks; showed library choice matters more than the inference method.

## 4. Other modalities
- **PPI:** STRING "experiments/physical" leaks (34% of PSP KS pairs have a score vs ~4% of controls; imports BioGRID "Biochemical Activity" = kinase assays). Use only systematic screens: BioPlex 3.0, HuRI (not Lit-BM), OpenCell, Buljan 2020 kinase AP-MS (300+ kinases), Human Cell Map. These are 7–8× enriched for KS pairs but cover 0.3–1.5% of pairs → precision, not coverage. Network propagation (RWR) increases coverage.
- **Localisation:** HPA v25 IF, OpenCell, Human Cell Map, COMPARTMENTS knowledge+experiments (not text mining). Missingness flags.
- **Co-expression/co-dependency:** GTEx/HPA RNA, CPTAC co-abundance, ProteomeHD, DepMap 26Q1 co-essentiality — leak-free, modest value.
- **Site structure/sequence:** AlphaFold DB pLDDT/RSA/disorder (structuremap); ESM-2 650M/3B, ESM-C 300M/600M embeddings (MIT); docking motifs from ELM *class regexes* (MAPK D-motif, RxL, PIF, polo-box custom); priming represented as (residue present at priming position, observed phosphorylated, predicted priming-kinase score).
- **Kinase side (zero-shot / sharing):** Cantley PSSMs as descriptors, kinase-domain PLM embeddings, KinCore/Manning group/family.
- **Structure-based KS (AF-Multimer):** infeasible at scale on 4×4090 (≥9,000 GPU-h for 400 kinases × 10k sites); KinoPlex (Nat Biotech 2026) as optional comparator.

## 5. Implications carried into DESIGN.md
1. Protein-cluster-grouped splits with orthologue/identical-window merging; hold-out test set touched once; leave-kinase(-family)-out tier.
2. Kinase-equal-weight metrics stratified by training-label count; per-site kinase-ranking metric.
3. Positive–unlabelled framing with two pre-declared negative sets; residue-type matching.
4. Fold-aware computation of every label-derived feature; leakage-prone resources excluded; provenance audit + "curated vs screen" ablation as a leakage alarm.
5. Baselines: Cantley as-published; Phosformer re-implemented & retrained; SELPHI2.0 retrained with fold-aware features; per-kinase PSSM; NetworKIN/KolossuS as-published where feasible.
6. Priming-aware motif features as a differentiator.
