> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# 03 — Large-scale phosphoproteomics data sources for kinase–substrate (KS) prediction

Author: literature/data-review agent. Date: 2026-10-26. Status: review, not yet acted on.
Scope: phosphoproteomics datasets that can be turned into features for site-level KS prediction (co-regulation, perturbation response, inhibitor-target-informed scores), how to access them, and how to use them without label leakage.
Everything marked **[verified]** was checked during this review (API call, HTTP HEAD, file listing or full-text read). **[unverified]** items come from memory or secondary sources and need checking before anyone relies on them.

---

## 0. TL;DR

1. The user compendium (`data/01-03`) is a harmonised log2 fold-change matrix: **1,937 contrasts from 178 PMIDs, 189,962 sites on 19,233 proteins, about 17.0 M non-missing values (about 4.6 % dense)**. It already contains **all of decryptM (790 contrasts)** and **Hijazi 2020 (183 contrasts)**. It is the main perturbation feature source. Its weaknesses are the extreme sparsity and the strong study/batch structure.
2. **CPTAC pan-cancer phosphoproteomes can be downloaded without authentication** **[verified]**. They come from the LinkedOmics public S3 mirror of the BCM-harmonised freeze v1.2: 10 tumour cohorts, about 30–65 MB per cohort, about 500 MB in total, plus normal-tissue files for most cohorts. This is the best steady-state co-regulation source. It is complementary to the compendium because it has about 1,100 tumours and real between-sample variance.
3. Mapping inhibitors to kinases is the critical side-table for perturbation features. We already hold the Hijazi matrix (`kinhib_inhib_comb__1.csv`, 420 kinases × 60 inhibitors). Klaeger 2017 Kinobeads (243 clinical drugs) and Reinecke 2024 (1,183 tool compounds) supply Kd_app values for the decryptM drugs and many more. Davis 2011 and Drug Target Commons add a second, recombinant-assay view.
4. The strongest new direct-KS dataset is **Bayer … Kuster, bioRxiv 2025** ("potency-coherence"). It has 17 M dose–response curves for 133 inhibitors with known affinities in 5 cell lines, and claims thousands of new KSRs for 96 kinases. The data do not appear to be in PRIDE yet (only ProteomicsDB). **Request this from the user or the authors.**
5. The benchmark/label-adjacent sets (Hernandez-Armenta 2017, Müller-Dott 2025 benchmarKIN, Sugiyama 2019 in vitro, Hijazi KS) are all openly downloadable. Each one must be classified up front as either **feature source** or **evaluation set**, never both.
6. **Leakage is the dominant risk.** Any feature computed from "known substrates of kinase K" falls in this class: kinase-activity profiles, PSSMs, correlation with the substrate centroid, CoPheeMap embeddings, PTMsigDB signatures, NetworKIN and the SELPHI2 outputs. Such features must be recomputed inside each CV fold from training-fold substrates only. Precomputed third-party versions cannot be used as features in evaluation.

---

## 1. The user-provided compendium (`data/`)

### 1.1 What is in it [verified by inspection; 03 was only line-counted]

| File | Content | Shape |
|---|---|---|
| `01_phosphoproteomics_conditions.csv` | per-contrast metadata: conditionID, Condition, Control, Timepoint, Sample (cell line/tissue), Labeling, Spectrometer, Enrichment, ID/quant software, PMID, `source` | 1,937 rows × 13 |
| `02_phosphosite_information.csv` | `peptide_UID` = `<UniProtAcc>_<pos>_M1` (all multiplicity 1), 15-mer `sequence_window` (centre = residue 8) | 189,962 sites (S 133,604 / T 42,354 / Y 14,004) on 19,233 accessions (3,996 isoform accessions, e.g. `P12345-2`) |
| `03_phosphosite_quantifications.csv` | long format `peptide_UID, conditionID, log2FC` | 17,032,346 values → **≈4.6 % of the 189,962 × 1,937 matrix observed** |
| `kinhib_inhib_comb__1.csv` | binary kinase × inhibitor target matrix for the Hijazi 2020 screen | 420 kinases × 60 inhibitors; 307 kinases hit by ≥1 inhibitor; per inhibitor from 0 (U73122, a PLC inhibitor) to 98 (PF3758309) kinases |

`source` column: **"This study" = 1,437 contrasts / 27 PMIDs** (re-processed by the compendium authors); **qPTM = 461 contrasts / 144 PMIDs**; **Ochoa2016 = 39 contrasts / 7 PMIDs**.

Largest contributing studies:

| PMID | Study | # contrasts | Notes |
|---|---|---|---|
| 36926954 | **decryptM**, Zecha et al. 2023 *Science* (doi:10.1126/science.ade3925) | 790 | dose series (e.g. Dasatinib 0.01 nM → µM) for 31 drugs in 13 lines. Includes kinase inhibitors (dasatinib, lapatinib, imatinib, gefitinib, afatinib, PD325901, selumetinib, MK2206, AZD8055, dactolisib, pictilisib, staurosporine, nintedanib, AZD4547, …) but also rituximab (306), trastuzumab/pertuzumab, proteasome, HDAC and KAT inhibitors and chemotherapies |
| 31959955 | **Hijazi et al. 2020** *Nat Biotechnol* (doi:10.1038/s41587-019-0391-9) | 183 | 61 inhibitors × HL60/MCF7/NTERA2, single dose |
| 32034161 | EGF ± MEK inhibitors (cobimetinib, PD0325901) | 66 | |
| 34663829, 37188730 | TNF time courses (15 s → 4 h) | 65, 24 | high temporal resolution, so useful for kinase-cascade ordering |
| 38047086 | BRAFi in CRC models | 42 | |
| 26060313 | (Ochoa2016 source) panel of kinase inhibitors (AKTi, CAMK2i, EGFRi …) | 22 | |
| others | ~170 studies with 1–36 contrasts: growth factors, stresses, cell cycle, infections, a few genetic perturbations (e.g. PLK1-as + 3-MB-PP1, 22199227; MARK3-KO; AMPK-DKO; FAM20C-KO; TBK1 shRNA; PKD1/2/3 KD; RSK1/2 RNAi + BI-D1870; IKKB K44M) | | |

### 1.2 What it most likely is
- It **extends the Ochoa et al. 2016 *Mol Syst Biol* "atlas of human kinase regulation" compendium** (PMID 27974457; ~435 conditions). It does so with the **qPTM** database (Yu et al. 2023 *NAR* 51(D1); 2,596 conditions across 4 species, of which the human phospho subset was taken) and with 27 recent studies re-processed "in house", among them decryptM and Hijazi.
- Ochoa 2020's compendium (6,801 runs, PXD012174) is an identification resource, not a fold-change compendium. Invergo et al. 2020 *Cell Syst* (signed kinase circuits) used the Ochoa 2016 fold-change compendium. The presence of qPTM and the Petsalaki authorship of SELPHI2.0 (paper1.pdf, Maier et al. 2025 MCP, doi:10.1016/j.mcpro.2025.100994) suggest that **this compendium is a Petsalaki-lab extension**. SELPHI2.0 itself used only Mertins (CPTAC BRCA) and Hijazi for co-regulation and fitted the Ochoa 436-condition compilation post hoc. So training on the full 1,937-contrast compendium is a modality the SELPHI2.0 baseline did not use. **[inference; ask the user for the provenance/citation and the processing pipeline]**

### 1.3 Value, identifier mapping, pitfalls
- **Value: high.** It is the only source here with many perturbations aimed at specific kinases, including dose information (decryptM). It supports both (a) co-regulation with training substrates and (b) "site responds to inhibitors of K" features.
- **Mapping:** UniProt accession + position. It must be reconciled with the PSP `SUB_ACC_ID`/`SUB_MOD_RSD` and the PSP site groups (`SITE_GRP_ID`). Isoform accessions (≈2 %) need sequence-window matching against canonical sequences. Mapping should be done on the **15-mer window**, not on position, because positions drift between UniProt releases.
- **Pitfalls.**
  - Rows are 95 % missing, and missingness is not random: low-abundance sites, Y sites only in pY-enriched studies, and study-specific detection.
  - Strong batch/study structure. Contrasts within a study share a control. decryptM alone is 41 % of the columns and rituximab is 16 %, so naive correlations will be dominated by B-cell rituximab dose series.
  - Fold changes come from different pipelines (SILAC/TMT/LFQ). TMT ratio compression shrinks effects.
  - qPTM-derived contrasts carry the depositors' processing (often only significant hits reported, so a zero/missing ambiguity).
  - Multiplicity M1 only: multiply phosphorylated peptides were dropped, so primed sites (GSK3/CK1 substrates) are under-represented.

---

## 2. CPTAC pan-cancer phosphoproteomics

### 2.1 Sources and access
| Route | What | Auth | Notes |
|---|---|---|---|
| **LinkedOmics public S3 mirror** of the BCM-harmonised freeze v1.2 **[verified]** — e.g. `https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/{COHORT}/{COHORT}_phospho_site_abundance_log2_reference_intensity_normalized_{Tumor,Normal}.txt` (index page: https://www.linkedomics.org/data_download/CPTAC-pancan-OV/ etc.) | site-level log2 MS1/reference-normalised abundance; matching `_proteomics_gene_abundance_…`, RNA, CNV, mutation, `_meta.txt` | **none** | Tumor file sizes [verified HEAD]: BRCA 50.6 MB, CCRCC 46.0, COAD 30.7, GBM 65.0, HNSCC 54.8, LSCC 62.0, LUAD 57.9, OV 29.1, PDAC 48.7, UCEC 57.2 MB (≈0.5 GB total). Normal files exist for most cohorts (OV 8.2 MB; BRCA none). The bucket listing is denied, so URLs must be constructed |
| PDC pan-cancer page https://pdc.cancer.gov/pdc/cptac-pancancer | `phosphoproteome_bcm_gencode_v34_harmonized_v1.zip`, `proteome_bcm_gencode_v34_harmonized_v1.zip`, `rna_bcm_v1.zip` (names from the CoPheeMap data statement) | open (JS portal) | same data as S3. The PDC also has per-study raw/CDAP outputs (PDC0000xx) |
| Zenodo 12566560 (Müller-Dott 2025, CC-BY-4.0) **[verified]** | median-centred, ≥30-sample-filtered phospho tables for BRCA, CCRCC, GBM, HNSCC, LSCC, LUAD, UCEC (24–52 MB each) + kinase "gold standard" | none | convenient and pre-filtered. Same underlying data |
| `cptac` python package (PayneLab, PyPI) | per-cohort dataframes from several pipelines (BCM, UMich, WashU, Broad) | none (downloads from Zenodo/Box) | handy, but adds an extra dependency layer. Prefer S3/Zenodo flat files |
| LinkedOmicsKB https://kb.linkedomics.org | precomputed associations for 126,547 phosphosites | none (web) | browsing only |
| HCC (Gao 2019) harmonised matrix | Zenodo 14553766 (cited by CoPheeMap) | none | 11th cohort **[unverified content]** |

Key papers: Li et al. 2023 *Cancer Cell* 41:1397 "Proteogenomic data and resources for pan-cancer analysis" (doi:10.1016/j.ccell.2023.06.009); Li et al. 2023 *Cell* "Pan-cancer proteogenomics connects oncogenic drivers to functional states" (doi:10.1016/j.cell.2023.07.014); Geffen et al. 2023 *Cell* "Pan-cancer analysis of post-translational modifications" **[DOI unverified]**.

### 2.2 Format and identifiers [verified from file head]
- Row ID: `ENSG…|ENSP…|S1053|PPTIRPNSPYDLCSR|1`, i.e. gene | protein isoform (GENCODE v34 / Ensembl 100) | residue+position | flanking sequence | multiplicity. Columns are case IDs. Values are log2 reference-normalised intensities. They are not ratios to a common reference within plex, which removes the TMT-plex structure but leaves cohort effects.
- **Mapping needed:** ENSP → UniProt (Ensembl 100 BioMart or UniProt ID-mapping). Sites on non-canonical ENSP need flanking-sequence matching to the UniProt canonical sequence. Plan for ~5–10 % loss **[estimate]**.
- Scale (CoPheeMap): 1,195 tumours across 11 cohorts; 158,796 sites in total; 77,442 quantified in ≥20 % of samples in at least one cohort.

### 2.3 Value and pitfalls
- **Value: high for co-regulation.** CoPheeMap (Ayati/Zhang lab, *Nat Commun* 2025, doi:10.1038/s41467-025-57993-2) shows that known same-kinase site pairs have higher across-tumour correlation. Its KS model built on CPTAC features reached AUROC 0.95 (under their own, permissive split).
- **Confound:** much of the co-phosphorylation signal is **host-protein co-expression**. Normalising site abundance by protein abundance "reduced the PCC differences between positive and negative pairs markedly" (CoPheeMap). Compute both raw and protein-residualised correlations, and keep the residualised one as the primary feature. The raw one should only be used if a fold-internal ablation shows it adds information beyond protein-protein co-expression (which a PPI/co-expression feature may already provide).
- Within-cohort correlation, followed by aggregation across cohorts (Fisher-z mean, or the count of cohorts with significant correlation), avoids tissue-of-origin confounding.
- Tumours carry little direct kinase-specific perturbation. The signal is steady-state variation, which favours kinases with high variance across tumours (RTKs, CDKs, MAPKs).
- License: CPTAC data are open access with citation required **[policy page not re-checked]**.

---

## 3. Kinase-inhibitor / kinase-perturbation datasets

### 3.1 Perturbation phosphoproteomics
| Dataset | Content | Access | License | In compendium? | KS value | Leakage notes |
|---|---|---|---|---|---|---|
| **Hijazi et al. 2020** *Nat Biotechnol* 38:493 | 61 kinase inhibitors (single dose) × HL60, MCF7, NTERA2; ~6,000 sites used to infer >1,500 kinase–kinase links; also derives KS associations | PRIDE **PXD015943** [verified, CC0]; processed fold-changes in Zenodo 12566560 `benchmark_data_hijazi.csv` (19 MB) + `benchmark_metadata_hijazi.csv` [verified] | CC0 / CC-BY | yes (183 contrasts) | high: inhibitor-target-informed scores | SELPHI2.0 uses Hijazi both as a co-regulation feature and as a "novel KS" test set, a circularity we must not copy. Decide: **feature source** (preferred) vs held-out test |
| **decryptM**, Zecha et al. 2023 *Science* 380:93 | 31 drugs, 13 lines, dose–response (≥10 doses), ~1.8 M curves; 47,502 regulated of 124,660 p-peptides | PRIDE **PXD037285** [verified, CC0]. Curve-fit zips per experiment (e.g. `10_Kinase_Inhibitors_Curves.zip` 1.0 GB, `3_EGFR_Inhibitors_Curves.zip` 1.1 GB, `BCRABL_Inhibitors_Curves.zip` 0.3 GB, `Dasatinib_Triplicates_Curves.zip` 0.3 GB) plus Kinobeads zips (AZD8055, dactolisib, pictilisib); ProteomicsDB | CC0 | yes as per-dose log2FC (790) | **very high:** per-site EC50 that can be compared with the drug's Kd_app on each candidate kinase ("potency coherence") | the compendium has the doses, but the curve fits (EC50, curve class) are only in the PRIDE zips. Download the `*_Curves.zip` files for the kinase-inhibitor experiments first (~3 GB) |
| **Bayer, Müller, …, Kuster, bioRxiv 2025.11.18.689017** "Chemical proteomics decrypts the kinases that shape the dynamic human phosphoproteome" | 133 kinase inhibitors with known targets/affinities × 5 cancer lines; 17 M peptidoform dose–response curves; potency-coherence analysis; "thousands of new KSRs for 96 kinases"; implemented in ProteomicsDB | not found in PRIDE (search of Kuster-lab projects, 2026-10) → probably ProteomicsDB only / embargoed | preprint CC-BY-ND | no | **highest potential:** quasi-direct KS evidence across 96 kinases | their KSR calls use motif plausibility (Kinase Library), so it is not independent of motif features. If used as labels or test, keep their motif filter in mind. Their re-evaluation of published KSRs is also a check on PSP label noise |
| Ochoa 2016 / qPTM-sourced inhibitor studies (e.g. PMID 26060313 inhibitor panel; 23532336 AML; 25002506 RSK; 22199227 PLK1-as) | small sets | in compendium | — | yes | medium | — |
| **LINCS P100** (Litichevskiy et al. 2018 *Cell Syst*) | ~96 targeted phosphopeptides × thousands of perturbagens incl. kinase inhibitors | LINCS / Panorama **[not re-checked]** | open | no | low (tiny site coverage) | — |
| **Frejno et al. 2020** *Nat Commun* (doi:10.1038/s41467-020-17336-9) | 125 cancer cell lines, ~55,000 sites, steady state; ATLANTiC portal | PRIDE **[accession not verified]** | CC-BY | no | medium: another cross-sample co-regulation matrix (cell lines) | — |
| **Koenig … Olsen, bioRxiv 2025.12.05.692571** "Human Cell Line Phosphoproteome Atlas" | 33 cell lines, nDIA, >200,000 class-I sites; spectral library; kinase activity scores | PRIDE **PXD068583** [project exists] | CC-BY-NC | no | medium-high for **site universe/coverage** (Priority 2); low-medium for co-regulation (n=33) | — |
| Genetic perturbations (KO/KD/analog-sensitive) | scattered single-kinase studies (PLK1-as, CK2 KO in C2C12, FAM20C KO, AMPK DKO, MARK3 KO, ERK2/PLK1/PAK4 KO in HCT116 2026 …) | mostly in compendium or PRIDE | various | partly | high per kinase but few kinases | CK2 KO showed only ~1/3 of motif-matching sites drop, a reminder that perturbation response ≠ direct substrate |

**No systematic kinome-wide KO/KD phosphoproteomics resource for human was found.** Analog-sensitive (AS) chemical-genetic datasets exist only for single kinases.

### 3.2 Inhibitor → kinase target maps (needed to turn perturbations into kinase-specific features)
| Resource | Content | Access | Notes |
|---|---|---|---|
| `kinhib_inhib_comb__1.csv` (provided) | binary 420 × 60, Hijazi inhibitors | local | binary >50 %-inhibition-style calls. Lacks potency |
| **Klaeger et al. 2017** *Science* 358:eaan4368 (doi:10.1126/science.aan4368) | Kinobeads, 243 clinical kinase drugs, Kd_app for ~250+ kinases in cell lysate | PRIDE **PXD005336** [verified exists]; Science supplementary tables; ProteomicsDB | covers most decryptM kinase drugs (dasatinib, lapatinib, imatinib, gefitinib, afatinib, nintedanib …). Lysate binding ≠ cellular inhibition. Kinobeads miss some kinases (poor bead binders) |
| **Reinecke et al. 2024** *Nat Chem Biol* 20:577 (doi:10.1038/s41589-023-01459-3) | 1,183 tool compounds (PKIS, PKIS2, KCGS, Roche), ~500,000 compound–target interactions | MassIVE **MSV000092248**; ProteomicsDB; source data | needed if Bayer 2025 data are obtained (their 133 inhibitors are probably drawn from Klaeger+Reinecke) |
| **Davis et al. 2011** *Nat Biotechnol* 29:1046 (doi:10.1038/nbt.1990) | 72 inhibitors × 442 kinases, Kd (KINOMEscan) | supplementary; also in DTC/TDC | recombinant binding |
| **Drug Target Commons** (Tang et al. 2018 *Cell Chem Biol*) https://drugtargetcommons.fimm.fi | crowd-curated bioactivities incl. Anastassiadis 2011, Davis 2011, Klaeger | bulk CSV (connection did not complete in this session) | use to fill inhibitors not in Kinobeads |

### 3.3 Kinase-activity benchmark datasets (perturbation → expected kinase regulation)
| Dataset | Content | Access [verified] | Use for us |
|---|---|---|---|
| **Hernandez-Armenta et al. 2017** *Bioinformatics* 33:1845 (doi:10.1093/bioinformatics/btx082) | 103 curated perturbation experiments with quantile-normalised fold changes; 27 kinases; 173 kinase-condition positives (100 up / 73 down), as re-used by Müller-Dott | Zenodo 12566560 `benchmark_data_hernandez.csv` (29.6 MB) + metadata (CC-BY-4.0); older decoupleR record 5645208 | **downstream evaluation** of the predictor's use for kinase-activity inference (a project use-case). Many experiments overlap the compendium (Ochoa 2016 sources), so use them for evaluation **or** as features, never both |
| **Müller-Dott et al. 2025** *Nat Commun* 16 (2025) (doi:10.1038/s41467-025-59779-y; benchmarKIN https://github.com/saezlab/benchmarKIN, GPL-3) | perturbation benchmark (Hernandez-Armenta + Hijazi + tyrosine-kinase set) **and** a CPTAC tumour benchmark (top/bottom 5 % kinase protein or activating-site level → 12,850 kinase–tumour pairs, 388 kinases, 7 cohorts) | Zenodo 12566560 [verified file list, CC-BY-4.0] | primary **external application benchmark**: kinase activity inferred with our predicted KS network vs PSP / Kinase Library / Phosformer / NetworKIN. They found NetworKIN-added targets help in the tumour benchmark. The tumour benchmark is defined from kinase protein levels, so it is not circular with KS labels, but CPTAC correlation features computed on the same tumours could inflate it |
| Sugiyama, Imamura & Ishihama 2019 *Sci Rep* 9:10503 (doi:10.1038/s41598-019-46385-4) | in vitro kinase reactions with 385 kinases on dephosphorylated HeLa lysate peptides; 175,574 KS candidates | jPOST/PRIDE **PXD011366**; supplementary tables (CC-BY) | an experimental but non-physiological KS set (in vitro, peptide/lysate level). Best as a **secondary/auxiliary label or independent test** for kinase-specificity components. Not literature-curated, so per the project rules it is not a primary label |
| iKiP-DB (Mari et al. 2022 *J Proteome Res*, doi:10.1021/acs.jproteome.2c00198) | in vitro kinase-to-phosphosite on lysates | ACS supplementary zip | same role as Sugiyama |

---

## 4. Other large compendia and site universes

| Source | Content / size | Access | License | Value | Leakage |
|---|---|---|---|---|---|
| **Ochoa et al. 2020** *Nat Biotechnol* 38:365 (doi:10.1038/s41587-019-0344-3) | joint MaxQuant re-analysis of 6,801 runs (112 datasets, 104 cell types/tissues); 119,809 sites at 1 % site FDR; **functional score for 116,258 sites** (59 features) | PRIDE **PXD012174** (raw + MQ output); functional scores in Supplementary Tables | CC-BY (OA) | high for coverage (Priority 2) and as a site-level prior (detectability, conservation, structure). SELPHI2.0 uses it | gold standard = 2,638 PSP "functional" sites. Features include "number of conditions regulated" and "matching to kinase motifs". The score is **study-bias-correlated with having a known kinase**, so its value comes partly from annotation bias. Ablate it, and report metrics with and without it |
| **Ochoa et al. 2016** *Mol Syst Biol* (atlas of kinase regulation, ~435 conditions) | subset of the compendium | supplementary | CC-BY | already contained | Ochoa 2016 kinase activities were computed from PSP substrates, so they are leaky if used directly |
| **qPTM** (Yu et al. 2023 *NAR*; successor of qPhos, Yu et al. 2019 *NAR*) | 2,596 conditions across 4 species, 660,030 sites in all PTM types | http://qptm.omicsbio.info (bulk-download path returned 404) | academic | partially contained (461 human contrasts). The remaining human conditions might add some value, but incremental | — |
| **PhosphoSitePlus** `Phosphorylation_site_dataset.gz` | ~240k human sites with LTP/HTP reference counts | https://www.phosphosite.org/staticDownloads — **gated by login (OAuth redirect) [verified]** | CC BY-NC-SA 3.0 | defines the candidate-site universe used by SELPHI2.0 (238,374 sites) and by most baselines; HTP/LTP counts are useful study-bias covariates | HTP/LTP counts strongly correlate with being annotated → bias-only features; use for stratified reporting, not as predictive features (or ablate) |
| Bekker-Jensen et al. 2017 *Cell Syst* / 2020 *Nat Commun* (DIA, library-free phosphoproteomics) | deep single-line phosphoproteomes (HeLa etc.; tens of thousands of class-I sites) | PRIDE **[accessions not re-verified]** | CC0/CC-BY | site universe only | — |
| Koenig … Olsen 2025 (see §3.1) | >200,000 class-I sites, 33 lines | PXD068583 | CC-BY-NC | best current **empirical site universe** | — |
| Human multi-tissue phosphoproteome atlas ("Sun et al. 2022/2025") | **not found**. The 2024 Astral paper (Lancaster et al. *Nat Commun* 2024, doi:10.1038/s41467-024-51274-0) is a **mouse** multi-tissue atlas (81,120 sites). CPTAC normal-adjacent tissue files are the best human-tissue phospho source available | — | — | low | ask the PI/user whether a specific paper was meant |
| **CoPheeMap / CoPheeKSA** (Nat Commun 2025, https://github.com/bzhanglab/CoPheeMap) | ML co-regulation network of 26,280 sites from CPTAC; KS predictions (24,015 for 104 S/T kinases) | GitHub + Source Data | OA | **useful as a method template**; also an external baseline | **its edges were trained on "pairs of sites known to be phosphorylated by the same kinase"** (161,832 positive pairs, PSP-derived), so the published network and embeddings encode our labels. **Do not use as a feature**; if used, retrain within each fold |
| CoPhosK (Ayati et al. 2019 *PLoS Comput Biol*), NetKSA (PSB 2023) | co-phosphorylation KS predictors | GitHub | OA | baselines/ideas | same KS-prior issue |
| PTMsigDB (Krug et al. 2019 *MCP*) | site-specific perturbation and kinase signatures | https://proteomics.broadapps.org/ptmsigdb/ | — | perturbation signatures are useful; **kinase signatures derive from PSP**, so leaky | — |
| NetworKIN 3.1 predictions | motif+STRING KS scores | in Zenodo 12566560 (77.5 MB) | — | baseline | trained on PSP/Phospho.ELM |
| SELPHI2.0 outputs | KS probabilities | paper1 supp / selphi2.com | — | baseline only | trained on PSP (2024) |
| Single-cell phosphoproteomics | too sparse (hundreds of sites) | — | — | not useful now | — |

---

## 5. Methods: turning phosphoproteomics into KS features

### 5.1 Feature families (all per (kinase K, site s) pair unless stated)
1. **Substrate-centroid co-regulation** (per dataset D: compendium, each CPTAC cohort, Hijazi, decryptM).
   - Compute corr(s, t) for each training substrate t of K and aggregate: mean of top-k, median, max, and fraction with |r| > τ. Alternatively, correlate s with the **profile of K's inferred activity**, where the activity per sample is the mean/z-score of K's training substrates (as in CoPheeMap and in SELPHI).
   - Report n_overlap (number of samples observed in both s and K's substrates) as a companion feature. When n_overlap < a threshold (e.g. 6), output NaN rather than a noisy r.
   - Use robust correlation (Spearman, or biweight midcorrelation) and Fisher-z shrinkage by n: z·√(n−3).
2. **Kinase-protein / activating-site co-variation (CPTAC):** correlate site s with K's protein abundance, and with K's known activation-loop phosphosites (an activity proxy that does not depend on known substrates; the activating-site list is from PSP "regulatory sites", so check for overlap with labels only when s itself is on K, i.e. autophosphorylation).
3. **Inhibitor-target-informed response** (Hijazi, decryptM, Bayer 2025 if obtained, and inhibitor contrasts in the compendium):
   - For each inhibitor i with target vector T_i (Kinobeads Kd_app or a binary matrix), compute the response of s (log2FC, or the EC50 from curves).
   - Features:
     - (a) mean down-regulation of s across inhibitors that hit K vs those that do not (signed t-statistic or AUC, weighted by Kd_app). This mirrors the Hijazi approach and needs **no KS labels**, so it is leakage-free.
     - (b) **potency coherence:** |log EC50_site − log Kd_app(K, drug)| minimised across drugs (decryptM/Bayer). This is direct evidence that s sits downstream of, or is a substrate of, K. Also label-free.
     - (c) inhibitor-specificity-deconvolved site → kinase assignment: a regression of the site-response matrix on the target matrix (NNLS/elastic-net), giving β_{s,K}. Also label-free.
   - Label-free features are the most valuable for understudied kinases, because they do not need training substrates of K. They do need K to be targeted by some profiled inhibitor; of the 420 kinases in the Hijazi matrix, 307 are hit by at least one inhibitor.
4. **Perturbation-class responses:** responses to stimuli with a known dominant kinase, e.g. insulin → AKT/mTOR, EGF → ERK, DNA damage → ATM/ATR, mitosis → CDK1/PLK1/AURK, TNF → IKK/p38. This needs a curated condition → kinase table (benchmarKIN metadata and Hernandez-Armenta give a starting point). It is label-free with respect to KS pairs.
5. **Site-level priors (kinase-independent):** detectability (number of studies/conditions observed), variance and regulation frequency in the compendium, Ochoa functional score. These help the "which sites are phosphorylated by anyone" component but carry study bias.

### 5.2 Missing values and batch effects
- **Never impute zeros into fold-change matrices.** Use pairwise-complete correlations with an n_overlap threshold, plus explicit missingness indicators.
  - Missingness itself is informative (e.g. a pY site only seen in pY-enriched studies), but it is also a technical artefact. Keep it as a separate feature so ablations can show its contribution.
  - Tree models (LightGBM/XGBoost) handle NaN natively. If missingness rates differ between positives and negatives, CoPheeMap's trick of "matching missing-value ratio between positive and negative data" is worth copying, because it prevents the model from learning "missing → negative".
- **Study weighting:** compute correlations **within study blocks** and aggregate, or down-weight large studies (e.g. one effective vote per study or per cell line). Otherwise decryptM-rituximab dose series (306 columns, mostly B-cell receptor signalling) dominate.
  - Remove within-study shared-control artefacts: in a dose series all columns share a DMSO reference, so ratios are correlated by construction. Consider summarising dose series into one EC50/AUC value per site × drug × cell line rather than treating doses as independent samples.
- **CPTAC:** within-cohort median centring, then protein residualisation (site ~ host protein, per cohort) and a correlation per cohort. Aggregate with Fisher-z mean and the count of significant cohorts.
- **Isoform/multiplicity:** collapse to (UniProt canonical, residue) via the 15-mer, and record which source multiplicities were merged.

### 5.3 Leakage: rules (proposed for DESIGN.md)
1. **Every feature that uses KS annotations must be computed inside each CV fold from training-fold KS pairs only.** In this project that covers substrate-centroid correlations, kinase-activity profiles, PSSMs/motif models trained on PSP, learned co-regulation networks (CoPheeMap-like), PTMsigDB kinase sets, and NetworKIN/SELPHI/GPS scores.
   - Test sites (and, under protein- or homology-grouped splits, all sites on test proteins) must be removed from the substrate sets before activity profiles are computed.
   - Autophosphorylation sites must be removed from K's own activity estimate when scoring K's own sites, as in benchmarKIN and CoPheeMap.
   - Leave-kinase-out benchmarks are the special case: the held-out kinase has **no** training substrates, so only label-free features (§5.1 items 2–4) and features transferred from similar kinases are allowed for it.
2. **Precomputed third-party KS predictions or networks are baselines, not features**, unless retrained per fold. SELPHI2.0 used PSP-derived features and evaluated with stratified K-fold on KS pairs (not grouped by site/protein), which will inflate its reported AUROC. Our retrained-baseline comparison must use our splits.
3. **Feature-source vs evaluation-set declaration.** Hijazi, Sugiyama, Hernandez-Armenta and benchmarKIN CPTAC must each be fixed as "feature" or "evaluation" before modelling starts. If Hijazi is a feature source, its derived KS calls cannot be used as a test. If the benchmarKIN tumour benchmark is used for downstream evaluation, CPTAC co-regulation features computed on those same tumours should be flagged (or computed on the non-benchmark cohorts: COAD, OV, PDAC).
4. **PSP in vivo labels can originate from the same phosphoproteomics papers** (PSP curates HTP inhibitor studies). Perturbation features can therefore "re-discover" the evidence behind a label, which flatters in vivo labels relative to in vitro ones. Mitigation: report metrics stratified by IN_VIVO_RXN/IN_VITRO_RXN. Where feasible, compare compendium PMIDs with the PSP evidence (PSP KS rows have no PMIDs; the PSP `Kinase_Substrate_Dataset` only has flags, so this needs the PSP web references or a sensitivity analysis that drops each study in turn).
5. **Study bias:** detectability, number of conditions observed and the Ochoa functional score correlate with annotation status. Report per-kinase metrics, and metrics stratified by site "study depth" (PSP LTP/HTP counts if obtainable), so that gains are not just annotation-bias recovery.
6. **Homology:** paralog sites (e.g. AKT1/2/3 substrates and paralogous substrate proteins) can be near-identical in co-regulation profiles because peptides are shared. Remove shared/ambiguous peptides (non-unique to one protein), and group homologous proteins in splits (already required by the brief).

---

## 6. Prioritised download plan

| Priority | Item | Size | Auth | Action |
|---|---|---|---|---|
| 1 | User compendium (local) | 1.1 GB | — | convert to a sparse site × contrast matrix (parquet), map to PSP site groups via the 15-mer (the data-profiling agent is doing this) |
| 1 | Zenodo **12566560** (benchmarKIN: Hernandez-Armenta + Hijazi fold changes, CPTAC 7-cohort tables, gold standards, NetworKIN) | ~0.5 GB | none, CC-BY-4.0 | download to `external/benchmarkin/` |
| 1 | CPTAC pan-cancer phospho + proteome + meta, 10 cohorts (Tumor + Normal) from the LinkedOmics S3 mirror | ~0.5 GB phospho + ~0.3 GB proteome **[estimate]** | none | download to `external/cptac_pancan_v1.2/`; ENSP→UniProt via Ensembl 100 |
| 1 | Klaeger 2017 Kinobeads Kd_app tables (Science supp. / ProteomicsDB) | <50 MB | none | needed for decryptM inhibitor → target mapping |
| 2 | decryptM `*_Curves.zip` for kinase-inhibitor experiments (10_Kinase_Inhibitors, 3_EGFR, BCRABL, HER2, Dasatinib_Triplicates, Selumetinib_MK2206, AZD4547 combos, SHP099) from PXD037285 | ~5 GB | none, CC0 | EC50/curve-class per site → potency-coherence features |
| 2 | Ochoa 2020 functional score (Supplementary Table) | <50 MB | none | site prior + universe |
| 2 | Reinecke 2024 target landscape (MSV000092248 / source data) and Davis 2011 / DTC | <200 MB | none | inhibitor-target matrix expansion |
| 2 | Sugiyama 2019 supplementary KS table | <100 MB | none | auxiliary in vitro set / independent test |
| 3 | Frejno 2020 (125 lines), Koenig/Olsen 2025 atlas (PXD068583, processed tables only) | variable | none | extra co-regulation / site universe |
| 3 | CoPheeMap code + source data | small | none | baseline / method template, not features |

### Requests to the user (gated or unknown)
1. **Provenance of the compendium.** Which paper/pipeline, how the qPTM contrasts were processed (were non-significant values retained?), and whether a site-level FDR or localisation filter was applied. Also, whether a **wide-format or multiplicity-resolved version** exists.
2. **PhosphoSitePlus `Phosphorylation_site_dataset.gz`** (and `Regulatory_sites.gz`) for the site universe and LTP/HTP counts. The download is login-gated [verified OAuth redirect].
3. **Bayer et al. 2025 (Kuster) potency-coherence data:** 133 kinase inhibitors, 17 M curves. These are likely available only via ProteomicsDB or by request. They are the single most valuable new dataset for direct KS evidence across 96 kinases.
4. Confirm which "tissue phosphoproteome atlas (Sun et al.)" was meant; none could be verified.

---

## 7. Key references (DOIs / URLs)
- Zecha J. et al. *Science* 2023, decryptM — doi:10.1126/science.ade3925; PRIDE PXD037285.
- Bayer F.P. … Kuster B. bioRxiv 2025 — doi:10.1101/2025.11.18.689017.
- Hijazi M. et al. *Nat Biotechnol* 2020 — doi:10.1038/s41587-019-0391-9; PRIDE PXD015943.
- Klaeger S. et al. *Science* 2017 — doi:10.1126/science.aan4368; PRIDE PXD005336.
- Reinecke M. et al. *Nat Chem Biol* 2024 — doi:10.1038/s41589-023-01459-3; MassIVE MSV000092248.
- Davis M.I. et al. *Nat Biotechnol* 2011 — doi:10.1038/nbt.1990. Drug Target Commons — https://drugtargetcommons.fimm.fi.
- Li Y. et al. *Cancer Cell* 2023 — doi:10.1016/j.ccell.2023.06.009; Li Y. et al. *Cell* 2023 — doi:10.1016/j.cell.2023.07.014; https://pdc.cancer.gov/pdc/cptac-pancancer; https://www.linkedomics.org/data_download/
- CoPheeMap/CoPheeKSA, *Nat Commun* 2025 — doi:10.1038/s41467-025-57993-2; https://github.com/bzhanglab/CoPheeMap.
- Hernandez-Armenta C. et al. *Bioinformatics* 2017 — doi:10.1093/bioinformatics/btx082.
- Müller-Dott S. et al. *Nat Commun* 2025 — doi:10.1038/s41467-025-59779-y; Zenodo doi:10.5281/zenodo.12566560; https://github.com/saezlab/benchmarKIN.
- Sugiyama N., Imamura H., Ishihama Y. *Sci Rep* 2019 — doi:10.1038/s41598-019-46385-4; PXD011366.
- Ochoa D. et al. *Nat Biotechnol* 2020 — doi:10.1038/s41587-019-0344-3; PRIDE PXD012174. Ochoa D. et al. *Mol Syst Biol* 2016 (PMID 27974457).
- Yu K. et al. qPTM, *NAR* 2023 — http://qptm.omicsbio.info.
- Frejno M. et al. *Nat Commun* 2020 — doi:10.1038/s41467-020-17336-9.
- Koenig C. … Olsen J.V. bioRxiv 2025 — doi:10.64898/2025.12.05.692571; PRIDE PXD068583.
- Lancaster N.M. et al. *Nat Commun* 2024 (Astral, mouse tissue atlas) — doi:10.1038/s41467-024-51274-0.
- Maier B.D. … Petsalaki E. SELPHI2.0, *MCP* 2025 — doi:10.1016/j.mcpro.2025.100994 (data/paper1.pdf).
- Krug K. et al. PTMsigDB, *MCP* 2019 — PMID 30563849.
