# Data sources

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

**No input data are included in this repository.** The only data files here are the released predictions and
metric tables, which are derived scores (`predictions/`, `results_summary/`). To re-run the pipeline you must obtain
every dataset below yourself, under its own licence, and place it in the layout given at the end of this file.

This file summarises the project's data register. The full register, with the SHA-256 of every file, access dates
and exact download URLs, is `docs/DATA_internal.md`; its paths are relative to the project directory. All
downloads were made on 2026-10-26.

Column "Redistribution" states what the licence allows, as far as the project could determine. It does not mean
that the file is in this repository; none is. "Obtain yourself" marks data that are licensed, restricted,
non-commercial or of unclear terms. Where the terms are not stated, ask the provider.

## 1. Data supplied by the human researcher (`data/`)

| Dataset | Used for | Version | Source | Licence | Redistribution | How to obtain |
|---|---|---|---|---|---|---|
| **PhosphoSitePlus** `Kinase_Substrate_Dataset` | all labels (training, CV, TEST), hard negatives | release header "051925" (2025-05-19) | Cell Signaling Technology, www.phosphosite.org (Hornbeck et al. 2012, doi:10.1093/nar/gkr1122; 2015, doi:10.1093/nar/gku1267; 2019, doi:10.1093/nar/gky1159) | CC BY-NC-SA 3.0 | non-commercial, share-alike, with attribution | **Obtain yourself under its licence**: register and download from www.phosphosite.org (Downloads) |
| PhosphoSitePlus `Phosphorylation_site_dataset` | site universe (observed sites) | same | same | CC BY-NC-SA 3.0 | as above | **Obtain yourself**, as above |
| PhosphoSitePlus `Regulatory_sites` | list of excluded studies (A12) | same | same | CC BY-NC-SA 3.0 | as above | **Obtain yourself**, as above |
| **Klaeger et al. 2017** Kinobeads, Tables S1–S3 | drug × kinase occupancy for the inhibitor feature (F6b); downstream drug annotation | publisher supplement | *Science*, doi:10.1126/science.aan4368 | publisher supplementary material | unclear; do not redistribute | **Obtain yourself** from the article's supplementary materials |
| **Perturbation phosphoproteomics compendium** (`01_phosphoproteomics_conditions.csv`, `02_phosphosite_information.csv`, `03_phosphosite_quantifications.csv`) | co-regulation (F6a), inhibitor deconvolution (F6b), site universe, downstream benchmark | as provided (1,937 contrasts, 178 PMIDs) | compiled from published studies (PMIDs per condition in file 01); the provenance of the compilation is not documented (`docs/DATAINEED.md` item 5) | not stated; derived from published studies | unclear; do not redistribute | **Obtain yourself**: ask the authors of the compilation, or rebuild it from the source studies listed by PMID |
| Hijazi et al. 2020 kinase × inhibitor matrix (`kinhib_inhib_comb__1.csv`) | inhibitor feature (F6b), kinase list check | as provided | Hijazi et al. 2020, PMID 31959955 | not stated | unclear | **Obtain yourself** from the publication [verify the exact supplementary file] |
| `paper1.pdf` (SELPHI2.0; Maier et al. 2025 *MCP*, doi:10.1016/j.mcpro.2025.100994) and `paper2.pdf` (Phosformer; Zhou et al. 2023 *Bioinformatics*, doi:10.1093/bioinformatics/btad046) | read by agents only; not used by code | — | publishers | publisher copyright | no | from the publishers |

## 2. Downloaded public and semi-public data (`external/`)

"Needed for" uses these tags: **model** = features of the released model; **labels** = extra curated labels
(external TEST stratum, ablations); **base** = baselines; **down** = downstream kinase-activity benchmark;
**audit** = leakage audits only; **unused** = downloaded and assessed, but the included code does not read it.

| Directory | Content | Needed for | Version | URL / DOI | Licence | Redistribution |
|---|---|---|---|---|---|---|
| `uniprot/` | human reference proteome UP000005640 (FASTA with isoforms; annotation, PTM and kinase tables; pkinfam) | model, labels | UniProt 2026_03 | https://rest.uniprot.org (queries in `scripts/download_uniprot.py`) | CC BY 4.0 | yes |
| `uniprot_topology/` | topology and cross-references, reviewed human | model | UniProt 2026_03 | https://rest.uniprot.org | CC BY 4.0 | yes |
| `alphafold/` | AFDB human proteome models (`UP000005640_9606_HUMAN_v6.tar`) and a PAE calibration subset (`pae/`, 399 proteins) | model | AFDB v6 | https://ftp.ebi.ac.uk/pub/databases/alphafold/ ; `scripts/download_afdb_pae.py` | CC BY 4.0 | yes |
| ESM-2 weights (Hugging Face cache, not under `external/`) | `facebook/esm2_t33_650M_UR50D` (revision 08e4846e…), `facebook/esm2_t36_3B_UR50D` (476b6399…) | model | see revisions | https://huggingface.co/facebook | MIT | yes |
| `kinase_library/` | kinase-library 1.8.0 sdist, extracted as data only (Cantley S/T and Y matrices, favourability tables) | model, base (B1) | 1.8.0 | https://pypi.org/project/kinase-library/1.8.0/ | CC BY-NC-SA 3.0 | non-commercial; **obtain yourself** |
| `cantley_johnson2023/` | S/T atlas supplement (MOESM1–6) | model, base | Johnson et al. 2023 *Nature*, doi:10.1038/s41586-022-05575-3 | Springer static content | publisher supplement; academic use | **obtain yourself** |
| `cantley_yaronbarir2024/` | Y atlas supplement (MOESM1–9; SH2 matrices) | model, base | Yaron-Barir et al. 2024 *Nature*, doi:10.1038/s41586-024-07407-y | Springer static content | publisher supplement; academic use | **obtain yourself** |
| `bioplex/` | BioPlex 3.0 293T and HCT116 networks, interactionList_v4a | model | BioPlex 3.0 | https://bioplex.hms.harvard.edu/data/ | CC BY 4.0 (check site) | likely yes |
| `huri/` | HuRI, Lit-BM, HI-union | model (HuRI only) | as of 2026-10-26 | https://interactome-atlas.org/data/ | CC BY 4.0 (check site) | likely yes |
| `opencell/` | OpenCell interactions, localisation, abundance | model | as of 2026-10-26 | https://opencell.czbiohub.org/data/datasets/ | CC BY 4.0 | yes |
| `buljan2020/` | Buljan et al. 2020 *Mol Cell* supplement (kinase AP-MS interactome) | model | doi:10.1016/j.molcel.2020.07.001 | Europe PMC PMC7427327 | CC BY-NC-ND (article) | non-commercial, no derivatives; **obtain yourself** |
| `humancellmap/` | Human Cell Map baits, preys, SAINT | model | "latest" | https://humancellmap.org/resources/downloads/ | CC BY 4.0 | yes |
| `hpa/` | HPA subcellular location, RNA tissue consensus, RNA cell line, proteinatlas.tsv | model, base | HPA v25 | https://www.proteinatlas.org/download/ | CC BY-SA 4.0 | yes, share-alike |
| `hein2025_zenodo_14000085/` | Hein et al. 2025 organelle IP-MS analysis tables | model | Zenodo 14000085 | https://zenodo.org/records/14000085 | CC BY 4.0 | yes |
| `compartments/` | COMPARTMENTS knowledge, experiments, integrated | model (experiments channel; knowledge channel ablation only) | as of 2026-10-26 | https://download.jensenlab.org/ | CC BY 4.0 | yes |
| `go/` | GO annotations (goa_human.gaf.gz), go-basic.obo | model (localisation) | as of 2026-10-26 | https://current.geneontology.org/ | CC BY 4.0 | yes |
| `depmap/` | DepMap 26Q1 Chronos gene effect | model (co-essentiality) | 26Q1 | https://api.figshare.com/v2/articles/31660582 | CC BY 4.0 | yes |
| `cptac_pancan_v1.2/` | CPTAC phospho and protein abundance, 10 cohorts (BCM harmonised freeze v1.2, LinkedOmics) | model (F6a, F6c) | v1.2 | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/ | open access under the CPTAC data use policy; cite | **obtain yourself** under the CPTAC policy |
| `cptac_pancancer_pdc/` | CPTAC pan-cancer harmonised files (Li et al. 2023 *Cancer Cell*) and clinical metadata | model (confounders) | harmonised v1 | https://pdc.cancer.gov (`scripts/dl_pdc_pancancer.py`) | open access under the NCI PDC / CPTAC policy; cite | **obtain yourself** under that policy |
| `reinecke2024/` | target landscape of 1,000 kinase inhibitors (supplement) | model (F6b) | Reinecke et al. 2024 *Nat Chem Biol*, doi:10.1038/s41589-023-01459-3 | Springer | publisher supplement | **obtain yourself** |
| `curvecurator_zenodo_8399823/` | Kinobeads dose–response example outputs | model (F6b) | Zenodo 8399823 | https://zenodo.org/records/8399823 | CC BY 4.0 | yes |
| `davis2011/` | Davis et al. 2011 kinase inhibitor Kd supplement | model (F6b) | doi:10.1038/nbt.1990 | Springer | publisher supplement | **obtain yourself** |
| `hgnc/` | HGNC complete set | model, labels (name resolution) | as of 2026-10-26 | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | CC0 | yes |
| `kinbase_manning2002/` | Manning et al. 2002 kinome tables S1/S2 | model (kinase groups) | 2002 | http://kinase.com/human/kinome/tables/ | supplementary data of a published paper | **obtain yourself** |
| `kinhub/` | KinHub kinase table (HTML) | labels (kinase list) | as of 2026-10-26 | http://www.kinhub.org/kinases.html | not stated | **obtain yourself** |
| `kincore/` | KinCore human kinase-domain alignment | labels (paralog masks), splits | as of 2026-10-26 | https://dunbrack.fccc.edu/kincore/ | not stated (academic); cite Modi & Dunbrack | **obtain yourself** |
| `omnipath/` | OmniPath enzyme–substrate table (academic licence set) | labels (external curated stratum; ablation) | web service, 2026-10-26 | https://omnipathdb.org/enzsub?...&license=academic | mixed per resource; includes non-commercial resources such as PSP | **obtain yourself** |
| `iptmnet/` | iPTMnet 6.2 | labels audit (study exclusion list) | 6.2 | https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ | CC BY-NC-SA 4.0 | non-commercial; **obtain yourself** |
| `ochoa2020/` | Ochoa et al. 2020 *Nat Biotechnol* supplement (functional score, phosphoproteome features) | base (B4 SELPHI2.0), site universe | doi:10.1038/s41587-019-0344-3 | Springer | publisher supplement | **obtain yourself** |
| `selphi2/` | SELPHI2.0 code (GitLab clone, read as reference, not executed) and Zenodo small files | base (B4) | GitLab commit f9822f90…; Zenodo 10.5281/zenodo.13821321 | https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0 ; https://zenodo.org/records/13821321 | code Apache-2.0; Zenodo CC BY 4.0 | yes |
| `gtex_v8/` | GTEx v8 gene median TPM | base (B4) | v8 (2017-06-05) | https://storage.googleapis.com/adult-gtex/ | GTEx open-access data; cite | yes (summary data) |
| `phosformer/` | Phosformer code (GitHub clone), weights (`hf_phosformer`), curated dataset and splits (`hf_dataset_curated`), Phosformer-ST space | base (B3) | GitHub commit 6e7c019b…; HF as of 2022-11-04 / 2023-01-25 / 2025-04-16 | https://github.com/esbgkannan/phosformer ; https://huggingface.co/waylandy/phosformer ; https://huggingface.co/datasets/waylandy/phosformer_curated | code and weights CC BY-NC-ND 4.0; dataset licence not stated (treated as non-commercial) | non-commercial, no derivatives; **obtain yourself** |
| `benchmarkin/` | benchmarKIN perturbation and tumour benchmarks, NetworKIN predictions | down; audit (study overlap) | Zenodo 10.5281/zenodo.12566560 (Müller-Dott et al. 2025 *Nat Commun*, doi:10.1038/s41467-025-59779-y) | https://zenodo.org/records/12566560 | CC BY 4.0 | yes |
| `benchmarkin_github/` | tyrosine perturbation set (`.rda`) | down | benchmarKIN main @ 2aa6e25c | https://github.com/saezlab/benchmarKIN | GPL-3 | yes, under GPL-3 |
| `string/` | STRING v12.0 aliases, info, links | ID mapping; audit | v12.0 | https://stringdb-downloads.org/download/ | CC BY 4.0 | yes |
| `biogrid/` | BioGRID 5.0.262 organism and PTM files | audit only | 5.0.262 | https://downloads.thebiogrid.org/ | MIT | yes |
| `intact/`, `sugiyama2019/`, `ikipdb_zenodo_*`, `kiba_zenodo_5105698/`, `kinasecanvas_zenodo_22941318/`, `kinoplex_*`, `kolossus/`, `orphan_top5_zenodo_16881580/`, `proteomehd_kustatscher2019/`, `zenodo_meta/` | assessed during design (see `docs/DATA_internal.md`, `docs/lit/`) | unused | see register | see register | see register | not needed |

## 3. Software tools used outside Python

- **MMseqs2** builds the substrate homology clusters (splits) and filters the Phosformer corpus. Set
  `KSP_MMSEQS` to the binary, or put `mmseqs` on `PATH`. The project used MMseqs2 commit `8cc5ce367b5638c4306c2d7cfc652dd099a4643f` (output of `mmseqs version`).

## 4. Expected local directory layout

Point `KSP_DATA_DIR` at a directory with this layout. If it is unset, the repository root is used; `data/`,
`external/`, `cache/`, `results/`, `sealed/` and `logs/` are in `.gitignore`. File names must match exactly.
`docs/DATA_internal.md` gives the SHA-256 of each file.

```
$KSP_DATA_DIR/
├── data/
│   ├── Kinase_Substrate_Dataset              # PhosphoSitePlus, release 051925 (uncompressed, as downloaded)
│   ├── Phosphorylation_site_dataset          # PhosphoSitePlus
│   ├── Regulatory_sites                      # PhosphoSitePlus
│   ├── Klaeger/aan4368_Table_S1.xlsx
│   ├── Klaeger/aan4368_Table_S2.xlsx
│   ├── Klaeger/Table_S3.xlsx
│   ├── 01_phosphoproteomics_conditions.csv
│   ├── 02_phosphosite_information.csv
│   ├── 03_phosphosite_quantifications.csv
│   └── kinhib_inhib_comb__1.csv
├── external/
│   ├── uniprot/{MANIFEST.json, proteome/UP000005640_all_isoforms.fasta.gz, proteome/UP000005640_reviewed_isoforms.fasta.gz,
│   │            human_tsv/UP000005640_annotations.tsv.gz, kinase_kw0418/human_kw0418_reviewed.tsv.gz,
│   │            kinase_detail/human_kw0418_detail.tsv.gz, ptm_tsv/UP000005640_reviewed_ptm.tsv.gz, pkinfam/pkinfam.txt}
│   ├── uniprot_topology/UP000005640_reviewed_topology.tsv.gz
│   ├── alphafold/{UP000005640_9606_HUMAN_v6.tar, README.txt, CHANGELOG.txt, download_metadata.json,
│   │              pae/AF-<acc>-F1-predicted_aligned_error_v6.json.gz, pae/MANIFEST.tsv}
│   ├── kinase_library/src/kinase_library-1.8.0/src/kinase_library/databases/{kinase_data,matrices,substrates}/
│   ├── cantley_johnson2023/41586_2022_5575_MOESM{1..6}_ESM.{pdf,xlsx}
│   ├── cantley_yaronbarir2024/41586_2024_7407_MOESM{1..9}_ESM.{pdf,xlsx}
│   ├── bioplex/{BioPlex_293T_Network_10K_Dec_2019.tsv, BioPlex_HCT116_Network_5.5K_Dec_2019.tsv, BioPlex_interactionList_v4a.tsv}
│   ├── huri/{HuRI.tsv, Lit-BM.tsv, HI-union.tsv}
│   ├── opencell/{opencell-library-metadata.csv, opencell-localization-annotations.csv,
│   │             opencell-protein-abundance.csv, opencell-protein-interactions.csv}
│   ├── buljan2020/supp_PMC7427327.zip
│   ├── humancellmap/{baits-latest.txt, preys-latest.txt, saint-latest.txt}
│   ├── hpa/{subcellular_location.tsv.zip, rna_tissue_consensus.tsv.zip, rna_celline.tsv.zip, proteinatlas.tsv.zip}
│   ├── hein2025_zenodo_14000085/Organelle_IP_analyses_and_figures-2024-10-28.zip
│   ├── compartments/human_compartment_{knowledge,experiments,integrated}_full.tsv
│   ├── go/{goa_human.gaf.gz, go-basic.obo}
│   ├── depmap/DepMap_26Q1_chronos_{gene,library}_effect.csv
│   ├── cptac_pancan_v1.2/<COHORT>/<COHORT>_{meta.txt, phospho_site_abundance_..._{Tumor,Normal}.txt,
│   │                      proteomics_gene_abundance_..._{Tumor,Normal}.txt}   # BRCA CCRCC COAD GBM HNSCC LSCC LUAD OV PDAC UCEC
│   ├── cptac_pancancer_pdc/{Clinical_meta_data_v1.zip, Phosphoproteome_*.zip|tsv.gz, Proteome_*.zip|tsv.gz,
│   │                        pancancer_file_manifest.json}
│   ├── reinecke2024/{41589_2023_1459_MOESM{1..8}_ESM.*, supp.zip}
│   ├── curvecurator_zenodo_8399823/Kinobeads.zip
│   ├── davis2011/41587_2011_BFnbt1990_MOESM{2..6}_ESM.{pdf,xls}
│   ├── hgnc/hgnc_complete_set.txt
│   ├── kinbase_manning2002/{TableS1.txt, TableS1.xls, TableS2.txt}
│   ├── kinhub/kinhub_kinases.html
│   ├── kincore/{alignment-files_Human-PK-alignment.aln, alignment-files_Human-PK-alignment.fasta,
│   │            alignment-files_Human-PK.fasta, fasta_with_labels_PK_labels_PDB.fasta}
│   ├── omnipath/omnipath_enzsub.tsv
│   ├── iptmnet/{ptm.txt, protein.txt, score.txt, readme.txt}
│   ├── ochoa2020/41587_2019_344_MOESM{3..7}_ESM.xlsx
│   ├── selphi2/{selphi_2.0/ (GitLab clone), zenodo_13821321/*.tsv|*.csv}
│   ├── gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
│   ├── phosformer/{github/ (clone), hf_phosformer/, hf_dataset_curated/{corpus,curated,splits}/, hf_space_phosformerST/}
│   ├── benchmarkin/*.csv|*.tsv|*.Rds          # 14 files, see docs/DATA_internal.md
│   ├── benchmarkin_github/{tyrosineData.rda, tyrosineMeta.rda}
│   ├── string/9606.protein.{aliases,info,links.detailed,physical.links.detailed}.v12.0.txt.gz
│   └── biogrid/{BIOGRID-ORGANISM-LATEST.tab3.zip, BIOGRID-PTMS-LATEST.ptm.zip}     # audit only
├── cache/      # written by the pipeline (indices, feature stores, predictions); tens of GB
├── results/    # written by the pipeline (splits, run manifests, metrics)
├── sealed/     # written once by scripts/build_splits.py (TEST labels); read only by ksp.eval_test and ksp.deploy
└── logs/
```

ESM-2 weights are read from the Hugging Face cache (`HF_HOME`). Run `scripts/download_uniprot.py`,
`scripts/download_afdb_pae.py`, `scripts/dl.sh`, `scripts/dl_w2.sh`, `scripts/dl_w2b.sh` and
`scripts/dl_pdc_pancancer.py` to fetch part of the public data. They write under `$KSP_DATA_DIR/external/`. Check
each file against the SHA-256 in `docs/DATA_internal.md`. Several sources change between releases (UniProt, GO,
HPA, OmniPath, STRING), so new downloads will not match the recorded checksums. The splits, the site universe and
all numbers will then differ.
