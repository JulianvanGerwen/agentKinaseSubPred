# DATA — datasets used (source, version, checksum, license)

Restricted/large data live in `data/` (provided) and `external/` (downloaded); neither is committed. Checksums are SHA-256.

## Provided by the user (data/)
| File | Source / version | SHA-256 | License / notes |
|---|---|---|---|
| Kinase_Substrate_Dataset | PhosphoSitePlus (CST), release header "051925" (2025-05-19) | bf057b8bb347881f8d0ac82d705adcc0be6f346f8c72882078fbd440a2ed3c4c | CC BY-NC-SA 3.0; attribution to PhosphoSitePlus required (Hornbeck et al. 2012/2015) |
| Phosphorylation_site_dataset | PhosphoSitePlus (CST), release header "051925" (2025-05-19); provided by user 2026-10-26 | 9d79ab45cbcf9f4813db1afcfe0520f3202e579b75e119c8455ef930443d49e6 | CC BY-NC-SA 3.0; attribution required |
| Regulatory_sites | PhosphoSitePlus (CST), release header "051925" (2025-05-19); provided by user 2026-10-26 | 1e654c316bee94dcf2713544796ed2dfbfdcb1ebd42b72db89166bf6d6c974c7 | CC BY-NC-SA 3.0; attribution required |
| 01_phosphoproteomics_conditions.csv | Perturbation phosphoproteomics compendium (condition metadata; PMIDs per condition) — provided by user | 1257eb673d1b6db532773f4d6a54a75f2960a7f4434ccb67b0239834170af690 | Provided by user; derived from published studies (see PMIDs) |
| 02_phosphosite_information.csv | Same compendium (peptide_UID → 15-mer window) | 56cab92d90ea1760302169674552161bee405d946ff32c8221935f483dac8bfc | as above |
| 03_phosphosite_quantifications.csv | Same compendium (long-format log2FC) | 8877144f25cc95b4625e9c55e05d33ffe2ef11af07b1297bf27791d13d3198d7 | as above |
| kinhib_inhib_comb__1.csv | Kinase × inhibitor target matrix for Hijazi et al. 2020 (PMID 31959955) | b0217e7962ba936f18a7e4d060d849d10d0b38d2f0e04d3aa730d0c0861ee545 | Provided by user |
| paper1.pdf | Provided paper | f33d8ab5fd2bd3cab05c2b7deceb11ed7b7770f330a918f8a065bf6800ef99af | Publisher copyright; not redistributed |
| paper2.pdf | Provided paper | 57b46d368a07f0d1ddc9693f6f901ad583fe11802c7ae9a80497449ffd12b0b0 | Publisher copyright; not redistributed |

## Downloaded (external/)
(added as acquired; each with URL, access date, version, SHA-256, license)

### UniProt (external/uniprot/; fetched by scripts/download_uniprot.py; manifest external/uniprot/MANIFEST.json)
| File | Description | URL | Release / date | Accessed (UTC) | SHA-256 | License |
|---|---|---|---|---|---|---|
| external/uniprot/proteome/UP000005640_all_isoforms.fasta.gz | UniProt human reference proteome UP000005640, all entries (20,416 Swiss-Prot + 127,104 TrEMBL), canonical + isoforms (169,651 sequences) | https://rest.uniprot.org/uniprotkb/stream?query=proteome%3AUP000005640&format=fasta&includeIsoform=true&compressed=true | 2026_03 (02-September-2026) | 2026-10-26T12:34:17 | 85c818d1a62dc2eed012c39068c2228b3126fe01366f54953bc56f72c8226fc3 | CC BY 4.0 (UniProt) |
| external/uniprot/proteome/UP000005640_reviewed_isoforms.fasta.gz | UP000005640 Swiss-Prot only, canonical + isoforms (42,547 sequences) | https://rest.uniprot.org/uniprotkb/stream?query=proteome%3AUP000005640+AND+reviewed%3Atrue&format=fasta&includeIsoform=true&compressed=true | 2026_03 (02-September-2026) | 2026-10-26T12:35:19 | 63011ee75fbfd0c92f35412fb660cfe920980b0c96276f3f941de91039e8ece1 | CC BY 4.0 (UniProt) |
| external/uniprot/human_tsv/UP000005640_annotations.tsv.gz | UP000005640 per-entry table (gene names, keywords, domains, families, EC, seq version) | https://rest.uniprot.org/uniprotkb/stream?query=proteome%3AUP000005640&format=tsv&fields=accession%2Cid%2Creviewed%2Cgene_primary%2Cgene_names%2Cprotein_name%2Corganism_id%2Clength%2Csequence_version%2Cdate_sequence_modified%2Ckeywordid%2Ckeyword%2Cft_domain%2Cprotein_families%2Cec&compressed=true | 2026_03 (02-September-2026) | 2026-10-26T12:38:26 | bb1ab565fc7b3cbee92cdf271f36980d5483a14fe172d8ccd5faa4f7d8c0c6dd | CC BY 4.0 (UniProt) |
| external/uniprot/kinase_kw0418/human_kw0418_reviewed.tsv.gz | Reviewed human entries with keyword KW-0418 Kinase (625; 459 with a 'Protein kinase' domain) | https://rest.uniprot.org/uniprotkb/stream?query=organism_id%3A9606+AND+reviewed%3Atrue+AND+keyword%3AKW-0418&format=tsv&fields=accession%2Cid%2Creviewed%2Cgene_primary%2Cgene_names%2Cprotein_name%2Corganism_id%2Clength%2Csequence_version%2Cdate_sequence_modified%2Ckeywordid%2Ckeyword%2Cft_domain%2Cprotein_families%2Cec&compressed=true | 2026_03 (02-September-2026) | 2026-10-26T12:38:29 | 509512f55f8d7ced99fc1d0bc7fbca30bdf2c76fbace1544f16408c99c24d465 | CC BY 4.0 (UniProt) |
| external/uniprot/pkinfam/pkinfam.txt | UniProt 'Human and mouse protein kinases: classification and index' (pkinfam; 512 human entries; release 2026_03 of 02-Sep-2026) | https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/pkinfam.txt | 2026_03 (02-Sep-2026) | 2026-10-26T12:38:58 | 3c6051ea748c91b1e02a60df2374eb14c56a01ccf97bbffe667a0f83b82d4e2d | CC BY 4.0 (UniProt) |

## Downloaded (external/) -- wave 1, access date 2026-10-26
All under /workspace/external/<dir>/. Per-file SHA-256 below each source; full list in logs/dl/sha256_wave1.txt (kinase_library extracted tree: logs/dl/sha256_kinase_library_extracted.txt).

### benchmarkin
- URL: https://zenodo.org/records/12566560
- Version/release: Zenodo 10.5281/zenodo.12566560 (2024-06-27; Muller-Dott et al.)
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 450.7 MB

| file | bytes | SHA-256 |
|---|---|---|
| benchmark_data_hernandez.csv | 29574322 | c1ec0bcd341aeab7c716162e1d52b3efa726ba2fbd510aea5b6b8b071d19ed17 |
| benchmark_data_hijazi.csv | 19015174 | 4d804f54f5891fc1c509fdf79bb69b0ada261366d549856673d095bc3cdd1c43 |
| benchmark_metadata_hernandez.csv | 9023 | d788452df69603aea0090ad1471871c5c92ea4303fc1c1f7b07f8a68aa8e71c8 |
| benchmark_metadata_hijazi.csv | 4594 | 7b1ea19411dfc9bccc7bbec9642d675bc3bd98467387501cb07fb758a012b84f |
| benchmarking_gold_standard_top5per.Rds | 1697996 | f4059a599e0ff2e967fe263a5057ce2daee31a8d562ce25e317107e1b5ed18bd |
| brca_original_medcent_30plus.tsv | 24225606 | 044f928f179df467f1750724b2f3d86ba5d5148e31081a80a1d854417bbebfb4 |
| ccrcc_original_medcent_30plus.tsv | 31109868 | b4c82fd51ff73f6a79f79ef5e7e30a95dc77b1bc5c254ab28792c331cf3fa79b |
| gbm_original_medcent_30plus.tsv | 52438617 | b51b8da6497647067f0b2b008955341bf33d65122bf994dbc5046248032ed1ef |
| hnscc_original_medcent_30plus.tsv | 42151714 | 09b8d055f5dd5a3636ebcc24e75d46958d6c50308184b7fe98795fbc2d931b17 |
| lscc_original__medcent_30plus.tsv | 45979460 | 3877ed4b9f002520dfb38cc2ea0a9f0b4bf295c12cddf86042eae1940b65c274 |
| lscc_original_medcent_30plus.tsv | 45979460 | 3877ed4b9f002520dfb38cc2ea0a9f0b4bf295c12cddf86042eae1940b65c274 |
| luad_original_medcent_30plus.tsv | 38972187 | a72763006ef62c86885d83eb3fff1d2d28fe000e98df57371dd88b6b0dd63f30 |
| networkin_human_predictions_3.1.tsv | 77469548 | ec5fef0cd34e828cfd2723cd543b501f238d42c1c7a38e8cad4bb4a227bbda03 |
| ucec_original_medcent_30plus.tsv | 42109467 | f1be3de5a2d927e32457270350a79a1e92107586c94e1c2058168988e907238f |

### biogrid
- URL: https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/ (ORGANISM tab3 zip = per-species tab3 files incl. Homo_sapiens; PTMS ptm zip)
- Version/release: BioGRID 5.0.262
- Accessed: 2026-10-26; License: MIT; total size: 248.1 MB

| file | bytes | SHA-256 |
|---|---|---|
| BIOGRID-ORGANISM-LATEST.tab3.zip | 188866951 | 096a9e4d871ad5d698bd503b2b9a2409297b88f5b55bc479c710f0c700429820 |
| BIOGRID-PTMS-LATEST.ptm.zip | 59206317 | 6db862dd3eaeed866314281b6fa516b4ba8206cbde8a6b4a5688ca2c4e4f419c |

### bioplex
- URL: https://bioplex.hms.harvard.edu/data/
- Version/release: BioPlex 3.0 (293T 10K Dec 2019, HCT116 5.5K Dec 2019) + interactionList_v4a (combined)
- Accessed: 2026-10-26; License: CC BY 4.0 (check site); total size: 26.1 MB

| file | bytes | SHA-256 |
|---|---|---|
| BioPlex_293T_Network_10K_Dec_2019.tsv | 12866060 | 5a98def087fc77c7bf07aac79e2553fd64ea1a27a7d0c4617870476f3f120833 |
| BioPlex_HCT116_Network_5.5K_Dec_2019.tsv | 7714364 | b37be6de6fe15e948151ab9480cfcb88a6795780cf7a4059b386b07f3c1e760d |
| BioPlex_interactionList_v4a.tsv | 5485753 | 8d0c4cbbcd4efc4bf994aac080ce4a8468bc076c298d13bd4e0269958129c2fd |

### cantley_johnson2023
- URL: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/ (MOESM1-6)
- Version/release: Johnson et al. 2023 Nature suppl. (MOESM3 T1 kinase info, MOESM4 T2 PSSMs, MOESM5 T3 site scores/percentiles [370MB], MOESM6 T4)
- Accessed: 2026-10-26; License: Publisher supplementary data; academic use, cite; total size: 429.7 MB

| file | bytes | SHA-256 |
|---|---|---|
| 41586_2022_5575_MOESM1_ESM.pdf | 56583127 | 3b0f8c34f031c62041c2d33fddff6abd0e65395e2b3a5d9975ce0a83b26a5c03 |
| 41586_2022_5575_MOESM2_ESM.pdf | 1290569 | d4f014418d82ed066364d8dc862cd06c02677dfc82bed1165ff0256b1a8856da |
| 41586_2022_5575_MOESM3_ESM.xlsx | 43791 | 4a0f4d39f748925879499f7f6c55e0843d30fc16401f3bcd9f2c05e49b33cd92 |
| 41586_2022_5575_MOESM4_ESM.xlsx | 1604728 | 9f690bf4e8e9710a8354b1e5371653ad8dd1e318d16bdf2c15e5a92f9fefe059 |
| 41586_2022_5575_MOESM5_ESM.xlsx | 370189209 | bb2c71e3daaff7d2570693a7fd168ce606cb84ab5ed9c2606437a8f2b5b4b9db |
| 41586_2022_5575_MOESM6_ESM.xlsx | 15550 | 6713094372286cfbd3c696327af03381d83229da2945ab93b2d270e261428331 |

### cantley_yaronbarir2024
- URL: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/ (MOESM1-9)
- Version/release: Yaron-Barir et al. 2024 Nature suppl. (MOESM4 T2 Tyr PSSMs, MOESM5 T3 annotation, MOESM7 SH2 matrices)
- Accessed: 2026-10-26; License: Publisher supplementary data; academic use, cite; total size: 34.7 MB

| file | bytes | SHA-256 |
|---|---|---|
| 41586_2024_7407_MOESM1_ESM.pdf | 15572092 | 97f20399e835332712245857a705e29a13067f85cb7ee0698842fccae9ddd440 |
| 41586_2024_7407_MOESM2_ESM.pdf | 1364793 | bf73c0bceed1c103c317b64622554a2f17d974a7265a13b863b86a8c4794321d |
| 41586_2024_7407_MOESM3_ESM.xlsx | 23167 | 6e1a14acc2f84378e2953e2bb8ff0789fcbb7eff491ab912563a66fb6c061669 |
| 41586_2024_7407_MOESM4_ESM.xlsx | 552780 | 13ff77e49f737e787cc297c4718b672f9b75cc152846147b404870f5a60f4709 |
| 41586_2024_7407_MOESM5_ESM.xlsx | 16776397 | 03080c0c9eb42fb953075dc118439e79e06d0796019caf259a009b77c3913343 |
| 41586_2024_7407_MOESM6_ESM.xlsx | 34625 | eca6ed81115366c5886c5a7041916bddc4b1d552a0ee645b974f8f17cf7e3a87 |
| 41586_2024_7407_MOESM7_ESM.xlsx | 214167 | b7f38b452720d0df763a9ec4c7186ab75f6ac0d5a43025202cb2b3f1207b3fbe |
| 41586_2024_7407_MOESM8_ESM.xlsx | 84025 | 841a359b6a8b674311fb545bf05ee46a8f4b7b2a97118f976a1f7e3466ac03c1 |
| 41586_2024_7407_MOESM9_ESM.xlsx | 118465 | b2817c7f1f39c96a020a7e89178bea25bc1f50e7de400c79d8788544d220a16b |

### compartments
- URL: https://download.jensenlab.org/human_compartment_{knowledge,experiments,integrated}_full.tsv
- Version/release: COMPARTMENTS current as of 2026-10-26
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 218.0 MB

| file | bytes | SHA-256 |
|---|---|---|
| human_compartment_experiments_full.tsv | 8631519 | 6ce5142afdebe226ff2979ff5dd7c5b98c8ce97cb318eb4ec1a80c3c2509e526 |
| human_compartment_integrated_full.tsv | 143799824 | 436e350e86940fdc43a6dc8655235ccd1243c70c201c4aeecfbafdf30b8770d6 |
| human_compartment_knowledge_full.tsv | 65587367 | d35f2d14aba6020e4b5084fa4502d4c3880d4f234ce328a8edab14bf22899e49 |

### cptac_pancan_v1.2
- URL: https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/<COHORT>/ (LinkedOmics mirror)
- Version/release: BCM-harmonised freeze v1.2, 10 cohorts
- Accessed: 2026-10-26; License: Open (CPTAC data-use policy); cite; total size: 1017.6 MB

| file | bytes | SHA-256 |
|---|---|---|
| BRCA/BRCA_meta.txt | 7369 | 16d00a3022fbafa9d29e74c0ee689ae63279911f04aeef2d21344ed51a0a5cb7 |
| BRCA/BRCA_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 50633530 | 64b2f52e9148747a0380b088892cd5d31cce6ac044c97e8bef88d599e084e190 |
| BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 21398351 | e584331bc77a440af50692c7874b736b3d135d01d96ce95d66af18cc64009b21 |
| CCRCC/CCRCC_meta.txt | 10979 | c332231f3cbd204094236dfd302977e70bbdddcf0049e862b60dd10244d04dae |
| CCRCC/CCRCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 36327869 | 0fb82e502a316db8ebd6950f536cc5ef6a6a0573d9d21e435f13b797f6cb43f0 |
| CCRCC/CCRCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 46018835 | e373b43f03df3ddb5c648f5af4b57e4caa61e767c911dfd5aa1b3ded1b772d98 |
| CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 13026061 | 1b00d684cd99d870bee199128485c22d6cb079e0a01de3c044b7a4fdb1735618 |
| CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16629110 | 090077c25f0eeb59b737912f39dec98901db949e1f143bbc481d5cd731550f60 |
| COAD/COAD_meta.txt | 11672 | 5242103f57fb5d919ab7e228dd36c97b669915d6bee88ae1e5a43c3aa270295e |
| COAD/COAD_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 31131857 | 84bddc6e7ac01c7e240959a23a104684cabe4a98b9062fd6f5811cfe98083451 |
| COAD/COAD_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 30662724 | dfa526d0ec65755195269148edd5c4f87079a32645ef8cdb78fac7cd07167c4b |
| COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 12515956 | a76a0cc0ab3986635e8745028eb0bcc64470035ab737cc0937ec7be114a09177 |
| COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 12293004 | fa2dde5f9165acdf6320d28e01dbad24b735ede320f5fab4a6ceab6c15aa0580 |
| GBM/GBM_meta.txt | 6809 | 6617ef2f9010acec1eba09a765ca3e35e4ba2ea7ba75ba9e8734ed9ad94feecc |
| GBM/GBM_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 65019693 | 6dfbe3a23e70fabb428b81b30ec3493a96cbb09ca7cf9c48e0d093f777f6c467 |
| GBM/GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 18033471 | d7393915eb8ec10c1c48bea03d78cd9ee2a3477f7bb4b93f808ffade03c2f01b |
| HNSCC/HNSCC_meta.txt | 12087 | 3c9e32ff8b51b4a22f967ca753185e01d79f69d860f9530f5587314bd06c80ac |
| HNSCC/HNSCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 32994426 | b23d8b000259af9ca39b6672cb32970ad67bcc6170e452d5ac5ec3cbd3de6dcf |
| HNSCC/HNSCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 54793177 | c793e056d97bb39dd11e0490f3ff44d8af9b2b5e368c609ca99912b38603e861 |
| HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 10542735 | 0a067886f9b553d3de24d2dbab21e18037c44a50a52b9d60004d22dba3442adc |
| HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 18184418 | 6d52dac50c55a4aa30e0aabaf3047a3a099875ea075b65d8f0bc83ac45a5c4b3 |
| LSCC/LSCC_meta.txt | 12051 | 6c334bbd8c9c4acfc94709a7fea5fd0181efa0881cee8ce1793085a5281a9d99 |
| LSCC/LSCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 56979342 | f1c91d8b7d09bb8115e70dcda7546d8265d2404a1c5a740e76c45ba4dca02ad1 |
| LSCC/LSCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 62042879 | e14fe8506376ed8f755b6066873756971013d6371092d1b7e9e7644f1e370077 |
| LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 18920447 | ea3f8e5a1e42f708af87f6cdfc86d3fa63a96bfdc4f6f3316ba54380ea579aaf |
| LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 20631271 | fd9f9cc0aed40942fd33d750d1d207bcdf9ea231fcddeff3cc8ac8fdf3e1384d |
| LUAD/LUAD_meta.txt | 12232 | f9facf713af130ad097c56794f1fbc9107334e02ad9bf0df3958b40d68b6af10 |
| LUAD/LUAD_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 53242041 | eec7e5570cd35ec8bdfa4a6d5cf44a6df9309d6f824023e86460d27fd2fb89b4 |
| LUAD/LUAD_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 57949691 | 5dbdb71093de80b5f8c4e8893c754ea7b3db065c461d5f10e424aa70559d1868 |
| LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 18491662 | 3cfc62b159276837dea644f457ebf6b212964320e0b3d8ac87168d82eb73747c |
| LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 20132337 | 746900f54072198481207ea13f8a7ba14dcabd52345bace0bef7438120c03845 |
| OV/OV_meta.txt | 6029 | ab6b3bb70460c682492a98f7620cbbeda61baab6e4d8bc1e25d19c0580390123 |
| OV/OV_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 8222958 | 16d2980e042d75e619f39f4672cff519911f8425af2ccf76b53a70b434928000 |
| OV/OV_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 29097272 | 509b6daee400feff78896eda5b9ce3ddba13ffb45e94f4165711834b0a0a943e |
| OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 3007610 | 6fc249b27c9284b598a94860c1e36f4fd45843f65b14c4ccdfbc0153725bb07a |
| OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 12506387 | 274c22512ab70e2a07197ca6e51d268038c90e903bac823c1b27b7950af4ae6f |
| PDAC/PDAC_meta.txt | 11358 | ce53f0acf92cbbbf165c8c27425c27b7b6de301c2e95279eaf124c37a2316d33 |
| PDAC/PDAC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 22119041 | 090f38557f48f8b615b077658385077b1a9e5625001547431db404ffee8d515e |
| PDAC/PDAC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 48706052 | be6834314f24daaebab9f56ad5736e314af87d7868115963e25b9e82aaf1a89f |
| PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 7177097 | beb1bc33ee52d10f0f7a6fa0bfc0de97916471c37d80c5cb655ee0d0a2452ce9 |
| PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16800925 | e1be4bed67bfdac6fff7b12487b8af9e0e2a0e93901039f5f9124e180f24f6bd |
| UCEC/UCEC_meta.txt | 13972 | 4058a8420129dd37f4d9035e1d87beb70a800ce9cb42754f74e7ad3d35f3020b |
| UCEC/UCEC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt | 14102105 | 45551e71b6eaed7ae8279b7c68d12d6d46bce0822f9580be0cd6164c4566df7f |
| UCEC/UCEC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt | 57212205 | 0cb6fe7b534b00018cffbc77a64d014f72571baa4b49ad4fc4009cb060c0c70f |
| UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 3354332 | 225df3e7827fbeec499deba3d773828935dbe722559730a902c66f9dabe1d229 |
| UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16560392 | aea5d956b3496530faf86052886de19e857e96294a2dbcb3adb9702457865d36 |

### cptac_pancancer_pdc
- URL: https://pdc.cancer.gov/pdc/cptac-pancancer via /pdcapi/file/signedURLFromUuid/<file_id> (open files)
- Version/release: CPTAC pan-cancer harmonized v1 (Li et al. 2023 Cancer Cell)
- Accessed: 2026-10-26; License: Open access (NCI PDC / CPTAC data-use policy); cite; total size: 4378.4 MB

| file | bytes | SHA-256 |
|---|---|---|
| Clinical_meta_data_v1.zip | 773971 | eb1b0397249f3f38d9594b60c2bf23a3e21b9495ae4a43d5c563f38f5e43c890 |
| Phosphoproteome_BCM_GENCODE_v34_harmonized_v1.zip | 289057948 | 57c3f2988defb02c2a911d37dc7586c20e3adf288e1eba2d8cb1e28510679d00 |
| Phosphoproteome_Broad_Institute_harmonized_v1.tsv.gz | 99914130 | 675da0584fb56c97d16339e7c597a45a2516caaa490d437aa51064960686014a |
| Phosphoproteome_UMich_GENCODE34_v1.zip | 1723125169 | 714e61ea805e2393fe40052d4513964152978fab0219fe4d16ce9c15c8ed987e |
| Phosphoproteome_UMich_SinaiPreprocessed_GENECODE34_v1.zip | 1439685853 | 1def96d7015c1a0b63c81fc6f48d6f17cf8f2538f9ba01ec2d123a99e2f6bf2a |
| Proteome_BCM_GENCODE_v34_harmonized_v1.zip | 121231401 | 4b14a6fe378bbd3d115b0cc7e89a6170c70ab6723d7f13654072bd186f827324 |
| Proteome_Broad_Institute_harmonized_v1.tsv.gz | 42879328 | 6ba578543d31d9e3c6ad9c02acd91e25b0502cc4d96f3bfacb1a5c5da767aa1b |
| Proteome_UMich_GENCODE34_v1.zip | 388514753 | db80eb170bcf56c94f27a2fe799cb009f77685cce43dac5a677cfd1258cfca4f |
| Proteome_UMich_SinaiPreprocessed_GENECODE34_v1.zip | 273250428 | 91dffd60966f16c7e8fbb4941afa8f3af918c4542cf643b1c29082b2e1838bf9 |
| pancancer_file_manifest.json | 12404 | 9f9ab922b4ba2f229ee7ab811e62aaa18ed70319c58d9821bad7ab7ada483e2f |

### depmap
- URL: https://api.figshare.com/v2/articles/31660582 (files 67214582, 67214579); depmap.org portal blocked by Cloudflare challenge
- Version/release: DepMap Public 26Q1 Chronos parameters: gene_effect.csv (= CRISPR gene effect), library_effect.csv
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 416.0 MB

| file | bytes | SHA-256 |
|---|---|---|
| DepMap_26Q1_chronos_gene_effect.csv | 413544937 | 5ea09d72b7789eb0fcf0dc81847cc71ef2b7dd8500216b0837085c3155b21607 |
| DepMap_26Q1_chronos_library_effect.csv | 2420581 | 5324c3e77292eb7e3280542770462489f4bcee8452a6b237259adcdc5a2f806f |

### go
- URL: https://current.geneontology.org/annotations/goa_human.gaf.gz ; https://purl.obolibrary.org/obo/go/go-basic.obo
- Version/release: GO current as of 2026-10-26
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 47.3 MB

| file | bytes | SHA-256 |
|---|---|---|
| go-basic.obo | 32227785 | b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc |
| goa_human.gaf.gz | 15041996 | db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30 |

### hpa
- URL: https://www.proteinatlas.org/download/ (tsv/subcellular_location, rna_tissue_consensus, rna_celline; proteinatlas.tsv.zip)
- Version/release: HPA current release as of 2026-10-26 (version not stated in file)
- Accessed: 2026-10-26; License: CC BY-SA 4.0; total size: 218.9 MB

| file | bytes | SHA-256 |
|---|---|---|
| proteinatlas.tsv.zip | 7460135 | ad401f8519ecdee5b67b0ad3bf28175d402d1ea0ad6c542b90fcab2bee39b470 |
| rna_celline.tsv.zip | 205861141 | e5b11779308c2720cfe89d9cb19932ad4c0c3fea0fa37e8c61a14d501995852a |
| rna_tissue_consensus.tsv.zip | 5293680 | b001903da23f802c8dfbd1b6b3a8c20b831fedd4932f4bb721b06622d83ae60b |
| subcellular_location.tsv.zip | 252331 | a9466882bd3328416bcde77c32a6bcb2a17d70183af90e98d032f6f5f781d504 |

### humancellmap
- URL: https://humancellmap.org/resources/downloads/
- Version/release: Human Cell Map "latest" (2026-10-26)
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 81.9 MB

| file | bytes | SHA-256 |
|---|---|---|
| baits-latest.txt | 21408 | c3bb253a0f59e61e76206585bb4ce79b4450196bfe889e4264c3de1fe3fbc053 |
| preys-latest.txt | 441953 | e38668873c7811336c8c50769366f36a5be68ea845440c61bacc6416ddd1f18c |
| saint-latest.txt | 81457216 | 0db37da51b94faa478ad41509f8fdc649517759c92954946669df89ab59c29da |

### huri
- URL: https://interactome-atlas.org/data/{HuRI,Lit-BM,HI-union}.tsv (ENSG pairs; www. host has a TLS cert mismatch, apex host used)
- Version/release: HuRI as of 2026-10-26
- Accessed: 2026-10-26; License: CC BY 4.0 (check site); total size: 4.2 MB

| file | bytes | SHA-256 |
|---|---|---|
| HI-union.tsv | 2112198 | 01ae99dac8deacfd34be74f37574062506068914b85ad2b83320028610c3c3fd |
| HuRI.tsv | 1681536 | 4cdf94366ff6e23e71ef9ec81922e85b2b3a1fb08e938f7b00d86c027cf1ecd6 |
| Lit-BM.tsv | 430112 | f39c9c4b8a372286c3eca77d0ff651d837a643d32235eae866daae082663391d |

### intact
- URL: https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
- Version/release: IntAct release of 2026-01-14 (MITAB 2.7, all species, 11 GB unzipped)
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 1350.8 MB

| file | bytes | SHA-256 |
|---|---|---|
| intact.zip | 1350789771 | 8409c87570d9dc9d3bc3676615e29146b2d8b89bfdbf2d30efa5b46484dc7a79 |

### iptmnet
- URL: https://research.bioinformatics.udel.edu/iptmnet_data/files/current/
- Version/release: iPTMnet 6.2
- Accessed: 2026-10-26; License: CC BY-NC-SA 4.0; total size: 79.6 MB

| file | bytes | SHA-256 |
|---|---|---|
| protein.txt | 9324847 | 381b90d7a9d5b1653c33bc3076db62ba9ffb3130c89f453276f310b2413e9309 |
| ptm.txt | 44116546 | d5b2ed7c7138bc77f6f712feb038e18ddc4c791e96823b286a680e1982eff20f |
| readme.txt | 1892 | 03526d7031f5270804b26605a67325a6e0a979d30fc62fbf6b896e0e2760ef2c |
| score.txt | 26110370 | e7dfa494099ab9d8cdff54c727aba65251067a041b931846d8656c622563fe63 |

### kinase_library
- URL: https://pypi.org/project/kinase-library/1.8.0/ (sdist; extracted as data only into src/, never installed/executed)
- Version/release: kinase-library 1.8.0
- Accessed: 2026-10-26; License: CC BY-NC-SA 3.0; total size: 70.0 MB

| file | bytes | SHA-256 |
|---|---|---|
| kinase_library-1.8.0.tar.gz | 70027403 | 0fa1929b1c2c561231608abeace01f803afc012b7ff4c5affcb2623b9884de6d |

### kinbase_manning2002
- URL: http://kinase.com/human/kinome/tables/ (Manning et al. 2002 Science SOM Table S1/S2)
- Version/release: Manning 2002
- Accessed: 2026-10-26; License: Supplementary data of published paper; cite; total size: 5.4 MB

| file | bytes | SHA-256 |
|---|---|---|
| TableS1.txt | 2607250 | 9788a469f28f53bb8ea94146209733ab718b3670287edb69ea2342c69542080a |
| TableS1.xls | 2733056 | a5673056c47aae3e985dea867648d5d1c40e56e29470738f6121a3961bf8d46b |
| TableS2.txt | 81184 | eed45ad222904371731a052aa830e3306d0eae898277b83ae7213a2e2526c17f |

### kinhub
- URL: http://www.kinhub.org/kinases.html (HTML table, to be parsed)
- Version/release: KinHub as of 2026-10-26
- Accessed: 2026-10-26; License: Not stated; academic use; total size: 0.2 MB

| file | bytes | SHA-256 |
|---|---|---|
| kinhub_kinases.html | 155599 | c81227b4ac53c732b508f3c4f1d3d55fe491ac3b873516f7f4eeea42dc09b1c6 |

### ochoa2020
- URL: https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0344-3/MediaObjects/ (MOESM3-7; MOESM5 = functional_score, 116,258 sites; MOESM4 = annotated_phosphoproteome features)
- Version/release: Ochoa et al. 2020 Nat Biotechnol suppl.
- Accessed: 2026-10-26; License: Publisher supplementary data; cite; total size: 59.1 MB

| file | bytes | SHA-256 |
|---|---|---|
| 41587_2019_344_MOESM3_ESM.xlsx | 190938 | 7557c0d661b221eacf5d0a88265ebbd774e40e3bbaa1989ea6d14ab4e4f15748 |
| 41587_2019_344_MOESM4_ESM.xlsx | 54465047 | bab244a8e48a5e52bf605efda9df6f1f5ef6df8b1cbfd9081e71e7de8ae32eb5 |
| 41587_2019_344_MOESM5_ESM.xlsx | 3556910 | a0b5ac3b096094726ec633f8e9787dc053087e03b1286ebec2e3fdaabb8d9df7 |
| 41587_2019_344_MOESM6_ESM.xlsx | 521922 | 5cbce7c21d1aed2fcc198d68f13b462edc9adb9245f1ba1d9df2f27218f6afae |
| 41587_2019_344_MOESM7_ESM.xlsx | 411932 | 96a44f46f9a4d881f30bea9afa5f3dfc53b95c34e41b8fdb73867c03c223da43 |

### omnipath
- URL: https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references&format=tsv&license=academic
- Version/release: OmniPath web service as of 2026-10-26
- Accessed: 2026-10-26; License: Mixed per-resource; academic licence set (includes non-commercial resources e.g. PSP); total size: 5.0 MB

| file | bytes | SHA-256 |
|---|---|---|
| omnipath_enzsub.tsv | 4974142 | 57fa74706ef9c54ad23fe12e7390ba11a0694a8168534416d183fedc070e4f5d |

### opencell
- URL: https://opencell.czbiohub.org/data/datasets/
- Version/release: OpenCell as of 2026-10-26
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 7.3 MB

| file | bytes | SHA-256 |
|---|---|---|
| opencell-library-metadata.csv | 647376 | 92da4a93cee640a50c9745b1c7da31e9316f2014400dd0feacbe4f541b7c82e6 |
| opencell-localization-annotations.csv | 61612 | 14cf3536d9d3699460dd34e166d001a34e5430fb551488538b9468107f9f10b6 |
| opencell-protein-abundance.csv | 2046131 | 1de0a35480601a6fce0bb46c4a348bb70358d02586c1c905c48f5d837901cdbc |
| opencell-protein-interactions.csv | 4568928 | a9422ca12b8b5fce52439f247a7439fbc693cb34b4ca45ced597aa02d267e376 |

### string
- URL: https://stringdb-downloads.org/download/ (protein.*.v12.0, 9606)
- Version/release: STRING v12.0
- Accessed: 2026-10-26; License: CC BY 4.0; total size: 172.5 MB

| file | bytes | SHA-256 |
|---|---|---|
| 9606.protein.aliases.v12.0.txt.gz | 19777800 | b65f730b993ed0c1bd72edf4565d3d425db42861101b29699704810e8f125680 |
| 9606.protein.info.v12.0.txt.gz | 1970090 | 144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c |
| 9606.protein.links.detailed.v12.0.txt.gz | 139634808 | 28a2ea3b22c964e691c0f2f934c5e3e397a4b0a12c833d8787fa02b7847b91b6 |
| 9606.protein.physical.links.detailed.v12.0.txt.gz | 11086354 | db1de86dc7a6e34eba78ca6e6747f9201195b0af1edb607c19bf651b0eddede5 |

