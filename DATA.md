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
