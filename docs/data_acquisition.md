# Data acquisition, wave 1 (access date 2026-10-26)

Per-file URL, version, licence, size and SHA-256 are in /workspace/DATA.md (section "Downloaded (external/)"). Raw hashes: /workspace/logs/dl/sha256_wave1.txt. Scripts: /workspace/scripts/dl.sh (generic curl), dl_pdc_pancancer.py (PDC pan-cancer), make_data_rows.py (DATA.md rows). Nothing downloaded was executed; only the kinase-library sdist was untarred (as data) to external/kinase_library/src/. Excel files were inspected via stdlib zip/XML only.

## Obtained (path under /workspace/external/; key identifiers; quick counts)
| Item | Dir | Contents / IDs / rows |
|---|---|---|
| 1 Cantley Ser/Thr (Johnson 2023) | cantley_johnson2023/ | MOESM3 (T1 kinase info), MOESM4 (T2: raw/norm/norm-scaled PSSMs, 303 S/T kinases), MOESM5 (T3: precomputed phosphosite scores/percentiles, 370 MB xlsx, ~89,785 rows x ~630 cols, Ochoa sites), MOESM6 (T4), 2 PDFs |
| 1 Cantley Tyr (Yaron-Barir 2024) | cantley_yaronbarir2024/ | MOESM4 (T2 Tyr PSSMs, 93 kinases), MOESM5 (T3 annotation), MOESM7 SH2 matrices, MOESM8 C. elegans, MOESM3/6/9 other tables; MOESM10-12 returned 403 (probably do not exist) |
| 1 kinase-library pkg | kinase_library/ | v1.8.0 sdist + extracted src: per-kinase TSV matrices (ser_thr 311 x {raw,norm,log2,densitometry}; tyrosine 93), combined matrices TSVs, kinome_information.tsv (group/family), st_favorability.tsv, Ochoa phosphoproteome, precomputed scored phosphoproteome parquet (ser_thr 63 MB, tyrosine 2 MB), PSP KS count file (07_2021). No separate "background distribution" file found; percentiles come from the scored parquets (verify when parsing) |
| 2 STRING v12 | string/ | physical detailed (1,477,610 rows), full detailed (13,715,404), info (19,699), aliases (3,889,207); ENSP IDs, map via aliases to UniProt |
| 3 BioGRID 5.0.262 | biogrid/ | ORGANISM tab3 zip (98 species files; use Homo_sapiens), PTMS zip (ptmtab 1.2 GB unzipped + ptmrel) |
| 4 IntAct | intact/ | intact.zip 1.35 GB (MITAB, all species, 11 GB unzipped; filter taxid 9606) |
| 5 BioPlex 3.0, HuRI | bioplex/, huri/ | 293T 118k, HCT116 71k, combined v4a 56k rows; HuRI 52.5k, Lit-BM 13.4k, HI-union 64k pairs (Ensembl ENSG) |
| 6 HPA | hpa/ | subcellular_location, rna_tissue_consensus, rna_celline (1.07 GB unzipped), proteinatlas.tsv.zip (43 MB unzipped) |
| 7 OpenCell | opencell/ | interactions (29.9k rows), localization annotations (1.3k), library metadata, abundance |
| 8 Human Cell Map | humancellmap/ | preys (4.4k), baits (193), SAINT (223k rows). Not fetched: not-significant, crapome matrix, fasta, supplementary-prediction zips (listed at humancellmap.org/templates/download.html) |
| 9 COMPARTMENTS | compartments/ | knowledge (918k), experiments (92k), integrated (2.3M) full files; no textmining |
| 10 KS annotation | omnipath/, iptmnet/ | OmniPath enzsub (41.5k rows, academic licence set, includes PSP-derived records: not independent of PSP labels); iPTMnet 6.2 ptm (438k), score (851k), protein (68k) |
| 11 GO | go/ | goa_human.gaf.gz (906k lines), go-basic.obo |
| 12 DepMap | depmap/ | DepMap Public 26Q1 Chronos gene_effect.csv (413 MB; 1,208 rows incl. header) + library_effect via figshare 31660582. NOTE: this is the Chronos-parameters deposit; I assume gene_effect.csv equals the portal CRISPRGeneEffect (not verified; no model metadata file) |
| 13 Kinase classification | kinhub/, kinbase_manning2002/ | KinHub HTML (to parse; UniProt links), Manning 2002 Table S1 (Group/Family/Subfamily, 638 kinases) and S2. Also uniprot/pkinfam from other agent |
| 14 CPTAC | cptac_pancan_v1.2/, cptac_pancancer_pdc/ | (a) LinkedOmics S3 BCM freeze v1.2: 10 cohorts x phospho site + proteome gene (Tumor and, where exist, Normal) + meta; BRCA/GBM have no Normal files (403). e.g. BRCA tumour phospho 49.9k lines, OV 37.3k. (b) PDC pan-cancer harmonized v1: Phosphoproteome BCM/UMich/UMich-Sinai/Broad and Proteome BCM/UMich/UMich-Sinai/Broad, Clinical_meta (4.1 GB; UMich zips contain .DS_Store/__MACOSX junk) + manifest of all 30 pan-cancer files (RNA, CNV, mutation, acetylome etc. not fetched) |
| 15 Ochoa 2020 | ochoa2020/ | MOESM5 functional_score (116,258 rows), MOESM4 annotated_phosphoproteome + feature tables (116,258 x 89 cols), MOESM3/6/7 |
| PI add: benchmarKIN | benchmarkin/ | Zenodo 12566560 (all 14 files incl. NetworKIN 3.1 predictions, Hernandez/Hijazi benchmarks, 7 CPTAC tables, gold standard .Rds -- do not load Rds with untrusted-code risk; use pyreadr/rdata only) |

## Requests to user / blocked
- **Phospho.ELM**: dump is behind an academic-licence click-through (non-commercial, no redistribution to third parties; http://phospho.elm.eu.org/dumps/Phospho.Elm_AcademicLicense.pdf). The file URL responds, but I did not accept a licence on your behalf. Please confirm acceptance; then download http://phospho.elm.eu.org/dumps/phosphoELM_all_latest.dump.tgz.
- **SIGNOR**: signor.uniroma2.it returned proxy 503 / connection failure from this host (tried http/https, downloads.php, API getHumanData.php). Retry later or supply the file manually. SIGNOR-derived KS edges partly present in OmniPath enzsub.
- **DepMap portal (depmap.org)**: Cloudflare Turnstile challenge; not bypassed. Used the figshare 26Q1 Chronos deposit instead (see caveat). Model.csv / lineage metadata not obtained.
- **Klaeger 2017 Kinobeads tables**: science.org supplementary returned 403 (paywall/bot block). PRIDE PXD005336 holds only per-drug supplement CSVs/PDFs (partly), raw files huge. Please supply Science Tables S1-S6, or authorise a ProteomicsDB route. Not fetched.
- Johnson 2023 background/score-distribution files were not found as a separate item; check scored parquets and MOESM5.
- Not attempted (not in wave-1 scope): decryptM curves (per PI), cptac python package, Reinecke 2024, Drug Target Commons.

## Notes / assumptions
- HuRI host www.interactome-atlas.org has a TLS certificate mismatch; the apex host interactome-atlas.org (valid cert) was used, verification not disabled.
- BioGRID "latest" resolved to 5.0.262; IntAct "current" = 2026-01-14 build; STRING v12.0 file names as in brief.
- HPA/COMPARTMENTS/GO/OmniPath/HCM are rolling releases, so SHA-256 pins the exact state obtained.
- Disk: ~9.5 GB in external/ for this wave (largest: PDC 4.1 GB, IntAct 1.3 GB, S3 CPTAC 1.0 GB).
