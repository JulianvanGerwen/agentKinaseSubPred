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

# Data acquisition, wave 2 (access date 2026-10-26)

Per-file URL/version/licence/size/SHA-256: /workspace/DATA.md section "wave 2" (raw hashes logs/dl/sha256_wave2.txt; git commits logs/dl/wave2_git_commits.txt). Scripts: scripts/dl_w2.sh, dl_w2b.sh, make_data_rows_w2.py. Nothing downloaded was executed; git clones/sdist are data only. Useful routes that worked: Europe PMC `.../rest/<PMCID>/supplementaryFiles` (open-access papers only), `media.springernature.com/original/springer-static/esm/art%3A<doi>/MediaObjects/<file>` (links scraped from the nature.com article page), Zenodo via `zenodo.org/api/records?q=recid:<id>` (metadata) and `.../api/records/<id>/files/<name>/content` (data; Zenodo rate-limits bursts with 403 "unusual traffic", wait several minutes).

## Obtained (dirs under /workspace/external/)
| # | Item | Dir | Notes |
|---|---|---|---|
| 1 | Buljan 2020 Mol Cell | buljan2020/ | supp zip (PMC7427327): mmc3 = Table S2 final AP-MS interactome (Bait_id, Protein_id, wdscore, gfpratio, spectral counts; ~7.3k rows), mmc4 S3 co-purified pairs (1.7k), mmc7 S6 predicted regulatory phosphosites (NetPhorest), mmc2 S1 kinase classification, mmc8-10 disease/neighbourhood/exp info |
| 2 | Davis 2011 | davis2011/ | MOESM2 PDF, MOESM3-6 .xls (OLE2). **xlrd not installed** in either env, so not parsed; ask PI to approve `xlrd` install (or I convert). Also KIBA bundle (kiba_zenodo_5105698/KiBA.xlsx, 73 MB, Davis/Metz/Klaeger-comparison) |
| 2 | Reinecke 2024 Nat Chem Biol | reinecke2024/ | MOESM2-6 xlsx (kinobeads target landscape of ~1,000 inhibitors; MOESM6 14 MB probably the main matrix), SI PDFs |
| 2 | Klaeger 2017 (partial substitute) | curvecurator_zenodo_8399823/ | Kinobeads.zip: CurveCurator reprocessing of kinobeads dose-response (all_binder_data_lfq_annotated.txt, curves.txt). Not the original Science tables |
| 3 | Sugiyama 2019 Sci Rep | sugiyama2019/ | Table S2 (MOESM3): 198,536 kinase-substrate-site rows with SIDIC/PTMscore confirmation flags, substrate by UniProt entry name + position; S1 peptides (25k), S3 motifs, S4/S5 per-kinase S/T/Y counts (386 kinase preps). iKiP-DB (Zenodo 21469258: Table7 GMT, Table6 activation-loop sites, Table5 C-score; Zenodo 20037245 kinase scores) also fetched in ikipdb_zenodo_*/ |
| 4 | Phosformer | phosformer/ | hf_phosformer (185 MB weights), hf_space_phosformerST (Phosformer-ST, ESM2-150M multitask, 602 MB + code), hf_dataset_curated (corpus kinases.fa 109 MB, substrates.fa, curated 11-mer KS tsv with PSP/PhosphoELM/PhosphoNetworks source flags, splits train 715 MB / val / test_pos_neg / shifted / neg_only), github/ (code, reference_human_kinases.csv). **Licence CC BY-NC-ND 4.0 (non-commercial, no derivatives)**; training data derives from PSP/Phospho.ELM so splits are NOT leakage-safe vs our PSP-based benchmark |
| 5 | SELPHI2.0 | selphi2/ | selphi_2.0 GitLab clone (Apache-2.0, 4.4 MB, README + src; commit f9822f9), Zenodo 13821321 small files (prediction_matrix_high_conf_curated.csv 87 MB, positive-set, kinase-citation tables, kinase-kinase rels). Big tables not fetched (8.3 GB matrix, 10.7 GB feature table, two 4 GB PSSM score tables). Also orphan_top5_zenodo_16881580 (top-5 kinase predictions for orphan phosphosites, 12 MB, author unverified), kinasecanvas_zenodo_22941318 (KinaseCanvas 1.5 GB, new 2026 KS resource; contains curated edges -> leakage risk) |
| 6 | KolossuS | kolossus/ | PyPI 1.0.6 sdist + GitHub clone + HF weights (aparekh2/kolossus; MIT): large (100 MB, expects ESM2-15B embeddings) and small (15.5 MB, ESM2-650M). ESM2 backbones not fetched. Not installed. **GPS 6.0**: web server http://gps.biocuckoo.cn/ with a DOWNLOAD page and GitHub link (source code); not fetched |
| 7 | KinCore | kincore/ | Human-PK.fasta (497 kinase-domain seqs, header `GROUP_GENE/start-end ENTRY GENE ACC`), Human-PK-alignment.{aln,fasta}, PK_labels_PDB.fasta (13,161 PDB chains labelled) |
| 9 | Hein 2025 Cell | hein2025_zenodo_14000085/ | analysis repo zip (133 MB): cluster_annotation_Dec6.csv, enrichment/volcano tables (per-protein), external refs. The paper's final localisation supplementary table (cell.com) is 403 |
| 10 | HGNC | hgnc/ | hgnc_complete_set.txt (17 MB, CC0) |
| 11 | ProteomeHD | proteomehd_kustatscher2019/ | Nat Methods supplement (MOESM3/4 .txt = likely co-regulation scores, MOESM5-9 xlsx); full ProteomeHD matrix at proteomehd.net not tried (host unreachable earlier test not completed) |
| 12 | KinoPlex Nat Biotech 2026 | kinoplex_natbiotech2026/, kinoplex_zenodo_20633591/ | Supplement MOESM4/5/7/8 xlsx (17-30 MB: prediction/atlas tables; not yet inspected), plus Zenodo structural atlas feather (all canonical STY x structural features, 1.4 GB). ST bootstrapped-quantile table (3.5 GB, Cantley-PSSM percentiles) and CSV duplicate not fetched |

## Blocked / not obtained
- **Klaeger 2017 Science tables**: science.org 403; PMC6542668 author-manuscript supplements are behind a reCAPTCHA/proof-of-work page (not bypassed); ProteomicsDB OData exists (`/proteomicsdb/logic/drugSearchEndpoint.xsodata`, `getCurveInformationByProteinID.xsjs?protein_id=&drug_name=&assay_type=PDB:xxxx&assay_variable=dose`) but I did not identify the kinobeads assay-type id within the time box. Substitutes: CurveCurator Kinobeads.zip, Reinecke 2024 (superset of Klaeger drugs), KIBA. Request: supply Science Tables S1-S6 if the original is needed.
- **Drug Target Commons**: drugtargetcommons.fimm.fi unreachable (no connection); Zenodo mirror 18450948 (DTC_DrugTargetInteractions.RDa, 87 MB, CC BY 4.0) returned 403 repeatedly (rate limit); retry later with `scripts/dl_w2b.sh`-style get.
- **decryptM (Zecha 2023)**: Science supplement 403; PMC7615311 supplements behind the same challenge; decryptM 2.0 Zenodo 17533475 is access-restricted. Not obtained.
- **SIGNOR**: signor.uniroma2.it still unreachable (000/503). No Zenodo dump, not in Pathway Commons v12. Only via OmniPath (already have).
- **Phospho.ELM** still awaiting PI licence decision (wave 1).
- ProteomicsDB, Cell Press supplements: bot-gated.

## Notes
- Disk: external/ now ~18 GB (wave 2 adds ~3.6 GB; largest KinaseCanvas 1.5 GB, KinoPlex feather 1.4 GB).
- Excel/OLE files not opened beyond header peeks (openpyxl via `python -I`).
