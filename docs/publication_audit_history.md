# Publication audit: full git history

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Independent read-only audit, run 2026-10-28 10:18:32 to 10:33:21 CET (about 15 min, from `date`).
Scope: all 108 commits on the single branch `master` (no tags, no stash), 3,030 unique blobs and 3,091 distinct paths.
3,075 files are in HEAD. 16 paths exist only in history.
Method:
- Every unique blob was streamed with `git cat-file --batch`, with gzip blobs decompressed, and regex-scanned.
- Restricted inputs in `data/` and `external/` were fingerprinted in memory and compared against every blob:
  - PSP kinase–site pairs, 15-mer windows, Regulatory_sites free text and the licence header;
  - Klaeger, Cantley and kinase-library numbers;
  - CPTAC values and sample IDs;
  - compendium log2FC values, peptide UIDs and condition metadata;
  - text of the publisher PDFs, and Phosformer, SELPHI2 and kinase_library code (10-word shingles).
- No repository state was changed. Only this file was written.
- Nothing from `sealed/` was opened except the tracked `sealed/kinase_label_counts.tsv` blob.
- No restricted rows are reproduced here.

## Prioritised findings

| # | Sev. | What | Where (commit / path) | HEAD or history-only | Recommended fix (least destructive) |
|---|---|---|---|---|---|
| 1 | **HIGH, conditional** | Condition-level metadata of the user-provided phosphoproteomics compendium (`data/01`) copied verbatim as column subsets. All 1,937 conditions: conditionID, Condition, Control, Sample, Labeling, Timepoint, PMID, source. 1,437 conditions are labelled `source = "This study"`. Provenance and licence are unknown (DATAINEED.md item 5). No quantitative log2FC values or peptide UIDs are present (0 matches). | Since commit #3 `3c3fc9c`. `results/data_profile/{compendium_per_condition_coverage, kinhib_hijazi_conditions, kinhib_condition_matches, compendium_conditions_by_pmid}.tsv`; `results/features/phospho/{conditions_meta, f6b_conditions, f6b_conditions_clean_bench, f6b_conditions_f6_clean, f6b_unmapped_or_dropped_drug_conditions, label_source_studies, benchmarkin_studies}.tsv`; plus `f6b_conditions_lam10/lam100.tsv` (history-only). | HEAD (11 paths) + history (2) | **The human must confirm the compendium is published.** If it is published or cleared: (iv) nothing needed, cite it in DATA.md. If it is unpublished or confidential: **(iii) history rewrite of these 13 paths is REQUIRED** before any public push. `git rm` at HEAD is not enough. |
| 2 | MEDIUM | Internal NFS server IP and lab share export path (`<internal-file-server>:<share>` (an RFC 1918 address and a lab share path)). This is the only internal-infrastructure string found. | `docs/reviews/r1_ai_failure_modes.md` line 39. Added in commit #6 `2109b65`; the same blob is in 103 commits. | HEAD + history | Edit at HEAD (ii). Removing it from history needs (iii) `git filter-repo --replace-text` on this one string, which changes the hashes of commits #6–108. This is **recommended but optional** (RFC1918 address, not routable). It is **required** if ETH or lab IT policy forbids publishing internal addresses. |
| 3 | MEDIUM (licence) | **(d) Sealed TEST labels at pair level.** 12 files `results/test/721dc1977366/*_site_level.tsv`, 2,205 rows each, the same TEST-positive (kinase_acc, site_id) set ranked per method. 2,144 of the pairs are PSP Kinase_Substrate pairs (category c). | Commit #85 `37a3e6a` | HEAD + history | Licence-permitted (PSP CC BY-NC-SA 3.0) with attribution + NC + SA. Fix: add a top-level LICENSE/NOTICE at HEAD (ii). A rewrite is **not required**: `results/splits/test_proteins.tsv` plus PSP rebuilds this set anyway, and the TEST set was scored once and is now open. (iii) is only needed if the human wants the TEST set to stay unpublished. `release_repo/.gitignore` already blocks `*_site_level.tsv`, but the project root does not. |
| 4 | MEDIUM (licence) | **(c) PSP-derived pair tables:** `config/bio_controls.tsv` (44 literature pairs; 42 in PSP; it keeps the `in_psp`/`in_ext`/`in_test_split` columns that were stripped in `release_repo`, and flags 8 TEST pairs). `results/strata/dev_hub_sites.tsv` (88 DEV hub sites with their kinase lists). | #10 `8020e49` (2 versions); #63 `07cff92` | HEAD + history | Licence-permitted: add attribution (ii). Optionally copy the stripped `release_repo/config/bio_controls.tsv` over the root copy at HEAD (ii). No rewrite needed. |
| 5 | MEDIUM (licence) | **(b) Model predictions** `results/release/site_top5.tsv.gz` (1.75 M rows) and `kinase_top200.tsv.gz` (100.6 k rows), also present under `release_repo/predictions/`. They are trained on PSP. They contain 8,506 and 3,767 (kinase, site) pairs that also appear in PSP, but as predictions, not labels. | #93 `4249203`; release_repo in #108 | HEAD + history | Licence notice: `LICENSE-PREDICTIONS.md` lives only in `release_repo/`, which will be removed. Move it next to `results/release/` or to the root (ii). No rewrite. |
| 6 | LOW | **(b)/(d) Sealed per-kinase label counts** `sealed/kinase_label_counts.tsv` (504 rows of aggregate counts including TEST). It was force-added although `sealed/` is gitignored; an earlier copy at `results/labels/kinase_label_counts.tsv` is history-only. | #7 `6225822`, moved in #11 `e863776` | HEAD + history | (iv) Aggregates only. |
| 7 | LOW | **(b) Derived aggregates** of PSP, compendium, Klaeger, CPTAC, kinase library and benchmarKIN: about 47 files in `results/data_profile/`; `results/labels/*` (kinase list built from UniProt/KinHub, synonyms, counts); `results/splits/*` (protein → cluster/fold, per-protein label counts); `results/sites/*`, `results/leakage/*`, `results/strata/*`, `results/features/**`, `results/downstream/**`, `results/models/**`; per-kinase metrics in `results/test/`. | many | HEAD + history | (iv) |
| 8 | LOW | Cluster host name `<gpu-node>` (an HPC GPU node), plus `pid` and `argv`, in every run manifest. | 402 `results/runs/*/manifest.json` since #16 `bee7b42`; 93 commits | HEAD only (all 402 blobs are in HEAD) | (iv), or (i) going forward: drop `host`/`pid` in the manifest writer (`src/ksp`). A rewrite is not worth the hash change. |
| 9 | LOW | Claude Code scratchpad paths `<scratchpad>/...`. They carry a numeric id and the local session directory UUID. | 22 files (10 `docs/audit/`, 7 `docs/reviews/`, 3 `results/runs/`, `docs/lit/`, `manuscript/BRIEF_COMMON.md`); 28 matches; since #3 | HEAD | (iv), or (ii) cosmetic. |
| 10 | LOW | Absolute paths. `/workspace/...`: 1,291 matches in 282 blobs (255 in HEAD, 201 of them in `results/runs/`). `/workspace/.venv/...`: 123 matches in 56 blobs, from tracebacks. Three paths under `/home` (a shared Python env, the MMseqs2 install and the Hugging Face cache): 28 matches in 19 blobs. These are shared environment paths, not a personal home directory. | since #1 | HEAD + history | (iv) |
| 11 | INFO | **Author identity** `<author name> <author e-mail>` (the researcher's institutional address), as both author and committer, on all 108 commits. The unix user name equals the author name. The unix group name does not appear anywhere. | all commits (metadata) | all | Expected attribution: (iv). Hiding it needs (iii) `filter-repo --mailmap` on every commit, which changes every hash. `.mailmap` alone does not change the objects. |
| 12 | INFO | **`Claude-Session: <session URL>`** trailer in **93 of 108** commits (the same URL each time). Commits #1, 11, 13, 18, 21, 27–29, 32, 34, 36, 37, 40, 41 and 45 lack it. It also appears in no blob. | commit messages | all | **Not a credential.** The link opens only for the owning account and anyone the transcript is shared with. It identifies one session and reveals nothing else. Keeping it matches the goal of making the agentic work visible: (iv). Removing it would need a message rewrite of 93 commits (iii), which is optional. |
| 13 | INFO | Size. Blobs over 1 MB, all in HEAD: `site_top5.tsv.gz` 34.75 MB (one blob at 2 paths), `kinase_top200.tsv.gz` 1.95 MB (2 paths), `results/strata/cv_strata_site_per_kinase.tsv` 1.63 MB, `results/downstream/validation_per_kinase.tsv` 1.44 MB, `final_r2_nophospho_zscore+ulm_per_kinase.tsv` 1.33 MB, `r2_clean_per_kinase.tsv` 1.22 MB. | — | HEAD | (iv). Under GitHub's 50/100 MB limits. Optionally host the predictions on Zenodo. |
| 14 | none | **Secrets.** 0 hits in all history for: `AKIA…`, `ghp_`, `github_pat_`, `sk-`, `sk-ant-`, `hf_…`, `xox?-`, `BEGIN … PRIVATE KEY`, `password=/:`, `passwd`, `token=/:`, `api_key=`, `secret=`, `Authorization:`, `Bearer …`, and signed-URL parameters. There are no `.env`, `.netrc`, ssh-key, `.pgpass` or AWS-credential files (`.env` appears only as a `.gitignore` entry and in the old audit text). High-entropy token scan: 0. Email addresses in blobs: 0. The author name and the institution domain in blobs: 0. Commit messages hold no paths, hosts or IPs. | — | — | (iv) |
| 15 | none | **Verbatim restricted data.** 0 PSP rows: no PSP-style lowercase-centre windows, and 0 Regulatory_sites free-text shingles. The PSP licence and attribution text is quoted only in the `src/ksp/io_psp.py` docstring, which is appropriate. Klaeger, CPTAC (values and sample IDs), compendium log2FC, kinase-library matrices, and Cantley Johnson 2023 / Yaron-Barir 2024 tables: 0 exact copies (at most 2 coincidental 8-decimal numbers per file; 0 matrix rows; 0 matrix headers). No publisher PDF, xlsx, parquet, npy or pt file has ever been committed. The only PDFs and PNGs are the project's own Matplotlib figures, whose metadata says Matplotlib only. Publisher PDF text shows up only as citations and short quotes (at most 18 shingles, in `manuscript/`). The Phosformer (CC BY-NC-ND) baseline is a clean reimplementation: 15 shingles in `src/ksp/baselines/phosformer/model.py`, all from config lists. `results/runs/*/` holds no prediction or label tables: only per_kinase, importance, train_stats, strata, bootstrap and manifest files. `logs/` holds only checksum lists, short `.out` logs and shell chains. | — | — | (iv) |

Sizes:
- Pack on disk: 12.73 MiB packed plus 42.78 MiB loose. A full push pack is about 53.0 MB, and a HEAD-only pack is the same size, so history adds almost nothing.
- HEAD working tree: 152.3 MB in 3,075 files (`results/` 104.2 MB, `release_repo/` 42.7 MB; the release_repo predictions share blobs with `results/release`).
- All unique blobs in history: 113.5 MB uncompressed.

**Cost of any rewrite (iii).** Hashes are pinned in tracked files:
- `FREEZE.json` `harness_commit` = `974d520` (commit #38);
- 406 manifests carry `git_hash` (31 distinct commits);
- `results/test/LEDGER.md`.

Rewriting from commit #3 or #6 would make every one of those references dangling. If you rewrite, publish filter-repo's `commit-map` (old→new) as a tracked file, so the pre-registration provenance can still be checked.

## Verdict

**Publish after HEAD cleanup.** That means:
- remove `release_repo/`;
- add a root LICENSE/NOTICE covering the PSP-derived content (CC BY-NC-SA 3.0 attribution, NC, SA), including next to `results/release/`;
- optionally strip the extra columns from `config/bio_controls.tsv`, and add `*_site_level.tsv` and `sealed/` guards to the root `.gitignore`;
- optionally scrub the IP from `docs/reviews/r1_ai_failure_modes.md` at HEAD.

**History rewrite is required only for X = the 13 compendium condition-metadata tables (finding 1), and only if the human cannot confirm that the compendium is published or cleared for release.**

A rewrite is optional and recommended for the NFS IP/share string (finding 2). It is optional for the author email (finding 11) and the session trailers (finding 12). Neither of those is a credential.
