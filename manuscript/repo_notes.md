# Release repository: privacy and confidentiality scan, and notes for the PI

Repository: `release_repo/` (no git; the PI runs `git init`). Scan run last, 2026-10-27, after all edits.
Runtime of the whole repo-assembly task: about 40 min wall clock (`date` start 1793098056, end 1793100471), plus two
8-minute test runs against the project data.

## 1. Scan (Python regex over every file; `.gz` files decompressed)

| Pattern | Hits | Resolution |
|---|---|---|
| unix user name, group name (and variants of the surname) | 0 | — |
| numeric uid | 2 files | false positives: the digits occur inside score values in `predictions/*.tsv.gz` (e.g. `2.588982`); left as is |
| email domain of the researcher; any email address | 0 real | the generic email regex matched only metric names like `Coverage@0.5`, `FPR@0.5` |
| hostname of the compute node; cluster names | 0 | — |
| `/workspace` | 0 (after scrub) | ~60 occurrences in docs and code removed: in code via `ksp.paths` (see §3); in docs replaced by relative paths |
| `/home` | 1 | false positive: URL `https://gtexportal.org/home/downloadsPage` in `docs/lit/04_other_modalities.md`. `/home/plm_env`, `/home/tools/mmseqs`, `/home/.cache/huggingface` were replaced by `<plm-env>`, `mmseqs`, `$HF_HOME` |
| `/tmp` | 0 (after scrub) | 8 reviewer-scratchpad paths in `docs/reviews/*` contained the uid and the session directory; replaced by `<agent-scratchpad>` |
| `session_`, session URL, session UUID | 0 | — |
| API keys, tokens, passwords, private keys, Bearer headers (generic and known formats: sk-, ghp_, hf_, AKIA, xox, PEM) | 0 | — |
| IPv4 addresses | 0 (after scrub) | 1 found: `docs/reviews/r1_ai_failure_modes.md` named the NFS server IP and a share path containing a lab name; replaced by `<nfs-server>:/<share>/...` |
| lab or institution names | literature only | the lab name above was the only non-citation occurrence; remaining hits are author names in published citations (`docs/lit/`) |

Not changed, but worth knowing: `ksp.manifest` writes `socket.gethostname()` into every run manifest
(`results/runs/*/manifest.json`). Those files are not in the repository and `results/` is in `.gitignore`, so the
code was left as it is (no computation changes). The PI should not commit run manifests from the project.

`PROMPT.md` is byte-identical to the project's `PROMPT.md` (SHA-256 `0391da01…`). It contains no personal data: no
names of the researcher, no emails, paths, hosts or tokens. It names one published lab ("Evangelia Petsalaki's lab")
as a pointer to the SELPHI2.0 paper; that is a literature reference, not personal data.

## 2. Size

Total 40,787,309 bytes (41 MB), 366 files; no file is larger than 100 MB.

| Bytes | File |
|---|---|
| 34,750,123 | `predictions/site_top5.tsv.gz` |
| 1,952,136 | `predictions/kinase_top200.tsv.gz` |
| 78,621 | `docs/REPORT_internal.md` |
| 73,592 | `results_summary/test_721dc1977366/B0_pop_per_kinase.tsv` |
| 68,531 | `docs/lit/01_methods_benchmarks.md` |
| 60,186 | `results_summary/test_721dc1977366/B2_t3_per_kinase.tsv` |
| 59,694 | `results_summary/test_721dc1977366/B2_per_kinase.tsv` |
| 58,972 | `results_summary/test_721dc1977366/B4_per_kinase.tsv` |
| 58,944 | `src/ksp/eval_test.py` |
| 58,913 | `results_summary/test_721dc1977366/final_per_kinase.tsv` |

## 3. Files copied from restricted sources

**None.** Nothing was copied from `data/`, `external/`, `cache/` or `sealed/`, and no parquet or numpy file is
included. The only PSP-derived files are the released predictions and the metric tables, which the brief allows.
Their checksums match `RELEASE_MANIFEST.json`, and they carry the CC BY-NC-SA 3.0 statement.
The 12 TEST `*_site_level.tsv` files were left out because they list labelled TEST kinase–site pairs.

Borderline items that the PI may want to review (kept, because they are small curation tables or well-known facts,
not database content):
- `config/bio_controls.tsv`: 45 textbook kinase–site pairs with literature citations. It has `in_psp`, `in_ext` and
  `in_test_split` flags (one bit per pair, derived from PSP and the split).
- `tests/test_labels_universe.py`: asserts 7 well-known pairs (e.g. AKT1→GSK3B S9) and one PSP site-group ID.
- `config/kinase_curation.tsv`: one comment names a single PSP Y label that was dropped as out of class.

## 4. Discrepancies and issues for the PI

1. **The released `site_top5` ranking can put a non-top kinase in rank 1.** `assignment_prob` maps the top kinase's
   softmax mass through the isotonic curve and spreads the rest over the other kinases in proportion. When the curve
   lowers the top mass, another kinase can receive more. In the released file, at 25,532 of 350,736 sites (S 17,826,
   T 7,663, Y 43) the rank-1 kinase is not the kinase with the highest `prob` among the five listed. 19 Y sites
   (INSR/IGF1R paralog pairs) reach `assignment_prob` 0.55. That contradicts the manifest and REPORT §9 wording "no
   Y-site assignment reaches 0.5"; the wording holds only for the top kinase, which is what `coverage.tsv` counts
   (rank-1 ≥0.5 in the file: 29,905; Coverage@0.5: 29,886; the difference is exactly these 19 Y sites). The files
   were not changed. The caveat is documented in `predictions/README.md` and `docs/RELEASE_NOTES.md`. The PI may
   want to state it in the manuscript, or regenerate the file ranked by the top-kinase rule. It probably also
   applies to the TEST `assignment_prob` file, but the A17 reliability check uses only the top kinase.
2. REPORT.md was updated by the PI during this task (B4 M2 upper CI 0.157 → 0.156; `superiority.tsv` gives
   0.15649, so 0.156 is right). The copy `docs/REPORT_internal.md` and `README.md` use 0.156.
3. The test suite, run in the original project before this task, called `git diff --quiet HEAD -- FREEZE.json`
   (read-only) inside `test_real_paths_untouched`. In the release verification that test was deselected, so no further
   git command was run.

## 5. Code changes and test results

- Paths: new `src/ksp/paths.py` (`KSP_DATA_DIR`, `KSP_MMSEQS`, `KSP_PYTHON`). 81 files had their
  `ROOT = Path(__file__)...` replaced; the shell wrappers, download scripts and absolute defaults were fixed as well.
  `grep` and a Python scan for `/workspace`, `/home`, `/tmp` and the user, host and IP patterns over `src/`,
  `scripts/`, `tests/` and `config/` find none.
- No computation changed. With `KSP_DATA_DIR` set to the project directory, all module-level data paths are
  identical to the original code.
- Tests, clean copy, no data: **51 passed, 50 skipped, 0 failed**. The same in a fresh venv with
  `pip install -e .` (uv, `--no-deps`, dependencies borrowed from the project venv; no network).
- Tests, clean copy, with the project data: **97 passed, 3 skipped, 0 failed, 1 deselected** (499 s).
- Excluded from the package: `ksp.models.r5` and `ksp.models.regime`, and the exploratory scripts (R5, regimes,
  post-TEST, diagnostics, PSP profiling, data-register generators, machine-specific queue chains).

## 6. What is missing for full reproduction

- All input data: PhosphoSitePlus (licensed; the user must download it), Klaeger 2017, the perturbation compendium
  (provenance undocumented; it must come from the researcher or its source studies), CPTAC (data use policy), the
  Cantley kinase-library matrices and the Phosformer weights (non-commercial), and about 12 GB of public downloads.
  Several sources change between releases, so new downloads will not match the recorded checksums, and the split
  and numbers will differ.
- Compute: a CUDA GPU for ESM-2 embeddings (650M and 3B), R3 and the Phosformer pre-training (B3 took about
  3,150 s to train and 6,480 s to score, on a shared node with 4× RTX 4090); about 20 GB peak RAM for the release
  step, ≤60 GB per job; MMseqs2.
- The intermediate caches (feature stores, cross-fitted stores, CV and TEST predictions, the 4.1 GB release score
  table) are not included. They can be recovered only by re-running the pipeline; their SHA-256 pins are in the
  project's run manifests and freeze entries, which are not included either.
- The sealed TEST set is spent; the TEST batch cannot be re-run honestly on the same split.
- Git history and run manifests are not included. The manuscript and supplement still have to be copied into
  `manuscript/`.
- Author names are placeholders in `LICENSE`, `CITATION.cff` and `pyproject.toml`.
