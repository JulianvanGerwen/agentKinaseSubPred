# Release notes: v1.0.0 (2026-10-28)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

v1.0.0 is the first public release. It is the main project repository, with its full git history, plus the changes
listed here. On 2026-10-27 an agent had built a separate release copy (`release_repo/`, without history). For v1.0.0
those changes were merged into the main repository and the copy was removed; it remains in the git history. No
computation changed, and no feature, model, TEST or release step was re-run.

## Where the files of the release copy went

| In the release copy | In v1.0.0 |
|---|---|
| `predictions/` | `results/release/` (the location that `RELEASE_MANIFEST.json` names; no second copy) |
| `results_summary/` | `results/test/`, `results/cv/` and `results/release/coverage.tsv` (the originals) |
| `docs/REPORT_internal.md`, `docs/DATA_internal.md` | `REPORT.md`, `DATA.md` |
| `docs/DESIGN.md`, `LITERATURE.md`, `RESEARCH_TREE.md`, `LESSONS.md`, `DATAINEED.md` | the same names at the root |
| `manuscript/` (paper draft, supplement, review, figures) | not public at v1.0.0 (kept outside the repository; earlier commits still contain it). The supplement is `docs/SUPPLEMENT_failed_approaches.md` |
| `LICENSE`, `NOTICE`, `CITATION.cff`, `DATA_SOURCES.md`, `pyproject.toml`, `requirements.txt`, `examples/`, `scripts/README.md`, `docs/README.md`, `docs/AGENT_SETUP.md` | the same paths at the root, with paths updated |

## Code changes (path handling only)

1. `src/ksp/paths.py` defines `REPO_ROOT` (the code checkout) and `ROOT` (the project data and working directory:
   `$KSP_DATA_DIR`, default `REPO_ROOT`). It also defines `MMSEQS` (`$KSP_MMSEQS`, or `mmseqs` on `PATH`) and
   `PYTHON` (`$KSP_PYTHON`, or the running interpreter).
2. In the modules, scripts and tests, `ROOT = Path(__file__).resolve().parents[N]` was replaced by an import of
   `ksp.paths.ROOT`. Uses of `ROOT / "src"`, `"scripts"` and `"config"` go through `REPO_ROOT`. This includes the
   modules and scripts that were not in the release copy: `ksp.models.r5`, `ksp.models.regime` and the
   exploratory scripts (`scripts/README.md`).
3. Other replacements: the absolute MMseqs2 path (`ksp.splits`, `scripts/baselines/phosformer_corpus_filter.py`);
   interpreter paths; absolute defaults in `scripts/audit_ppi_ks_leakage.py`, `scripts/dl_pdc_pancancer.py`,
   `scripts/make_data_rows*.py` and `scripts/models/regime_queue.py`; `cd` to the absolute project path in the shell
   wrappers, which now change into `$KSP_DATA_DIR` (or the repository root) and call Python as `$KSP_PYTHON`.
4. The git look-ups for code identity (`ksp.manifest`, the frozen-harness checks in `ksp.eval_test`) run in
   `REPO_ROOT`.

**Default.** With `KSP_DATA_DIR` unset, every module-level path resolves to the same location as before: the 17
module `ROOT` values, `io_psp.DEFAULT_PATH`, `uniprot.DEFAULT_DIR` and `kinases.CURATION` were printed and all point
into the repository root. One default did change: the MMseqs2 binary was an absolute path on the study machine and
is now `$KSP_MMSEQS` or `mmseqs` on `PATH`. It is used only to rebuild the frozen splits.

**Frozen harness.** `eval.py`, `labels.py`, `splits.py`, `crossfit.py` and `canary.py` were frozen at commit
`974d520` (A28). In v1.0.0, `labels.py` and `splits.py` differ from that commit in their path lines only. The new
test `test_frozen_files_changed_only_paths` checks this against git. As a result, `ksp.eval_test` refuses real
mode at v1.0.0. This is intended, because the single TEST batch is spent. The tag `v0.9-test-batch` holds the code
that ran it.

## Data and documentation changes

- `config/bio_controls.tsv`: the columns `in_psp`, `in_ext`, `in_U_ms` and `in_test_split` were removed (per-pair
  PSP and sealed-TEST membership). `scripts/build_splits.py` re-adds `in_test_split` when it runs.
- The 11 condition tables copied from the provided compendium are listed in `.gitignore` and are to be untracked
  (`docs/PUBLISHING.md` §3).
- `manuscript/` is listed in `.gitignore` and is to be untracked: the paper is published separately. Its supplement
  (S1 Text, the negative results) is copied to `docs/SUPPLEMENT_failed_approaches.md`. The research briefing in
  `briefing/` is the entry point for readers; its links to the paper now read "the accompanying paper (link to be
  added)". Its figures and HTML were rebuilt without `manuscript/` and are byte-identical.
- Documents: the AI-exercise flag was added at the top of every document except `PROMPT.md`. Absolute paths,
  scratchpad paths and one internal file-server address were replaced by placeholders. `results/runs/` and
  `results/test/` were not edited.
- New: `README.md`, `.gitignore` (rewritten), `docs/PUBLISHING.md`, `docs/GITHUB_RELEASE_v1.0.0.md`,
  `scripts/publication/check_publishable.py`, `scripts/publication/rewrite_history.sh`, this file, and
  `docs/publication_prep_report.md`.

## Test changes

- `tests/conftest.py` registers a `requires_data(*paths)` marker (from the release copy). Tests that need the index,
  labels or splits under `cache/` are skipped with the missing paths in the message.
- The module-level skips of `test_eval_harness.py`, `test_eval_test_batch.py` and `test_splits.py` now also require
  the `cache/` index and DEV labels. The frozen split files in `results/splits/` are tracked, so a fresh clone has
  them without the data.
- `test_real_policy` expects real mode to be refused when the frozen files differ from `974d520`, and to pass
  when they do not. `test_real_paths_untouched` is skipped outside a git checkout.
- `tests/test_no_test_access.py` allows `scripts/publication/`: the scanners match the file names of restricted
  inputs and read no labels.
- `tests/test_knn.py` skips when torch or Biopython is missing (from the release copy).

## Verification

| Check | Result |
|---|---|
| full test suite in the project directory, data present (`.venv/bin/python -I -m pytest -q`) | 100 passed, 2 skipped (two tests need `cache/labels/ks_pairs.parquet`, deleted after sealing), 0 failed |
| fresh copy without `data/`, `external/`, `cache/`, `sealed/`, `results/runs/` and the compendium tables; new venv; `pip install -e ".[test]"` from PyPI (core dependencies only) | 43 passed, 51 skipped (each with its reason), 0 failed |
| code of that fresh copy with `KSP_DATA_DIR` = the project directory | 99 passed, 3 skipped, 0 failed |
| `examples/read_predictions.py` | SHA-256 of both prediction files OK |
| no absolute `/workspace`, `/home` or `/tmp` paths in `src/`, `scripts/`, `tests/`, `config/` | none, except the patterns inside `scripts/publication/check_publishable.py` |
| `scripts/publication/check_publishable.py` | 0 errors; warnings for the unedited run manifests and other records (`docs/publication_prep_report.md` §3) |

## Known issue in the released predictions

In `results/release/site_top5.tsv.gz` the rank is by `assignment_prob`. When the isotonic map lowers the top
kinase's mass, the redistribution formula can give another kinase a higher `assignment_prob`. This happens at
25,532 of 350,736 sites. It also lets 19 Y sites (INSR/IGF1R paralog pairs) reach 0.55, although the manifest says
that no Y site reaches 0.5; that statement holds for the top kinase only. The released files were not changed. See
`results/release/README.md`.
