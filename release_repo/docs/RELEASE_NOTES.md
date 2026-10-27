# Release notes: MVP repository (2026-10-27)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

An agent assembled this repository from the project directory after the study had finished. Git history is not
included. This note lists what was copied, what was changed, and how the result was checked.

## Contents

- Copied without change: `PROMPT.md` (byte-identical); `predictions/*` (byte-identical, SHA-256 checked);
  `config/*.tsv`; the metric tables in `results_summary/` (only the absolute output path in the TEST ledger was made
  relative).
- Copied with paths made configurable: `src/ksp/` (without `models/r5` and `models/regime`), the pipeline scripts
  listed in `scripts/README.md`, and `tests/`.
- Copied with absolute paths scrubbed and the AI-exercise flag added: the documents in `docs/`.
- New: `README.md`, `DATA_SOURCES.md`, `LICENSE`, `NOTICE`, `CITATION.cff`, `.gitignore`, `pyproject.toml`,
  `requirements.txt` (pinned to the versions of the project environment), `src/ksp/paths.py`,
  `examples/read_predictions.py`, `predictions/README.md`, `predictions/LICENSE-PREDICTIONS.md`,
  `results_summary/README.md`, `scripts/README.md`, `docs/README.md`, `docs/AGENT_SETUP.md`, this file and
  `manuscript/README.md`.

Left out on purpose:
- all data, caches, run outputs, logs, the sealed TEST labels, `FREEZE.json`, parquet and numpy files;
- the TEST `*_site_level.tsv` tables, which list labelled TEST pairs;
- from `config/bio_controls.tsv`, the columns `in_psp`, `in_ext`, `in_U_ms` and `in_test_split` (per-pair PSP
  membership and sealed-TEST membership; removed in the final audit). The 44 literature pairs and citations remain;
- scripts and modules of exploratory or negative branches (`scripts/README.md`, "Not included");
- the audit and reviewer helper scripts that were kept with the audit reports.

## Code changes (path handling only)

The project resolved every path from the location of the source files, so the code and the data had to sit in one
directory. Some paths were absolute. For this release:

1. A new module, `src/ksp/paths.py`, defines `REPO_ROOT` (the code checkout) and `ROOT` (the project data and working
   directory: `$KSP_DATA_DIR`, default `REPO_ROOT`). It also defines `MMSEQS` (`$KSP_MMSEQS`, or `mmseqs` on `PATH`)
   and `PYTHON` (`$KSP_PYTHON`, or the running interpreter).
2. In 81 modules, scripts and tests, the line `ROOT = Path(__file__).resolve().parents[N]` was replaced by an import
   of `ksp.paths.ROOT`. Uses of `ROOT / "src"`, `"scripts"`, `"config"` and `".venv"` now go through `REPO_ROOT`.
3. Other replacements:
   - the hard-coded MMseqs2 binary path (`ksp.splits`, `scripts/baselines/phosformer_corpus_filter.py`);
   - the interpreter path in `scripts/baselines/phosformer_cv_run.py`;
   - absolute defaults in `scripts/audit_ppi_ks_leakage.py` and `scripts/dl_pdc_pancancer.py`;
   - the defaults of `KSP_DATA` and `KSP_UNIPROT` in `ksp.io_psp` and `ksp.uniprot`, which now follow `ROOT`;
   - `cd` to the absolute project path in the shell wrappers. They now change into `$KSP_DATA_DIR` (or the repository
     root), call scripts by their repository path, and call Python as `$KSP_PYTHON`;
   - the download directory in `scripts/dl*.sh`;
   - interpreter paths in usage docstrings.
4. The git look-ups for code identity (`ksp.manifest.git_hash`, the frozen-harness checks in `ksp.eval_test`) now
   run in `REPO_ROOT`. The FREEZE.json checks still run in the directory that holds FREEZE.json.

No computation was changed. Check: with `KSP_DATA_DIR` pointing at the original project directory, every
module-level path constant of every `ksp` module was compared with the original code. All data paths were
identical. Only code-relative paths differed (`__file__` and the three `config/*.tsv` paths, which point to
identical copies).

Consequence: the harness files frozen at commit 974d520 (`eval.py`, `labels.py`, `splits.py`, `crossfit.py`,
`canary.py`) now differ from that commit. So `ksp.eval_test.check_real_policy` cannot pass on this copy. That is
intended: the single TEST batch is spent.

## Test changes

- `tests/conftest.py` registers a `requires_data(*paths)` marker. A test with the marker is skipped, with the
  missing paths in the message, when any of them is absent under `KSP_DATA_DIR`. Ten tests that read the project
  index or labels carry the marker.
- Existing skips that said "split not frozen" now say that project data are needed, and point to
  `DATA_SOURCES.md`. `tests/test_labels_universe.py` used to print "SKIP" and pass; it now calls `pytest.skip`.
- `test_eval_test_batch.py::test_real_policy` checks that the frozen harness files are unchanged since git commit
  974d520. It runs only when `KSP_ORIGINAL_REPO=1` is set, because this release has neither that commit nor the
  unchanged files.

## Verification

| Check | Result |
|---|---|
| SHA-256 of `kinase_top200.tsv.gz` and `site_top5.tsv.gz` against `RELEASE_MANIFEST.json` | match; byte sizes match |
| `PROMPT.md` against the original | byte-identical (SHA-256 `0391da01…`) |
| syntax check of all 171 Python files; `bash -n` on all shell scripts | no errors |
| unit tests in a clean copy, no project data (`KSP_DATA_DIR` set to an empty directory, and unset) | **51 passed, 50 skipped, 0 failed** |
| the same, in a fresh virtual environment with `pip install -e .` (uv, `--no-deps`; dependencies taken from the project environment) | 51 passed, 50 skipped, 0 failed; `import ksp` resolves to the checkout |
| unit tests in a clean copy, with `KSP_DATA_DIR` = the original project directory | **97 passed, 3 skipped, 0 failed** (499 s). Skipped: 2 tests in `test_labels_universe.py` that need `cache/labels/ks_pairs.parquet`, which the project deletes after sealing; `test_real_policy` (needs the original git history). Not run: `test_real_paths_untouched`, which calls `git diff` on the project's FREEZE.json (git is reserved for the PI). In the original project the same suite gave 101 passed; the difference is the 2 tests that used to pass silently when their input was missing, plus the 2 tests above. |
| `examples/read_predictions.py` | runs; checksums OK |
| final audit: fresh copy, new virtual environment, `pip install -e ".[test]"` with dependencies resolved from PyPI (core dependencies only, no torch or Biopython) | 42 passed, 51 skipped, 0 failed. `tests/test_knn.py` now skips as a module when torch or Biopython is missing; before this fix its 7 tests failed with ModuleNotFoundError in a core-only install |
| final audit: `manuscript/figures/make_figures.py` from the repository root, inputs from `results_summary/` and `predictions/` | all five figures regenerate |

Not verified: no feature, model, TEST or release step was re-run from this repository. The pipeline needs the
licensed data and a GPU.

## Known issue in the released predictions

In `predictions/site_top5.tsv.gz` the rank is by `assignment_prob`. When the isotonic map lowers the top kinase's
mass, the redistribution formula can give another kinase a higher `assignment_prob`. This happens at 25,532 of
350,736 sites. It also lets 19 Y sites (INSR/IGF1R paralog pairs) reach 0.55, although the manifest says that no Y
site reaches 0.5; that statement holds for the top kinase only. The released files were not changed. See
`predictions/README.md`.
