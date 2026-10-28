# Report: preparation of v1.0.0 for publication

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Agent task, started 2026-10-28 10:34:27 CET and finished 11:14:41 CET (`date`), about 40 min. No git command was run. Nothing
under `sealed/` was opened, and no TEST label file was read. `FREEZE.json`, `results/test/**`, the prediction files
in `results/release/` and the run manifests in `results/runs/` were not edited.

## 1. Tests

| Run | Result |
|---|---|
| Project directory, data present, `KSP_DATA_DIR` unset (`.venv/bin/python -I -m pytest -q`), final code and tests | **100 passed, 2 skipped, 0 failed** (496 s). Skipped: 2 tests in `test_labels_universe.py` that need `cache/labels/ks_pairs.parquet`, which the project deletes after sealing (they used to pass silently) |
| Same, intermediate run before the two guard-test fixes below | 98 passed, 2 skipped, 2 failed in 505 s. The failures were `test_no_test_access.py`: the new scanner contains the names of restricted inputs. Fixed by allowing `scripts/publication/` in that guard |
| Fresh copy without `data/`, `external/`, `cache/`, `sealed/` (except the tracked counts file), `results/runs/` and the 11 compendium tables; new uv venv; `pip install -e ".[test]"` from PyPI (core dependencies, no torch) | **43 passed, 51 skipped, 0 failed** (3.5 s). Every skip gives its reason (needs project data under `KSP_DATA_DIR`; torch for `test_knn.py`; git history for the frozen-file test) |
| First fresh-copy run, before the skip fixes | 49 passed, 10 failed, 15 errors: `results/splits/SEAL.json` is tracked, so the module-level skips did not fire although `cache/` was absent. Fixed by also requiring the `cache/` index and DEV labels |
| Code from the fresh copy, data from the project (`KSP_DATA_DIR` = project directory) | **99 passed, 3 skipped, 0 failed** (506 s): the code reads and writes data only through `KSP_DATA_DIR`. Skipped: the 2 tests above and the frozen-file test (the copy has no git history) |

The new test `test_frozen_files_changed_only_paths` passed in the project directory: the frozen harness files differ
from commit `974d520` only in path lines. `test_real_policy` now expects real mode to be refused for that reason.
`examples/read_predictions.py` runs in the fresh copy and reports SHA-256 OK for both prediction files.

## 2. Grep for absolute paths in `src/`, `scripts/`, `tests/`, `config/`

`grep -rnE "/workspace|/home/|/tmp/|/Users/"` finds only the pattern definitions and docstring of
`scripts/publication/check_publishable.py` (5 lines). `src/ksp/paths.py` needs no absolute default: `REPO_ROOT` is
derived from the file location and `ROOT` defaults to it. With `KSP_DATA_DIR` unset, the 17 module-level `ROOT`
values, `io_psp.DEFAULT_PATH`, `uniprot.DEFAULT_DIR` and `kinases.CURATION` all resolve into the repository root,
as before. The only changed default is the MMseqs2 binary (`$KSP_MMSEQS` or `mmseqs` on `PATH`; it is not on `PATH`
here, so rebuilding the splits needs `KSP_MMSEQS`). All 232 Python files parse, and `bash -n` passes on all 23
shell scripts.

## 3. `check_publishable.py`

`python -I scripts/publication/check_publishable.py` (walks the tree and applies `.gitignore`; git not used):

```
scanned 2691 files, 109.9 MB
ERROR: 0
WARN: 1036
   402  host field in run manifest           (results/runs/*/manifest.json; kept unedited)
   402  user/group/host name                 (the same host names in those manifests)
   222  absolute /workspace path (per file)  201 results/runs, 11 docs/audit scripts, 6 logs/regime chains,
                                             1 results/test/LEDGER.md, 3 notes that name the pattern
     9  scratchpad-path                      (3 results/runs manifests)
     1  file>5MB                             results/release/site_top5.tsv.gz, 34.8 MB
```

Fixed before this result (first run: 55 errors):
- 30 scratchpad paths in 22 files (`docs/audit/`, `docs/reviews/`, `docs/lit/`, `briefing/BRIEF.md`,
  `manuscript/BRIEF_COMMON.md`) → `<scratchpad>`;
- the internal file-server address and share path in `docs/reviews/r1_ai_failure_modes.md` →
  `<internal-file-server>:<share>`;
- in `docs/publication_audit_history.md`: the address, the host name, the author name and e-mail, the session URL
  and three `/home` paths → placeholders;
- 10 shared-environment `/home/...` paths in 6 documents → `<plm-env>`, `<tools>`, `<mmseqs>`, `$HF_HOME` (3 more in
  `manuscript/repo_notes.md` were reworded);
- absolute `/workspace/...` paths in 25 documents → relative paths (not in `results/`, not in `PROMPT.md`).

The `.gitignore` matcher in the script gives the same result as `pathspec` (current PyPI release) on all files of the tree
and on the HEAD path list (checked before `manuscript/` was excluded: 2,711 files).

## 4. Size

- Files that a push would add at HEAD (tree walk with `.gitignore`): 2,691 files, 109.9 MB (`results/runs/` is about
  two thirds of the files). `release_repo/` (43 MB) was removed; `manuscript/` is excluded.
- Files larger than 5 MB: only `results/release/site_top5.tsv.gz` (34.8 MB; under GitHub's 50 MB warning).
- Not published: `data/` 1.1 GB, `external/` 19 GB, `cache/` about 800 GB, `.venv/` 5.9 GB.

## 5. What was added, removed and changed

Added:
- root: `LICENSE`, `NOTICE`, `CITATION.cff`, `DATA_SOURCES.md`, `pyproject.toml` (version 1.0.0), new `README.md`;
- `src/ksp/paths.py`; `examples/read_predictions.py` (reads `results/release/`);
- `results/release/README.md`, `results/release/LICENSE-PREDICTIONS.md`;
- `scripts/README.md`, `scripts/publication/check_publishable.py`, `scripts/publication/rewrite_history.sh`;
- `docs/README.md`, `docs/AGENT_SETUP.md`, `docs/RELEASE_NOTES.md`, `docs/PUBLISHING.md`,
  `docs/GITHUB_RELEASE_v1.0.0.md`, this report; `manuscript/README.md`.

Removed: `release_repo/` (384 files; kept in the history).

To be untracked by the PI (still on disk, now ignored): the 11 compendium condition tables (`docs/PUBLISHING.md` §3).

Changed:
- code (path handling only): 24 files in `src/ksp/`, 96 in `scripts/`, 15 in `tests/` (skip conditions, the
  frozen-file test, the guard allowlist), `config/bio_controls.tsv` (4 columns removed);
- `requirements.txt`: now the pinned versions of the study environment (from the release copy; includes
  `rdata==1.1.0`);
- `.gitignore`: rewritten in 8 commented sections;
- documents: the AI flag added to 67 Markdown files; path scrubs as in section 3; short notes on the removed
  compendium tables in `docs/features_phospho.md` and `docs/data_profile.md`; historical notes at the top of
  `manuscript/repo_notes.md` and `manuscript/repo_audit.md`; one link in `REPORT.md` (`release_repo/predictions/...`
  → `results/release/README.md`).

Tracked files that newly match `.gitignore` (36 of the HEAD paths; they stay tracked unless removed): the 21 tracked
files in `manuscript/` and the 11 compendium tables (both intended, to be untracked), and `logs/r2/cv5_nophospho.out`, `logs/r2/cv5_r2_full.out`, `logs/r2/cv5_r2_g2_main.out` (rule `/logs/**/*.out`;
harmless). `sealed/kinase_label_counts.tsv` matched the old rules too. Nothing in `results/release/`,
`results/test/`, `briefing/` or `docs/` is matched; the 12 TEST `*_site_level.tsv` tables are
re-included explicitly.

## 6. Change of plan during the task: briefing in, manuscript out

On request of the PI during this task:
- **`manuscript/` is excluded** from the public repository: a labelled section in `.gitignore` ("manuscript kept
  outside the public repo"); `git rm -r -q --cached manuscript` is in `docs/PUBLISHING.md` §3. It stays on disk.
  Earlier commits still contain it (decision 7 in `docs/PUBLISHING.md`; `rewrite_history.sh --manuscript`).
- **The supplement** (S1 Text, the negative results) is copied to `docs/SUPPLEMENT_failed_approaches.md`, flagged,
  with a note that it belongs to the accompanying paper. The scan finds nothing in it. `README.md` links to it.
- **`briefing/` is included** with all its files. `BRIEF.md`, `CHANGES.md` and `review_briefing.md` contain nothing
  sensitive (scan: 0 findings), so they were kept; the one scratchpad path in `BRIEF.md` had already been scrubbed.
  References to `manuscript/...` files in the briefing, its number map and its process notes were replaced by "the
  accompanying paper (link to be added)", "the paper draft", or `docs/SUPPLEMENT_failed_approaches.md`.
- **The briefing builds without `manuscript/`.** `briefing/figures/make_figures.py` and `briefing/build_html.py`
  read only `results/` and `briefing/`. In a copy of the repository without `manuscript/`, all 8 figures came out
  byte-identical, and `briefing.html` rebuilt from the edited `briefing.md` is byte-identical to the one here. The
  HTML differs from the previous build only in the two edited paragraphs; the embedded images are unchanged.
- **`README.md` puts the briefing first**: after the flag, a "Start here" block with the briefing's title, links to
  the HTML and Markdown, a 3-sentence summary and Figures 1 and 5 of the briefing (`fig1_pipeline.png`,
  `fig4_results.png`). `CITATION.cff` and `docs/GITHUB_RELEASE_v1.0.0.md` start with the briefing. A GitHub "About"
  line (343 characters) and suggested topics are in `docs/PUBLISHING.md` §4.
- The full test suite was run before these changes; they touched no code and no test. The fresh-copy suite was run
  again without `manuscript/` (43 passed, 51 skipped).

## 7. Risks and open points

- **The compendium decision** (`docs/PUBLISHING.md` §2.1). If it is not cleared, the history rewrite is required,
  and it breaks the hash pins in `FREEZE.json`, 406 manifests and the ledger.
- **Frozen harness at HEAD.** `labels.py` and `splits.py` differ from `974d520` in path lines, so real mode refuses
  at v1.0.0. The checks pass at the tag `v0.9-test-batch` (`c3e2a14`). The PI should confirm the tag hashes with
  `git log` before tagging.
- **Line numbers.** The flag adds 2 lines at the top of 67 documents. References of the form `FILE.md:line` to these
  documents are now off by 2. The number maps cite sections, not lines, so they are not affected.
- **Remaining absolute paths** in `results/runs/` (201 files), `results/test/LEDGER.md`, `logs/regime/*.sh` and the
  audit scripts in `docs/audit/` were left as the record. They are shared environment paths, not personal ones.
- **Other agents were working** in `briefing/` and `manuscript/` during this task (files changed after 10:34 that this
  task did not write). If the briefing agent edits `briefing.md` again, it must keep the new wording about the
  paper and rebuild `briefing.html`. The AI-exercise flag was added to the three process notes in `briefing/`.
- `manuscript/BRIEF_COMMON.md` no longer contains the scratchpad path of the clarity skill (now `<scratchpad>/...`).
  Agents launched later need the path in their own brief.
- `scripts/publication/rewrite_history.sh` was syntax-checked and its replacement pattern was tested on the
  affected line, but the script was not run (git-filter-repo is not installed, and it must not run here).
