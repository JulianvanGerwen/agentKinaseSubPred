# Final audit of the public repository (`release_repo/`)

Auditor: independent agent (did not build the repository). Started 2026-10-27 12:51:59 CET, finished about 13:00 CET.
No git was run. Nothing under `sealed/`, `data/` or `external/` was read. Edits were made inside `release_repo/` only.

## 1. Privacy, confidentiality and licence scan (all 384 files)

Searched every text file for the unix user and group names, the hostname, `ethz`, `@` email patterns, IPv4
addresses, absolute paths (`/workspace`, `/home`, `/tmp`, `/cluster`, `/scratch`, `/Users`, `/mnt`), session URLs and
ids, scratchpad ids, keys and tokens (`sk-`, `ghp_`, `hf_…`, `api_key`, `password`, `Bearer`), file-server shares and
HPC names. I checked PDF and PNG metadata (Matplotlib only) and the gzipped predictions (no paths). Result: **no hits**.
The only matches are false positives: a GTEx URL containing `/home/`, "protein/cluster" in prose, "Euler-Rodrigues"
in a docstring, and the word "scratchpad" or the placeholder `<agent-scratchpad>`.

Restricted data:
- No PSP, Klaeger or CPTAC tables. `docs/DATA_internal.md` lists file names, checksums and licences only.
- `results_summary/` holds per-kinase, strata, CI and reliability metrics only. I checked the columns of every TSV and
  the keys of `summary.json`; none lists sites. No `*_site_level.tsv` file is present.
- The tests and docs name individual illustrative sites (≤11 per file), which is acceptable.
- The predictions carry the CC BY-NC-SA 3.0 notice (`predictions/LICENSE-PREDICTIONS.md`, README, NOTICE, CITATION.cff, manifest).

**Fixes:**
1. `config/bio_controls.tsv`: I removed the columns `in_psp`, `in_ext`, `in_U_ms` and `in_test_split`. Together they
   formed a 44-row table of PSP membership per pair, and they flagged 8 substrate proteins as members of the sealed
   TEST split. The literature pairs, windows and citations are kept, and the comment header explains the removal.
   `scripts/build_splits.py` re-adds `in_test_split` when it runs locally, and no code reads the other columns.
   Recorded in `docs/RELEASE_NOTES.md`.
2. `.gitignore`: I added `*_site_level.tsv` as a guard against committing label-level TEST tables.
3. I added the AI-exercise flag to `manuscript/CLARITY_CHANGES.md`, `INDEPENDENT_REVIEW_OF_DRAFT.md` and
   `references.md`, which lacked it. `PROMPT.md` was left byte-identical on purpose, and `CITATION.cff` names the
   exercise in its title.

## 2. `docs/REPORT_internal.md` refreshed

The file is now the 4-line release header followed by the current `REPORT.md`. One scrub was needed:
`` `release_repo/predictions/README.md` `` → `` `predictions/README.md` in this repository ``. The diff from the old
copy has 7 changed lines, all of them the intended corrections:
- F6b and full-model M2 changed from 0.154 / 0.208 to 0.153 / 0.207.
- The wording on Tyr hubs, the in-vivo-only labels and label noise (§10 item 14).
- The rank-1 caveat for `assignment_prob`.

B4 M2 CI [0.091, 0.156] was already correct in both copies. No sensitive string entered the file (re-scanned).

## 3. `manuscript/`

- I rewrote `README.md`. It now lists the files and their roles, says that the manuscript is Markdown for PLOS Comp
  Biol and that conversion to Word or LaTeX is left to the researcher, marks the review as a review of an earlier
  draft, and maps the `results/...` paths to the repository folders.
- The manuscript files contain no paths or names that need scrubbing. Their content matches `manuscript/` in the
  project byte for byte, apart from the added flags. The review is `review_manuscript.md` under a new name.
- `figures/make_figures.py` now reads `results_summary/test_721dc1977366/` and `predictions/coverage.tsv` and takes an
  optional output directory. All required inputs were already in the repository, so I copied nothing.
- I ran it from the root of a fresh copy into a temporary directory. **All 5 PNGs are pixel-identical** to the
  committed ones. Fig 1 and Fig 5 are diagrams with no data inputs.

**For the PI (not edited, because it is manuscript content):** `manuscript.md:204` ("in-vivo-only labels are
enriched for indirect relations") and `supplement.md:63` ("since the latter are enriched for indirect events") are
stronger than the corrected REPORT v3 §10 item 14 ("may reflect indirect events; … cannot separate this from label
count"). Line 63 gives the motivation for the experiment, so it is defensible. Line 204 should be softened in
`manuscript/` and copied again.

## 4. Fresh-clone test

I copied the repository to a temporary directory and made a new uv venv. `pip install -e ".[test]"` resolved the
dependencies from PyPI (core dependencies only), and `import ksp` resolves to the copy.
- **pytest, before the fix:** 7 failed, all in `tests/test_knn.py` (ModuleNotFoundError for torch or Bio, which are
  optional extras).
- **Fix:** the module now calls `importorskip` for torch and Bio, and its wrong `parents[2]` path was corrected.
- **pytest, after the fix:** 42 passed, 51 skipped, 0 failed. Every skip gives its reason. With the project
  environment (torch present), `test_knn.py` gives 9 passed.
- `examples/read_predictions.py` runs in about 1 s and reports SHA-256 OK.
- I checked `sha256sum` independently: both `.gz` files match `RELEASE_MANIFEST.json` and the project
  `results/release/`, and the byte sizes match.
- The rank-1 count printed by the example (29,905 at ≥0.5) is 19 above Coverage@0.5 (29,886), as
  `predictions/README.md` explains.
- **Links:** the docs contain no Markdown hyperlinks, only backticked paths. All backticked paths in the six README
  files resolve, except module names and the original file names `REPORT.md` and `DATA.md`, which are explained.
- I added the result rows to the `docs/RELEASE_NOTES.md` table.

## 5. Honesty of README.md

All of these match REPORT v3 and the manuscript:
- TEST M1 0.794, M2 0.363, M3 0.363, with their CIs and deltas.
- Tyr not trusted (0.595 vs floor 0.585; 10 of 23).
- No ≥0.8 tier; Y capped at 0.445; ≥0.5 tier 0.615 claimed vs 0.603 observed.
- M2 rise is composition (CV 0.348). Cantley fell by 0.025, a weak-baseline draw.
- Baselines are our retrains; SELPHI2.0 defaults on 66.6% of its TEST pairs.
- About 90 CV configurations; one batch spent early; post-TEST work exploratory.
- The rank-1 caveat for `assignment_prob` is stated.

The AI flag is at the top, and `cmp PROMPT.md /workspace/PROMPT.md` shows the files are identical.

**Fixed:** "manuscript submitted to" → "manuscript prepared for" (it has not been submitted). The layout line said
`manuscript/` was a placeholder, which is no longer true.

**`.gitignore`:** checked with `pathspec` (gitignore semantics), without git.
- None of the 384 repository files is ignored, including `predictions/`, `results_summary/`, `manuscript/figures/`
  and `make_figures.py`.
- These are ignored: `data/`, `external/`, `cache/`, `sealed/`, `results/`, `logs/`, `FREEZE.json`, `*.parquet`,
  `*.npy`, `*.pt`, `*.safetensors`, `*.zip`, `.venv/`, `.env`, `__pycache__`, and now `*_site_level.tsv`.

## 6. Final state

Size 43 MB (42,746,683 bytes), 384 files. Files >5 MB: `predictions/site_top5.tsv.gz` (34.8 MB; under GitHub's 50 MB
warning, no LFS needed).

```
.gitignore  CITATION.cff  DATA_SOURCES.md  LICENSE  NOTICE  PROMPT.md  README.md  pyproject.toml  requirements.txt
config/            bio_controls.tsv  complex_subunits.tsv  kinase_curation.tsv  kinase_homonyms.tsv
docs/              README.md  REPORT_internal.md  RELEASE_NOTES.md  DESIGN.md  DATA_internal.md  LITERATURE.md ...
                   (35 .md files) + audit/  lit/  retro/  reviews/
examples/          read_predictions.py
manuscript/        README.md  manuscript.md  supplement.md  references.md  NUMBER_MAP.tsv  SUPP_NUMBER_MAP.tsv
                   CLARITY_CHANGES.md  INDEPENDENT_REVIEW_OF_DRAFT.md  figures/
predictions/       site_top5.tsv.gz  kinase_top200.tsv.gz  coverage.tsv  RELEASE_MANIFEST.json  README.md  LICENSE-PREDICTIONS.md
results_summary/   README.md  TEST_LEDGER.md  coverage.tsv  cv/  test_721dc1977366/
scripts/           README.md  build_*.py  dl*.sh  ... + baselines/  downstream/  features/  models/
src/               ksp/
tests/             conftest.py  test_*.py (17)
```

## Publishing checklist (human)

1. Choose the repository name and account; update `CITATION.cff` (authors, repository URL) and replace the author
   placeholders in `LICENSE`, `pyproject.toml`, `manuscript/*.md`.
2. Confirm the licences: MIT for code, CC BY-NC-SA 3.0 for predictions and `results_summary/`, and CC BY 4.0 for
   docs. Check that the PSP terms allow redistributing derived predictions; also check the Cantley matrices and
   the Phosformer weights (NOTICE).
3. `git init`, then set `user.name` and `user.email` to the identity you want to publish under (the repository
   carries none), `git add -A`, and check `git status` against this tree (384 files). Make the first commit, then push.
4. Before pushing, re-run `python -m pytest -q` and `python examples/read_predictions.py` once in a fresh venv.
5. Optional: archive a tagged release on Zenodo for a DOI and add it to `CITATION.cff` and the manuscript's
   data-availability statement.
6. Optional: soften `manuscript.md:204` (see §3) before converting the manuscript.

## Remaining risks

- The predictions are a PSP-derived product. The NC-SA terms are stated, but their redistribution relies on the
  PSP licence; the publishing researcher should confirm it.
- The docs reproduce the internal record verbatim and contain many project-relative paths to files that are not
  shipped. This is documented, but readers will find dead references.
- Some docs give label-derived examples in prose (for example the GLO1 Y136 hub, labelled to 15 TKs). They are
  illustrative sites, not tables, but they do describe PSP content.
- The frozen-harness checks in `ksp.eval_test` do not pass on this copy, because the paths changed. This is
  documented.
