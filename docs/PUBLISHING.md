# Publishing this repository (for the researcher)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

This note is for the researcher who publishes the repository. It says what was prepared for v1.0.0, which decisions
are yours, and the exact commands. The agents did not run git. The history audit is `docs/publication_audit_history.md`;
the report of this preparation, with test counts and scan output, is `docs/publication_prep_report.md`.

## 1. What was done for v1.0.0

- **Release merged into the main repository.** The separate copy `release_repo/` was removed. Its changes were
  applied here, and the git history keeps the copy itself:
  - `src/ksp/paths.py` makes every path configurable (`KSP_DATA_DIR`, `KSP_MMSEQS`, `KSP_PYTHON`). With
    `KSP_DATA_DIR` unset, every path resolves to the same location as before. The modules that existed only here
    (`ksp.models.r5`, `ksp.models.regime`, and the exploratory scripts) were changed in the same way.
  - `config/bio_controls.tsv` lost the columns `in_psp`, `in_ext`, `in_U_ms` and `in_test_split` (per-pair PSP and
    TEST membership).
  - Promoted to the root: `LICENSE`, `NOTICE`, `CITATION.cff`, `DATA_SOURCES.md`, `pyproject.toml` (with pinned
    versions in `requirements.txt`), `examples/`, `scripts/README.md`, `docs/README.md`, `docs/AGENT_SETUP.md`,
    `docs/RELEASE_NOTES.md`, `results/release/README.md` and
    `results/release/LICENSE-PREDICTIONS.md`. The predictions stay at `results/release/`, which is the location that
    `RELEASE_MANIFEST.json` names; there is no second copy.
- **New `README.md`.** It starts with the AI-exercise flag and then the research briefing (`briefing/`), which is
  the entry point for outside readers. It also gives the results, the limitations and a guide to files.
- **`manuscript/` is kept out of the public repository** (listed in `.gitignore`; to be untracked). The paper will be
  published separately. Links to it in `README.md`, `docs/` and `briefing/` now read "the accompanying paper (link
  to be added)". Its supplement (S1 Text, the negative results) is copied to `docs/SUPPLEMENT_failed_approaches.md`.
  `briefing/figures/make_figures.py` and `briefing/build_html.py` read nothing in `manuscript/`; rebuilt in a copy
  without it, all 8 figures and `briefing.html` are byte-identical to the files here.
- **Hygiene.** Scratchpad paths, absolute paths in documents, and the internal file-server address were replaced
  by placeholders. The AI-exercise flag was added to every document except `PROMPT.md` (kept verbatim). Run
  manifests in `results/runs/` and everything in `results/test/` were not edited.
- **`.gitignore`** was rewritten in commented sections. `scripts/publication/check_publishable.py` scans what a
  push would publish (section 4).
- **Tests.** One test of the frozen harness now checks that the frozen files changed only in path handling. Tests
  that need data skip cleanly in a fresh clone.

## 2. Decisions that are yours

1. **Is the phosphoproteomics compendium in `data/` published or cleared for release?** Its provenance and licence
   are unknown (`DATAINEED.md` item 5). Eleven tracked tables copy its condition metadata (section 3).
   - If it is published or cleared: keep the tables tracked (delete the section-2 block of `.gitignore`) and cite it
     in `DATA.md`.
   - If it is not: untrack them at HEAD (section 3), **and** rewrite the history (section 5), because 13 paths are
     in earlier commits. `git rm --cached` alone is not enough.
2. **The internal file-server address** was in `docs/reviews/r1_ai_failure_modes.md` in 103 commits. It is a
   private (RFC 1918) address with a lab share path, so it cannot be reached from outside. It is removed at HEAD.
   Removing it from the history needs the rewrite (`--nfs`). That is recommended but optional, unless your
   institution's IT policy forbids publishing internal addresses.
3. **Author name and e-mail, and the `Claude-Session` trailers, in the history.** All 108 commits carry your name
   and institutional address as author and committer. 93 commit messages carry the same `Claude-Session` URL. The
   URL is not a credential; it opens only for the owning account and for people the transcript is shared with.
   Keeping both matches the goal of showing the agentic record. Changing them needs `--mailmap` and
   `--strip-session`, which rewrite every commit.
4. **Host names in run manifests.** 402 `results/runs/*/manifest.json` files contain the cluster GPU node name, a
   process id and the command line. They were left unedited because they are the run record. Rewriting them is not
   worth the hash change. To stop this for new runs, drop `host` and `pid` in `ksp.manifest` (not done: it would
   change a module that the run records describe).
5. **Licence of PSP-derived material.** The predictions and the derived tables in `results/` are released under
   CC BY-NC-SA 3.0, the PhosphoSitePlus licence (`results/release/LICENSE-PREDICTIONS.md`, `NOTICE`). Please confirm
   that the PSP terms allow redistributing derived predictions and tables. Also check the Cantley kinase-library
   matrices (CC BY-NC-SA 3.0) and the Phosformer weights (CC BY-NC-ND 4.0; not included).
6. **The sealed TEST labels are now public.** `results/test/721dc1977366/*_site_level.tsv` (12 files) list the
   TEST-positive kinase–site pairs, and `results/splits/test_proteins.tsv` plus PSP rebuild the same set. The TEST
   set was scored once and is spent, so this costs nothing scientifically. A rewrite is needed only if you want
   the TEST set to stay unpublished.
7. **`manuscript/` in the history.** It is untracked at v1.0.0, but earlier commits still contain the paper draft,
   its supplement, the review of the draft and the working notes of the publication phase. Keep them as part of
   the record, or remove them from all commits with the rewrite (`--manuscript`, section 5).
8. **Authors and repository URL.** Replace "[authors to be added by the publishing researcher]" in `LICENSE`,
   `CITATION.cff`, `pyproject.toml`, `NOTICE` and `briefing/briefing.md`, and add the repository URL to `CITATION.cff`.
   Replace "the accompanying paper (link to be added)" in `README.md`, `briefing/`, `docs/SUPPLEMENT_failed_approaches.md`
   and `docs/GITHUB_RELEASE_v1.0.0.md` once the paper is out.

## 3. Files to untrack at v1.0.0

`manuscript/` (all of it; it stays on disk):

```
git rm -r -q --cached manuscript
```

The compendium condition tables follow.

### 3.1 The compendium condition tables

These 11 tracked files copy condition-level metadata of the compendium (conditionID, Condition, Control, Sample,
Labeling, Timepoint, PMID, source). They contain no log2FC values or peptide IDs. They are still on disk and are
listed in `.gitignore` (section 2). No test needs them, and a fresh clone without them passes its tests. Two
pipeline modules read `f6b_conditions.tsv` when they run with data (`ksp.downstream.ldo`,
`ksp.features.phospho_inhib.condition_ids_for_drugs`), so regenerate it before you run those steps.

| File | Regenerate with |
|---|---|
| `results/data_profile/compendium_per_condition_coverage.tsv` | `scripts/profile_compendium.py` |
| `results/data_profile/compendium_conditions_by_pmid.tsv` | `scripts/profile_compendium.py` |
| `results/data_profile/kinhib_hijazi_conditions.tsv` | `scripts/profile_compendium.py` |
| `results/data_profile/kinhib_condition_matches.tsv` | `scripts/profile_compendium.py` |
| `results/features/phospho/conditions_meta.tsv` | `scripts/features/phospho_build_blocks.py` |
| `results/features/phospho/f6b_conditions.tsv` | `scripts/features/phospho_f6b.py` |
| `results/features/phospho/f6b_conditions_clean_bench.tsv` | `scripts/features/phospho_f6b.py --variant clean_bench ...` (arguments in `docs/features_phospho.md`) |
| `results/features/phospho/f6b_conditions_f6_clean.tsv` | `scripts/features/phospho_f6b.py --variant f6_clean ...` (arguments in `docs/features_phospho.md`) |
| `results/features/phospho/f6b_unmapped_or_dropped_drug_conditions.tsv` | no committed script; the drug-like conditions of `f6b_conditions.tsv` that were unmapped or dropped |
| `results/features/phospho/label_source_studies.tsv` | `scripts/features/phospho_study_lists.py` |
| `results/features/phospho/benchmarkin_studies.tsv` | `scripts/features/phospho_study_lists.py` |

Two more paths exist only in the history: `results/features/phospho/f6b_conditions_lam10.tsv` and
`f6b_conditions_lam100.tsv`.

To untrack the 11 files at HEAD (they stay on disk):

```
git rm --cached results/data_profile/{compendium_per_condition_coverage,kinhib_hijazi_conditions,kinhib_condition_matches,compendium_conditions_by_pmid}.tsv \
  results/features/phospho/{conditions_meta,f6b_conditions,f6b_conditions_clean_bench,f6b_conditions_f6_clean,f6b_unmapped_or_dropped_drug_conditions,label_source_studies,benchmarkin_studies}.tsv
```

## 4. Commit, tag and push

Run these in the project directory. Check every step before you go on.

```
# 1. Scan what would be published (uses git ls-files; exit code 1 means findings)
python -I scripts/publication/check_publishable.py --git

# 2. Stage the v1.0.0 changes
git rm -r -q --cached release_repo                 # already deleted on disk
git rm -r -q --cached manuscript                   # the paper is published separately (section 3)
git rm --cached <the 11 paths of section 3.1>      # only if the compendium is not cleared
git add -A
git status --short | less                          # review: new, changed and deleted files
python -I scripts/publication/check_publishable.py --git

# 3. Run the tests once more (about 9 minutes with the project data)
.venv/bin/python -I -m pytest -q

# 4. Commit
git commit -m "prepare for publication / v1.0.0"

# 5. Tags (annotated). Check the hashes first:  git log --oneline | grep -E "5fc7870|c3e2a14|974d520"
git tag -a v0.9-freeze 5fc7870 -m "FREEZE.json committed before the single TEST batch (A28 harness 974d520)"
git tag -a v0.9-test-batch c3e2a14 -m "Code that ran the single TEST batch (git_head in results/test/LEDGER.md)"
git tag -a v1.0.0 -m "v1.0.0: release merged into the main repository; configurable paths; licences"

# 6. Push to a new, empty GitHub repository
git remote add origin git@github.com:<account>/<repository>.git
git push -u origin master
git push origin v0.9-freeze v0.9-test-batch v1.0.0

# 7. GitHub release with the release notes and the prediction files as assets
gh release create v1.0.0 --title "ksp v1.0.0" --notes-file docs/GITHUB_RELEASE_v1.0.0.md \
  results/release/site_top5.tsv.gz results/release/kinase_top200.tsv.gz results/release/RELEASE_MANIFEST.json
```

On GitHub, set the repository description ("About", at most 350 characters) to:

> Kinase–substrate site predictor for the human proteome, built and tested in 24 h by AI agents (Claude) under a human-written brief. Leakage-controlled benchmark with a sealed test; released predictions (CC BY-NC-SA 3.0, from PhosphoSitePlus); full agentic record in the git history. Start with briefing/briefing.html. Check results before use.

Suggested topics: `kinase`, `phosphorylation`, `kinase-substrate-prediction`, `phosphoproteomics`, `bioinformatics`,
`computational-biology`, `machine-learning`, `positive-unlabeled-learning`, `protein-language-models`, `benchmark`,
`reproducibility`, `ai-agents`, `llm-agents`. Set the website field to the GitHub Pages URL of
`briefing/briefing.html` if you enable Pages, or leave it empty.

The largest file is `results/release/site_top5.tsv.gz` (34.8 MB), under GitHub's 50 MB warning limit, so Git LFS is
not needed. A full push is about 53 MB. Optionally, archive the tag on Zenodo for a DOI and add it to
`CITATION.cff` and the paper's data-availability statement.

## 5. Optional: rewrite the history

`scripts/publication/rewrite_history.sh` wraps `git filter-repo`. It can:
- remove the 13 compendium paths from all commits (`--paths`);
- remove `manuscript/` from all commits (`--manuscript`);
- replace the internal file-server address and share path with `<internal-file-server>:<share>` (`--nfs`);
- rewrite author identities from a mailmap file that you write (`--mailmap FILE`);
- remove the `Claude-Session` trailers (`--strip-session`).

**Warnings.**
- **It changes commit hashes.** Every commit from the first rewritten one onwards gets a new hash. The hashes that
  tracked files pin will then dangle: `FREEZE.json` `harness_commit` (`974d520`), `git_hash` in 406 run manifests,
  and the `git_head` in `results/test/LEDGER.md`. The tags in section 4 must use the new hashes. `ksp.eval_test`
  can no longer verify the A28 harness commit.
- **Publish the commit map.** The script copies filter-repo's old → new map to `docs/history_commit_map.tsv`.
  Commit it, so that readers can still check the pre-registration against the old hashes.
- **`git-filter-repo` is not installed** in the project environment: `pip install git-filter-repo`.
- **Run it only in a fresh clone, never in the project directory.** The script refuses to run where `data/`,
  `external/` or `cache/` exist or the tree is not clean, and filter-repo itself refuses a repository that is not
  a fresh clone.

```
git clone --no-local <project-dir> ksp-rewrite && cd ksp-rewrite
pip install git-filter-repo
bash scripts/publication/rewrite_history.sh --paths --nfs --yes     # add --mailmap FILE / --strip-session if wanted
git log --oneline | head; grep -c . docs/history_commit_map.tsv
git add docs/history_commit_map.tsv && git commit -m "Add commit map of the history rewrite"
# then tag the new hashes (look them up in docs/history_commit_map.tsv) and push to a NEW remote
```

The agents did not run this script.
