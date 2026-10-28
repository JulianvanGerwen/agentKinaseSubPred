# Evaluating agentic AI for kinase–substrate prediction

**(Human-written text starts)**

This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (https://github.com/evocellnet/ai-hackathon/). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.

The goal of this specific project was to create a predictive model for substrates of human protein kinases. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 24 hours on a computational cluster with 128 Gb RAM and 4 RTX4090 GPUs.

Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It is the working directory used by the agent. It also includes code for implementation of the predictive model and a AI-written briefing ([`briefing/briefing.pdf`](briefing/briefing.pdf) or [`briefing/briefing.html`](briefing/briefing.html)) detailing the final method, the path followed to reach it, and benchmarking against other methods.

**(Human-written text ends)**

## Start here: the research briefing

**[Which kinase phosphorylates which site? A fair test, a modest predictor, and what AI agents did when they ran the study](briefing/briefing.html)**
(also as Markdown: [`briefing/briefing.md`](briefing/briefing.md); GitHub shows the HTML file as source, so download it
or open the Markdown version)

## What this repository is

This is the working directory of a 24-hour study that AI agents ran from a human-written brief (`PROMPT.md`,
verbatim), without the input data and the paper draft. The goal was a predictor of which human protein kinase directly phosphorylates which human
S, T or Y site, tested against established methods without data leakage.

The repository is both the software release (v1.0.0) and the record of how the agents worked:
- **The git history is the agentic record.** The 108 commits before this release are kept in full. Each commit shows what the agents did and when. The pre-registration is pinned to commit
  hashes in it (see "Reproduce the frozen state").
- The reports, design amendments, reviews, audits and retros are the agents' own documents, kept as written. For
  v1.0.0 the AI-exercise flag and a few release notes were added, and absolute paths were replaced by placeholders
  (`docs/README.md`).
- The research briefing (`briefing/`) and an accompanying paper (link to be added) were written from these
  records in the publication phase. The paper's supplement, which records the approaches that failed, is in
  `docs/SUPPLEMENT_failed_approaches.md`.

For each pair of a human protein kinase and a human phosphosite, the model scores how likely it is that the kinase
phosphorylates the site directly. Labels are experimentally supported kinase–substrate pairs from PhosphoSitePlus
(PSP). The problem is treated as positive–unlabelled learning. The final method is a late fusion (z-scored logits)
of two models:
- **R2**, a global LightGBM pair model. Its features are motif, structure and ESM-2 embeddings, protein interaction
  screens, localisation, co-expression, phosphoproteomics and kinase descriptors (5-seed ensemble).
- **R3**, a neural kinase × site model that is pre-trained by distillation of the Cantley kinase atlas.

The evaluation holds out substrate-protein clusters (11,609 homology clusters, 20.3% sealed as TEST) and cross-fits
every label-derived feature by cluster. Hard negatives are sites labelled to other, non-paralogous kinases. Every
kinase counts equally in the endpoints. Superiority tests against retrained baselines were pre-registered, and the
sealed TEST set was read in a single batch.

## Headline results

These numbers are from `REPORT.md` §1 (v3, independently reviewed in `docs/review_report_v3.md`).

TEST, single batch (FREEZE 721dc197), with 95% cluster-bootstrap CIs. Each Δ is the final model minus the baseline in
the matched-universe primary test (B = 10,000, Holm correction over 9 tests).

| Endpoint | Final | Cantley (B1) | Phosformer (B3, our retrain) | SELPHI2.0 (B4, our retrain) |
|---|---|---|---|---|
| M1: kinase-macro AUROC against hard negatives | **0.794** [0.774, 0.811] | 0.714; Δ +0.075 [0.046, 0.099] | 0.708; Δ +0.086 [0.059, 0.105] | 0.705; Δ +0.089 [0.063, 0.117] |
| M2: site-centric Hit@5 | **0.363** [0.324, 0.397] | 0.189; Δ +0.179 [0.147, 0.213] | 0.163; Δ +0.199 [0.151, 0.229] | 0.233; Δ +0.130 [0.091, 0.156] |
| M3: top-1% recall | **0.363** [0.309, 0.401] | 0.161; Δ +0.205 [0.165, 0.235] | 0.104; Δ +0.259 [0.213, 0.291] | 0.304; Δ +0.059 [0.010, 0.120] |

All 9 primary tests meet the pre-registered "beats" rule: 8 at the bootstrap floor (p < 0.0002), and B4 on M3 at
p = 0.018. On pooled 5-fold CV (203 kinases) the final model reaches M1 0.786, against 0.717 for Cantley. Sources:
`results/test/LEDGER.md`, `results/test/721dc1977366/` and `results/cv/LEADERBOARD.md`.

## Limitations

Read these before you use the numbers or the predictions. The full list is in `REPORT.md` §10 and in the audit
reports in `docs/audit/`.

1. **Tyrosine kinases: specificity is not demonstrated, and Tyr results should not be trusted.** On TEST, Y kinases
   reach M1 0.595, barely above the site-popularity floor (0.585). The final model beats that floor on only 10 of
   23 Y kinases.
2. **There is no reliable high-confidence tier.** On labelled TEST S/T sites the ≥0.5 assignment tier is
   calibrated (0.615 claimed vs 0.603 observed, n = 345). There is no ≥0.8 tier, and on TEST no Y site can exceed
   0.445. Use the "which kinase?" output as a ranked shortlist, not as a call.
3. **The rise of M2 from CV to TEST comes from the set of kinases, not from a better model.** On the same 117
   kinases, CV M2 is 0.348.
4. **The margin over Cantley on TEST is partly a weak-baseline draw.** On the same kinases Cantley fell by 0.025
   from CV to TEST, mostly on Y kinases.
5. **M2 is still weak for understudied kinases.** TEST holds only 13 kinases with 5–19 labels.
6. **The baselines are our retrains on this benchmark**, not the published models. The Phosformer fidelity check
   passed. SELPHI2.0 uses its default values for missing site features on 66.6% of its TEST pairs. The B4 win on M3
   is modest.
7. **Priming.** The advantage on primed sites holds only against Cantley without priming awareness. The
   priming-aware Cantley arm was not part of the frozen TEST comparison.
8. **One sealed TEST set, spent early.** About 90 configurations were tried on CV. The work after TEST
   (`docs/posttest_*.md`) was exploratory and did not change the released model.
9. **Labels are positive–unlabelled and curated.** Probabilities estimate P(curated PSP-type positive | score) and
   are lower bounds. The precision of predictions for unlabelled sites cannot be checked.
10. **Deployment model.** The released predictions come from the frozen configuration retrained on DEV+TEST labels.
    Its calibrators were fitted on DEV out-of-fold predictions and were not tested again.

**Caveat on `assignment_prob` rank 1.** In `results/release/site_top5.tsv.gz` the rank is by `assignment_prob`. At
25,532 of 350,736 sites (7.3%) this puts a kinase other than the model's top-scoring kinase at rank 1, and it lets
19 Y sites reach 0.55, although no Y site's top kinase reaches 0.5. For the model's own ordering, sort the five rows
by `prob` or `score_site`. The released files were not changed; details are in `results/release/README.md`.

## Guide to files

Main reports (repository root):

| File | What it is |
|---|---|
| `REPORT.md` | the final report (v3): methods, all results, limitations, follow-ups. Its numbers are the reference values |
| `RESEARCH_TREE.md` | the research tree: every branch tried, with status, decision notes and dated summaries |
| `DESIGN.md` | the design, frozen at v2, with amendments A1–A34 |
| `LITERATURE.md` | the literature synthesis (full reviews in `docs/lit/`) |
| `LESSONS.md` | process rules that the PI kept for running the project, each with its evidence |
| `DATA.md` | the data register: every file used, with source, version, licence and SHA-256 |
| `DATAINEED.md` | data requests from the agents to the human researcher, and their outcome |
| `PROMPT.md` | the original human-written brief, verbatim |
| `DATA_SOURCES.md` | how to obtain each dataset, and the directory layout the code expects |
| `FREEZE.json` | the pre-registered TEST configuration (do not edit; it is checked against git) |

Directories:

| Path | What it holds |
|---|---|
| `docs/DELEGATION.md` | which agent owned which files; `docs/AGENT_BRIEF.md` and `docs/AGENT_SETUP.md` describe the agent set-up |
| `docs/audit/` | independent audits of CV rounds 1 and 2, the TEST batch and the post-TEST work |
| `docs/retro/` | milestone retros (M1–M3) and the final outside-view retro |
| `docs/reviews/` | design reviews (rounds 1 and 2) |
| `docs/` (other) | module, feature, model and baseline notes; index in `docs/README.md` |
| `briefing/` | the research briefing (Markdown and HTML), its figures and the scripts that build them (`figures/make_figures.py`, `build_html.py`), its number map, and the brief, review and change log behind it |
| `docs/SUPPLEMENT_failed_approaches.md` | approaches that were tried and did not work, and other negative results (S1 Text of the accompanying paper) |
| `results/release/` | the released predictions, their manifest, licence and README |
| `results/test/` | the single TEST batch: ledger and metric tables |
| `results/cv/`, `results/runs/` | CV leaderboard and registry; one manifest and metric set per run |
| `results/` (other) | derived tables of the pipeline (labels, splits, strata, features, models, downstream) |
| `src/ksp/` | the Python package: labels, splits, features, models (R2, R3, and the negative R5 and regime branches), baselines, evaluation harness |
| `scripts/` | pipeline entry points in stage order, plus exploratory scripts (`scripts/README.md`) |
| `config/` | small curation tables (kinase curation, homonyms, complex subunits, biological controls) |
| `tests/` | unit tests; tests that need project data are skipped with the reason |
| `examples/` | `read_predictions.py`, a runnable example that needs only pandas |
| `logs/` | checksum lists of the downloads and short logs of the run |

Inputs are not in the repository: `data/` (provided data, including PSP), `external/` (downloads), `cache/`
(intermediates, tens of GB) and `sealed/` (TEST labels) are in `.gitignore`.

## Quick start

### Use the released predictions

You need Python 3 and pandas. Nothing else is required.

```
python examples/read_predictions.py                      # checks SHA-256, prints a summary and two look-ups
python examples/read_predictions.py --site P04637_S15    # top-5 kinases for TP53 S15
python examples/read_predictions.py --kinase CDK1        # top-20 sites for CDK1
```

Or directly:

```python
import pandas as pd
site5 = pd.read_csv("results/release/site_top5.tsv.gz", sep="\t")        # 350,736 sites x top-5 kinases
kin200 = pd.read_csv("results/release/kinase_top200.tsv.gz", sep="\t")   # 503 kinases x top-200 sites
print(site5[site5.site_id == "P04637_S15"])
```

Column definitions, score meanings and the caveats for `assignment_prob` are in `results/release/README.md`.

### Reproduce the pipeline (conceptually)

The code cannot re-run the pipeline without data that you must obtain yourself. You need:
- your own PhosphoSitePlus download under its licence, and the other datasets listed in `DATA_SOURCES.md`;
- about 60 GB of RAM, 8 or more cores, and a CUDA GPU for ESM-2 embeddings, R3 and the Phosformer baseline. The
  study used a shared node with 4 GPUs for 24 hours;
- MMseqs2 for the homology clusters (set `KSP_MMSEQS` to the binary, or put `mmseqs` on `PATH`).

Set up the code:

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt          # exact versions of the study environment (CUDA 12.4 torch build)
pip install -e .                         # or: pip install -e ".[gpu,full,test]"
python -m pytest -q                      # tests that need project data are skipped, with the reason
```

All paths go through `src/ksp/paths.py`:
- `KSP_DATA_DIR` is the project data and working directory. It holds `data/`, `external/`, `cache/`, `results/`,
  `logs/`, `sealed/` and `FREEZE.json`, in the layout given in `DATA_SOURCES.md` §4. **If it is unset, the
  repository root is used. This is the layout in which the study was run, and every data path resolves exactly as
  it did then.** The one exception is the MMseqs2 binary, which was an absolute path on the study machine.
- `KSP_MMSEQS` (MMseqs2 binary) and `KSP_PYTHON` (interpreter for the shell wrappers) are optional.

The order of the stages is in `scripts/README.md` and `REPORT.md` §12. In short: build the kinase list, site
universe, labels and frozen splits; build the feature stores; run 5-fold CV for R2 and R3 and fuse them; build the
TEST deliverables and the FREEZE; run the single TEST batch; retrain on DEV+TEST and write the release files. The
splits depend on the PSP release and on UniProt 2026_03. With other releases you will get a different split, so the
numbers will differ.

### Reproduce the frozen state exactly

The pre-registration is tied to git commits, so use the tags:

| Tag | Commit | State |
|---|---|---|
| `v0.9-freeze` | `5fc7870` | `FREEZE.json` committed, before the TEST batch ran |
| `v0.9-test-batch` | `c3e2a14` | the code that ran the single TEST batch (the `git_head` recorded in `results/test/LEDGER.md`) |
| `v1.0.0` | this release | the same computation, with configurable paths, licences and publication files |

`FREEZE.json` names the harness commit `974d520` (amendment A28). `ksp.eval_test` checks that the frozen harness
files (`eval.py`, `labels.py`, `splits.py`, `crossfit.py`, `canary.py`) are unchanged since that commit. At
`v0.9-test-batch` they are unchanged. At `v1.0.0` they differ only in how they locate files. The test
`test_frozen_files_changed_only_paths` checks this, and the real-mode TEST run now refuses to start. The TEST batch
is spent and must not be repeated on the same split. The 406 run manifests in `results/runs/` record the commit of
each run in `git_hash`.

## What was verified for v1.0.0

- With data present, the full test suite passes; in a fresh copy without data, the tests that need data are
  skipped with a message. Counts are in `docs/RELEASE_NOTES.md`.
- The SHA-256 of both prediction files matches `results/release/RELEASE_MANIFEST.json`.
- No feature, model, TEST or release step was re-run for v1.0.0.

## Licences

- **Code** (`src/`, `scripts/`, `tests/`, `examples/`): MIT, see `LICENSE`.
- **Predictions and derived tables** (`results/`): derived from PhosphoSitePlus labels, so they are released under
  **CC BY-NC-SA 3.0** (non-commercial, share-alike, attribution to PhosphoSitePlus, www.phosphosite.org). See
  `results/release/LICENSE-PREDICTIONS.md`. Some evaluation tables list labelled kinase–site pairs, for example the
  TEST tables `results/test/721dc1977366/*_site_level.tsv`.
- **Documentation** (`docs/`, `briefing/`, the reports): CC BY 4.0, unless a file says otherwise.
- **Third-party data and models** are not included. Their licences are listed in `NOTICE` and `DATA_SOURCES.md`.
  They include the non-commercial terms of the Cantley kinase-library matrices and the Phosformer weights.

## Citation

See `CITATION.cff`. Authors: [authors to be added by the publishing researcher].
