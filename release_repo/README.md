# ksp: kinase–substrate site prediction (MVP release)

> **This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.**

This repository accompanies a manuscript prepared for *PLOS Computational Biology* (see `manuscript/`). It is a
minimum viable release. It contains the released predictions, the code of the final method, the evaluation
harness and the retrained baselines, the unit tests, and the design and audit record. It does not contain any
input data. The full pipeline cannot be re-run without data that you must obtain yourself under their licences
(`DATA_SOURCES.md`), and parts of it need a GPU.

## What the project is

For each pair of a human protein kinase and a human phosphosite (S, T or Y), the model scores how likely the
kinase phosphorylates the site directly. Training and test labels are experimentally supported kinase–substrate
pairs from PhosphoSitePlus (PSP). The problem is treated as positive–unlabelled learning.

The evaluation design:
- Substrate-protein clusters (homology components) are held out; 20.3% of 11,609 clusters were sealed as TEST.
- Every feature derived from labels is cross-fitted by substrate cluster.
- Hard negatives are sites labelled to other, non-paralogous kinases.
- Endpoints weight each kinase equally (macro averages), with strata for understudied, Tyr, priming and atypical
  kinases.
- Superiority tests against retrained baselines were pre-registered, and the sealed TEST set was read in a single
  batch.

The final method is a late fusion (z-scored logits) of two models:
- **R2**, a global LightGBM pair model. Its features are motif, structure and ESM-2 embeddings, protein interaction
  screens, localisation, co-expression, phosphoproteomics and kinase descriptors (5-seed ensemble).
- **R3**, a neural kinase × site model that is pre-trained by distillation of the Cantley kinase atlas.

## Headline results (from `docs/REPORT_internal.md` §1)

TEST, single batch (FREEZE 721dc197), 95% cluster-bootstrap CIs. Each Δ is the final model minus the baseline in the
matched-universe primary test (B = 10,000, Holm correction over 9 tests).

| Endpoint | Final | Cantley (B1) | Phosformer (B3, our retrain) | SELPHI2.0 (B4, our retrain) |
|---|---|---|---|---|
| M1: kinase-macro AUROC against hard negatives | **0.794** [0.774, 0.811] | 0.714; Δ +0.075 [0.046, 0.099] | 0.708; Δ +0.086 [0.059, 0.105] | 0.705; Δ +0.089 [0.063, 0.117] |
| M2: site-centric Hit@5 | **0.363** [0.324, 0.397] | 0.189; Δ +0.179 [0.147, 0.213] | 0.163; Δ +0.199 [0.151, 0.229] | 0.233; Δ +0.130 [0.091, 0.156] |
| M3: top-1% recall | **0.363** [0.309, 0.401] | 0.161; Δ +0.205 [0.165, 0.235] | 0.104; Δ +0.259 [0.213, 0.291] | 0.304; Δ +0.059 [0.010, 0.120] |

All 9 primary tests meet the pre-registered "beats" rule (8 at the bootstrap floor p < 0.0002; B4 on M3 at
p = 0.018). On CV (pooled 5-fold, 203 kinases) the final model reaches M1 0.786 against 0.717 for Cantley.
Sources: `results_summary/test_721dc1977366/` and `results_summary/TEST_LEDGER.md`.

## Limitations (read these before using the numbers or the predictions)

1. **Tyrosine kinases: specificity is not demonstrated (not trusted).** On TEST, Y kinases reach M1 0.595, barely
   above the site-popularity floor (0.585). The final model beats that floor on only 10 of 23 Y kinases.
2. **There is no reliable high-confidence tier.** On labelled TEST S/T sites the ≥0.5 assignment tier is
   calibrated (0.615 claimed vs 0.603 observed, n = 345). There is no ≥0.8 tier, and on TEST no Y site can exceed
   0.445. Use the "which kinase?" output as a ranked shortlist, not as a call.
3. **The rise of M2 from CV to TEST is composition, not improvement.** On the same 117 kinases, CV M2 is 0.348.
4. **The margin over Cantley on TEST is partly a weak-baseline draw.** On the same kinases Cantley fell by 0.025
   from CV to TEST, mostly on Y kinases.
5. **M2 is still weak for understudied kinases**, and TEST holds only 13 kinases with 5–19 labels.
6. **The baselines are our retrains on this benchmark**, not the published models. The Phosformer fidelity check
   passed. SELPHI2.0 uses its default values for missing site features on 66.6% of its TEST pairs. The B4 M3 win
   is modest.
7. **Priming.** The advantage on primed sites holds only against Cantley without priming awareness. The
   priming-aware Cantley arm was not part of the frozen TEST comparison.
8. **One sealed TEST set, spent early.** About 90 configurations were tried on CV. Work after TEST
   (`docs/posttest_*.md`) was exploratory and did not change the released model.
9. **Labels are PU and curated.** Probabilities estimate P(curated PSP-type positive | score) and are lower bounds.
   Coverage on unlabelled sites cannot be checked for precision.
10. **Deployment model.** The released predictions come from the frozen configuration retrained on DEV+TEST
    labels. Its calibrators were fitted on DEV out-of-fold predictions and were not re-tested.

The full list is in `docs/REPORT_internal.md` §10 and in the audit reports under `docs/audit/`.

## Repository layout

```
README.md                 this file
PROMPT.md                 the human-written brief, verbatim
DATA_SOURCES.md           every dataset: source, version, licence, how to obtain it, expected local layout
LICENSE                   MIT, for the code
NOTICE                    third-party data and model attributions and licences
CITATION.cff              citation metadata (authors to be added)
pyproject.toml            package metadata (pip install -e .)
requirements.txt          versions used for the released results
predictions/              released predictions (CC BY-NC-SA 3.0, derived from PhosphoSitePlus) + manifest
results_summary/          small metric tables: TEST CIs, superiority tests, strata, CV leaderboard, coverage
examples/                 read_predictions.py: a runnable example that needs only pandas
src/ksp/                  the package: labels, splits, features, models (R2, R3), baselines, evaluation harness
scripts/                  pipeline entry points (see scripts/README.md)
config/                   small curation tables (kinase curation, homonyms, complex subunits, biological controls)
tests/                    unit tests; tests that need project data skip with a message
docs/                     design record, literature, research tree, lessons, audits, reviews, retros, model docs
manuscript/               manuscript, supplement, references, figures (with make_figures.py) and number maps
```

## Quick start

### (a) Use the released predictions

You need Python 3 and pandas. Nothing else is required.

```
python examples/read_predictions.py                      # checks SHA-256, prints a summary and two look-ups
python examples/read_predictions.py --site P04637_S15    # top-5 kinases for TP53 S15
python examples/read_predictions.py --kinase CDK1        # top-20 sites for CDK1
```

Or directly:

```python
import pandas as pd
site5 = pd.read_csv("predictions/site_top5.tsv.gz", sep="\t")        # 350,736 sites x top-5 kinases
kin200 = pd.read_csv("predictions/kinase_top200.tsv.gz", sep="\t")   # 503 kinases x top-200 sites
print(site5[site5.site_id == "P04637_S15"])
```

Column definitions, score meanings and the caveats for `assignment_prob` are in `predictions/README.md`. One
caveat matters for every user: in `site_top5`, the rank is by `assignment_prob`. At about 7% of sites this puts a
kinase other than the model's top-scoring kinase in rank 1. For the model's own ordering, sort the five rows by
`prob` or `score_site` (details in `predictions/README.md`).

### (b) Reproduce the pipeline (conceptually)

You cannot re-run the pipeline from this repository alone. You need:
- your own PhosphoSitePlus download under its licence, and the other datasets in `DATA_SOURCES.md`;
- about 60 GB of RAM, 8 or more cores, and a CUDA GPU for ESM-2 embeddings, R3 and the Phosformer baseline;
- MMseqs2 for the homology clusters (`KSP_MMSEQS`, or `mmseqs` on `PATH`).

Set up the code:

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt          # or: pip install -e ".[gpu,full,test]"
pip install -e .
python -m pytest -q                      # tests that need project data are skipped, with the reason
```

All paths go through `src/ksp/paths.py`. Put the inputs under one project directory with the layout in
`DATA_SOURCES.md` and point `KSP_DATA_DIR` at it. If `KSP_DATA_DIR` is unset, the repository root is used, which is
the layout in which the study was run. Outputs (`cache/`, `results/`, `logs/`) are written to the same directory.

The order of the stages is given in `scripts/README.md` and `docs/REPORT_internal.md` §12. In short: build the
kinase list, site universe, labels and frozen splits; build the feature stores; run 5-fold CV for R2 and R3 and
fuse them; build the TEST deliverables and the FREEZE; run the single TEST batch; retrain on DEV+TEST and write the
release files. The splits depend on the PSP release and on UniProt 2026_03. With another release you will get a
different split, so the numbers will differ.

## What has been verified, and what has not

Verified for this release:
- The SHA-256 of `predictions/kinase_top200.tsv.gz` and `predictions/site_top5.tsv.gz` matches
  `predictions/RELEASE_MANIFEST.json`.
- The unit tests that need no project data pass in a clean copy of this repository. The tests that need data are
  skipped with a message (counts are in `docs/RELEASE_NOTES.md`).
- The path changes are configuration only. With `KSP_DATA_DIR` pointing at the original project directory, every
  module-level data path resolves to the same file as before. The test suite then gives 97 passed and 3 skipped
  (two tests need a label table that the original project deletes after sealing, and one needs the original git
  history); one test that calls git was not run. Details are in `docs/RELEASE_NOTES.md`.

Not verified for this release:
- No feature, model, TEST or release step was re-run from this repository.
- The five harness files frozen at commit 974d520 (`eval.py`, `labels.py`, `splits.py`, `crossfit.py`,
  `canary.py`) differ from that commit, but only in how paths are resolved. The checks in `ksp.eval_test`
  against the frozen code therefore do not pass on this copy. The TEST batch is spent and cannot be re-run anyway.
- Scripts for exploratory and negative branches (R5 cross-encoder, regime sweep, post-TEST experiments,
  diagnostics) are not included; their reports are in `docs/`.

## Licences

- **Code** (`src/`, `scripts/`, `tests/`, `examples/`): MIT, see `LICENSE`.
- **Predictions** (`predictions/`) and the metric tables in `results_summary/`: derived from PhosphoSitePlus labels,
  so they are released under **CC BY-NC-SA 3.0** (non-commercial, share-alike, attribution to PhosphoSitePlus,
  www.phosphosite.org). See `predictions/LICENSE-PREDICTIONS.md`.
- **Documentation** (`docs/`, this README): CC BY 4.0, unless a file says otherwise.
- **Third-party data and models** are not included. Their licences, including the non-commercial terms of the
  Cantley kinase-library matrices and the Phosformer weights, are listed in `NOTICE` and `DATA_SOURCES.md`.

## Citation

See `CITATION.cff`. Authors: [authors to be added by the publishing researcher].
