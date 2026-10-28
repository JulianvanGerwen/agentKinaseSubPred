> **This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.**

# ksp v1.0.0: kinase–substrate site predictions from a 24-hour AI-agent study

**Start here:** the research briefing, *Which kinase phosphorylates which site? A fair test, a modest predictor, and
what AI agents did when they ran the study* (`briefing/briefing.html`, also `briefing/briefing.md`). It explains the
study for scientists outside the field, with its results and its limits.

This release scores, for each human protein kinase and human phosphosite (S, T or Y), how likely it is that the
kinase phosphorylates the site directly. The model is a late fusion of a LightGBM pair model (R2) and a neural
kinase × site model pre-trained on the Cantley kinase atlas (R3). The repository holds the code, the evaluation
harness, the retrained baselines, the unit tests and the complete record of the study. Its git history (108 commits
before this release) is the record of what the agents did.

## What is in this release

- **Predictions** (also attached as release assets):
  - `site_top5.tsv.gz`: the top 5 kinases for each of 350,736 observed sites (1,753,680 rows);
  - `kinase_top200.tsv.gz`: the top 200 sites for each of 503 kinases (100,600 rows);
  - `RELEASE_MANIFEST.json`: provenance and checksums. `coverage.tsv` is in the repository.
- **Code**: the `ksp` package (`pip install -e .`), the pipeline scripts, the tests and `examples/read_predictions.py`.
- **Record**: `REPORT.md`, `DESIGN.md` (amendments A1–A34), `RESEARCH_TREE.md`, `LITERATURE.md`, `LESSONS.md`,
  `DATA.md`, `PROMPT.md` (the original brief), and the audits, reviews and retros in `docs/`. The negative
  results (S1 Text of the accompanying paper, link to be added) are in `docs/SUPPLEMENT_failed_approaches.md`.
- No input data. PhosphoSitePlus and the other datasets must be obtained separately (`DATA_SOURCES.md`).

Changes from the last pre-release state: the paths are now configurable (`KSP_DATA_DIR`; unset = the original
layout), and licences, a citation file and publication checks were added. No computation changed. See
`docs/RELEASE_NOTES.md`.

## How to use the predictions

```
python examples/read_predictions.py                      # checks SHA-256, prints a summary
python examples/read_predictions.py --site P04637_S15    # top-5 kinases for TP53 S15
python examples/read_predictions.py --kinase CDK1        # top-20 sites for CDK1
```

```python
import pandas as pd
site5 = pd.read_csv("results/release/site_top5.tsv.gz", sep="\t")
print(site5[site5.site_id == "P04637_S15"])
```

`site_id` is `<UniProt accession>_<residue><position>` on the UniProt 2026_03 canonical sequence. The columns are
described in `results/release/README.md`.

**Rank-1 caveat.** In `site_top5` the rank is by `assignment_prob`. At 25,532 of 350,736 sites (7.3%) this puts a
kinase other than the model's top-scoring kinase at rank 1. For the model's own ordering, sort the five rows by
`prob` or `score_site`.

## Checksums (SHA-256)

From `RELEASE_MANIFEST.json`:

| File | Bytes | SHA-256 |
|---|---|---|
| `kinase_top200.tsv.gz` | 1,952,136 | `98f348ed330784c8dbe52f9d327a14226975a74ab07b1e7612ebb6cf58d64567` |
| `site_top5.tsv.gz` | 34,750,123 | `61e20f51b4eb043365ef492a0f20c31168f30bc6bacd637d40154fd0a8c3729f` |
| full score table (not included, 131,144,410 pairs) | 4,123,189,379 | `cceb32c73be9c092923e8d123bf8dd74a676afe2bea845e002f33f296de839f5` |

Computed for this release: `RELEASE_MANIFEST.json` `3d606802af8e73be0d35d32387b5c6f24e6db8273aa81a007db53882e7bc885e`,
`coverage.tsv` `596751540aef68eb9f221a7f43421c26c5629ca44f5a509c83dcfb879d259062`.

## Results (single sealed TEST batch, FREEZE 721dc197)

M1 (kinase-macro AUROC against hard negatives) 0.794 [0.774, 0.811]; M2 (site-centric Hit@5) 0.363 [0.324, 0.397];
M3 (top-1% recall) 0.363 [0.309, 0.401]. All 9 pre-registered superiority tests against Cantley (B1) and our
retrains of Phosformer (B3) and SELPHI2.0 (B4) meet the "beats" rule. Full tables: `README.md` and `REPORT.md` §1.

## Limitations

- **Tyrosine kinases: specificity is not demonstrated.** Y-kinase M1 on TEST is 0.595, barely above the
  site-popularity floor (0.585). Do not trust Tyr predictions.
- **There is no reliable high-confidence tier.** No site reaches an assignment probability of 0.8. Use the output
  as a ranked shortlist, not as a call.
- **The rise of M2 from CV to TEST comes from the set of kinases**; on the same kinases CV M2 is 0.348.
- **The margin over Cantley is partly a weak-baseline draw**, and the baselines are our retrains, not the
  published models.
- **One sealed TEST set, spent early**, after about 90 CV configurations; the work after TEST did not change the
  model.
- **Probabilities are lower bounds** under positive–unlabelled labels. The deployment model's calibrators were not
  tested again after retraining on DEV+TEST.

## Licences

- Code: MIT (`LICENSE`).
- Predictions and derived tables: **CC BY-NC-SA 3.0**, because they are derived from PhosphoSitePlus labels.
  Non-commercial use only, share-alike, attribution to PhosphoSitePlus (www.phosphosite.org). See
  `results/release/LICENSE-PREDICTIONS.md`.
- Documentation: CC BY 4.0 unless a file says otherwise. Third-party terms: `NOTICE`.

## Citation

`CITATION.cff`. Authors: [authors to be added by the publishing researcher].
