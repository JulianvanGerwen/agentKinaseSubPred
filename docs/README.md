# Design and research record (`docs/`)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

These files are the working record of the study, as the agents wrote them. For v1.0.0 the AI-exercise flag was added
at the top, and absolute paths, scratchpad paths and one internal file-server address were replaced by relative
paths or placeholders (`<scratchpad>`, `<plm-env>`, `<project root>`, `<internal-file-server>:<share>`). Short notes
were added to `features_phospho.md` and `data_profile.md` about the condition tables that are no longer tracked.
Nothing else was changed. Many `cache/` files that they cite are not in the repository, because `cache/` holds large or
restricted intermediates. The main reports (`REPORT.md`, `DESIGN.md`, `LITERATURE.md`, `RESEARCH_TREE.md`,
`LESSONS.md`, `DATA.md`, `DATAINEED.md`) are at the repository root.

| File | What it is |
|---|---|
| `DELEGATION.md`, `AGENT_BRIEF.md`, `AGENT_SETUP.md` | who owned what, the shared brief for sub-agents, and a short summary of the agent roles |
| `reviews/` | design reviews, rounds 1 and 2, and the reviewer brief |
| `audit/` | independent audits: CV rounds 1 and 2, the TEST batch, post-TEST work, with the auditors' scripts |
| `retro/` | milestone retros (M1–M3) and the final outside-view retro |
| `review_report_v3.md` | the independent review of `REPORT.md` v3 |
| `report_parts/` | the sections from which `REPORT.md` was assembled |
| `lit/` | the four full literature reviews behind `LITERATURE.md` |
| `data_acquisition.md`, `data_profile.md` | acquisition and profiling notes for the data in `DATA.md` |
| `labels_universe.md`, `core_infra.md`, `G4_RUNBOOK.md` | labels, site universe, splits, evaluation harness, TEST procedure |
| `features_*.md`, `feature_provenance.md` | feature families, and their provenance and leakage register |
| `model_r2.md`, `model_r3.md`, `model_r5.md` | the two final model components and the negative R5 cross-encoder |
| `baseline_phosformer.md`, `baseline_selphi2.md` | the retrained baselines B3 and B4 |
| `regimes_labels.md`, `strata_cv.md`, `downstream.md` | training regimes and label variants, CV strata, downstream benchmark |
| `posttest_cv.md`, `posttest_tyr.md` | exploratory work after TEST (it did not change the model) |
| `RELEASE_NOTES.md` | what changed for v1.0.0, and how it was checked |
| `PUBLISHING.md`, `GITHUB_RELEASE_v1.0.0.md` | publication steps and open decisions for the researcher; text of the GitHub release |
| `publication_audit_history.md`, `publication_prep_report.md` | audit of the full git history before publication; report of the v1.0.0 preparation |
