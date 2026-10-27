# Design and research record

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

These files are copies of the project's working record. Their content was not edited, with three exceptions: the
flag above was added, absolute paths were replaced by relative paths or placeholders (`<agent-scratchpad>`,
`<plm-env>`, `$HF_HOME`), and the files listed below were renamed. Paths inside the files refer to the original
project directory. Many `cache/`, `results/` and `logs/` files that they cite are not in this repository; the
metric tables are in `results_summary/`.

Renamed when copied: `REPORT.md` → `docs/REPORT_internal.md`; `DATA.md` → `docs/DATA_internal.md`; the top-level
`DESIGN.md`, `LITERATURE.md`, `RESEARCH_TREE.md`, `LESSONS.md` and `DATAINEED.md` → `docs/`.

| File | What it is |
|---|---|
| `REPORT_internal.md` | the final internal report (v3, reviewed in `review_report_v3.md`). Its numbers are the reference values. |
| `DESIGN.md` | the design, frozen at v2, with amendments A1–A34 (§16–17) |
| `LITERATURE.md`, `lit/` | literature synthesis and the four full reviews |
| `RESEARCH_TREE.md` | the research tree, with decision notes and status summaries |
| `LESSONS.md`, `retro/` | process rules, milestone retros (M1–M3) and the final outside-view retro |
| `DELEGATION.md`, `AGENT_BRIEF.md`, `AGENT_SETUP.md` | who owned what, the shared agent brief, and a summary of the agent roles |
| `reviews/` | design reviews, rounds 1 and 2, and the reviewer brief |
| `audit/` | independent audits: CV rounds 1 and 2, the TEST batch, post-TEST work |
| `DATA_internal.md`, `DATAINEED.md`, `data_acquisition.md`, `data_profile.md` | data register with checksums; data requests; acquisition and profiling notes |
| `labels_universe.md`, `core_infra.md`, `G4_RUNBOOK.md` | labels, site universe, splits, evaluation harness, TEST procedure |
| `features_*.md`, `feature_provenance.md` | feature families and their provenance and leakage register |
| `model_r2.md`, `model_r3.md`, `model_r5.md` | the two final model components and the negative R5 cross-encoder |
| `baseline_phosformer.md`, `baseline_selphi2.md` | the retrained baselines B3 and B4 |
| `regimes_labels.md`, `strata_cv.md`, `downstream.md` | training regimes and label variants, CV strata, downstream benchmark |
| `posttest_cv.md`, `posttest_tyr.md` | exploratory work after TEST (did not change the model) |
| `RELEASE_NOTES.md` | what was changed for this release, and how it was checked |
