# Result summaries (metrics only)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

These tables contain metrics and counts only. They contain no kinase–site labels. They are derived from
PhosphoSitePlus labels, so they are released under CC BY-NC-SA 3.0, like the predictions
(`predictions/LICENSE-PREDICTIONS.md`).

| Path | Content |
|---|---|
| `TEST_LEDGER.md` | ledger of the single sealed TEST batch (FREEZE 721dc197, attempt 1, COMPLETED). The absolute output path in the original ledger was replaced by the relative path `results/test/721dc1977366`; nothing else was changed. |
| `test_721dc1977366/summary.md`, `summary.json` | standalone TEST metrics for every configuration; primary and secondary tests |
| `test_721dc1977366/superiority.tsv` | the 9 primary tests (paired cluster bootstrap, B = 10,000, matched universe, Holm over 9) |
| `test_721dc1977366/secondary.tsv` | secondary comparisons (Holm over 17, separately) |
| `test_721dc1977366/<config>_ci.tsv`, `<config>_ci_kinases.tsv` | 95% CIs, resampling substrate clusters or kinases (B = 1,000) |
| `test_721dc1977366/<config>_per_kinase.tsv` | per-kinase metrics and per-kinase counts of TEST positives |
| `test_721dc1977366/<config>_strata.tsv`, `<config>_ext_tiers.tsv` | strata (descriptive; no CIs) and external curated tiers |
| `test_721dc1977366/<config>_t3.json`, `<config>_t3_per_kinase.tsv` | T3: unseen kinase family and unseen protein |
| `test_721dc1977366/final_reliability.tsv` | A17 assignment reliability on labelled TEST sites (final model only; the other configurations' files are empty) |
| `test_721dc1977366/preflight.json` | row counts and coverage checks of the submitted prediction files |
| `test_721dc1977366/<config>.json` | full metric record of each configuration |
| `cv/LEADERBOARD.md` | CV leaderboard of every registered configuration (pooled 5-fold DEV OOF) |
| `coverage.tsv` | coverage of the deployment model over U_obs (same file as `predictions/coverage.tsv`) |

Configurations: `final` (R2+R3 fusion), `R2`, `R3` (components), `B1_i` and `B1_ii` (Cantley as published and
universe-recalibrated), `B2` (learned PSSM), `B3` (Phosformer, our retrain), `B4` (SELPHI2.0, our retrain), `B7`
(simple logistic model), `B0_prior`, `B0_pop`, `B0_miss` (controls).

**Left out on purpose:** the 12 `<config>_site_level.tsv` files of the TEST batch. They list the labelled TEST
kinase–site pairs, which are PhosphoSitePlus labels. Use `docs/REPORT_internal.md` §7 to interpret the numbers. Its
caveats apply: Tyr not trusted; no ≥0.8 tier; the M2 rise on TEST is composition; the Cantley margin is partly a
weak-baseline draw.
