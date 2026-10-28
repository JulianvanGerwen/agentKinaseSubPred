# Agent set-up

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

This is a short factual description of how the work was organised. It is taken from `docs/DELEGATION.md`,
`docs/AGENT_BRIEF.md`, `docs/reviews/REVIEW_BRIEF.md`, the audit reports and the retros (`docs/retro/`). The
delegation table does not record which model (Opus 5.5 or Sonnet 5.5) each sub-agent ran on.

## Brief and time box

- The human researcher wrote the brief (`PROMPT.md`). They supplied the data in `data/`, including PhosphoSitePlus
  access, and answered data requests through `DATAINEED.md`.
- The run had a 24-hour budget: T0 = 2026-10-26 13:20 CET, deadline 2026-10-27 13:20 CET (`docs/retro/final_retro.md`).
- Hardware: one shared node with 64 cores, 4× RTX 4090 GPUs, and a 128 GiB memory limit (cgroup). Agents were
  capped at ≤60 GB resident memory and given explicit core and GPU allocations.

## Roles

**Lead agent (the "PI").**
- It kept the research record: `DESIGN.md` with its amendments A1–A34, `RESEARCH_TREE.md`, `LESSONS.md`,
  `DATAINEED.md`, the data register (`DATA.md`) and `docs/DELEGATION.md`. It also owned the
  canonical index module (`src/ksp/index.py`).
- It delegated tasks, reviewed and integrated the results, and decided on adoption against pre-registered rules.
- It was the only agent that ran git. It wrote `REPORT.md`.

**Sub-agents with a write scope.** Each implementation agent could write only in its assigned paths. Shared modules
were changed only by their owner. The tasks in `docs/DELEGATION.md` were:
- labels and site universe;
- structure and ESM features;
- core infrastructure (splits, evaluation harness, canaries, run manifests, TEST harness, tests);
- the feature families: motif, context, phosphoproteomics, similarity to training substrates, and Tyr;
- the models: R2 (LightGBM), R3 (neural) and R5 (cross-encoder; negative);
- the baselines: simple baselines, B3 Phosformer and B4 SELPHI2.0;
- strata, training regimes and label variants;
- the downstream benchmark.

Every agent received the same rules (`docs/AGENT_BRIEF.md`):
- labels only from experimentally verified annotations;
- splits by protein and homology, never by kinase–substrate pair;
- never change metrics, labels or tests to improve results;
- fixed seeds;
- no restricted data in git;
- downloaded files treated as untrusted;
- summaries back to the PI, with details in files.

**Literature agents.** Four agents wrote the literature reviews in `docs/lit/`. The PI condensed them into
`LITERATURE.md`.

**Design reviewers.** Independent reviewer agents critiqued `DESIGN.md` under `docs/reviews/REVIEW_BRIEF.md`. Each
proposed concrete changes within scope, or stated that none were needed. Round 1 had 7 reviewers: external PI,
prior work, statistics, data leakage, algorithms and compute, AI-agent failure modes, and kinase biology. Round 2
had 4: external PI, statistics, data leakage, and AI-agent failure modes. A third round was not needed.

**Auditors.** Results auditors were agents that had produced none of the audited work. They read everything and
changed nothing, used no git, and never touched `sealed/`. Each audit is in `docs/audit/`:
- `audit_cv_r1.md`: R3, simple baselines and early R2, on 2026-10-26;
- `audit_cv_r2.md`: the final R2 candidate and the fusion, before FREEZE;
- `audit_test.md`: an independent reproduction of the single TEST batch;
- `audit_posttest.md`: the exploratory work after TEST.

An independent reviewer also checked `REPORT.md` v3 against the committed artefacts (`docs/review_report_v3.md`).

**Retro agents.** At each milestone a fresh agent was asked what to stop doing and what worked
(`docs/retro/M1_retro.md` to `M3_retro.md`). A fresh agent with no prior involvement wrote the final outside-view
retro (`docs/retro/final_retro.md`).

## Safeguards against agent failure modes

- The TEST labels were sealed at split time (`sealed/`, SHA-256 in `SEAL.json`). Only `ksp.eval_test` and, after
  TEST was completed, `ksp.deploy` may read them. `tests/test_no_test_access.py` enforces this.
- TEST ran as one pre-registered batch: a FREEZE file pinned every deliverable by SHA-256, and a ledger recorded the
  STARTED and COMPLETED events.
- Every label-derived feature goes through `ksp.crossfit`. Canary pseudo-kinases with random labels serve as leak
  detectors (`ksp.canary`).
- Adoption rules were fixed before each experiment. Negative results were accepted: R5, A32a and A32b, and the three
  post-TEST experiments.

## What went wrong (from the final retro)

- The delivered model was finished about 6 h into the 24 h budget. Later CV runs did not beat it.
- TEST was spent 5.3 h before its planned gate. No second sealed shard was kept, so the post-TEST experiments
  could not be adopted.
- Review of the TEST harness code came late (fixes landed 34 min before FREEZE). The priming-aware Cantley arm was
  known to be stronger on primed sites but was left out of FREEZE.
- The regime sweep used 38% of CV compute and produced no adoptions. Two out-of-memory kills happened because
  memory was budgeted against the node instead of the cgroup.
- The status column of `docs/DELEGATION.md` was not kept up to date.

See `docs/retro/final_retro.md` §3–5 for the evidence and the proposed rules.
