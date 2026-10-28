# LESSONS — process rules for running this project efficiently

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

(One-line rule + evidence. Pruned at each milestone; last prune M3 2026-10-27 00:06 — see docs/retro/M1_retro.md, M2_retro.md, M3_retro.md.)

1. Reviewers get read-only data access and must give a measured number per change (r1: 6/7 found quantified leaks).
2. Key by UniProt accession, never symbol (174/404 Cantley names differ).
3. Treat label-derived features and "unbiased" resources as leaky until cross-fitted by substrate cluster or audited (68% co-labelled; STRING 34% vs 4%; iKiP 61% of PSP).
4. Null controls must match positives on every label-free axis (canary v1 AUROC 0.76 by site density alone).
5. Freeze the eval harness; version every number; log each post-freeze amendment and its reason (A24 made all docs stale).
6. Pooled-OOF metrics need fold-invariant scores; report per-fold (B0-prior 0.40 vs 0.5; M4 0.821 vs 0.764).
7. Write the config log, selection rule and adoption threshold before running (A32 8107905, R5 gate); count configs, not rows.
8. Compare against the strongest comparator with paired CIs, and label contaminated references "reference only" (Phosformer 0.796).
9. An auditor or the PI recomputes headline numbers from saved predictions; agent claims are untrusted (audit r1 caught 4 majors).
10. Verify artefacts by mtime/hash before claiming "rebuilt" (R2 folds 0/3 stale).
11. Fan out implementers only against frozen interfaces; a measured baseline comes before more review (M1: 84 changes, no baseline at T+1.5h).
12. One owner per path, PI-only git, status per event, timestamps from `date`; delete status columns that aren't maintained (drifted at M1, M2, M3).
13. Audit PROMPT.md asks against the leaderboard at every milestone (regimes unrun at 20:00), and budget memory against the job cgroup (128 GiB, not the node's 375 GB): ≤60 GB per agent (OOM kills at 17:50 and ~23:00).
14. After FREEZE, no protocol amendments and no reselection after seeing TEST; deploy only frozen configs (M3 retro: A23–A33 landed ~1/h before the freeze).
