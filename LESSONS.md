# LESSONS — process rules for running this project efficiently
(One-line rule + evidence. Pruned at each milestone; last prune M2 2026-10-26 17:35 — see docs/retro/M1_retro.md, M2_retro.md.)

1. Give reviewers read-only data access and require a measured number per change — r1: 6/7 found quantified leaks (68% same-protein co-labelled pairs; >50% paralog-contaminated SRC hard negatives).
2. Key by UniProt accession, never symbol — 174/404 Cantley names ≠ gene symbols; "PDK1" = PDPK1 vs PDHK1.
3. Benchmark compute before planning — co-regulation 3.6 s on GPU vs the guessed bottleneck; LightGBM 400k×160 took 15 min.
4. Label-derived features are leaky until cross-fitted by substrate cluster — 68% co-labelled same-protein sites.
5. Audit "unbiased" resources for KS-label leakage — STRING experiments 34% vs 4% controls; iKiP "combined" GMT contains 61% of PSP.
6. Fan out implementers only against frozen interfaces, and get a measured baseline before more review — r2 interface changes (crossfit, store) had to be broadcast to 6 running agents.
7. An auditor (or PI) recomputes headline numbers from saved predictions — audit r1 reproduced everything and caught 4 major issues.
8. Freeze the eval harness and stamp every number with its harness version — A24 landed mid-audit and made all docs stale (eval frozen at commit after A24, 17:35).
9. Pooled-OOF metrics need fold-invariant scores; always report per-fold values — B0-prior pooled M1 0.40 instead of 0.5.
10. Log every config tried and the selection rule before looking; report DEV optimism — ~20 fold-0 configs for R2, 0.0017 picks for R3.
11. Verify artefacts by mtime/hash before claiming they were rebuilt — R2 folds 0/3 canary stores were stale despite the doc saying otherwise.
12. One owner per path, PI-only git, status updated per event, timestamps from `date` — PI notes drifted ~40 min twice.
13. Null controls must match the positives' distribution on every label-free axis — canary v1 (uniform proteins) was separable by protein site density alone (AUROC 0.76).
14. Verify file identities before recording them — paper1/paper2 were swapped in LITERATURE.md until an agent read the PDFs.
15. The node is shared with other users' jobs — budget memory per agent (≤60 GB) and check `free -g` before parallel launches — a re-eval worker was OOM-killed at ~100 GB total (17:50).
