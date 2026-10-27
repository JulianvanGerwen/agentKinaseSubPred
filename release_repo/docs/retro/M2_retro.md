> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# M2 retro (T+4.7h, 18:00)

## 1. Stop doing
- Post-freeze harness/metric amendments reacting to results: A23-A27 landed in ~1 h (A24 changed eval.py at 17:01 mid-audit; every doc number is pre-A24). Batch them; stamp each number with harness version.
- Unverified bookkeeping: R2 docs claim "all xfit stores rebuilt with v2" but folds 0/3 are stale (audit M-2); R2/R3 docs carry raw pooled M4 now known inflated (0.821 vs 0.764 in 1-4 bin).
- Chasing marginal DEV gains: ESM-3B windows no gain; F9 SH2 gain is a NaN-indicator artefact (+0.019 [-0.018,0.057]); R3 selection within noise (0.0017); ~20 fold-0 configs tried, counted as rows not configs. Each adds DEV optimism before a single TEST shot.
- Treating contaminated references as results (released Phosformer 0.796 on PSP-trained weights): label "reference only" everywhere.
- Re-narrating status in RESEARCH_TREE (14.5 KB, duplicates DESIGN/docs); DELEGATION Status is again stale (e.g. R2-main, auditor-1 "running").

## 2. What worked
- Independent audit recomputed from saved predictions: all numbers reproduced, and it caught A24 artefact, stale stores, pooled-only canary gate.
- G2 met at 16:27, 4.4 h ahead of 20:50; cuts and gates pre-stated.
- Paired comparisons vs the stronger comparator, with CIs (R2-B1 +0.023 [0.009,0.034]).
- Honest negative findings: SELPHI2 below Cantley on matched pairs; its 0.97 is from random negatives, not split leakage; hub-site/macro-dilution explains Tyr gap.
- Canary v2 replaced a biased sampler before it misled gating.

## 3. Biggest risk (next ~19 h)
Selection optimism and harness drift converging on the single-shot TEST (G4 05:20): many configs tuned on DEV, canary gate within 0.01-0.03 of threshold, M2 (~0.15-0.22) and Tyr unsolved, and Phosformer finetune (ETA 22:30) leaving little slack before G3 02:20. Mitigation: freeze the eval harness now, rebuild stale stores, define the final config-selection rule and number of configs in writing before G3, cap exploration at ~21:00.

## 4. Proposed LESSONS.md (12)
1. Reviewers get read-only data access and must give a measured number per change (r1: 6/7 found quantified leaks).
2. Key by UniProt accession, never symbol (174/404 Cantley names differ).
3. Benchmark compute before planning (GPU 3.6 s vs LightGBM 15 min).
4. Label-derived features are leaky until cross-fitted by substrate cluster (68% co-labelled).
5. Audit "unbiased" resources for KS leakage (STRING 34% vs 4%).
6. Fan out implementers only against frozen interfaces; get a baseline before more review (M1).
7. Auditor/PI recomputes headline numbers from saved predictions; agent claims untrusted (audit: all reproduced, caught 4 majors).
8. Freeze the eval harness; any change is versioned and every number states its version (A24 made all docs stale).
9. Pooled-OOF metrics need fold/kinase-invariant scores; check per-fold values (B0-prior 0.40 vs 0.5; M4 0.821 inflated).
10. Log every configuration tried and the selection rule before looking; report DEV optimism (R3 +0.0017 pick).
11. Verify artifacts by mtime/hash before claiming "rebuilt" (R2 folds 0/3 stale).
12. One owner per path, PI-only git; update status per event; timestamps from `date` (M1 drift 40 min).
