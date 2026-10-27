> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# M1 retro (T+1.5h)

## 1. Stop doing
- Over-specifying before any measured result: DESIGN.md is 28 KB / 15 sections; r1 produced 7 reviews (~2,000 lines) and 84 change requests (~80 accepted). No baseline has run at T+1.5h, and G1 (16:50) needs one. Cap design at what gates G1-G2 need; defer the rest (DESIGN §0 already lists deferred items).
- Review round 2 running in parallel with 6 implementers: any accepted r2 change may force rework of code written against v1. Gate implementers on frozen interfaces only (splits, eval API).
- Stale status bookkeeping: RESEARCH_TREE marks "[~] 2 profiling" and "[~] 3 DESIGN" although their output files exist; DELEGATION.md Status column is blank for 5 of 7 rows. Update on completion or drop the column.
- Reviewer sprawl: 7 reviewers overlap (prior_work 407 lines, algorithms 429). Use <=4 with distinct remits and a length cap.
- Unsupported lessons: "<1 h M0->M1" is not evidence that parallelism helped (no counterfactual).
- Prose duplication: facts repeated in RESEARCH_TREE, DESIGN, docs/lit and LITERATURE; link instead.

## 2. What worked
- Reviewers with read-only data access found decisive issues: 68% same-protein co-labelled pairs (leave-site-out leakage), >50% paralog hard negatives in SRC family, STRING leak (34% vs ~4% controls).
- Parallel literature/profiling/acquisition (12.6 GB, 2 waves) with a leakage audit under results/leakage.
- Hard gates G1-G6 with pre-stated scope cuts; sealed TEST, FREEZE.json, canaries; PI recomputes metrics rather than trusting agents.
- File-ownership table + "only PI runs git" avoids NFS/concurrency collisions.
- Measuring compute before planning (co-regulation 3.6 s on GPU vs LightGBM 15 min).

## 3. Proposed LESSONS.md
1. Give reviewers read-only data access; require a measured number per requested change (r1: 6/7 did; leaks found were invisible by reading).
2. Key everything by UniProt accession (Cantley names: 174/404 differ from gene symbols).
3. Benchmark before planning compute (3.6 s GPU vs 15 min LightGBM surprised the design).
4. Treat any label-derived feature as leaky until cross-fitted by protein cluster (68% co-labelled).
5. Audit "unbiased" resources for KS-label leakage before use (STRING 34% vs 4%).
6. Run lit/profiling/acquisition in parallel, but only fan out implementers after interfaces are frozen.
7. Cap design effort: first measured baseline before G1 beats further review rounds (no baseline at T+1.5h).
8. PI recomputes headline metrics from predictions; agent-reported numbers are untrusted.
9. One owner per path, PI-only git; update DELEGATION/RESEARCH_TREE status at each milestone.
10. Max 4 reviewers with distinct remits and a length cap; reviews end with a ranked top-5 changes.
