# LESSONS — process rules for running this project efficiently
(One-line rule + evidence. Pruned at each milestone; last prune M1 2026-10-26 14:55, see docs/retro/M1_retro.md.)

1. Give reviewers read-only data access and require a measured number per change — r1: 6/7 reviewers found quantified leaks (68% same-protein co-labelled KS pairs; >50% SRC-family hard negatives are paralog substrates) that reading alone missed.
2. Key everything by UniProt accession, never symbol — 174/404 Cantley names ≠ gene symbols; Cantley "PDK1" = PDPK1, gene PDK1 = PDHK1.
3. Benchmark compute before planning around it — design guessed co-regulation was the bottleneck (actual 3.6 s GPU); LightGBM 400k×160 took 15 min.
4. Treat every label-derived feature as leaky until cross-fitted by substrate cluster — 68% co-labelled same-protein sites make leave-site-out leak.
5. Audit "unbiased" resources for KS-label leakage before use — STRING experiments channel covers 34% of PSP KS pairs vs ~4% controls (imports kinase assays).
6. Parallelise literature/profiling/acquisition from the start; fan out implementers only against frozen interfaces (index, splits, eval API) — avoids rework when design changes.
7. Get a first measured baseline before more review rounds — at T+1.5 h no baseline existed while 84 review changes were processed (M1 retro).
8. PI recomputes headline metrics from saved predictions with ksp.eval; agent-reported numbers are untrusted.
9. One owner per path, PI-only git (NFS), and update tree/delegation status at each notification, not in batches (M1 retro: stale statuses).
10. ≤4 reviewers per round with distinct remits, length cap, ranked top-5 changes — r1 produced ~26k words.
11. Take timestamps from `date`, never estimate — PI notes drifted ~40 min ahead of the real clock by T+1h (DESIGN/TREE timestamps corrected 14:35).
12. Verify file identities before writing them into records — LITERATURE.md had paper1/paper2 swapped until the Phosformer agent read the PDFs (fixed 15:22).
