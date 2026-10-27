> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# M3 retro (T+10.8h, 00:10; FREEZE 5fc7870)

## Stop
1. **Post-hoc protocol amendments.** A23-A27 (b9a2a3f, efcf167), A29 (faf6372 20:45), A31/A32 (6b65e7e, 8107905 22:16), A33 (6b9440c 00:03): about one amendment per hour, each after seeing results. After FREEZE (5fc7870) allow none. If TEST disagrees with CV, report it; do not re-pick the model, the metric or the comparator.
2. **Reselecting the deployment model.** Retrain DEV+TEST with the frozen configs, seeds and fusion weights only. About 90 configs were tried (audit r2), and all negatives (R5 5a93287, knn_sim, F9, A32a/b 2609ef5, L2+ext) were within noise. Run no new variants and no "one more ablation".
3. **TEST review sprawl.** Report the 9 Holm primaries first and label everything else (secondary strata, T3, calibration) exploratory. Don't hunt subsets where we win. Weak spots (Tyr 0.63, M2 5-19 labels 0.142 vs 0.162, primed sites 0.893 vs 0.945, unseen drugs -0.04 to -0.06) must be in the report as written, not softened.
4. **Uncalibrated coverage claims.** Don't ship the assignment tier >=0.8 (overconfident, 0.87 stated vs 0.54 observed; Tyr 0/41) or A17 numbers without the A31 isotonic fix. U_obs predictions should be a ranked list plus calibrated probabilities, with the 66.6% SELPHI-imputed baseline pairs flagged.
5. **Stale bookkeeping (third retro in a row).** DELEGATION.md still lists R2-main, core-infra(G4), strata and regimes as "running" though all are done. The tree is 142 lines, and its status block repeats DESIGN and docs. Stale values keep slipping in (B1 0.719 vs 0.717, fixed only in cd70ce7). Update the status column per event or delete it, and keep the status block to ≤10 lines.
6. **Discovering PROMPT asks late.** The label-variant and regime comparisons were only delegated at 19:55 (3face5b) and closed at 00:03. Check the report outline against PROMPT.md now, before writing.

## Keep
1. **Independent audits that recompute from saved predictions** (audit r1 6941a74, r2 6b65e7e): everything reproduced, and they caught the A24 pooling artefact, stale stores, the overconfident tier, and a single-batch guard that could be bypassed. Add an audit of the TEST batch and the report numbers.
2. **Pre-declared rules and gates.** G2 met 16:27 against a 20:50 deadline. R5 failed its gate and was dropped. The A32 rules were written (8107905) before the experiments ran (2609ef5). The CV-to-TEST optimism was bounded (≲0.005).
3. **Sealed TEST, FREEZE.json (sha 721dc197), hash-pinned deliverables and a STARTED/COMPLETED ledger.** Retrained baselines were pinned (23:24).
4. **Fidelity checks and honest negatives.** Phosformer reproduces 0.888 vs 0.885 on its own split (65f1cbb), so the 0.679 on our benchmark is credible. SELPHI2 sits below Cantley on matched pairs, and the contaminated 0.793/0.796 is labelled reference-only.
5. **Paired comparisons with CIs against the stronger comparator** (ΔM1 +0.063 [0.048, 0.074]).
6. **File ownership, PI-only git, read-only data access for reviewers**: no collisions in 150 commits.

## Proposed LESSONS (13)
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
13. Audit PROMPT.md asks against the leaderboard at every milestone (regimes unrun at 20:00), and check `free -g` and cap ≤60 GB per agent on the shared node (OOM 17:50).
