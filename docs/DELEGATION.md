# DELEGATION — who owns which files (DESIGN §11.8)
Agents write only in their assigned paths; shared modules are changed only by their owner. Only the PI runs git.

| Agent (task) | Owns (write) | GPU / cores | Status |
|---|---|---|---|
| PI | DESIGN.md, RESEARCH_TREE.md, LESSONS.md, DATAINEED.md, DATA.md (shared append), src/ksp/index.py, docs/DELEGATION.md | — | — |
| labels/universe | src/ksp/{ids,kinases,sites,labels}.py, scripts/build_{kinase_list,site_universe,labels}.py | 0 / 16 | done (handed to core-infra) |
| struct+ESM | src/ksp/features/structure*.py, src/ksp/features/esm*.py, cache/structure, cache/esm2_650m | GPU 2,3 / 16 | done 14:40 |
| core-infra | src/ksp/{splits,eval,canary,manifest,eval_test}.py, src/ksp/labels.py (load_dev, KP, masks, residue-class overrides), config/kinase_curation.tsv, config/bio_controls.tsv, sealed/, results/splits, tests/ | 0 / 12 | done 15:05 |
| feat-motif | src/ksp/features/motif*.py, cache/features/motif, results/features/motif | 0 / 12 | done |
| feat-context | src/ksp/features/{ppi,loc,coexpr}*.py, cache/features/context, results/features/context | 0 / 10 | done |
| feat-phospho | src/ksp/features/phospho*.py, cache/features/phospho, results/features/phospho | GPU 1 / 10 | done |
| R3-neural | src/ksp/models/r3/, cache/r3, results/models/r3, scripts/models/r3_*.py | GPU 2,3 / 8 | done |
| B4-selphi2 | src/ksp/baselines/selphi2/, cache/selphi2, results/baselines/selphi2, scripts/baselines/selphi2_*.py | 0 / 12 | done |
| B3-phosformer | src/ksp/baselines/phosformer/, cache/phosformer, results/baselines/phosformer | GPU 0–3 / 6 | done 23:30 |
| simple-baselines | src/ksp/baselines/simple.py, results/baselines/simple | 0/8 | done 17:00 |
| R2-main | src/ksp/models/r2/, cache/r2, results/models/r2 | 0/24 | done 23:03 (TEST/T3 delivered) |
| auditor-1 | docs/audit/ | 0/8 | done 18:10 |
| feat-knn | src/ksp/features/knn*.py, cache/features/knn | GPU 2 / 6 | done |
| feat-tyr | src/ksp/features/tyr*.py, cache/features/tyr | 0/6 | done |
| R5-crossenc | src/ksp/models/r5/, scripts/models/r5_*.py, results/models/r5 | GPU 1–3 | done 19:53 (negative) |
| core-infra (G4) | src/ksp/eval_test.py, tests/, docs/G4_RUNBOOK.md | 0/8 | done 22:00; deploy loader running (00:12) |
| strata | src/ksp/strata.py, tests/test_strata.py, results/strata, docs/strata_cv.md | 0/8 | done 23:05 |
| regimes | scripts/models/regime_*.py, src/ksp/models/regime/, results/models/regime, docs/regimes_labels.md | 0/24 | done 00:03 |
