# DELEGATION — who owns which files (DESIGN §11.8)
Agents write only in their assigned paths; shared modules are changed only by their owner. Only the PI runs git.

| Agent (task) | Owns (write) | GPU / cores | Status |
|---|---|---|---|
| PI | DESIGN.md, RESEARCH_TREE.md, LESSONS.md, DATAINEED.md, DATA.md (shared append), src/ksp/index.py, docs/DELEGATION.md | — | — |
| labels/universe | src/ksp/{ids,kinases,sites,labels}.py, scripts/build_{kinase_list,site_universe,labels}.py | 0 / 16 | done (handed to core-infra) |
| struct+ESM | src/ksp/features/structure*.py, src/ksp/features/esm*.py, cache/structure, cache/esm2_650m | GPU 2,3 / 16 | running |
| core-infra | src/ksp/{splits,eval,canary,manifest,eval_test}.py, src/ksp/labels.py (load_dev, KP, masks, residue-class overrides), config/kinase_curation.tsv, config/bio_controls.tsv, sealed/, results/splits, tests/ | 0 / 12 | running (14:40) |
| feat-motif | src/ksp/features/motif*.py, cache/features/motif, results/features/motif | 0 / 12 | running (14:40) |
| feat-context | src/ksp/features/{ppi,loc,coexpr}*.py, cache/features/context, results/features/context | 0 / 10 | running (14:40) |
| feat-phospho | src/ksp/features/phospho*.py, cache/features/phospho, results/features/phospho | GPU 1 / 10 | running (14:40) |
| B3-phosformer | src/ksp/baselines/phosformer/, cache/phosformer, results/baselines/phosformer | GPU 0 / 6 | running (14:40) |
