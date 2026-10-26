# Research Tree — Human kinase–substrate prediction

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-26 — M0: project start
Project initialised: git repo, PROMPT.md (verbatim), DATA.md, LESSONS.md, env (.venv, requirements.txt). Provided data: PhosphoSitePlus Kinase_Substrate_Dataset (23,904 rows; 14,724 human kinase–human substrate rows; 421 human kinases; 10,224 human sites; highly skewed: top kinase CDK1 694 rows, 146 kinases with <5 rows), a perturbation phosphoproteomics compendium (1,937 conditions × 189,962 phosphopeptides, long-format log2FC), Hijazi et al. inhibitor→kinase target matrix, and two papers (Phosformer; Petsalaki-lab integrative method).
Hardware: 64 CPU cores, 375 GB RAM, 4× RTX 4090 (24 GB). Internet available. Next: parallel literature review (methods/benchmarks, KS annotation sources, phosphoproteomics sources, other modalities) + data profiling, then DESIGN.md and up to 3 review rounds.

## Time plan (24 h budget: 2026-10-26 13:20 → 2026-10-27 13:20)
- T+0–3h: literature, data profiling, data acquisition (parallel); DESIGN.md + review rounds
- T+3–6h: labels, splits, feature pipelines, baselines (Cantley PSSM, Phosformer, Petsalaki method) retrained on our splits
- T+6–18h: model iterations (training regime, modalities, unorthodox kinases), ablations, verification
- T+18–22h: freeze; final test-set evaluation once; independent review of results
- T+22–24h: REPORT.md, final commit
User notes: PSP full downloads will be provided by user; CPTAC/others: try to fetch ourselves, ask user only if slow.

## Tree
- [~] 0. Setup
  - [x] 0.1 git, PROMPT.md, .gitignore (data/ never committed)
  - [x] 0.2 Python env (.venv from requirements.txt; torch 2.6 cu124, 4 GPUs visible)
- [~] 1. Literature (LITERATURE.md) — 4 parallel agents -> docs/lit/0{1..4}_*.md
  - [ ] 1.1 Prior methods & benchmarks (Cantley atlas, Phosformer, Petsalaki MCP 2025, others)
  - [ ] 1.2 Kinase–substrate annotation sources & quality
  - [ ] 1.3 Phosphoproteomics sources (compendium, CPTAC, inhibitor/decryptM, kinobeads)
  - [ ] 1.4 Other modalities (PPI, localisation, co-expression, structure, PLMs)
- [~] 2. Data profiling of provided data (agent -> docs/data_profile.md, src/ksp/io_psp.py)
- [ ] 3. DESIGN.md + review rounds (≤3)
