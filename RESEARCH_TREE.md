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

### 2026-10-26 ~14:15 — M1: design v1 + implementation start
Decision note: accepted ~80 of 84 round-1 changes (map in DESIGN §15). Biggest design changes: (1) hard gates G1–G6 with scope cuts; (2) leave-substrate-cluster-out cross-fitting for every label-derived feature (DL-1: 68% of KS pairs have a same-protein co-labelled site → leave-site-out leaks); (3) paralog-masked hard negatives (SRC family: >50% of hard negatives were paralog substrates); (4) endpoints M1 (macro E-hard AUROC), M2 (filtered site-centric Hit@5), M3 (top-1% recall on U_ms), M4/M5 for understudied kinases; (5) sealed TEST + FREEZE.json + canary kinases; (6) only systematic-screen PPI (STRING leaks).
Status: literature (4 reviews), data profiling, two acquisition waves (~12.6 GB), labels (L2 = 18,638 pairs, 401 kinases) and site universe (350,736 sites) are done. Round-2 review (4 reviewers) running in parallel with 6 implementation agents: core infra (splits/seal/eval/canaries), motif features, context features, phospho features, Phosformer reimplementation (GPU 0), structure+ESM (GPUs 2–3). DATAINEED.md created for the user (7 open requests; none blocking).

## Tree
- [~] 0. Setup
  - [x] 0.1 git, PROMPT.md, .gitignore (data/ never committed)
  - [x] 0.2 Python env (.venv from requirements.txt; torch 2.6 cu124, 4 GPUs visible)
- [~] 1. Literature (LITERATURE.md) — 4 parallel agents -> docs/lit/0{1..4}_*.md
  - [x] 1.1 Prior methods & benchmarks -> docs/lit/01_methods_benchmarks.md
    - [!] Phosformer: random pair split; 51% of test-positive peptides are training positives under another kinase → reported 0.885 AUROC is inflated. Code CC BY-NC-ND, no train script → reimplement architecture for retraining.
    - [!] SELPHI2.0 (Petsalaki): RF, 45 features (~5 kinase-specific); CV over pairs; coreg features leak (KSEA activities from same PSP substrates). Apache-2.0. Retrain with fold-aware features.
    - [!] Cantley PSSMs (303 S/T + 78 Y) via kinase-library package; percentiles vs Ochoa reference, not PSP-trained; priming supported (pS/pT/pY rows).
    - [!] KolossuS (ESM-2 15B, bioRxiv 2025): PSP AUROC 0.755 S/T / 0.838 Y on phylogeny-aware split — realistic bar for sequence-only. LinkPhinder 0.97 random→0.53 independent.
    - Recommended: protein-cluster splits w/ orthologs merged; tiers unseen-protein / unseen-kinase-cluster / both; macro per-kinase AUPRC + per-site Hit@k.
  - [x] 1.2 KS annotation sources -> docs/lit/02_ks_annotation_sources.md
    - [!] PSP is core. 40% in-vitro-only / 25% in-vivo-only / 35% both. In-vivo-only enriched for indirect (TK→S/T: 69% in-vivo-only). Map sites by 15-mer to canonical, never by number (893 isoform rows).
    - [!] Extra curated labels: SIGNOR (unreachable directly: connection fails; use OmniPath copy), UniProt ECO:0000269 (~960 new pairs), OmniPath curated non-PSP rows (~5.8k). Sugiyama 2019 in vitro → features only (Müller-Dott 2025: hurts when added as labels).
    - Decision: treat as PU learning; negatives from observed phosphosites matched on residue; pre-register easy and hard negative sets.
  - [x] 1.3 Phosphoproteomics sources -> docs/lit/03_phosphoproteomics_sources.md
    - [!] Compendium = 1,937 contrasts/178 PMIDs/189,962 sites, 4.6% dense; contains decryptM (790 contrasts, 41% of columns) + Hijazi (183). Study-structured → need study-aware aggregation.
    - [!] CPTAC 10 cohorts public on LinkedOmics S3 (no login); benchmarKIN Zenodo 12566560.
    - [!] Leakage: CoPheeMap trained on PSP; SELPHI2.0 splits on pairs and uses Hijazi as both feature and test; Ochoa functional score trained on PSP → ablate.
  - [x] 1.4 Other modalities -> docs/lit/04_other_modalities.md; leakage audit results/leakage/
    - [!] STRING 'experiments/physical' is NOT leak-free: 34% of PSP KS pairs have a score vs ~4% controls (imports BioGRID Biochemical Activity = kinase assays). Use only unbiased screens (BioPlex, HuRI, OpenCell, Buljan2020 kinase AP-MS, HumanCellMap) for PPI; ~7-8x enriched but 0.3–1.5% coverage.
    - Localisation: HPA IF, OpenCell, COMPARTMENTS knowledge+experiments (not textmining). Structure: AFDB pLDDT/RSA via structuremap. PLMs: ESM-2 650M/3B (MIT) downloadable; ESM-C via plm_env transformers 5.17.
    - Structure-based KS (AF-multimer) infeasible at scale; optional pilot only.
- [~] 2. Data profiling of provided data (agent -> docs/data_profile.md, src/ksp/io_psp.py)
- [~] 3. DESIGN.md + review rounds (≤3)
  - [x] 3.0 LITERATURE.md synthesis; DESIGN.md v0
  - [x] 3.1 Review round 1 (7 reviewers) → docs/reviews/r1_*.md; DESIGN v1 (§15 change map)
  - [x] 3.2 Review round 2 (4 reviewers) → docs/reviews/r2_*.md; DESIGN v2 FROZEN (§16 A1–A22). Round 3 skipped (no design-level critical left).
    - [!] Coverage metric as planned was vacuous (random scores put 98.5% of S/T sites in some kinase's top 1%) → redefined as site-level assignment probability vs random baseline (A17).
    - [!] iKiP-DB 'combined' GMT contains 61% of PSP pairs → forbidden (A11). Sugiyama S2 clean.
    - [!] 18% of human-native positives lack MS evidence vs 0.4% of unlabelled → U_ms-only candidates (A6).
    - [!] Cache row order ≠ index order for 94% of sites → ksp.store with index-hash checks (A2).
- [x] 3.3 Klaeger 2017 Kinobeads tables provided by user (data/Klaeger) → F6b primary drug→target source
- [x] 4. Foundations (design-independent, parallel with review)
  - [x] 4.1 ids.py, kinase list (503 incl.), site universe (350,736), labels (L2 18,638 pairs/401 kinases), external curated (3,788 non-PSP pairs) → docs/labels_universe.md
  - [x] 4.1b canonical feature index cache/index (src/ksp/index.py)
  - [x] 4.4 core infra → docs/core_infra.md, results/splits/split_report.md: eval (21M pairs/fold in 27 s), crossfit (toy leak 1.0 → 0.5), 13 canaries, FREEZE-gated eval_test, 63 tests pass. TEST external 215 pairs (158 tier A); temporal 263 positives.
  - [x] 4.5 F1/F2a motif features → docs/features_motif.md
    - [!] Cantley engine reproduces kinase-library 1.8.0 (100% within 1e-3); 396/503 kinases have own matrix. Sanity (fold 0, random negs, 94 kinases): Cantley 0.785, learned PSSM 0.766 (+phospho channel 0.782) — complementary.
    - [!] Priming-aware scoring helps GSK3B 0.67→0.79, GSK3A 0.77→0.87, CK1 family; hurts AKT1/PKC/CHEK1 slightly → keep as separate feature.
  - [x] 4.6 F3/F4/F5 context features → docs/features_context.md (43 matrices)
    - PPI screens graph 16,178 proteins/279,677 edges; 3.6% of KS pairs have a direct edge (4.5× vs shuffled). RWR AUROC 0.67 vs random but 0.54 vs shuffled → mostly 'studied-ness'.
    - COMPARTMENTS knowledge strongest (0.74) but imports GO/UniProt curation → ablation only. FAM20C sites 87.5% lumenal.
  - [x] 4.7 F6 phospho features → docs/features_phospho.md
    - Sanity macro AUROC: F6a coreg compendium 0.63, CPTAC 0.61; F6b inhibitor 0.61 (384 kinases); F6c 0.53–0.56. Labelled sites 2–2.5× more often covered → B0-miss floor mandatory.
    - Note: ridge strength picked on DEV sanity AUROC (5 values) — mild DEV optimism, logged. Some drug potencies hand-curated (documented in results/features/phospho).
  - [~] 4.8 B3 Phosformer reimplementation + MLM pretraining (GPU 0)
  - [x] 4.2 Structure (AFDB v6, 99.45% S/T/Y) + ESM-2 650M & 3B embeddings (2.8 / 7.9 GPU-min) → docs/features_structure_esm.md
  - [x] 4.4a G1 split freeze (11,609 clusters; TEST 20.0%; K_eval 119, K_any 256, K_cv 203; 13 canaries) → results/splits/
  - [x] 4.3 Data acquisition wave 2 (Buljan2020, Kinobeads/Davis/DTC, Sugiyama, Phosformer weights, SELPHI2 code, KinCore, HGNC, SIGNOR retry)
- [~] 5. Models & baselines on CV (G2 20:50)
  - [~] 5.1 Simple baselines B0-prior/pop/miss, B1 Cantley arms, B2 PSSM, B7
  - [~] 5.2 R2 LightGBM main model + staged search + ablation curve
  - [~] 5.3 R3 neural (GPUs 2–3)
  - [~] 5.4 B3 Phosformer reimplementation (GPU 0)
  - [~] 5.5 B4 SELPHI2.0 retrain
