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
  - [~] 3.1 Review round 1: external PI, statistics, data/leakage, algorithms/compute, AI failure modes, kinase biology, prior work → docs/reviews/r1_*.md
- [~] 4. Foundations (design-independent, parallel with review)
  - [~] 4.1 ids.py, kinase list, site universe, labels L1–L3 + external curated (agent)
  - [~] 4.2 AlphaFold structure features + ESM-2 650M per-residue embeddings for all S/T/Y (agent; GPUs 2–3)
  - [~] 4.3 Data acquisition wave 2 (Buljan2020, Kinobeads/Davis/DTC, Sugiyama, Phosformer weights, SELPHI2 code, KinCore, HGNC, SIGNOR retry)
