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

### 2026-10-26 17:28 — M2: G1+G2 passed; first trusted CV results
Decision note: design frozen at v2 after 2 review rounds (round 3 skipped). Post-freeze fixes (DESIGN §17 A23–A27) are correctness or conservatism changes applied to all methods: within-fold rank normalisation for pooled CV (A24), the stronger Cantley comparator per endpoint (A23), canary v2 drawn from real labelled sites (the v1 sampler was biased), matched-kinase M2, and a hub-site stratum (reporting only). A label change for Tyr activation-loop pairs was rejected. Independent audit: all numbers reproduce, no leakage.
Status: on 5-fold DEV CV (A24 metrics; K_cv=203 kinases), R2 LightGBM using all feature families reaches M1 0.773 and M3 0.312. Comparators: R3 neural (Cantley-distilled) M1 0.754; Cantley as-published 0.717; simple-strong B7 0.726; SELPHI2.0 faithful retrain 0.708 on its 34%-coverage matched universe (below Cantley 0.728 on the same pairs); released Phosformer weights 0.796 (contaminated: trained on PSP). M2 (site-centric Hit@5) is weak for all methods (~0.15–0.22) and is the main open problem, together with Tyr kinases (motif ~0.53–0.60; the gap is partly from hub sites and weak atlas matrices). Running: R2 seeds/ablations/late fusion with R3/calibration, kNN-similarity feature, Phosformer faithful retrain (pretraining ETA 19:15), leaderboard re-evaluation. Next gate: G3 at 02:20.

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
  - [x] 5.1 Simple baselines → results/baselines/simple/summary.md (audit pending)
    - CV pooled OOF (K_cv=203): B1 Cantley published M1 0.717 / M2 0.157 / M3 0.133 (official B1; arm (ii) better on M2 → A23); B7 LR 0.726/0.196/0.177; B2 learned PSSM 0.660; floors B0-pop 0.552 (E-all AUROC 0.782!), B0-miss 0.538.
    - [!] Tyr kinases: motif baselines only ~0.53 E-hard AUROC (65 kinases) — Y specificity is the hard part.
    - [!] Pooled-OOF artefact: fold-varying offsets distort per-kinase AUROC (B0-prior → 0.40) → A24 within-fold rank normalisation.
  - [~] 5.2 R2 LightGBM main model + staged search + ablation curve
    - [!] Ablation (A24, pooled CV M1): F1 0.723 → +F0 0.733 → +F2a 0.738 → +F3 0.748 → +F4 0.753 → +F6a 0.755 → +F6bc (G2) 0.754 → +F5+F7 (full) 0.773. Alarms: +COMPARTMENTS-knowledge 0.776 (no big jump → no strong leak); F6-clean 0.773 (= full → no label-source-study circularity). Y kinases: R2 full 0.644 vs Cantley 0.534.
    - G2 MET 16:27: R2 (F1+F0+F2a+F3+F4+F6) M1 0.746 [0.728–0.760], M3 0.267; vs B1 ΔM1 +0.023 [0.009–0.034], ΔM3 +0.130. Full model (+F5+F7): M1 0.763, M2 0.208, M3 0.308, M4 0.782; canary v2 pass. (Pre-A24; audit pending.)
  - [x] 5.3 R3 neural → docs/model_r3.md (pending independent review)
    - CV pooled OOF: M1 0.729 [0.710–0.747] vs B1 Cantley 0.719; M2 Hit@5 0.176 vs 0.171; M3 0.247 vs 0.133; M4 0.728. LKO T2 (unseen kinase family): M1 0.664, M3 0.203 → descriptors carry zero-shot signal. 2.4 min/fold.
    - [!] Canary v1 sampler biased (sites on low-density proteins; label-free score gets 0.76) → canary v2 from real labelled sites, gate on E-hard (PI decision).
    - Positive canary 0.825 (<0.9): R3 lacks PPI input; planted motif alone would give 0.96.
    - [!] R3 round 2: Cantley-atlas distillation pretraining (+0.2 aux loss) → CV M1 0.740 [0.725–0.757] (+0.019 vs base, p=0.002), M2 0.223, M3 0.282; LKO T2 0.686 (vs 0.650); canary v2 gates pass (neg 0.501, pos 0.929). ESM-2 3B windows: no gain (not adopted). Final OOF: results/models/r3/FINAL_OOF.json.
  - [~] 5.4 B3 Phosformer reimplementation (GPU 0) → docs/baseline_phosformer.md
    - Released weights reproduce paper (AUROC 0.8855 on their split). On our CV (contaminated: trained on PSP incl. our positives) M1=0.793, M4=0.801, Hit@5=0.248 — reference only.
    - Faithful MLM pretraining (corpus filtered: 642 TEST substrates + 458 homologs + 152 kinase domains removed) running, ETA 19:15; fine-tune CV ETA 22:30.
  - [x] 5.5 B4 SELPHI2.0 retrain → docs/baseline_selphi2.md (audit pending)
    - [!] Faithful retrain (39/44 features regenerated, 2 substituted, label-derived via crossfit): matched-universe M1 0.708 vs Cantley 0.728 on same pairs (−0.020); M3 0.316 (+0.187, driven by site-popularity features: dropping 4 leak-suspect features → M1 0.724, M3 0.265). Coverage 34% of pairs (Ochoa-table sites only). 109 min on 12 cores.
    - [!] Their 0.97 AUROC comes from random-pair negatives, not split leakage: same metric ≈0.947 on our cluster split.
  - [x] 5.6 knn_sim crossfit family → docs/features_knn.md: fold-0 sanity BLOSUM log-ratio 0.750, ESM centre 0.744 vs PSSM 0.765 (r≈0.25 with PSSM; rank-avg with PSSM 0.789). Better for ≥100-label kinases. Stores for 5 CV folds; final TEST rows use inner-fold fit (A1, PI decision).
  - [x] 5.7 Tyrosine-kinase specificity → docs/features_tyr.md
    - [!] Y gap real but partly artefactual: macro dilution (42% of Y kinases ≤3 positives); E-hard Y pool dominated by promiscuous hub sites (GLO1 Y136 → 15 TKs, STAT5A Y694 → 9); Cantley Tyr atlas within-family discrimination ≈0.51 and some matrices (FGFR1, BTK, Eph, VEGFR, TAM) anti-correlated with own labels.
    - F9 SH2-docking (Yaron-Barir SH2 matrices, 24/28 SH2 kinases): modest; most GBDT gain is a NaN-indicator artefact (+0.019 [−0.018, 0.057] for kinases with SH2 matrix). Recommended small set → R2 (ablation).
    - Decision A26: hub-site stratum added (reporting only); relabelling A-loop cross-kinase labels rejected.
  - [x] 5.8 Independent audit → docs/audit/audit_cv_r1.md: all numbers reproduce exactly; no leakage. A24 values: R2 full M1 0.773 / M3 0.312; R3 final M1 0.754. M4 & 1–4-label claims NOT trusted pre-A24 (pooling inflation). R2 canary stores for folds 0/3 stale → rebuild. Matched-kinase M2 rule added (A27).
  - [~] 5.9 R5 fine-tuned ESM cross-encoder (GPUs 1–3; hard stop 00:30)
  - [~] 5.10 Downstream benchmarKIN harness
