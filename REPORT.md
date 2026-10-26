# REPORT: Predicting direct kinase–substrate sites across the human proteome

**Status: DRAFT for PI review (2026-10-26, about 23:30 CET).** The single TEST batch has not run yet. Every TEST number is a marked placeholder, `[[TEST: …]]`. Each CV number cites the committed file it comes from.

**Conventions.** All CV metrics use the frozen harness `ksp.eval`, commit 974d520, with A24 within-(fold, kinase) rank normalisation, unless a row says otherwise. Where a document gives both a pre-A24 and an A24 value, this report uses the A24/audited value and says so. Square brackets are 95% cluster-bootstrap CIs.

---

## 1. Summary

**Task.** For each (human protein kinase, human phosphosite S/T/Y) pair, score the probability that the kinase phosphorylates the site directly. Labels are experimentally supported PhosphoSitePlus (PSP) pairs. The problem is treated as positive–unlabelled (PU) learning.

**Evaluation design.**
- Substrate-protein clusters are held out (homology components): 11,609 clusters, 20.3% sealed as TEST.
- Every label-derived feature is cross-fitted by substrate cluster.
- Hard negatives are paralog-masked: sites labelled to *other*, non-paralogous kinases.
- Kinase-equal macro endpoints, with understudied, Tyr, priming and atypical strata.
- Pre-registered superiority tests against retrained baselines; a single sealed TEST batch.

**Best method.** Late fusion (z-scored logits) of two models:
- **R2**, a global LightGBM pair model over motif, structure/ESM, PPI screens, localisation, co-expression, phosphoproteomics and kinase-descriptor features (5-seed ensemble);
- **R3**, a neural kinase × site model with Cantley-atlas distillation pretraining.

**Headline CV results** (pooled 5-fold DEV OOF, K_cv = 203 kinases; `docs/model_r2.md` §6, audited in `docs/audit/audit_cv_r2.md`):

| Endpoint | Final (R2+R3 fusion) | R2 alone | R3 alone | Cantley (B1) | Paired Δ final − B1 (matched universe) |
|---|---|---|---|---|---|
| M1: macro E-hard AUROC | **0.786** (folds 0.767–0.794) | 0.776 | 0.754 | 0.717 | +0.063 [0.048, 0.074] |
| M2: site-centric Hit@5 | **0.263** (nested 0.259) | 0.200 | 0.223 | 0.171 (arm ii) | +0.105 [0.082, 0.125] |
| M3: top-1% recall on E-all | **0.361** | 0.312 | 0.294 | 0.133 | +0.227 [0.189, 0.251] |
| M4-CV (339 kinases) | 0.774 | 0.767 | 0.731 | 0.690 | +0.071 [0.053, 0.080] |

**TEST (placeholders).**

| | Final | Cantley (as published) | B3 Phosformer (retrained) | B4 SELPHI2.0 (retrained) |
|---|---|---|---|---|
| M1 | [[TEST: final M1 + CI]] | [[TEST: B1(i) M1]] | [[TEST: B3 M1]] | [[TEST: B4 M1]] |
| M2 | [[TEST: final M2 + CI]] | [[TEST: B1(ii) M2]] | [[TEST: B3 M2]] | [[TEST: B4 M2]] |
| M3 | [[TEST: final M3 + CI]] | [[TEST: B1(i) M3]] | [[TEST: B3 M3]] | [[TEST: B4 M3]] |

**Runtime of the final model.**
- Training on all of DEV: 1,137 s (R2 418 s + R3 719 s).
- Scoring all 25.9 M TEST candidate pairs: 317 s.
- Source: `results/freeze_entries/final.json`.

**The main caveats, stated up front:**
1. **Tyrosine kinases remain weak.** CV M1 is 0.63 for Y kinases against 0.81 for S/T kinases.
2. **No reliable high-confidence site assignment.** After the A31 recalibration, the model has no ≥0.8 assignment tier and no Y site reaches 0.5.
3. **M2 is weak for understudied kinases.** Hit@5 is 0.14 for kinases with 5–19 labels, below Cantley's 0.16.
4. **Downstream, the library is at parity with PSP, not better.** It trails PSP by 0.04–0.06 on drugs it has not seen.
5. **Both retrained baselines fall below label-free Cantley on M1.** Phosformer (B3) and SELPHI2.0 (B4) do.
6. **About 90 configurations were tried on CV.** The audited CV→TEST optimism for M1 is ≲0.005.

---

## 2. Problem and evaluation design

### 2.1 Labels and PU framing (`DESIGN.md` §2–3, §5; `docs/labels_universe.md`)

| Item | Value | Source |
|---|---|---|
| Kinases included | 503 (ST 391 / Y 92 / STY 20; 38 pseudokinases flagged) | `docs/labels_universe.md` |
| Training labels L2 | human-native PSP pairs plus mouse/rat pairs mapped to human by 15-mer: 18,638 pairs, 401 kinases | `RESEARCH_TREE.md` 4.1 |
| After U_ms restriction (A6) | L2 used 14,747; evaluation pairs 10,961 (956 pairs with no index row and 2,934 non-U_ms pairs dropped) | `results/splits/split_report.md` |
| Site universe | U_obs 350,736 sites (PSP ∪ compendium ∪ CPTAC ∪ Ochoa); U_ms = sites with ≥1 MS observation | `DESIGN.md` §2 |
| Test positives | human-native PSP only, in residue class, not autophosphorylation; 2,205 on TEST | `split_report.md` |
| Evidence mix (PSP) | 40% in vitro only / 25% in vivo only / 35% both | `LITERATURE.md` §2 |

**PU framing.** There are no true negatives. For kinase k, the negatives are U_ms sites of k's residue class, minus three sets:
- every known positive of k from any source or species, plus any site with an identical 15-mer;
- the paralog mask Par(k): ≥60% kinase-domain identity or the same Manning subfamily;
- 15-mers within Hamming distance 2 of a known positive.

Two negative sets are pre-registered:
- **E-hard:** sites labelled to some kinase outside k ∪ Par(k). This is the primary set.
- **E-all:** every U_ms site, scored exhaustively.

**Why Par(k) is needed.** Without it, more than 50% of SRC-family hard negatives were paralog substrates (round-1 review).

### 2.2 Splits and leakage controls (`DESIGN.md` §4, §7, A1–A6; `results/splits/split_report.md`)

**Substrate clusters.** mmseqs2 all-vs-all search (≥40% identity, ≥50% coverage), with connected components unioned with identical-15-mer edges. That gives 11,609 clusters; the largest holds 2.35% of labels.

**TEST.** 2,353 clusters (20.3%), 3,957 proteins and 69,430 U_ms sites, chosen by simulated annealing (seed 0) to give per-kinase test fractions of 15–25%. TEST labels are sealed and their SHA-256 is recorded.

**DEV CV.** 5-fold, grouped by cluster, under the same constraints.

**Leave-kinase-out (LKO).** 120 kinase clusters (≥50% KinCore domain identity ∪ Manning family), split into 5 folds. The tiers are:
- T1: seen kinase, unseen protein;
- T2: unseen family, seen proteins (CV);
- T3: unseen family and unseen protein (TEST).

**Cross-fitting (DL-1).** 68% of KS pairs have a co-labelled site on the same protein, so leave-site-out evaluation leaks. Every label-derived feature therefore goes through one driver, `ksp.crossfit`, with inner 5-fold cross-fitting by substrate cluster. The driver raises `LeakageError` on any protein or cluster overlap. Eval rows use a randomly assigned inner-fold fit (A1).

**Canaries.** 12 negative canary kinases (v2: real DEV-labelled sites, clustered by protein) and 1 positive canary with a planted motif and PPI edge. The gate is negative ≤0.57 on E-hard (pooled) and positive ≥0.90.

### 2.3 Evaluable kinase sets (frozen at G1)

| Set | Definition | n |
|---|---|---|
| K_eval | ≥3 TEST positives on ≥2 clusters | 119 |
| K_any | ≥1 TEST positive | 256 |
| K_cv | ≥5 DEV positives on ≥2 clusters | 203 |

Source: `results/splits/split_report.md`. DESIGN §13 anticipated about 160 for K_eval and A19 about 135; the frozen value is 119.

### 2.4 Endpoints (`DESIGN.md` §6, A17–A20, A24, A27, A31)

| Endpoint | Definition | Score used |
|---|---|---|
| **M1** (primary) | macro over K_eval of per-kinase AUROC, E-hard, autophosphorylation excluded | `score` |
| **M2** (co-primary) | site-centric filtered Hit@5: rank of each true kinase among kinases of matching class, other true kinases and unlabelled Par(k) removed; macro within kinase, then over K_eval | `score_site`, the raw fused score (A31a) |
| **M3** (co-primary) | macro recall of TEST positives in the kinase's top 1% of E-all | `score` |
| **M4** | M1 over K_any; a single positive counts as its percentile | `score` |
| **M5** | CV M1 over K_cv, stratified by training-label bin | `score` |
| **A17 coverage** | sites whose top kinase has assignment probability ≥τ; reliability checked on labelled sites | `assignment_prob` |

**Strata:** Manning group; residue class; evidence class; priming kinases and primed sites; PIKK; FAM20C; training-label bins; hub sites (A26, reporting only); new since 2021; external curated.

**Floors:**
- B0-prior: the kinase's training count;
- B0-pop: site popularity (diagnostic only);
- B0-miss: LightGBM on feature missingness only.

### 2.5 Statistics, sealing and the single-batch rule

**Bootstrap.** Cluster bootstrap over substrate clusters; B = 1,000 on CV, B = 10,000 for TEST primaries. A kinase-resampling sensitivity analysis is run in the batch.

**Primary family: exactly 9 tests, Holm-corrected.** "Beats" means Holm p < 0.05 and the paired CI lower bound > 0. The tests are:

| Comparator | Endpoints |
|---|---|
| B1 arm (i), Cantley as published | M1, M3 |
| B1 arm (ii), universe-recalibrated | M2 |
| B3, Phosformer retrained | M1, M2, M3 |
| B4, SELPHI2.0 retrained | M1, M2, M3 |

The comparator rules follow A19 and A23, which uses the stronger Cantley arm for each endpoint.

**Matched-universe rule.** Head-to-heads use the pairs both methods score (PRIOR-1); M2 competitors are restricted to kinases both methods score (A27).

**Sealing (A3, A28, A29).**
- Only `src/ksp/eval_test.py` reads `sealed/`.
- FREEZE.json pins every prediction file by SHA-256.
- A STARTED/COMPLETED ledger is written around the read; `--real` and B = 10,000 are required.
- The harness is frozen at 974d520.
- Any COMPLETED block refuses further batches.

**Current state.**
- Every TEST deliverable was pinned by 23:24 (`results/freeze_entries/`, `scripts/build_freeze.py`).
- The draft FREEZE validates with all 9 primaries.
- FREEZE.json is not yet committed.

**Selection discipline.**
- Exploration on CV closed at 01:00 (A28).
- Selection rule: best M1-CV, with M2-CV as tie-break.
- About 90 configurations were evaluated on CV, about 55 of them selectable.

---

## 3. Data and features

### 3.1 Modalities (`docs/features_*.md`; `docs/feature_provenance.md`)

| Family | Content | Label-derived? | Coverage | Doc |
|---|---|---|---|---|
| F0 site | AlphaFold v6 pLDDT, exposure, IDR, RSA, SS; ESM-2 650M per-residue embedding (PCA-64) | no | 99.45% of S/T/Y | `features_structure_esm.md` |
| F1 Cantley | atlas score; published and universe percentiles; primed hard/soft; nearest-atlas imputation with a distance feature (396/503 kinases have their own matrix) | no | 396 own + imputed | `features_motif.md` |
| F2a learned PSSM | ±7 log-odds with an empirical-Bayes family/group prior | **yes**, cross-fitted | 503 | `features_motif.md` |
| F3 PPI | systematic screens only (BioPlex 3.0, HuRI, OpenCell, Buljan 2020, HCM): 279,677 edges; direct edges, RWR; F3-label neighbourhood | F3-label **yes** | edges on 3.6% of KS pairs | `features_context.md` |
| F4 localisation | HPA IF v25, OpenCell, HCM, Hein 2025 cosines; UniProt topology (FAM20C sites 87.5% lumenal) | no | 92% of KS pairs | `features_context.md` |
| F5 co-expression | HPA tissue RNA, DepMap 26Q1 | no | ≥94% | `features_context.md` |
| F6a co-regulation | compendium study blocks and CPTAC (protein-residualised), top-3/median Fisher z with training substrates | **yes**, cross-fitted | about 300 kinases × 60k sites per fold | `features_phospho.md` |
| F6b inhibitor | ridge deconvolution of 692 inhibitor conditions on a drug × kinase occupancy matrix (Klaeger, Reinecke, CurveCurator, Davis, Hijazi) | no | 384 kinases × 75k sites | `features_phospho.md` |
| F6c CPTAC kinase state | correlation with kinase abundance and A-loop activity, 10 cohorts | no | 427 / 99 kinases | `features_phospho.md` |
| F7 kinase descriptors | Manning group/family, ESM-2 kinase domain, Cantley score-profile PCA | no | 503 | `model_r2.md` |
| F8 Sugiyama 2019 | planned | — | **never built**: the R2 group is empty (`F8_unavailable` fallback) | `model_r2.md` §2 |
| F9 Tyr / SH2 docking | Yaron-Barir SH2 matrices on other pY sites, ITAM, RTK × membrane | no | 24/28 SH2 kinases | `features_tyr.md` (not adopted) |
| knn_sim | similarity to training substrates (ESM, BLOSUM) | **yes**, cross-fitted | U_ms | `features_knn.md` (not adopted) |

### 3.2 Leakage exclusions and the evidence for each

| Excluded or restricted | Evidence | Source |
|---|---|---|
| STRING, every channel incl. "experiments" | 34% of PSP KS pairs have a score vs ≈4% of controls; it imports BioGRID "Biochemical Activity" (kinase assays) | `LITERATURE.md` §4, `results/leakage/` |
| iKiP-DB "combined" GMT and activity tables | contain 61% of PSP pairs | `RESEARCH_TREE.md` 3.2; A11 |
| COMPARTMENTS knowledge channel (ablation only) | imports GO/UniProt curation; strongest single feature (0.74 vs random) but only +0.003 M1 as an alarm (0.776 vs 0.773) | `features_context.md`; `model_r2.md` §5 |
| Ochoa functional score, PSP Regulatory_sites, LT/MS counts, third-party predictors, CoPheeMap, KinaseCanvas | trained on, or derived from, PSP. In SELPHI2.0, dropping 4 such features raised M1 by +0.016 and lowered M3 by −0.051 (a popularity effect) | `DESIGN.md` §7; `baseline_selphi2.md` §0 |
| Centre-site observation flags and counts (A6) | 18% of human-native positives lack MS evidence vs 0.4% of unlabelled sites; positives and negatives are both restricted to U_ms | `RESEARCH_TREE.md` 3.2 |
| HPA reliability, `loc_n_sources`, Buljan publication counts (A11) | study-bias proxies; `loc_n_sources` alone reaches AUROC 0.656 | `features_context.md` §6 |
| Sugiyama 2019 as labels | in vitro lysate screen; reported to hurt activity inference when used as labels (Müller-Dott 2025) | `LITERATURE.md` §2 |
| Uncross-fitted label-derived features | knn_sim leaky fit 0.93–1.00 vs 0.59–0.76 cross-fitted; F3-label in-sample vs cross-fitted 0.583 vs 0.552 (vs SHUF) | `features_knn.md`; `features_context.md` |
| Phosformer pretraining corpus | 642 TEST substrates, 458 TEST homologs and 152 near-identical kinase domains removed | `baseline_phosformer.md` |
| Label-source study circularity (A12) | F6-clean, with 18 label-source PMIDs removed: M1 0.773 = full model | `model_r2.md` §5 |
| Missingness as signal | B0-miss inside R2: M1 0.522 | `model_r2.md` §5 |

**Rejected after seeing results.** A26(b) proposed treating cross-kinase activation-loop Tyr labels as mask-only. It was rejected because it would change labels after looking.

### 3.3 Provenance

- Every feature family has a row in `docs/feature_provenance.md`: source, whether label-derived, leakage risk, coverage.
- Every dataset is listed with source, version, SHA-256 and licence in `DATA.md`.
- Feature matrices load only through `ksp.store`, which carries an index hash and raises on mismatch. 94% of sites had cache row order different from index order (A2).

---

## 4. Models and training regimes

### 4.1 R2: global LightGBM pair model (`docs/model_r2.md`)

**Final configuration v3.**
- Binary objective; 127 leaves; lr 0.03; inner early stopping on real kinases.
- Sample weights n_k^−0.5.
- 1:50 negatives with 50% hard (A15 fill).
- Feature groups F1+F0+F2a+F3+F4+F6abc+F5+F7.
- 5-seed ensemble by mean logit.

**Staged search on fold 0** (19 configurations; noise floor about ±0.005):
- tree size;
- α ∈ {0, 0.5, 1};
- negative mix (1:20, 1:50, 1:200 × 25%, 50%, 75% hard);
- objectives.

| Objective (fold 0) | M1 | M2 |
|---|---|---|
| binary (chosen) | 0.754 | 0.223 |
| lambdarank by kinase | 0.754 | 0.086 |
| site-grouped lambdarank | 0.730 | 0.189 |

Source: `model_r2.md` §4. The ranking objectives did not beat the binary loss.

| | M1 | per-fold mean | M2 | M3 | M4 |
|---|---|---|---|---|---|
| seeds 0–4, mean ± sd | 0.7711 ± 0.0028 | 0.7669 | 0.187 | 0.300 | 0.761 |
| 5-seed ensemble | **0.776** | 0.771 | 0.200 | 0.312 | 0.767 |

Source: `model_r2.md` §6; audit r2 reproduced 0.7760.

### 4.2 R3: neural model with Cantley distillation (`docs/model_r3.md`; A24 values from `docs/audit/audit_cv_r1.md`)

**Architecture.**
- Site encoder: frozen ESM-2 650M ±7 window in full-protein context, plus a phospho channel for flanking U_ms sites, through a 2-layer transformer.
- Kinase encoder: ESM domain embedding, Cantley score-profile PCA and Manning group/family, with descriptor dropout of 0.2.
- Head: bilinear plus MLP.
- Loss: masked BCE plus 0.5 × filtered InfoNCE over kinases.
- Pretraining: 3,000 steps distilling label-free Cantley percentiles, then an auxiliary distillation loss of 0.2 during fine-tuning.
- 5 seeds.

| | M1 (A24) | M2 | M3 (A24) | LKO T2 M1 |
|---|---|---|---|---|
| R3 base | 0.735 | 0.183 | 0.247 | 0.650 |
| R3 final (distilled) | **0.754** | 0.223 | 0.294 | 0.686 |

**About the R3 numbers.** `docs/model_r3.md` reports 0.740 / 0.282. Those are pre-A24 values; the audited A24 values are shown here. The distillation gain is +0.019 M1 [0.010, 0.033] under the pre-A24 harness, and it holds out of the selection fold.

**Distillation also fixes the positive canary**, from 0.830 to 0.929.

**Not adopted:** ESM-2 3B windows (+0.005 alone, −0.015 combined with cpt_aux).

### 4.3 Late fusion (A13)

**Mechanism.** score = w·z(R2 ensemble) + (1−w)·z(R3), with z-scoring per fold. w = 0.6 is chosen on M1, and w_site = 0.5 on M2.

| w_R2 | 0 | 0.3 | 0.5 | **0.6** | 0.7 | 1 |
|---|---|---|---|---|---|---|
| M1 | 0.754 | 0.780 | 0.786 | **0.786** | 0.785 | 0.776 |
| M2 | 0.223 | 0.258 | **0.263** | 0.257 | 0.246 | 0.200 |
| M3 | 0.294 | 0.353 | 0.360 | 0.361 | 0.351 | 0.312 |

Source: `model_r2.md` §6.

- **Gain over the R2 ensemble:** ΔM1 +0.010 [0.005, 0.016], ΔM2 +0.063 [0.051, 0.078] (`results/models/r2/compare__r2_fusion_r3__vs__r2_final_ens5.tsv`).
- **Nested leave-one-fold-out selection** picks w = 0.6 in every fold, so M1 optimism is 0.000. Nested M2 is 0.259 against 0.263 (`docs/audit/audit_cv_r2.md` §3).

### 4.4 Regime and label-variant comparison: PLACEHOLDER

`docs/regimes_labels.md` is present but **uncommitted**, and its results sections are still template placeholders (RESULTS_TABLE, VERDICT_SECTION, …). This section will be replaced when the regimes agent finishes.

The pre-declared adoption rule is ΔM1 ≥ +0.005, an M2 drop ≤ 0.01, and the canary passing.

**Preliminary rows already in the committed `results/cv/REGISTRY.tsv`** (seed 0 unless stated; pooled CV; base = R2 v3 seed 0, M1 0.768, M2 0.192). These are not reviewed and carry no verdict:

| Variant or regime | M1 | M2 |
|---|---|---|
| L1 (human-only), features refit | 0.771 | 0.163 |
| L2_vitro, refit | 0.765 | 0.145 |
| L2_vivo-only, refit | 0.739 | 0.103 |
| L2 + external curated, refit (s0 / s1 / s2) | 0.774 / 0.769 / 0.771 | 0.218 / 0.200 / 0.198 |
| Partial pooling, ridge residual, λ = 300 (λ chosen on fold 0) | 0.768 | 0.199 |
| Per-family residual boosting, λ = 300 | 0.768 | 0.195 |
| Per-kinase LightGBM (N ≥ 20, otherwise global) | 0.756 | 0.174 |
| Per-Manning-group LightGBM | 0.761 | 0.170 |

**Reading, provisional.**
- Per-kinase and per-group models are worse than the global model.
- Partial pooling converges to the global model.
- In-vivo-only labels are much weaker.
- For L2+ext, the seed-0 gain (+0.006) does not hold across seeds. Seed-matched means over seeds 0–2 are 0.772 against 0.770 for L2 (computed from REGISTRY rows `r2_final_s{0,1,2}` and `regime_L2ext_refit_s{0,1,2}`).

All TEST deliverables were produced on L2.

[[REGIMES: final table, verdict, runtime from docs/regimes_labels.md once committed]]

### 4.5 Calibration and site assignment (A16, A17, A31; `model_r2.md` §6.1; `audit_cv_r2.md` R2-M1)

**Per-kinase calibration.** Cross-fitted tiered Platt scaling:
- per kinase for kinases with ≥30 OOF positives (61 kinases);
- family slope with a per-kinase intercept (261 kinases);
- global slope (17 kinases).

Calibration-in-the-large holds out of sample (expected vs observed positives by bin):

| Training-label bin | Expected / observed |
|---|---|
| 1–4 | 173 / 147 |
| 5–19 | 863 / 848 |
| 20–99 | 2,558 / 2,361 |
| ≥100 | 5,739 / 5,400 |

These probabilities are lower bounds under PU labelling.

**Site assignment, as first built.** p(k|s) is a softmax over the calibrated logits, T ≈ 0.93. On labelled DEV sites it was overconfident:

| p_top tier | n | claimed | observed |
|---|---|---|---|
| [0.5, 0.8) | 794 | 0.62 | 0.61 |
| [0.8, 1] | 167 | 0.87 | **0.54** |
| [0.8, 1], ST sites | 126 | 0.85 | 0.72 |
| [0.8, 1], **Y sites** | 41 | 0.91 | **0.00** |

On Y sites the Platt slope was negative in 4 of 5 folds.

**A31(b) fix.** A CV-only isotonic recalibration of p_top, separately for ST and Y sites, fitted on DEV OOF and cross-fitted by fold. After the fix:
- ST [0.5, 0.8): 0.62 claimed / 0.61 observed, n = 1,272.
- Only 2 sites reach ≥0.8.
- **No Y site reaches 0.5.**
- DEV coverage@0.5 = 21,115 sites, all ST; coverage@0.8 = 9.

**Conclusion.** The model has no reliable ≥0.8 tier and no reliable Y-site assignment at ≥0.5.

- The earlier coverage figures (14,235 at 0.5 and 1,674 at 0.8) are superseded.
- TEST `assignment_prob` uses maps fitted on all DEV OOF. The isotonic maxima are ST 0.756 and Y 0.445.
- Overall top-1 accuracy is 0.365, against 0.282 for Cantley.

---

## 5. Baselines and fair retraining

### 5.1 Cantley (B1) arms (`results/baselines/simple/summary.md`; re-evaluated values)

| Arm | M1 | M2 | M3 | M4-CV | Kinases scored (of 339) |
|---|---|---|---|---|---|
| (i) as published, Ochoa reference percentile | 0.717 [0.701, 0.731] | 0.157 | 0.133 | 0.690 | 318 |
| (ii) universe-recalibrated percentile | 0.717 | **0.171** | 0.133 | 0.690 | 318 |
| (iii) priming-aware, hard | 0.714 | 0.157 | 0.116 | 0.688 | 318 |
| imputed nearest-atlas (supplementary; not eligible) | 0.722 | 0.149 | 0.132 | 0.701 | 339 |
| B2 learned PSSM | 0.660 | 0.118 | 0.108 | 0.639 | 339 |
| B7 simple-strong logistic regression | 0.726 | 0.196 | 0.177 | 0.709 | 339 |
| B0-pop / B0-miss / B0-prior | 0.552 / 0.538 / 0.500 | | | | |

- Arms (i) and (ii) tie on M1 (difference 9e-7). Arm (ii) is better on M2 by +0.0135 [0.007, 0.020], hence A23.
- **TEST (A30).** The 10 STY kinases with a single-residue matrix are left unscored, so 386 of 503 kinases are scored. This can only help the baseline.
- `model_r2.md` and `RESEARCH_TREE.md` 5.14 quote B1 M1 = 0.719. That is the stale pre-15:31 REGISTRY value (audit r1, m-6). This report uses the re-evaluated 0.717.

### 5.2 B3 Phosformer (`docs/baseline_phosformer.md`)

| Arm | M1 | M2 | M3 | M4-CV | Note |
|---|---|---|---|---|---|
| B3-pub, released weights | 0.796 (A24; 0.793 in the doc) | 0.248 | n/a | 0.804 | **contaminated**: trained on PSP including our positives; positives + E-hard only; contrast only |
| **B3a retrained** (our MLM pretraining on the TEST-filtered corpus, then fine-tuned on our folds) | **0.679** (fold mean 0.690) | 0.132 | 0.104 | 0.657 | 339 kinases, E-all complete |
| Cantley B1 (i), same CV | 0.717 | 0.157 | 0.133 | 0.690 | |

**How B3a was built.**
- Architecture as in the paper: 6 layers, 12 heads, d = 768, focal loss, kinase domain + 11-mer input.
- 50-epoch MLM on the TEST-filtered corpus.
- Fine-tuned at 1:16 negatives. The PI fixed the ratio at the paper's value, replacing the pre-registered selection between 1:16 and 1:50.
- Seed 0 only.

**Result.** Retrained fairly on cluster splits, Phosformer is **below the label-free Cantley atlas on every endpoint**.

**Fidelity check: passed.** Our MLM checkpoint and fine-tuning code, trained on Phosformer's own published train split (60% uniform subsample, one epoch), give AUROC **0.888** on their test split. The paper reports 0.885 and the released weights give 0.8855 through our scorer. Per-kinase mean AUROC is 0.858 against the paper's 0.860 (`results/baselines/phosformer/fidelity_check.json`). B3a's shortfall on our benchmark therefore reflects the split, labels and negatives, not the implementation.

### 5.3 B4 SELPHI2.0 (`docs/baseline_selphi2.md`)

**Faithful port.**
- 44 model columns: 39 regenerated, 3 regenerated with a method substitution, 2 substituted; label-derived columns cross-fitted.
- The same random forest and full 72-configuration grid.
- 10 bags instead of 100; 10× masked random negatives.

| Arm | Universe | M1 | M2 | M3 |
|---|---|---|---|---|
| faithful | matched U_B4 (Ochoa sites, 34% of pairs) | 0.708 [0.685, 0.728] | 0.176 | 0.316 |
| noleak (4 suspect features dropped) | matched U_B4 | 0.724 | 0.190 | 0.265 |
| Cantley B1, same pairs | U_B4 | 0.728 [0.711, 0.749] | 0.182 | 0.123 |
| faithful | full universe (unscored pairs = worst) | 0.587 | 0.112 | 0.297 |

**Paired comparison, SELPHI − Cantley on shared pairs** (unadjusted):
- M1 −0.020 [−0.048, −0.001];
- M3 +0.187 [0.137, 0.224]. The M3 gain comes from site-popularity features.

**Why SELPHI reports 0.97.** The inflation comes from random-pair negatives, not from the split: the same "positives vs random pairs" metric gives 0.947 on our cluster split.

**Imputation caveat for TEST.**
- Only 31.3% of TEST sites have Ochoa features. 66.6% of TEST pairs use SELPHI's own NA defaults.
- On CV fold 0, imputed scoring gives M1 0.683, against 0.696 on Ochoa sites only and 0.587 when unscored pairs are ranked worst.
- So the TEST B4 is exhaustive and close to SELPHI's own behaviour, but it is evaluated on a universe it was not trained on.

---

## 6. CV results, with runtime

### 6.1 Leaderboard: main methods (pooled 5-fold OOF, A24, K_cv = 203)

| Method | M1 | fold-mean M1 | M2 | M3 | M4-CV | M1 ST / Y (K_any) | Neg/pos canary | Source |
|---|---|---|---|---|---|---|---|---|
| **Final: R2+R3 fusion** | **0.786** | 0.781 | **0.263** | **0.361** | 0.774 | — | 0.539 / 1.000 | `model_r2.md` §6; audit r2 |
| R2, 5-seed ensemble | 0.776 | 0.771 | 0.200 | 0.312 | 0.767 | — | 0.548 / 0.999 | `model_r2.md` §6 |
| R2 full (seed 0, pre-v3) | 0.773 | 0.768 | 0.207 | 0.312 | 0.770 | 0.797 / 0.644 | 0.544 / 0.998 | LEADERBOARD #5 |
| R2 clean (no benchmark PMIDs, no CPTAC) | 0.771 | — | 0.189 | 0.300 | 0.768 | — | 0.551 / 0.999 | `model_r2.md` §6.2 |
| R2 no-phospho (no F6) | 0.771 | 0.767 | 0.187 | 0.296 | 0.761 | — | 0.549 / 0.999 | REGISTRY |
| R3 final (distilled) | 0.754 | 0.749 | 0.223 | 0.294 | 0.731 | 0.761 / 0.586 | 0.524 / 0.929 | LEADERBOARD #11 |
| B7 simple-strong | 0.726 | 0.731 | 0.196 | 0.179 | 0.708 | 0.744 / 0.540 | — | LEADERBOARD #19 |
| B1 Cantley (i) / (ii) | 0.717 | 0.721 | 0.157 / 0.171 | 0.136 (full-cov, 196 kinases) | 0.690 | 0.732 / 0.534 | — | LEADERBOARD #26–27 |
| B4 SELPHI2.0 (matched) | 0.708 | — | 0.176 | 0.316 | 0.714 | — | not run | `baseline_selphi2.md` |
| B3a Phosformer retrained | 0.679 | 0.690 | 0.132 | 0.104 | 0.657 | — | not run | `baseline_phosformer.md` |
| B2 learned PSSM | 0.660 | 0.673 | 0.118 | 0.109 | 0.638 | 0.659 / 0.529 | — | LEADERBOARD #38 |
| B0-pop / B0-miss | 0.552 / 0.539 | | 0.020 / 0.040 | 0.142 / 0.053 | | | — | LEADERBOARD |

**`results/cv/LEADERBOARD.md` is stale.** It was generated at 18:54, so it lacks the ensemble, the fusion, B3a, R2 LKO, clean, no-phospho, F9 and A32 runs. Those rows come from the cited docs and REGISTRY.

### 6.2 Ablation curve (R2, seed 0, pooled, A24; `model_r2.md` §5; `results/models/r2/cv_table.tsv`)

| Step | M1 [CI] | M2 | M3 |
|---|---|---|---|
| F1 Cantley | 0.723 [0.707, 0.737] | 0.136 | 0.125 |
| + F0 site | 0.733 | 0.169 | 0.236 |
| + F2a learned PSSM | 0.738 | 0.094 | 0.239 |
| + F3 PPI | 0.748 | 0.143 | 0.260 |
| + F4 localisation | 0.753 | 0.141 | 0.269 |
| + F6a co-regulation | 0.755 | 0.149 | 0.274 |
| + F6b inhibitor | 0.755 | 0.154 | 0.281 |
| + F6c CPTAC (= G2 model) | 0.754 | 0.146 | 0.272 |
| + F5 co-expression | 0.759 | 0.160 | 0.283 |
| + F7 kinase descriptors (= full) | **0.773** [0.755, 0.785] | 0.208 | 0.312 |

**Reading.**
- F0 drives M3 (0.125 → 0.236).
- F3 adds +0.010 M1 and F7 adds +0.014 M1.
- **Phosphoproteomics adds about nothing to CV M1** (F6 in total about +0.002; the no-phospho variant scores 0.771). This is a negative result for the modality the brief emphasised.

### 6.3 Experiments that failed their adoption rule

| Experiment | Result | Rule | Source |
|---|---|---|---|
| **R5**: fine-tuned ESM late-interaction (ColBERT) cross-encoder | fold-0 M1 0.725 / M2 0.167 vs R3 0.741 / 0.221; removing the late-interaction term *helps* (0.727); the trunk overfits (early stop at step 100–600) | gate M1 ≥ 0.751 or M2 ≥ 0.241: **not met** | `model_r5.md` |
| Full cross-encoder | 3.6–14.3 GPU-h per fold for inference alone | infeasible | `model_r5.md` §1 |
| **knn_sim** | +0.003 mean M1 over 5 seeds; **ΔM2 −0.018 [−0.030, −0.006]** | ΔM1 ≥ 0.005: **not met** | `model_r2.md` §6 |
| **F9** Tyr/SH2 | +0.0017 M1; most of the GBDT gain in the sanity harness was a NaN-indicator artefact | **not met** | `model_r2.md` §6; `features_tyr.md` §3.2 |
| **A32(a)** priming-offset features | M1 +0.004, M2 +0.010; primed-positives Δ −0.006 [−0.022, 0.008] | ≥ +0.02 on primed: **not met** | `strata_cv.md` §9 |
| **A32(b)** prior-debiased score_site | β = 0 best pooled and in every nested fold; β = 0.1 costs 0.027 M2, and the 5–19 bin moves only 0.142 → 0.146 | **not met** | `strata_cv.md` §9 |
| R3 + ESM-2 3B windows | +0.005 alone, −0.015 combined | **not met** | `model_r3.md` §2 |

### 6.4 Strata (`docs/strata_cv.md`; R2 = `r2_final_knn_s0`, **not** the final fusion; M1 macro E-hard over K_cv)

| Stratum | n kinases | R2 [CI] | R3 | B1 (i) | B1 primed | B7 | B0-pop |
|---|---|---|---|---|---|---|---|
| all | 203 | 0.774 [0.756, 0.786] | 0.754 | 0.717 | 0.714 | 0.726 | 0.552 |
| ST kinases | 161 | 0.808 | 0.783 | 0.751 | 0.747 | 0.763 | 0.538 |
| **Y kinases** | 38 | **0.629** [0.580, 0.659] | 0.624 | 0.594 | 0.594 | 0.569 | **0.592** |
| Y, hub sites excluded | 38 | 0.627 | 0.623 | 0.586 | 0.589 | 0.560 | 0.583 |
| 1–4 labels (M5, K_any) | 79 | 0.741 [0.690, 0.787] | 0.647 | 0.624 | 0.628 | 0.651 | 0.565 |
| 5–19 labels | 83 | 0.758 | 0.748 | 0.689 | 0.683 | 0.705 | 0.571 |
| ≥100 labels | 31 | 0.815 | 0.786 | 0.766 | 0.768 | 0.775 | 0.543 |
| in vivo only | 172 | 0.766 | 0.768 | 0.707 | 0.704 | 0.720 | 0.583 |
| in vitro only | 181 | 0.781 | 0.750 | 0.729 | 0.720 | 0.734 | 0.529 |

**M2 (Hit@5) by training bin.** The cross-kinase ranking is dominated by well-studied kinases:

| Bin | R2 (knn run) | R3 | B1 (ii) | Final fusion |
|---|---|---|---|---|
| 5–19 | **0.099** | 0.092 | **0.162** | 0.142 |
| 20–99 | 0.164 | 0.216 | 0.171 | — |
| ≥100 | 0.410 | 0.592 | 0.193 | 0.621 |

The B0-prior floor alone reaches 0.367 in the ≥100 bin. Final-fusion values are from `strata_cv.md` §9 at β = 0.

**The Tyr gap** (`features_tyr.md` §1; `strata_cv.md` §5):
- Y kinases sit barely above the B0-pop floor (0.592) for every method.
- About 1/4 of the raw gap is macro dilution: 42% of Y kinases have ≤3 positives.
- The rest has two sources:
  - an E-hard pool dominated by promiscuous pTyr hubs (GLO1 Y136 is labelled to 15 TKs) and probable indirect in-vivo events;
  - Tyr atlas matrices that are anti-correlated with their own labels for FGFR1, BTK, Eph, VEGFR and TAM.
- Excluding hub sites does not change R2 (0.627). F9 did not help.
- **The gap is unresolved.**

### 6.5 Unorthodox kinases (`strata_cv.md` §4–5, §9)

| Group | R2 | Best comparator | Note |
|---|---|---|---|
| Priming kinases (9: GSK3, CK1, CK2α) | 0.837 | B7 0.812; B1 primed 0.792 | primed positives score above unprimed within kinase (+0.078; 6/7 kinases) |
| **Primed positives** (231) | 0.870 (knn run); **0.893** (r2_final_s0) | **B1 primed 0.945** [0.924, 0.962] | Cantley's priming-aware arm wins; the A32(a) fix failed |
| Unprimed GSK3B positives | 0.684 | — | hard for every method |
| PIKK (4; small) | 0.880; S/T-Q PIKKs 0.926 | B7 0.876; B1 0.833 | MTOR weakest (0.742) |
| FAM20C (1 kinase, 6 positives) | 0.980 | B1 0.974 | anecdotal |

### 6.6 Leave-kinase-out T2 (unseen family, seen proteins; `model_r2.md` §6.3; `audit_cv_r2.md` R2-m2)

| Model | M1 | M3 |
|---|---|---|
| R2 v3, seed 0 | **0.751** | 0.295 |
| R2+R3 fusion | 0.749 | 0.337 |
| R3 final (A24 off) | 0.686 | 0.226 |

T2 M1 split by whether the held-out site is labelled to another kinase:

| Held-out positives | share | T2 M1 | CV control |
|---|---|---|---|
| labelled to another kinase | 45% | 0.798 | 0.799 |
| **not labelled to any other kinase** | 55% | **0.703** | 0.757 |

**Reading.**
- On truly novel sites, the cost of an unseen family is about 0.05 M1, not the 0.02 average.
- T2 is "unseen labels", not kinase-naive: label-free Cantley percentiles of held-out kinases enter R3 pretraining.
- LKO M2 is not meaningful and is not reported.

### 6.7 Canaries, per fold (`audit_cv_r2.md`; `audit_cv_r1.md`)

| Model | Pooled negative (gated, ≤0.57) | Per fold f0–f4 | Per-fold mean | Positive |
|---|---|---|---|---|
| Final fusion | 0.539 PASS | 0.551 / **0.575** / 0.481 / 0.520 / **0.676** | 0.561 | 1.000 |
| R3 final | 0.501 (raw) / 0.524 (A24) | 0.498 / 0.522 / 0.456 / 0.530 / **0.667** | 0.535 | 0.929 |
| R2 5 seeds (pooled) | 0.525–0.550 | — | — | 0.999 |

**Two of five folds exceed 0.57.** The A31(c) ruling keeps the gate pooled.

- Every model, including label-free controls (R2 F0+F7: 0.539 on fold 0), sits above 0.5.
- This fits a sampler-population difference (protein-first canary sampling against site-weighted E-hard), not leakage.
- **Not proven.** A site-uniform v3 canary was recommended and not built.

### 6.8 Runtime (wall-clock; shared node: 64 cores, 4× RTX 4090; contention inflates GPU times by about 2×)

| Method | CV | TEST: train on DEV | TEST: infer 25.9 M pairs | Source |
|---|---|---|---|---|
| **Final R2+R3 fusion** | sum of its components | **1,136.6 s** | **316.9 s** | `results/freeze_entries/final.json` |
| R2, 5 seeds | per fold: cross-fitted features ≈21 CPU-min, fit 37 s, OOF predict 28 s | 417.9 s | 298.9 s | `r2_alone.json`; `model_r2.md` §3 |
| R2, T3 (5 LKO folds × 5 seeds) | — | 2,520 s | 406 s | `model_r2.md` §6.4 |
| R3, 5 seeds (GPU) | 5.5–7.5 min per fold | 718.7 s | 3.1 s | `r3.json`; `model_r3.md` §4 |
| R3, T3 | — | 3,446 s | 9.4 s | `r3.json` |
| B1 Cantley (i) / (ii) | ≈1 s fit+predict | 0 | 4.4 / 1.0 s | `results/baselines/simple/freeze_entries/ALL.json` |
| B2 PSSM | 146 s | 3.3 s | 39.7 s | ALL.json; `summary.md` |
| B7 LR | 925 s | 89.1 s | 287.0 s | same |
| B0-miss | 138 s | 19.0 s | 30.4 s | same |
| B3 Phosformer | MLM pretraining ≈4.3 h on 1 GPU; fine-tune and E-all scoring per fold (not recorded in manifests) | 3,152 s (fine-tune only) | **6,479 s** | `freeze_entries/b3.json`; `baseline_phosformer.md` |
| B4 SELPHI2.0 (12–16 cores) | 109 min for 5 folds | 737 s | 2,043 s | `freeze_entries/b4.json`; `baseline_selphi2.md` |
| R5 (not adopted) | 11 min per fold (35M) / 21 min (150M); scoring 12 s per fold | — | — | `model_r5.md` §4 |
| R2 clean, U_obs library (131 M rows) | — | 274 s | 729 s | `model_r2.md` §6.2 |
| Feature builds | ESM-2 650M 2.8 min (0.09 GPU-h), 3B 7.9 min; AFDB structure 15.7 min (15 CPU); F1 151 s; context ≈2 min; F6a 17 s per full fit (GPU); knn dense pass 29 s | | | feature docs |
| TEST batch, B = 10k | rehearsal on fold 0: 21 min; TEST estimated 45–60 min | | | `docs/G4_RUNBOOK.md` |

**Cost–accuracy.** The final model is both the most accurate and among the cheapest to score. It infers about 20× faster than B3 and about 6× faster than B4, measured on the same candidate set.

---

## 7. TEST results: PLACEHOLDERS (single batch not yet run)

Every table below is filled from `results/test/<sha12>/` after `eval_test.py --real`. **No TEST number may be quoted before the COMPLETED ledger block exists.**

### 7.1 T1: seen kinase, unseen protein (K_eval = 119; K_any = 256)

| Config | M1 [CI] | M2 [CI] | M3 [CI] | M4 [CI] | Kinases scored |
|---|---|---|---|---|---|
| final (R2+R3) | [[TEST: final M1 + CI]] | [[TEST]] | [[TEST]] | [[TEST]] | 503 |
| R2 alone | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | 503 |
| R3 alone | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | 503 |
| B1 (i) / (ii) | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | 386 |
| B3 retrained | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | 503 |
| B4 retrained | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | 503 (66.6% of pairs imputed) |
| B7 / B2 | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | |
| B0-prior / pop / miss | [[TEST]] | [[TEST]] | [[TEST]] | [[TEST]] | |

Kinase-resampling sensitivity CIs: [[TEST: kinase-resampling CIs for final]].

### 7.2 The 9 primary tests (B = 10,000, paired, matched universe, Holm over 9)

| # | Comparison | Endpoint | Δ [95% CI] | raw p | Holm p | Verdict |
|---|---|---|---|---|---|---|
| 1 | final − B1 (i) | M1 | | | | |
| 2 | final − B1 (i) | M3 | | | | |
| 3 | final − B1 (ii) | M2 | | | | |
| 4 | final − B3 | M1 | | | | |
| 5 | final − B3 | M2 | | | | |
| 6 | final − B3 | M3 | | | | |
| 7 | final − B4 | M1 | | | | |
| 8 | final − B4 | M2 | | | | |
| 9 | final − B4 | M3 | | | | |

On the fold-0 rehearsal, final − B1(i) on M1 gave Holm p 0.018 at B = 1,000 and 0.050 at B = 10,000 (A29). Expect TEST p-values near the threshold to be sensitive to B.

### 7.3 Secondary comparisons (Holm separately; `scripts/build_freeze.py` SECONDARY)

| Comparison | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| final − R2 | | | | |
| final − R3 | | | | |
| final − B7 | | | | |
| final − B2 | | | | |
| final − B1 (i) | — | — | — | |

### 7.4 T3: unseen family, unseen protein

| Config | M1 | M3 | M4 | Missing kinases |
|---|---|---|---|---|
| final T3 (R2 5 seeds + R3, fused) | [[TEST: T3 M1]] | | | |
| R2 T3 | | | | |
| R3 T3 | | | | |
| B2 / B7 T3 | | | | |

### 7.5 TEST strata

| Stratum | final | B1 (i) | B3 | B4 |
|---|---|---|---|---|
| ST kinases | | | | |
| Y kinases | | | | |
| Y, hub sites excluded | | | | |
| training bins 0 / 1–4 / 5–19 / 20–99 / ≥100 (M4/M1) | | | | |
| evidence: in vivo only / in vitro only / both | | | | |
| priming kinases; primed positives | | | | |
| PIKK; FAM20C | | | | |
| Manning groups | | | | |

### 7.6 Temporal (new since 2021: 263 TEST positives) and external curated (215 pairs; tier A 158)

| Set | final | B1 (i) | B3 | B4 | Calibration (expected / observed) |
|---|---|---|---|---|---|
| new since 2021 | [[TEST: temporal]] | | | | |
| external tier A | [[TEST: external A]] | | | | |
| external tier B | | | | | |

### 7.7 A17 assignment reliability on labelled TEST sites

| p_top tier | ST: n / claimed / observed | Y: n / claimed / observed |
|---|---|---|
| [0, 0.2) | | |
| [0.2, 0.5) | | |
| [0.5, 0.8) | | |
| [0.8, 1] | expected about 0 sites (isotonic maximum ST 0.756) | expected 0 (maximum Y 0.445) |

ECE: [[TEST: ECE]]. TEST coverage@0.5 over all TEST candidate sites (label-free): 5,773 sites, all ST; coverage@0.8 = 0 (`model_r2.md` §6.4).

### 7.8 CV→TEST gap (leakage alarm)

[[TEST: ΔM1 TEST − CV for final, R2 and B1]].

**Prior expectation.** About 90 configurations were tried (about 55 selectable). The audited optimism estimate is ≲0.005 M1, before any split-shift effect. A gap well beyond that should trigger a leakage review.

---

## 8. Downstream: kinase-activity inference (benchmarKIN; `docs/downstream.md`)

**Harness reproduces Müller-Dott 2025.**

| Library | Our perturbation AUROC | Paper | Our tumour AUROC | Paper |
|---|---|---|---|---|
| PSP + z-score | 0.782 | 0.79 | 0.664 | 0.662 |
| PSP + KinLib | 0.650 | 0.65 | 0.583 | 0.58 |

Shuffled control: 0.518 ± 0.028.

**Library tested.** Clean R2: DEV-trained, with the 22 benchmark PMIDs and all CPTAC features removed (CV M1 0.771). Disclosed residual channel: site existence in U_ms draws on the compendium and CPTAC (A31e).

Matched-kinase comparisons (z-score, U_obs, ours − reference; unadjusted CIs; descriptive endpoint):

| Comparison | Perturbation, all | **Drug unseen** (114 experiments) | Tumour (7 CPTAC cohorts) |
|---|---|---|---|
| ours topPSP vs PSP_all | −0.017 [−0.054, 0.018] | −0.058 [−0.123, 0.003] | +0.016 [−0.018, 0.051] |
| ours top50 vs PSP_all | −0.013 [−0.043, 0.012] | **−0.050 [−0.096, −0.001]** | +0.011 [−0.027, 0.046] |
| ours top100 vs PSP_all | −0.019 [−0.053, 0.012] | −0.041 [−0.099, 0.017] | +0.029 [−0.010, 0.071] |
| ours top100 ∪ PSP vs PSP_all | +0.001 (n.s.) | −0.011 | **+0.034 [0.005, 0.068]** |
| ours top{PSP,50,100} − PSP vs PSP_all | **−0.045 / −0.052 / −0.068** | worse | +0.013 to +0.020 (n.s.) |
| ours top100 vs Cantley top100 | **+0.176 [0.091, 0.263]** | +0.095 [−0.019, 0.208] | **+0.114 [0.078, 0.149]** |

**Coverage of evaluable kinases:**
- perturbation: 50–65 kinases (top50 / top100), against 41 for PSP;
- tumour: 149–274 kinases, against 73 for PSP.

**Reading.**
- On matched kinases the library is **at parity with PSP, not better**.
- It is much better than Cantley.
- It makes more kinases evaluable.
- Adding it to PSP helps the tumour benchmark slightly.
- **Unseen-drug gap.** On experiments whose drug the model never saw (mostly Hijazi), it trails PSP by 0.04–0.06.
- Removing *all* phospho features (no-phospho variant) does not close the gap; it widens it to −0.08 at top100 (`model_r2.md` §6.2b).
- The drug effect is confounded with the Hijazi dataset. A true leave-drug-out retrain was not done.

**`model_r2.md` §6.2 overstates the tumour result.** It says the clean library "beats PSP on the tumour benchmark (0.692 vs 0.664)". That is a native, not like-for-like, contrast; matched, the difference is +0.016 and not significant.

---

## 9. Coverage: PLACEHOLDER (deployment model over U_obs)

**Not yet run.** The deployment model (retrained on DEV+TEST with the frozen configuration, §8.6) and its U_obs predictions are G5 items.

| | Coverage@0.5 | Coverage@0.8 | Sites in some kinase's top 1% |
|---|---|---|---|
| Deployment model, all U_obs (350,736) | [[COVERAGE: U_obs @0.5]] | [[COVERAGE: U_obs @0.8]] | [[COVERAGE]] |
| … of which sites with no PSP kinase | [[COVERAGE]] | [[COVERAGE]] | [[COVERAGE]] |
| Random-score baseline | [[COVERAGE]] | [[COVERAGE]] | |
| Cantley-only baseline | [[COVERAGE]] | [[COVERAGE]] | |

**Expectation from CV under A31.** DEV coverage@0.5 = 21,115 sites, all ST; coverage@0.8 = 9.

Three caveats apply to any coverage number:
- These are model-based estimates, conditional on the site being a target of a modelled kinase.
- Under PU labelling they are lower-bound precision statements.
- Counts are comparable across methods only when the methods share a candidate set (`model_r2.md` §6.1).

---

## 10. Limitations

1. **Tyr specificity is unsolved.** Y-kinase M1 is about 0.63 against about 0.81 for S/T; no Y assignment reaches 0.5; F9 failed. Part of the gap is endpoint composition (hubs, indirect in-vivo labels), which no feature can fix.
2. **Assignment probabilities were overconfident until A31.** After the fix there is no ≥0.8 tier, so "which kinase?" answers are ranked suggestions (M2 Hit@5 0.263), not confident calls.
3. **M2 favours well-studied kinases.** Hit@5 is 0.142 at 5–19 labels against 0.621 at ≥100; Cantley gets 0.162 at 5–19. The debiasing attempt (A32b) failed. This conflicts with the brief's "harmonious across kinases" goal.
4. **Phosphoproteomics adds almost nothing to CV M1** (F6 ≈ +0.002; no-phospho 0.771). Downstream, the phospho features may help on unseen drugs, but the evidence is indirect.
5. **Downstream parity, not superiority, with PSP**, plus an unseen-drug gap of 0.04–0.06.
6. **The canary gate passes only pooled.** Two of five folds exceed 0.57, and every model sits above 0.5. The sampler explanation is plausible but unproven.
7. **LKO generalisation is weaker on novel sites** (0.703 vs 0.757 control), and T2 is not kinase-naive.
8. **Baseline fidelity.**
   - The B3 training pipeline reproduces the published number (0.888 vs 0.885).
   - B3 has 1 seed and a PI-fixed ratio.
   - B4 has 10 of 100 bags, a substituted PSSM/coreg, and 66.6% of TEST pairs imputed.
   - B1 is exact, reproducing kinase-library 1.8.0 to 100% within 1e-3.
9. **Optimism.** About 90 CV configurations were tried. A24 was a post-freeze metric correction applied to all methods; it widened learned-vs-Cantley margins by about 0.01.
10. **Label noise and PU.** In-vivo-only labels are enriched for indirect events (S/T sites assigned to TKs are 69% in vivo only). All precision and probability statements are lower bounds.
11. **Small strata.** PIKK (4 kinases), FAM20C (6 positives) and STY (4 in K_cv) are anecdotal.
12. **Unbuilt items.** F8 (Sugiyama), docking motifs, F6d, U_all (unobserved S/T/Y) and the deployment model. `scripts/run_all.sh`, promised in DESIGN §11.8, does not exist.
13. **Data gaps** (DATAINEED): SIGNOR direct export, decryptM curve fits, CPTAC plex/ischemia confounders and PTB matrices were unavailable; workarounds are in use.

---

## 11. Follow-ups and remaining ideas

| Priority | Idea | Rationale / evidence |
|---|---|---|
| 1 | Run the TEST batch; fill §7; commit FREEZE.json and the ledger | G4 |
| 1 | Deployment model plus U_obs coverage (§9) | Priority 2 of the brief |
| 2 | Kinase-gated priming: use primed-hard Cantley only for priming kinases at their own offset, or a per-kinase blend | B1 primed beats R2 on primed sites (0.945 vs 0.893); A32(a) additive features failed |
| 2 | M2 for understudied kinases: per-kinase prior removal fitted on CV, or a blend toward Cantley below 20 labels | 5–19 bin Hit@5 0.142 vs Cantley 0.162; a global β failed |
| 2 | Tyr: nearest-atlas imputation for anti-correlated matrices (FGFR, Eph, VEGFR, TAM, BTK); a pre-registered A-loop/hub sensitivity analysis | `features_tyr.md` §1.5–1.6 |
| 2 | True leave-drug-out retrain (F6 with `exclude_drugs`) | separates the drug effect from the Hijazi-dataset effect |
| 3 | Site-uniform v3 negative canary | turns the pooled-only pass into a clean null |
| 3 | R5 revisited: LoRA or last-k-layer fine-tuning, a lower trunk learning rate, protein-grouped early stopping | the trunk overfit by step 100–600 |
| 3 | decryptM potency coherence (EC50 vs Kd); Bayer–Kuster 2025 dose–response data | closer to direct-substrate evidence (DATAINEED 4, 8) |
| 3 | Temporal retrain on 2021-only labels; L3 variant; KolossuS / AF-Multimer pilot | deferred in DESIGN §0 |
| 3 | Seeds 1–4 for B3; 100 bags for B4 | baseline strength |

---

## 12. Reproducibility

**Commit and harness.**
- Harness commit: `974d520` (A28), covering `eval.py`, `labels.py`, `splits.py`, `crossfit.py` and `canary.py`.
- Repository HEAD at drafting time: `be635d1`.
- TEST deliverables are pinned by SHA-256 in `results/freeze_entries/*.json` and `results/baselines/*/freeze_entr*`.

**Seeds.**
- Splits and annealing: seed 0.
- R2 and R3 final models: seeds 0–4.
- Ablations, B3, B4 and regimes: seed 0 (L2ext also seeds 1–2).
- Bootstrap: seed 0.
- downstream mcauroc: seed 42; tumour: seed 1.
- There is no single central seed list in code; seeds are recorded per run manifest.

**Data checksums.**
- `DATA.md`: every source, version, SHA-256 and licence.
- `logs/dl/sha256_wave1.txt`.
- `results/index_manifest.json`: index hashes for kinases, sites and proteins.
- `results/splits/SEAL.json`: SHA-256 of the sealed files.

**Commands.** These come from the cited docs; `PYTHONPATH=src` and `.venv/bin/python -I` are assumed throughout.

```
# foundations
scripts/build_kinase_list.py; scripts/build_site_universe.py; scripts/build_labels.py; scripts/build_splits.py
python -m pytest -q tests/                                   # core harness tests
# features
scripts/features/motif_cantley_build.py; scripts/features/context_{ppi,loc,coexpr}.py
scripts/features/phospho_{build_blocks,f6b,f6c}.py; scripts/export_structure_esm_features.py
scripts/models/r2_xfit.py --family f2a f3nb f6a f6a_clean --folds 0 1 2 3 4   # cross-fitted stores (--toy = leak test)
# CV
scripts/models/r2_cv5.sh ... ; scripts/models/r3_cv.py train --fold F --config w7_nce05_cpt_aux ; scripts/models/r2_fusion.py
# TEST deliverables
scripts/models/r2_final.py --name r2_final_test --target test --seeds 0 1 2 3 4 --set num_leaves=127 learning_rate=0.03 f2a_family=f2a_sep desc_v=2 es_real_only=true
scripts/models/r2_final.py --name r2_final_t3 --target t3 --seeds 0 1 2 3 4 --set ...
scripts/models/r3_cv.py test|assemble-test ; scripts/models/r2_fuse_test.py --name r2r3_final_test ... ; scripts/models/r2_assign_test.py
scripts/baselines/simple_test.py ; scripts/baselines/selphi2_final.py --mode test ; scripts/baselines/phosformer_test.py
src/ksp/eval_test.py --check-predictions cache/runs/<run_id>/test_predictions.parquet
# FREEZE and the single batch
scripts/build_freeze.py ; git commit FREEZE.json ; src/ksp/eval_test.py --real
# downstream
scripts/downstream/pred_topk.py --tag r2_clean ; scripts/downstream/run_final.py --tag r2_clean --methods zscore
```

**Run records.**
- Every run writes `results/runs/<run_id>/manifest.json`: git hash, config, seed, input SHA-256, runtime, fallbacks.
- Predictions go to `cache/runs/<run_id>/`. They are gitignored and recoverable only by re-running; their pins record their contents.
- Configurations tried are in `results/cv/REGISTRY.tsv`. The audit r2 notes that REGISTRY undercounts fusion grid points and the M2-arm choice.
