# REPORT: Predicting direct kinase–substrate sites across the human proteome

**Status: v3 (2026-10-27 09:46 CET).** The single TEST batch ran once (FREEZE 721dc197) and was independently reproduced (`docs/audit/audit_test.md`). Deployment model and coverage (§9) done. Post-TEST exploratory CV (§11.1) audited (`docs/audit/audit_posttest.md`). Each number cites the committed file it comes from.

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

**TEST results: single batch, FREEZE 721dc197, completed 2026-10-27 00:57, independently reproduced.**

Standalone values use missing = 'worst', with 95% cluster-bootstrap CIs (B = 1,000). Each Δ is the final model minus the baseline in the matched-universe primary test (B = 10,000, Holm over 9 tests).

| Endpoint | Final | Cantley B1 | B3 Phosformer (our retrain) | B4 SELPHI2.0 (our retrain) |
|---|---|---|---|---|
| M1 | **0.794** [0.774, 0.811] | 0.714 (arm i); Δ **+0.075** [0.046, 0.099] | 0.708; Δ **+0.086** [0.059, 0.105] | 0.705; Δ **+0.089** [0.063, 0.117] |
| M2 Hit@5 | **0.363** [0.324, 0.397] | 0.189 (arm ii); Δ **+0.179** [0.147, 0.213] | 0.163; Δ **+0.199** [0.151, 0.229] | 0.233; Δ **+0.130** [0.091, 0.157] |
| M3 top-1% recall | **0.363** [0.309, 0.401] | 0.161 (arm i); Δ **+0.205** [0.165, 0.235] | 0.104; Δ **+0.259** [0.213, 0.291] | 0.304; Δ +0.059 [0.010, 0.120] (modest) |

Sources: `results/test/LEDGER.md`; `results/test/721dc1977366/{final,B1_i,B1_ii,B3,B4}_ci.tsv`; `superiority.tsv`.

**All 9 primary tests meet the "beats" rule.**
- 8 of them have p < 0.0002 (bootstrap floor), with Holm p ≤ 0.0018.
- The B4 M3 test has p 0.018 (Holm 0.018).

**Unseen family and unseen protein (T3).** M1 is 0.752 (`results/test/LEDGER.md`).

**Runtime.** Training takes 1,137 s and scoring the 25.9 M TEST pairs 317 s (`results/test/LEDGER.md`).

**Read these numbers with the following in mind.** The independent review's caveats (`docs/audit/audit_test.md` §0, §7) govern how they may be stated.

- **M1 shows no CV→TEST optimism.** On the same 117 kinases, CV M1 is 0.792 and TEST is 0.793 (paired Δ +0.001 ± 0.012).
- **TEST M2 (0.363) is not an improvement over CV.** The rise is K_eval composition: on the same 117 kinases, CV M2 is 0.348.
- **The TEST margin over Cantley (+0.075) is partly a weak-baseline draw.** On CV the margin was +0.063 [0.048, 0.074]. On the same kinases Cantley fell by 0.025 from CV to TEST, mostly on Y kinases.

The main caveats:
1. **No demonstrated Tyr-kinase specificity on TEST.**
   - Y kinases score M1 0.595, barely above the site-popularity floor B0-pop (0.585), and the final model beats that floor on only 10 of 23 Y kinases.
   - Cantley is at chance (0.496).
   - On unseen families (T3), Y kinases score 0.531.
   - Y-kinase M1 dropped by 0.072 ± 0.033 from CV.
   - The review rates this NOT TRUSTED as a skill claim (`docs/audit/audit_test.md` §0 #10, §4).
2. **No reliable high-confidence site assignment.**
   - On labelled TEST S/T sites, the ≥0.5 tier is honest: 0.615 claimed vs 0.603 observed (n = 345).
   - There is no ≥0.8 tier, and no Y site can exceed 0.445.
   - Coverage is 5,773 TEST sites at 0.5 and 0 at 0.8 (`results/test/721dc1977366/summary.json`; audit §6).
3. **M2 is still weak for understudied kinases.** TEST includes only 13 kinases with 5–19 labels, and they are the best-labelled of that bin. The CV finding (Hit@5 0.142 against 0.162 for Cantley) is therefore not refuted (audit §2).
4. **B3 and B4 are our faithful retrains on this benchmark.**
   - The B3 fidelity check passed.
   - 66.6% of B4's TEST pairs use SELPHI's NA-default site features.
   - The B4 M3 win is modest: CI lower bound 0.010, with 52 kinases favouring the final model and 27 favouring B4.
5. **Priming.** The primed-site advantage (0.897 vs 0.704) holds only against non-priming-aware Cantley. On unprimed positives the two are near-tied (0.740 vs 0.725). On CV the priming-aware Cantley arm led (0.945); that arm was not in FREEZE.
6. **Strata carry no CIs in the batch, so all stratum statements are descriptive.** The gains on hub sites and external tier B are mostly popularity (B0-pop 0.788 and 0.724).
7. **Downstream, the library is at parity with PSP, not better, and trails PSP on unseen drugs** by 0.04–0.06 (§8).
8. **About 90 configurations were tried on CV.** TEST confirms the audited optimism bound of ≲0.005.
9. **Post-TEST CV exploration (A34, §11.1) found no clean fix** for caveats 1, 3 or 5:
   - Tyr: a Y-specific fusion weight and nested motif repair both fail.
   - Understudied-kinase M2: per-kinase intercepts gain only in a post-hoc variant, and only by taking top-5 slots from well-studied kinases.
   - Priming: the gated blend is a rank shift that a label-free priming-site flag reproduces.
   - The released model is unchanged.

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
- The FREEZE validated with all 9 primaries and was committed (5fc7870, sha 721dc197…) before the single TEST batch ran (00:10–00:57; §7, `results/test/LEDGER.md`).

**Selection discipline.**
- Model-selection exploration on CV closed at 01:00 (A28). Post-TEST exploratory CV analyses (§11.1, A34) did not change the model.
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

### 4.4 Training regime and label variants (`docs/regimes_labels.md`; A33)

**Setup.**
- Pooled 5-fold CV with the frozen R2 config and seed 0, unless stated.
- Every label-derived feature family (F2a, F3-label, F6a) was refitted for each label variant.
- Sanity check: L2 run through the regime pipeline reproduces R2 seed 0 exactly (M1 0.768153).
- Each Δ is against the same-seed R2 run.

| Config | M1 | ΔM1 | M2 | ΔM2 | M3 | M1 ST / Y | M5 bins 1–4 / 5–19 / 20–99 / ≥100 |
|---|---|---|---|---|---|---|---|
| **L2 (human + mouse/rat), global (default)** | 0.768 | — | 0.192 | — | 0.298 | 0.789 / 0.626 | 0.731 / 0.759 / 0.772 / 0.806 |
| L1 human-only | 0.771 | +0.003 | 0.163 | **−0.029** | 0.295 | 0.794 / 0.587 | 0.713 / 0.761 / 0.772 / 0.807 |
| L2_vitro (in vitro evidence only) | 0.765 | −0.003 | 0.145 | −0.047 | 0.276 | 0.791 / 0.586 | 0.726 / 0.744 / 0.773 / 0.807 |
| L2_vivo (in vivo only) | 0.739 | −0.029 | 0.103 | −0.090 | 0.225 | 0.754 / 0.569 | 0.684 / 0.716 / 0.741 / 0.788 |
| L2 + external curated, seeds 0 / 1 / 2 | 0.774 / 0.769 / 0.771 | +0.006 / −0.004 / +0.002 | 0.218 / 0.200 / 0.198 | +0.026 / +0.018 / +0.002 | 0.309 | 0.797 / 0.585 | 0.718 / 0.760 / 0.774 / 0.806 |
| Per-kinase LightGBM (N ≥ 20, else global) | 0.756 | −0.012 | 0.174 | −0.019 | 0.259 | 0.781 / 0.622 | 0.731 / 0.759 / 0.749 / 0.793 |
| Per-kinase logistic | 0.755 | −0.013 | 0.171 | −0.021 | 0.261 | 0.782 / 0.617 | |
| Partial pooling, ridge residual (λ = 300, chosen on fold 0) | 0.768 | +0.000 | 0.199 | +0.007 | 0.298 | 0.790 / 0.626 | ≈ global |
| Partial pooling, boosted residual (λ = 300) | 0.768 | +0.000 | 0.195 | +0.002 | 0.298 | 0.789 / 0.625 | ≈ global |
| Per-Manning-group LightGBM | 0.761 | −0.007 | 0.170 | −0.022 | 0.298 | 0.776 / 0.625 | 0.728 / 0.735 / 0.772 / 0.803 |

**Findings.**
- **One global model is best on CV, within the tested regimes.**
  - Per-kinase models lose 0.012 M1. On the 78 kinases that reached N ≥ 20 in every fold, the global model still wins (0.805 vs 0.779).
  - Partial pooling selects λ = 300, which is effectively the global model.
  - Per-group models are in between.
- **Pooling mouse and rat labels helps the site-centric endpoint** (M2 +0.029 over human-only) at no M1 cost.
- **In vitro vs in vivo.**
  - In-vivo-only labels are clearly worse on every endpoint, consistent with their enrichment for indirect events.
  - In-vitro-only training raises AUROC on in-vitro-only test positives slightly (0.781 vs 0.772) but lowers M2 sharply.
- **External curated labels (L2+ext).** The seed-0 M1 gain does not replicate (3-seed mean +0.001). They were pre-registered as ablation-only, to keep the external TEST stratum independent, so they were not adopted (A33).
- **Side finding.** Per-kinase ridge intercepts (λ = 3–10) lift R2-alone M2 from 0.192 to 0.232 for 0.001–0.003 M1. CV had closed, so this was not explored (§11).
- **Runtime.** A 5-fold variant run takes 6–29 min. Refitting the label-derived features for all variants took 272 job-minutes.
- **Configurations tried.** 24 configs, logged in REGISTRY as 37 rows.

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
- `model_r2.md` quotes B1 M1 = 0.719 (RESEARCH_TREE 5.14 quotes 0.717 and flags 0.719 as stale). 0.719 is the stale pre-15:31 REGISTRY value (audit r1, m-6). This report uses the re-evaluated 0.717.

### 5.2 B3 Phosformer (`docs/baseline_phosformer.md`)

| Arm | M1 | M2 | M3 | M4-CV | Note |
|---|---|---|---|---|---|
| B3-pub, released weights | 0.796 (A24; 0.793 in the doc) | 0.248 | n/a | 0.804 | **contaminated**: trained on PSP including our positives; positives + E-hard only; contrast only |
| **B3a retrained** (our MLM pretraining on the TEST-filtered corpus, then fine-tuned on our folds) | **0.677** (fold mean 0.689; leaderboard re-evaluation — the doc gives 0.679) | 0.132 | 0.104 | 0.655 (LEADERBOARD; `baseline_phosformer.md` gives 0.657) | 339 kinases, E-all complete |
| Cantley B1 (i), same CV | 0.717 | 0.157 | 0.133 | 0.690 | |

**How B3a was built.**
- Architecture as in the paper: 6 layers, 12 heads, d = 768, focal loss, kinase domain + 11-mer input.
- 50-epoch MLM on the TEST-filtered corpus.
- Fine-tuned at 1:16 negatives. The PI fixed the ratio at the paper's value, replacing the pre-registered selection between 1:16 and 1:50.
- Seed 0 only.

**Result.** Retrained fairly on cluster splits, Phosformer is **below the label-free Cantley atlas on every endpoint**.

**Fidelity check: passed.** Our MLM checkpoint and fine-tuning code, trained on Phosformer's own published train split (60% uniform subsample, one epoch), give AUROC **0.888** on their test split. The paper reports 0.885 and the released weights give 0.8855 through our scorer. Per-kinase mean AUROC is 0.858 against the paper's 0.860 (`results/baselines/phosformer/fidelity_check.json`). B3a's shortfall on our benchmark therefore is consistent with the split, labels and negatives rather than the implementation; the fidelity check removes the obvious implementation errors but cannot exclude all of them.

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
| R2 full (seed 0, pre-v3) | 0.773 | 0.768 | 0.207 | 0.312 | 0.770 | 0.797 / 0.644 | 0.544 / 0.998 | LEADERBOARD #15 |
| R2 clean (no benchmark PMIDs, no CPTAC) | 0.771 | — | 0.189 | 0.300 | 0.768 | — | 0.551 / 0.999 | `model_r2.md` §6.2 |
| R2 no-phospho (no F6) | 0.771 | 0.767 | 0.187 | 0.296 | 0.761 | — | 0.549 / 0.999 | REGISTRY |
| R3 final (distilled) | 0.754 | 0.749 | 0.223 | 0.294 | 0.731 | 0.761 / 0.586 | 0.524 / 0.929 | LEADERBOARD #42 |
| B7 simple-strong | 0.726 | 0.731 | 0.196 | 0.179 | 0.708 | 0.744 / 0.540 | — | LEADERBOARD #53 |
| B1 Cantley (i) / (ii) | 0.717 | 0.721 | 0.157 / 0.171 | 0.136 (full-cov, 196 kinases) | 0.690 | 0.732 / 0.534 | — | LEADERBOARD #60/#63 |
| B4 SELPHI2.0 (matched) | 0.708 | — | 0.176 | 0.316 | 0.714 | — | not run | `baseline_selphi2.md` |
| B3a Phosformer retrained | 0.677 | 0.689 | 0.132 | 0.104 | 0.655 | — | not run | LEADERBOARD.md (A24 re-eval) |
| B2 learned PSSM | 0.660 | 0.673 | 0.118 | 0.109 | 0.638 | 0.659 / 0.529 | — | LEADERBOARD #74/#75 |
| B0-pop / B0-miss | 0.552 / 0.539 | | 0.020 / 0.040 | 0.142 / 0.053 | | | — | LEADERBOARD |

`results/cv/LEADERBOARD.md` was regenerated at 2026-10-27 01:23 (135 runs) and includes the runs cited here (it does not list every historical run).

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

**Cost–accuracy.** The final model is the most accurate and cheaper to score than B3 and B4. It infers about 20× faster than B3 and about 6× faster than B4, measured on the same candidate set.

---

## 7. TEST results (single batch, FREEZE 721dc197; attempt 1; COMPLETED 2026-10-27 00:57:29)

**About the batch.**
- One sealed access. The batch took 2,838 s for 12 configurations at B = 10,000, with harness 974d520 (`results/test/LEDGER.md`).
- Outputs are in `results/test/721dc1977366/`.
- **Independent review:** `docs/audit/audit_test.md`.
  - Every endpoint, stratum and reliability table reproduces exactly (max |diff| = 0).
  - All 17 prediction pins match FREEZE.
  - TEST and DEV share 0 proteins and 0 clusters.
  - The final score columns are an exact label-free function of the pinned R2 and R3 files.
- **Wording rules.** The review's verdicts and caveats bind the wording below.
- **Audit-diagnostic values.** Values marked "(audit diagnostic)" are the review's re-cuts of batch quantities: kinase-resampling, matched subsets, corrected truth sets. They are not TEST endpoints.

### 7.1 T1: seen kinase, unseen protein (K_eval = 119, K_any = 256)

Standalone metrics use missing = 'worst' and are not head-to-heads; see §7.2 for those. CIs are cluster bootstrap, B = 1,000.

| Config | M1 [95% CI] | M2 Hit@5 [CI] | M3 [CI] | M4 | Kinases scored | Train s | Infer s |
|---|---|---|---|---|---|---|---|
| **final (R2+R3)** | **0.794** [0.774, 0.811] | **0.363** [0.324, 0.397] (score_site) | **0.363** [0.309, 0.401] | **0.799** | 503 | 1,136.6 | 316.9 |
| R2 alone | 0.778 [0.758, 0.795] | 0.264 [0.224, 0.293] | 0.312 [0.261, 0.348] | 0.781 | 503 | 417.9 | 298.9 |
| R3 alone | 0.775 [0.752, 0.793] | 0.322 [0.290, 0.362] | 0.292 [0.242, 0.330] | 0.769 | 503 | 718.7 | 3.1 |
| B7 simple-strong | 0.725 [0.700, 0.751] | 0.243 | 0.200 | 0.717 | 503 | 89.1 | 287.0 |
| B1 (i) Cantley as published | 0.714 [0.685, 0.743] | 0.163 | 0.161 [0.129, 0.194] | 0.699 | 386 | 0 | 4.4 |
| B1 (ii) universe percentile | 0.714 | 0.189 [0.155, 0.214] | 0.161 | 0.699 | 386 | 0 | 1.0 |
| B3 Phosformer (our retrain) | 0.708 [0.686, 0.733] | 0.163 [0.145, 0.190] | 0.104 [0.077, 0.129] | 0.688 | 503 | 3,152.4 | 6,478.8 |
| B4 SELPHI2.0 (our retrain; 66.6% of pairs imputed) | 0.705 [0.673, 0.733] | 0.233 [0.204, 0.266] | 0.304 [0.233, 0.351] | 0.701 | 503 | 737.0 | 2,043.4 |
| B2 learned PSSM | 0.685 [0.647, 0.717] | 0.155 | 0.138 | 0.664 | 503 | 3.3 | 39.7 |
| B0-pop (site popularity) | 0.578 [0.543, 0.608] | 0.020 | 0.127 | 0.588 | 503 | 0 | 0.4 |
| B0-miss | 0.565 | 0.037 | 0.066 | 0.579 | 503 | 19.0 | 30.4 |
| B0-prior | 0.500 | 0.109 | 0.010 | 0.500 | 503 | 0 | 0.8 |

Sources: `results/test/LEDGER.md` and `results/test/721dc1977366/<config>_ci.tsv`.

**Coverage and runtime.**
- B1 leaves 1,199,434 of 13,230,299 required pairs unscored (117 kinases).
- Runtimes are the FREEZE manifest values; B3's training time excludes its about 4.3 h of MLM pretraining.

**Kinase-resampling CIs for final** (`final_ci_kinases.tsv`): M1 [0.767, 0.821], M2 [0.310, 0.414], M3 [0.318, 0.407], M4 [0.775, 0.820].

### 7.2 The 9 primary tests (paired cluster bootstrap B = 10,000, matched universe, Holm over 9)

| # | Comparison | Endpoint | Δ [95% CI] | p (bootstrap) | Holm p | Beats | Kinase-resampling CI (audit diagnostic) | Kinases favouring A / B (audit) |
|---|---|---|---|---|---|---|---|---|
| 1 | final − B1 (i) | M1 | +0.075 [0.046, 0.099] | < 0.0002 (bootstrap floor) | 0.0018 | yes | [0.057, 0.094] | 97 / 19 |
| 2 | final − B1 (i) | M3 | +0.205 [0.165, 0.235] | < 0.0002 (floor) | 0.0018 | yes | [0.167, 0.243] | 78 / 4 |
| 3 | final − B1 (ii) | M2 | +0.179 [0.147, 0.213] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 4 | final − B3 | M1 | +0.086 [0.059, 0.105] | < 0.0002 (floor) | 0.0018 | yes | [0.065, 0.110] | 97 / 22 |
| 5 | final − B3 | M2 | +0.199 [0.151, 0.229] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 6 | final − B3 | M3 | +0.259 [0.213, 0.291] | < 0.0002 (floor) | 0.0018 | yes | — | — |
| 7 | final − B4 | M1 | +0.089 [0.063, 0.117] | < 0.0002 (floor) | 0.0018 | yes | [0.069, 0.111] | 98 / 20 |
| 8 | final − B4 | M2 | +0.130 [0.091, 0.157] | < 0.0002 (floor) | 0.0018 | yes | [0.089, 0.172] | 70 / 16 |
| 9 | final − B4 | M3 | **+0.059 [0.010, 0.120], modest** | 0.0178 | 0.0178 | yes | [0.021, 0.097] | 52 / 27 (40 ties) |

Sources: `results/test/721dc1977366/superiority.tsv`; `docs/audit/audit_test.md` §3. The audit re-ran tests 1, 2 and 7–9 at B = 10,000 and reproduced the batch to the last digit.

**Matched universe.**
- Against B1: 386 of 503 kinases, 74.1% of pairs, and 116 K_eval kinases for M1/M3. TTK, MAP2K1 and LIMK2 are not scored by B1.
- M2 against B1 is matched-kinase (A27); there the final model's M2 is 0.372.
- B3 and B4 are exhaustive (503 kinases).

**Required wording.**
- B1 is Cantley arms (i) and (ii). The priming-aware arm was not a TEST comparator.
- B3 and B4 are "our faithful retrained reimplementations under this benchmark's split, labels and negatives".
  - B3's pipeline reproduces the published AUROC on its own split (0.888 vs 0.885; `results/baselines/phosformer/fidelity_check.json`).
  - 66.6% of B4's TEST pairs use NA-default site features. The fold-0 CV check puts the imputation cost at about 0.013 M1, which does not explain the gap.
- The B4 M3 win (#9) is **modest**: it is the only test off the p floor, and its CI lower bound is 0.010.

### 7.3 Secondary comparisons (Holm over 17, separately; `results/test/721dc1977366/secondary.tsv`)

| Comparison | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| final − R2 | +0.016 [0.009, 0.022] | +0.099 [0.073, 0.128] | +0.051 [0.034, 0.071] | +0.017 [0.007, 0.025] |
| final − R3 | +0.019 [0.008, 0.030] | +0.041 [0.007, 0.060] (Holm p 0.013; weakest) | +0.071 [0.041, 0.097] | +0.030 [0.017, 0.040] |
| final − B7 | +0.069 [0.044, 0.088] | +0.120 [0.083, 0.146] | +0.163 [0.123, 0.198] | +0.082 [0.051, 0.098] |
| final − B2 | +0.109 [0.076, 0.140] | +0.208 [0.168, 0.239] | +0.225 [0.171, 0.263] | +0.135 [0.098, 0.157] |
| final − B1 (i) | — | — | — | +0.087 [0.057, 0.103] |

All 17 meet the "beats" rule; Holm p ≤ 0.0034 except final − R3 M2 (0.013).

**Fusion helps on TEST.** Both components are beaten on every endpoint, which replicates the CV choice of late fusion.

### 7.4 T3: unseen kinase family and unseen protein

Each kinase is scored by the model trained without its LKO family (`results/test/LEDGER.md`; `<config>_t3.json`).

| Config | T3 M1 | T3 M3 | T3 M4 | Missing kinases |
|---|---|---|---|---|
| **final** | **0.752** | 0.285 | 0.753 | 0 |
| R2 | 0.744 | 0.247 | 0.746 | 0 |
| B7 | 0.708 | 0.171 | 0.699 | 0 |
| R3 | 0.703 | 0.209 | 0.703 | 0 |
| B2 | 0.616 | 0.072 | 0.609 | 0 |

T3 M2 is not reported, because it mixes LKO models.

**Reading.**
- **Consistent with CV.** CV LKO-T2 for the fusion was 0.749 (`docs/model_r2.md` §6.3).
- **Relative to Cantley.** The final T3 M1 sits above label-free Cantley's T1 M1 (0.714). On B1-scored kinases (audit diagnostic) it is 0.755 vs 0.719, with 79 of 116 kinases ahead. There is **no paired T3 test**, so this is descriptive.
- **Y kinases.** T3 on Y kinases is 0.531 (23 kinases; `final_t3_per_kinase.tsv`), below the B0-pop floor of 0.585. There is no skill there.

### 7.5 TEST strata (M1 macro E-hard over K_eval ∩ stratum; standalone; no CIs in the batch, so descriptive)

| Stratum | Kinases / pairs | final | B1 (i) | B3 | B4 | B0-pop (floor) |
|---|---|---|---|---|---|---|
| all | 119 | 0.794 | 0.714 | 0.708 | 0.705 | 0.578 |
| ST kinases | 94 | 0.841 | 0.772 | 0.751 | 0.736 | 0.576 |
| **Y kinases** | 23 | **0.595** | 0.496 | 0.534 | 0.575 | **0.585** |
| STY | 2 | 0.865 | 0.500† | 0.695 | 0.733 | 0.604 |
| bin 5–19 | 13 | 0.759 | 0.646 | 0.610 | 0.694 | 0.639 |
| bin 20–99 | 73 | 0.799 | 0.710 | 0.706 | 0.701 | 0.577 |
| bin ≥100 | 31 | 0.795 | 0.748 | 0.757 | 0.703 | 0.542 |
| M4, bin 0 (K_any) | 16 | 0.807 | 0.602† | 0.508 | 0.676 | 0.677 |
| M4, bin 1–4 (K_any) | 19 | 0.757 | 0.606† | 0.712 | 0.636 | 0.548 |
| priming kinases | 8 | 0.786 | 0.726 | 0.740 | 0.635 | 0.449 |
| primed positives | 6 / 65 | 0.897 | 0.704 | 0.788 | 0.665 | 0.517 |
| unprimed positives (priming kinases) | 8 / 194 | 0.740 | 0.725 (near-tie) | 0.729 | 0.644 | 0.445 |
| PIKK (S/T-Q PIKK) | 4 (3) | 0.868 (0.903) | 0.796 (0.834) | 0.859 | 0.739 | 0.517 |
| FAM20C (M4, K_any) | 1 | 0.999 | 0.969 | 0.855 | 0.766 | 0.371 |
| hub sites (A26) | 85 / 272 | 0.863 | 0.712 | 0.731 | 0.802 | **0.788** |
| non-hub | 117 / 1,932 | 0.784 | 0.721 | 0.707 | 0.682 | 0.544 |
| in vivo only | 88 / 468 | 0.788 | 0.713 | 0.713 | 0.706 | 0.589 |
| in vitro only | 96 / 721 | 0.827 | 0.730 | 0.715 | 0.751 | 0.615 |
| both | 107 / 771 | 0.779 | 0.718 | 0.701 | 0.676 | 0.539 |

Source: `results/test/721dc1977366/{final,B1_i,B3,B4,B0_pop}_strata.tsv`.

† Standalone B1 counts its unscored kinases as AUROC 0.5. This deflates the STY row (TTK and MAP2K1, per A30) and the bin-0 and bin-1–4 rows. Restricted to B1-scored kinases (audit diagnostic, `audit_test.md` §5), the values are M4 bin 0: 0.793 vs 0.704, and M4 bin 1–4: 0.756 vs 0.634.

**M1 by Manning group.**

| Group | Kinases | final | B1 (i) |
|---|---|---|---|
| AGC | 24 | 0.816 | 0.741 |
| CAMK | 15 | 0.829 | 0.774 |
| CMGC | 23 | 0.880 | 0.823 |
| STE | 8 | 0.913 | 0.801 |
| Other | 15 | 0.804 | 0.721 |
| CK1 | 4 | 0.817 | 0.725 |
| Atypical | 4 | 0.868 | 0.796 |
| TKL | 3 | 0.819 | 0.638 |
| TK | 23 | 0.595 | 0.496 |

**M2 (Hit@5) by training bin, final vs B1 (ii).**

| Bin | final | B1 (ii) |
|---|---|---|
| 5–19 | 0.208 | 0.165 |
| 20–99 | 0.303 | 0.190 |
| ≥100 | 0.568 | 0.177 |

These are selected subsets: the 13 TEST kinases with 5–19 labels are the best-labelled of that bin. The CV M2 weakness for understudied kinases is not refuted (audit §2).

**Where the final model does not clearly beat Cantley, or the comparison is uninformative:**
- Tyr kinases: both sit at or near the popularity floor.
- Unprimed positives of priming kinases: a near-tie.
- STE M2: a near-tie (0.255 vs 0.196, arm ii).
- TKL M2: 3 kinases, uninformative.
- Primed sites against the *priming-aware* Cantley arm: not tested.

Hub-site gains are mostly popularity (B0-pop 0.788 on hub positives).

### 7.6 Temporal and external curated sets

| Set | K_eval kinases / pairs | final | R2 | R3 | B1 (i) | B3 | B4 | B0-pop |
|---|---|---|---|---|---|---|---|---|
| new since 2021 | 61 / 215 | **0.766** | — | — | 0.710 | 0.660 | 0.696 | 0.556 |
| pre-2021 | 118 / 1,745 | 0.800 | — | — | 0.719 | 0.721 | 0.711 | 0.580 |
| external tier A | 41 / 134 | **0.788** | 0.771 | 0.766 | 0.684 | 0.721 | 0.679 | 0.531 |
| external tier B (popularity-heavy) | 19 / 49 | 0.804 | 0.752 | 0.842 | 0.661 | 0.778 | 0.754 | **0.724** |
| external, all (M1 / M4) | 49 / 130 (K_any 75 / 183) | 0.800 / 0.800 | | | 0.699 / 0.691 | 0.734 / 0.710 | 0.696 / 0.726 | 0.577 / 0.589 |

Sources: `<config>_strata.tsv` (temporal, external all) and `<config>_ext_tiers.tsv`.

- **Discoveries since 2021:** the margin over Cantley is smaller (+0.056) but still positive. M3 is 0.351 against 0.178.
- **Count drift:** the batch evaluates 262 temporal pairs and 134 tier-A pairs. The earlier draft quoted 263 and 158, from the split report's pre-filter counts (audit §7 #9).

### 7.7 A17 assignment reliability on labelled TEST sites (final model only)

Overall: n = 1,458 sites, top-1 accuracy 0.353, ECE 0.026 (`final_reliability.tsv`).

| p_top tier | n | claimed | observed |
|---|---|---|---|
| [0, 0.2) | 261 | 0.12 | 0.14 |
| [0.2, 0.5) | 852 | 0.34 | 0.32 |
| [0.5, 0.6) | 236 | 0.563 | 0.564 |
| [0.6, 0.7) | 37 | 0.679 | 0.676 |
| [0.7, 0.8) | 72 | 0.755 | 0.694 |
| ≥ 0.8 | 0 | — | — |

By residue (audit diagnostic, `audit_test.md` §6):
- **ST:** ECE 0.021; the ≥0.5 tier is 0.615 claimed vs 0.603 observed (n = 345; Wilson [0.55, 0.65]).
- **Y:** ECE 0.069; no site can exceed 0.445; the [0.2, 0.5) tier is 0.41 claimed vs 0.34 observed (n = 195), mildly overconfident.

**Coverage over TEST candidate sites.** 5,773 at 0.5 (5,428 unlabelled), all ST; 0 at 0.8.

The A31(b) isotonic map is honest for S/T and replicates CV (0.62 / 0.61). There is no ≥0.8 tier. The ≥0.5 tier means about 60% top-1 precision on *labelled* sites; under PU selection, that is not a precision on unlabelled sites.

**Temporal and external ECE.** The batch reports ECE 0.122 (new since 2021) and 0.218 (external), with ≥0.5-tier accuracy of 0.36 and 0.21. **These are not miscalibration.** They are mostly a truth-set artefact:
- The batch scores each subset only against that subset's new or external kinases at the site, so a top call that is the site's older PSP kinase counts as wrong.
- With each site's full truth set (audit diagnostic), the ≥0.5 tier is 0.55 observed vs 0.63 claimed (n = 53, new since 2021) and 0.54 vs 0.61 (n = 28, external). ECE becomes 0.062 and 0.068.
- That is mild overconfidence at most, within the Wilson intervals (`audit_test.md` §6).

### 7.8 CV→TEST gap (leakage alarm)

**M1 on common kinases** (K_eval ∩ K_cv; `docs/audit/audit_test.md` §4):

| Config | CV (K_cv) | CV on common kinases | TEST on the same kinases | paired Δ ± SE |
|---|---|---|---|---|
| **final** | 0.786 | 0.792 | 0.793 | **+0.001 ± 0.012** |
| R2 | 0.776 | 0.786 | 0.778 | −0.007 ± 0.013 |
| R3 | 0.754 | 0.763 | 0.773 | +0.010 ± 0.012 |
| B1 (i) | 0.717 | 0.737 | 0.713 | −0.024 ± 0.013 |
| B3 | 0.677 | 0.715 | 0.709 | −0.007 ± 0.012 |

- **No optimism in the final model.** This matches the audited bound of ≲0.005 over about 90 CV configurations. M3 is equally stable: CV 0.358 on common kinases vs TEST 0.361.
- **Residue split, final model on the same kinases:**

  | | CV | TEST | paired Δ ± SE |
  |---|---|---|---|
  | ST (92 kinases) | 0.827 | 0.841 | +0.014 ± 0.012 |
  | **Y (23 kinases)** | 0.667 | 0.595 | **−0.072 ± 0.033** |

  SYK (0.85 → 0.56), IGF1R (0.82 → 0.38) and LYN (0.68 → 0.49) collapse; most Y kinases have only 3–4 TEST positives.

**M2 composition.** "On the same 117 kinases, CV M2 = 0.348 vs TEST 0.363 (paired Δ +0.014 ± 0.021). No leak signal."

| | final | R2 | R3 | B1 (ii) | B0-prior |
|---|---|---|---|---|---|
| CV M2 over K_cv (203) | 0.263 | 0.200 | 0.223 | 0.171 | 0.064 |
| CV M2 over K_eval ∩ K_cv (117) | **0.348** | 0.246 | 0.320 | 0.180 | 0.111 |
| TEST M2 over K_eval | **0.363** | 0.264 | 0.322 | 0.189 | 0.109 |

Source: `audit_test.md` §2.

- Of the +0.100 CV→TEST change, about +0.085 is kinase selection. K_eval holds 13 / 73 / 31 kinases in the 5–19 / 20–99 / ≥100 bins, against 83 / 89 / 31 for K_cv, and Hit@5 rises steeply with label count.
- The label-free B0-prior control moves the same way (0.064 → 0.109; CV on the same kinases 0.111).
- **TEST M2 0.363 must not be presented as an improvement over CV.**

**Margin over Cantley, CV vs TEST.**

| | Δ M1, final − B1 |
|---|---|
| CV, matched (`results/models/r2/compare__r2_fusion_r3__vs__simple_b1_universe.tsv`) | +0.063 [0.048, 0.074] |
| TEST, matched (`superiority.tsv`) | +0.075 [0.046, 0.099] |

- Part of the larger TEST margin is a weak-baseline draw. On the 115 kinases B1 scores in both, Cantley fell from 0.741 to 0.717 (−0.025 ± 0.013), almost entirely on Y kinases (0.607 → 0.496). The final model was flat (−0.003) (`audit_test.md` §3).

### 7.9 Not computed in the batch (cannot be back-filled without a second sealed access, A29)

Source: `docs/audit/audit_test.md` §5, §6, §7 #10. The review has no §8; this list collects its "not delivered" items.

1. **Y kinases with hub sites excluded.** The planned §7.5 row; only the hub/non-hub split of positives exists.
2. **A16 calibration on the external and temporal strata** (A29).
3. **A17 coverage baselines on TEST** (random scores and Cantley-only). Only the final model supplied `assignment_prob`. The CV versions are in `docs/model_r2.md` §6.1; they predate A31, so they are not like-for-like.
4. **A paired T3 test** against any baseline. The T3 comparisons above are descriptive.
5. **CIs for strata.** Every §7.5–7.6 stratum comparison is descriptive.
6. **The priming-aware Cantley arm** as a TEST comparator. It was not in FREEZE, so no TEST claim can be made against it.

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

## 9. Coverage: deployment model over U_obs (`results/release/coverage.tsv`, `RELEASE_MANIFEST.json`, `docs/model_r2.md` §7.1)

**What the deployment model is.**
- It uses the frozen final configuration (R2 5 seeds + R3 5 seeds, fused with the same weights), retrained on DEV+TEST labels (15,202 training pairs, L2) through `ksp.deploy.load_deploy_labels`. That loader refuses until the TEST ledger is COMPLETED.
- It scores all 350,736 U_obs sites × the 503 kinases of matching class: 131,144,410 pairs, in `cache/release/…/release_scores.parquet`, sha256 `cceb32c7…`, pinned in the manifest.
- Calibrators, the temperature and the ST/Y isotonic maps are the ones fitted on DEV OOF. No OOF predictions exist for DEV+TEST.

| Method | Coverage@0.5 (no-PSP sites) | Coverage@0.8 (no-PSP) | Sites in some kinase's top 1% (no-PSP) |
|---|---|---|---|
| **Final (deployment)** | **29,886** (26,762) | 0 (0) | 162,783 (152,377) |
| Cantley-only (B1 imputed-universe percentile, its own DEV-OOF calibration recipe) | 30,258 (29,220) | 1,828 (1,780) | 170,065 (162,611) |
| Random scores (same softmax, no isotonic map) | 6 (6) | 0 (0) | 329,349 (318,032) |

U_obs has 350,736 sites; 338,525 of them have no PSP kinase.

**Reading.**
- **Coverage@τ counts are not precision.** The denominator is the 350,736-site U_obs, so 29,886 sites at 0.5 is 8.5% of it, and 0 at 0.8. Final and Cantley-only reach about 30k sites at 0.5. The TEST reliability check (§7.7) applies to the frozen DEV-trained model; the DEV+TEST deployment model's calibrators were not re-tested. On labelled sites the final model's ≥0.5 tier is calibrated on TEST (S/T: 0.615 claimed vs 0.603 observed, §7.7), and its top-1 accuracy is 0.365 against 0.282 for Cantley (CV, §4.5).
- **Cantley's ≥0.8 tier is not evidence of higher precision.** It comes from the tail of an isotonic map fitted on very few labelled sites: only 152 Cantley sites reached ≥0.5 on CV.
- **"Top-1%" coverage is uninformative as a precision measure.** Random scores spread their top 1% over the most sites.
- **The model has no ≥0.8 tier and no Y-site assignment ≥0.5** (A31b). The "which kinase?" output should be used as a ranked shortlist; M2 Hit@5 on TEST is 0.363.
- **Coverage on unlabelled sites cannot be precision-verified.** All coverage numbers are model-based and conditional on the site being a substrate of a modelled kinase.

**Spot checks.**
- TP53 S15 → ATM, assignment probability 0.70 (correct).
- AKT1 S473 → PDPK1, AKT1, then MTOR third (MTOR/mTORC2 is the known kinase, ranked 3).
- EIF2S1 S52 → PKC family: a miss; the known kinases are EIF2AK1–4.

**Release files**, PSP-derived, so CC BY-NC-SA and non-commercial:
- `results/release/kinase_top200.tsv.gz`: 100,600 rows, the top-200 sites for each of the 503 kinases.
- `results/release/site_top5.tsv.gz`: 1,753,680 rows, the top-5 kinases for each of the 350,736 sites, with assignment_prob.
- `RELEASE_MANIFEST.json`: SHA-256 pins, licence and limitations.

**Runtime.**

| Step | Wall-clock |
|---|---|
| Deployment cross-fit stores | 208–322 s, in parallel |
| R2 training | 310 s |
| R2 inference (131 M rows) | 985 s |
| R3 training | 447 s |
| R3 inference | 12 s |
| Fusion, calibration, assignment and files | 169 s |

Peak RSS was 18.7 GB.

**Note.** In the deployment run, seed 0's early stopping hit its 53-round floor, against 235–543 rounds for seeds 1–4. That is the frozen procedure, and the 5-seed mean is dominated by the other seeds.

---

## 10. Limitations

1. **No demonstrated Tyr-kinase specificity.**
   - On TEST, Y kinases score M1 0.595 against the site-popularity floor B0-pop at 0.585. The final model beats that floor on only 10 of 23 Y kinases.
   - The review rates this **NOT TRUSTED** as a skill claim. Cantley is at chance (0.496).
   - T3 on Y kinases is 0.531.
   - Y-kinase M1 fell from CV to TEST: 0.667 → 0.595 (paired Δ −0.072 ± 0.033). SYK, IGF1R and LYN collapse.
   - No Y site can receive an assignment probability above 0.445.
   - On CV the Y gap was already large (0.63 vs 0.81). F9 did not help. Part of the gap is endpoint composition (hubs, indirect in-vivo labels).
   - Sources: `docs/audit/audit_test.md` §0 #10, §4; `results/test/721dc1977366/final_t3_per_kinase.tsv`; `docs/strata_cv.md`.
2. **No high-confidence assignment tier.**
   - After A31 there is no ≥0.8 tier. TEST coverage@0.8 is 0, and coverage@0.5 is 5,773 sites, all S/T.
   - The S/T ≥0.5 tier is honest on labelled TEST sites (0.615 claimed vs 0.603 observed, n = 345).
   - Y probabilities are mildly overconfident (0.41 vs 0.34 in [0.2, 0.5)).
   - These are precisions on labelled sites only (PU); they do not transfer to unlabelled sites.
   - Sources: `final_reliability.tsv`, `summary.json` (coverage counts); `audit_test.md` §6.
3. **M2 favours well-studied kinases.**
   - CV Hit@5 is 0.142 at 5–19 labels against Cantley's 0.162. A32(b) debiasing failed.
   - TEST cannot refute this: it holds only 13 such kinases, the best-labelled of the bin.
   - TEST M2 (0.363) is K_eval composition, not improvement; CV on the same 117 kinases is 0.348 (`audit_test.md` §2).
4. **The margin over Cantley is partly a weak-baseline draw.**
   - The TEST M1 margin (+0.075 [0.046, 0.099]) exceeds the CV margin (+0.063 [0.048, 0.074]).
   - That is because Cantley fell by 0.025 ± 0.013 on the same kinases, mostly on Y kinases. Quote both margins (`audit_test.md` §3).
5. **Priming.**
   - The TEST primed-site advantage (0.897 vs 0.704; 6 kinases, 65 pairs) is only against non-priming-aware Cantley.
   - On unprimed positives the final model and Cantley are near-tied (0.740 vs 0.725).
   - On CV, priming-aware Cantley beat the model on primed sites (0.945 vs 0.893), and it was not a TEST comparator.
6. **Baselines are our reimplementations.**
   - B3: the fidelity check passed (0.888 vs 0.885), but it has 1 seed and a PI-fixed 1:16 negative ratio.
   - B4: 10 of 100 bags, a substituted PSSM/coreg, and 66.6% of TEST pairs on NA-default site features.
   - The B4 M3 win is modest: Δ +0.059 [0.010, 0.120], Holm p 0.018.
   - The superiority claims are against "our faithful retrained reimplementations under this benchmark", not against the published systems as distributed.
7. **TEST strata carry no CIs, so they are descriptive.**
   - Small strata are anecdotal: PIKK (4 kinases), FAM20C (1 kinase), STY (2 in K_eval), TKL (3) and the 5–19 bin (13).
   - Hub-site and external-tier-B gains are mostly popularity (B0-pop 0.788 and 0.724).
8. **Items the single batch could not provide, and which cannot be added without a second sealed access:**
   - hub-excluded Y;
   - A16 calibration on the temporal and external strata;
   - TEST coverage baselines;
   - a paired T3 test;
   - stratum CIs.

   See §7.9.
9. **Phosphoproteomics adds almost nothing to CV M1** (F6 ≈ +0.002; the no-phospho variant scores 0.771). Downstream, these features may help on unseen drugs, but the evidence is indirect.
10. **Downstream parity, not superiority, with PSP**, plus an unseen-drug gap of 0.04–0.06 (`docs/downstream.md` §7).
11. **The canary gate passes only pooled.**
    - Two of five CV folds exceed 0.57 (0.575 and 0.676).
    - Every model sits above 0.5, including label-free controls. The sampler explanation is plausible but unproven.
12. **LKO generalisation is weaker on novel sites.**
    - On CV T2, sites with no other-kinase label score 0.703 against a 0.757 control.
    - T2 and T3 are "unseen labels", not kinase-naive: label-free Cantley percentiles enter R3 pretraining.
    - TEST T3 (0.752) has no paired test against a baseline.
13. **Optimism and post-freeze changes.**
    - About 90 configurations were tried on CV. TEST confirms no M1 optimism: paired Δ on common kinases +0.001 ± 0.012.
    - A24 was a post-freeze metric correction applied to all methods; it widened learned-vs-Cantley CV margins by about 0.01.
14. **Label noise and PU.** In-vivo-only labels are enriched for indirect events. All precision and probability statements are lower bounds.
15. **Unbuilt items.**
    - F8 (Sugiyama), docking motifs, F6d and U_all (unobserved S/T/Y).
    - The deployment model is done (§9).
    - `scripts/run_all.sh`, promised in DESIGN §11.8, does not exist.
16. **Data gaps** (DATAINEED): the SIGNOR direct export, decryptM curve fits, CPTAC plex/ischemia confounders and PTB matrices were unavailable; workarounds are in use.

---

## 11. Follow-ups and remaining ideas

The TEST batch, FREEZE and ledger are done (`results/test/LEDGER.md`, 2026-10-27 00:57), and so is the Phosformer fidelity check (`results/baselines/phosformer/fidelity_check.json`); both are removed from this list.

**Any further evaluation must use new data, a new split, or a logged `--override-completed` (A29).** The sealed TEST set has been spent.

| Priority | Idea | Rationale / evidence |
|---|---|---|
| — | ~~Deployment model (DEV+TEST) plus U_obs coverage~~ | Done (G5, §9, `results/release/`) |
| 1 | Tyr specificity: nearest-atlas imputation for anti-correlated matrices (FGFR, Eph, VEGFR, TAM, BTK); a pre-registered A-loop/hub sensitivity analysis; PTB-adaptor recruitment | TEST Y 0.595 vs popularity floor 0.585; T3 Y 0.531 (`audit_test.md` §4); `features_tyr.md` §1.5–1.6. **Post-TEST exploration done (A34, `docs/posttest_tyr.md`): both cheap score-level fixes fail the A34 bar on CV** — Y-specific fusion weight −0.0039 [−0.0097, +0.0010]; nested nearest-atlas motif repair +0.0042 [−0.0103, +0.0167]. Neither adopted. Remaining value is in the endpoint (hub/A-loop stratum), not more CV search |
| — | ~~Kinase-gated priming blend~~ | Done post-TEST (§11.1 (b)): a rank shift reproduced by a label-free priming-site flag; +0.001 on all positives. A real test needs flag-matched negatives |
| 2 | M2 for understudied kinases: a blend toward Cantley below 20 labels (not yet tried) | CV 5–19 bin 0.142 vs Cantley 0.162; TEST has only 13 such kinases (`audit_test.md` §2). The intercept route is done (§11.1 (a)): it trades M2 between kinases |
| — | ~~Per-kinase intercepts for the site-centric score~~ | Done post-TEST (§11.1 (a)): the macro-M2 gain needs a post-hoc balanced variant and costs ≥100-label kinases −0.124 |
| 2 | Y-site assignment calibration: Y probabilities are capped at 0.445 and mildly overconfident | `audit_test.md` §6 |
| 2 | True leave-drug-out retrain (F6 with `exclude_drugs`) | separates the drug effect from the Hijazi-dataset effect in the downstream gap |
| 2 | Prospective or external evaluation for the analyses the batch could not provide (hub-excluded Y, temporal/external calibration, paired T3, stratum CIs), for example on a future PSP release | §7.9; TEST cannot be reopened |
| 3 | Site-uniform v3 negative canary | turns the pooled-only gate into a clean null |
| 3 | R5 revisited: LoRA or last-k-layer fine-tuning, a lower trunk learning rate, protein-grouped early stopping | the trunk overfit by step 100–600 |
| 3 | decryptM potency coherence (EC50 vs Kd); Bayer–Kuster 2025 dose–response data | closer to direct-substrate evidence (DATAINEED 4, 8) |
| 3 | Temporal retrain on 2021-only labels; L3 variant; KolossuS / AF-Multimer pilot | deferred in DESIGN §0 |
| 3 | Seeds 1–4 for B3; 100 bags for B4 | baseline strength (only relevant for a future benchmark) |

---


### 11.1 Post-TEST exploratory results (A34)

These ran after the TEST batch was spent, on **DEV CV only**, with every tunable chosen inside each outer fold (nested). The released model is unchanged. Bootstrap CIs resample substrate clusters (B = 1,000). An independent audit re-ran all three scripts and reproduced every number exactly (`docs/audit/audit_posttest.md`). One caveat is shared by all three: they are second-stage fits on OOF predictions whose base models saw the held-out fold's *training* partner folds (ordinary stacking). Write-ups: `docs/posttest_cv.md` for (a) and (b), `docs/posttest_tyr.md` for (c).

**None of the three is a clean improvement, and none is adopted.**

| Target | Intervention (nested) | Target stratum Δ [95 % CI] | Cost elsewhere | Audit verdict |
|---|---|---|---|---|
| (a) M2, understudied kinases | per-kinase intercepts on `score_site`, **kinase-balanced variant added after seeing the plain variant's first grid point** | M2 5–19-label bin 0.142 → 0.227, +0.085 [0.048, 0.116]; pooled macro M2 0.263 → 0.299 | M2 ≥100-label bin 0.621 → 0.497 (−0.124; 28/31 kinases worse); site-averaged (micro) Hit@5 0.535 → 0.465; unlabelled kinases' share of top-5 slots 14 % → 2.3 %. **Pre-declared plain variant alone: 5–19 bin −0.013 [−0.021, −0.004], pooled −0.006** | confirmed with caveats |
| (b) primed sites | priming-gated blend with priming-aware Cantley (γ 0.2–0.3 per fold) | primed positives 0.893 → 0.911, +0.018 [0.013, 0.024] | unprimed −0.006 [−0.010, −0.001]; all priming-kinase positives +0.001 [−0.003, +0.005]; M3 −0.0017 (CI < 0); overall M1 ±0.000 | confirmed with caveats |
| (c) Tyr specificity | Y-specific fusion weight | Y M1 (38 K_cv kinases) 0.6387 → 0.6347, −0.0039 [−0.0097, +0.0010] | none (ST bit-identical) | confirmed (fails) |
| (c) Tyr specificity | motif repair: own / nearest / cross-family Cantley Tyr row | Y M1 → 0.6429, +0.0042 [−0.0103, +0.0167] | none (ST bit-identical) | confirmed (fails) |

**(a) moves M2 around rather than improving it.** The balanced intercepts take top-5 slots from well-studied kinases and from kinases with no labels at all, which M2 never scores, and give them to the 5–19-label kinases. That raises the kinase-macro M2 and lowers the per-site one. The gain exists only for the variant added after the first result was seen. The pre-declared variant makes the target bin worse. A per-kinase intercept is cheap to apply, so a user who wants rankings balanced across kinases can add one. It would need a pre-registered test on new labels before it could be claimed as an improvement.

**(b) is a rank shift, not new information.** The primed stratum is defined by the same property the priming-aware Cantley arm rewards: a phosphosite at the priming position. A label-free 0/1 indicator for "the priming position is an observed phosphosite", with no Cantley term, reproduces the gain (0.909–0.917). Against negatives that carry the same flag, the gain is +0.002. It would also fail the earlier A32 bar of ≥ +0.02.

**(c) both fail the A34 bar.** γ was selected at the top of its grid (0.5) in all five folds, so "fails" holds for γ ≤ 0.5. The Y-specific weight re-chose w = 0.6 in three folds (Δ = 0 by construction) and was negative in the other two. Motif repair gains +0.012 to +0.027 on the inner folds but only +0.0056 on the held-out folds (3 of 5 positive). With 38 kinases and ~1,150 Y positives, this stratum cannot support selecting even one parameter. The per-family pattern is mixed:
- **Gains:** 18 of 38 kinases, with FGFR +0.060, Abl +0.038 and Jak +0.024.
- **Losses:** 19 kinases. They include two of the families `features_tyr.md` §1.5 predicted to benefit, TAM/Axl −0.053 and Tec −0.031, as well as Ret −0.052, PDGFR −0.038 and Trk −0.035.
- **Untested:** Eph and VEGFR, which have no K_cv kinase.

B0-pop on the Y stratum is 0.592, so the fusion's lead over popularity on Y is marginal on CV as well as on TEST (§7). Tyr remains the main open problem.

## 12. Reproducibility

**Process notes (not elsewhere in this report).**
- **Design review.** DESIGN.md and LITERATURE.md went through two independent review rounds (`docs/reviews/r1_*.md`, `r2_*.md`), within the ≤3 allowed.
- **Cantley is not retrained.** B1 is the published kinase-library atlas used as-is (`ksp` reproduces kinase-library 1.8.0 within 1e-3); its matrices are not fitted on our labels, so it has no fold-dependent training.
- **Time.** Run timings are in `results/runs/*/manifest.json`; `docs/retro/final_retro.md` analyses where the 24 h went.

**Commit and harness.**
- Harness commit: `974d520` (A28), covering `eval.py`, `labels.py`, `splits.py`, `crossfit.py` and `canary.py`.
- Repository HEAD at first drafting: `be635d1` ("tree: all TEST deliverables pinned; report draft started"); the final commit is in `git log`.
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
