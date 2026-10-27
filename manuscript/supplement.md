# S1 Text. Approaches that were tried and did not work, and other negative results

> **AI-agent exercise.** This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Authors: [authors to be added by the publishing researcher]

## Scope and conventions

This text lists the methods, resources and ideas that were tried and not adopted, the evaluation defects that we found and fixed, and the failures in how the agent team ran the project. The main text describes the released model and its results; this text does not repeat them.

- **Metrics.** CV numbers are pooled 5-fold DEV out-of-fold values over K_cv = 203 kinases, from the frozen harness with A24 rank normalisation unless stated. M1 is macro E-hard AUROC, M2 filtered site-centric Hit@5, M3 top-1% recall on E-all. Each Δ is against the matched baseline named in the source. Square brackets are 95% cluster-bootstrap CIs.
- **Sources.** Each entry names the repository file that holds its numbers. `SUPP_NUMBER_MAP.tsv` maps every number to its file.
- **Firmness.** We use four labels.
  - *Firm*: the CI excludes the gain that was sought, or the result replicated.
  - *Not shown to help here*: the CI spans zero, or the test was small, single-seed or specific to this benchmark. The idea may still work elsewhere.
  - *Not tested*: planned or proposed but never run.
  - *Defect*: an error in our own evaluation or process, found and corrected.
- **Rules.** Unless an entry gives another rule, CV changes needed ΔM1 ≥ +0.005, an M2 drop of at most 0.01 and a passing canary gate.

## Summary table

| # | Idea | Outcome | Key evidence (CV unless stated) | Firmness |
|---|---|---|---|---|
| 1.1 | Human-only labels (L1) | not adopted | ΔM2 −0.029 [−0.046, −0.015] | Firm (M2) |
| 1.2 | In-vitro-only / in-vivo-only labels | not adopted | ΔM2 −0.047 / −0.090 | Firm |
| 1.3 | Add external curated labels (L2+ext) | not adopted | ΔM1 +0.0058 on seed 0, mean +0.0013 over 3 seeds | Not shown to help |
| 1.4 | Activation-loop Tyr labels as mask-only | rejected before running | would change labels after seeing results | Not tested |
| 1.5 | STRING, iKiP-DB, PSP-derived scores | excluded as leaky | STRING scores 34% of PSP pairs vs ≈4% of controls; iKiP holds 61% of PSP pairs | Firm |
| 1.6 | Sugiyama 2019 as labels; F8 feature | not used; F8 never built | literature report; no feature store | Not tested |
| 2.1 | knn_sim similarity feature | not adopted | ΔM1 +0.003; ΔM2 −0.018 [−0.030, −0.006] | Firm (M2 harm) |
| 2.2 | F9 SH2-docking / Tyr features | not adopted | ΔM1 +0.0017; sanity gain mostly a NaN artefact | Not shown to help |
| 2.3 | A32(a) priming-offset features | not adopted | primed Δ −0.006 [−0.022, +0.008] vs rule ≥ +0.02 | Firm (vs rule) |
| 2.4 | Phosphoproteomics (F6) | kept, ≈0 CV gain | F6 ≈ +0.002 M1; no-phospho M1 0.771 | Not shown to help |
| 2.5 | ESM-2 3B windows in R3 | not adopted | +0.005 alone, −0.015 combined | Not shown to help |
| 2.6 | PPI random walk, study-bias proxies | weak or excluded | RWR 0.54 against shuffled controls | Firm |
| 3.1 | Per-kinase models | not adopted | ΔM1 −0.0119 [−0.0168, −0.0085] | Firm |
| 3.2 | Per-Manning-group models | not adopted | ΔM1 −0.007, ΔM2 −0.022 | Firm (single seed) |
| 3.3 | Partial pooling | not adopted | λ = 300 selected; ΔM1 +0.0002 | Firm (no gain) |
| 3.4 | Ranking objectives | not adopted | fold-0 M2 0.086 / 0.189 vs 0.223 | Not shown to help |
| 3.5 | R5 fine-tuned ESM late-interaction model | gate failed | fold-0 M1 0.725 vs gate 0.751 | Firm for this set-up |
| 3.6 | A32(b) prior-debiased site score | not adopted | β = 0 chosen in all 5 folds | Firm |
| 3.7 | B2 learned PSSM; B7 simple logistic regression | below or level with Cantley | B2 ΔM1 −0.056; B7 matched ΔM1 +0.005 (CI spans 0) | Firm |
| 3.8 | nnPU, hard-negative mining, KolossuS, AF-Multimer, U_all, docking motifs | not built | — | Not tested |
| 4.1 | Pooled-OOF metric (pre-A24) | defect, fixed | B0-prior M1 0.40 instead of 0.50; 47 runs re-scored | Defect |
| 4.2 | Canary v1 negative control | failed as a null | label-free score 0.76; v2 passes only pooled | Defect |
| 4.3 | Coverage metric as first designed | vacuous | random scores: 98.5% of S/T sites covered | Defect |
| 4.4 | Assignment probability ≥0.8 tier (A31) | overconfident, removed | 0.87 claimed vs 0.54 observed; Y 0/41 | Firm |
| 4.5 | Calibrated assignment probability as the M2 score | replaced | M2 0.205 vs 0.263 for the raw score | Firm (late change) |
| 4.6 | Cantley on single-residue kinases (A30) | left unscored on TEST | 10 kinases; 386 of 503 scored | Defect |
| 4.7 | Cross-fitted stores built from a stale canary table | defect, rebuilt | folds 0 and 3 | Defect |
| 5.1 | Per-kinase intercepts for M2 | reallocation, not a gain | pre-declared variant −0.013 [−0.021, −0.004] | Firm (pre-declared) |
| 5.2 | Priming-gated blend with Cantley | rank shift | all priming-kinase positives +0.001 [−0.003, +0.005] | Firm |
| 5.3 | Y-specific fusion weight | failed | −0.0039 [−0.0097, +0.0010] | Not shown to help |
| 5.4 | Tyr motif repair (nearest atlas) | failed | +0.0042 [−0.0103, +0.0167] | Not shown to help |
| 6 | Agent-team process | see §6 | TEST spent 5.3 h before its gate | Process |

---

## 1. Data, labels and resources

**1.1 Human-only labels (L1).** *Tried:* train on human-native PSP pairs only (9,234 pairs) instead of human plus mouse/rat pairs mapped to human (L2, 11,795 pairs). *Result:* ΔM1 +0.0030 [−0.0019, +0.0097], ΔM2 −0.0294 [−0.0462, −0.0148]. Y-kinase M1 fell from 0.626 to 0.587. *Why:* orthologue labels add positives for the cross-kinase ranking at no measurable M1 cost. *Firmness:* firm for M2; the M1 difference is inside seed noise (R2 seed SD 0.0028). Source: `docs/regimes_labels.md` §1–2.

**1.2 Evidence-restricted labels (L2_vitro, L2_vivo).** *Tried:* train only on pairs with in vitro evidence (8,987 pairs) or only on in-vivo-only pairs (2,808), since the latter are enriched for indirect events. *Result:* L2_vitro ΔM1 −0.0034 [−0.0102, +0.0059], ΔM2 −0.0470 [−0.0627, −0.0283], with a small gain on in-vitro-only test positives (0.7806 vs 0.7720). L2_vivo was worse everywhere (ΔM1 −0.0289, ΔM2 −0.0896) and gained nothing on in-vivo-only test positives. *Why:* L2_vivo is about three times smaller, so it mainly measures label count. *Firmness:* firm that neither beats L2 here. Source: `docs/regimes_labels.md` §1–2.

**1.3 External curated labels (L2+ext).** *Tried:* add 1,177 filtered OmniPath and UniProt pairs on DEV proteins. *Result:* seed 0 met the rule (ΔM1 +0.0058 [−0.0006, +0.0121]; ΔM2 +0.0257 [+0.0124, +0.0351]). The pre-registered replicates did not (seed 1 −0.0040 [−0.0097, +0.0025], seed 2 +0.0022; 3-seed mean +0.0013). M2 was positive in all three seeds (+0.026 / +0.018 / +0.002); Y-kinase M1 was lower in all three. External pairs had also been pre-registered as ablation-only, to keep the external TEST stratum independent (A33). *Firmness:* not shown to help M1; a possible small M2 gain was not pursued. Source: `docs/regimes_labels.md` §2, §4; `docs/DESIGN.md` A33.

**1.4 Activation-loop Tyr labels as mask-only (A26(b)).** The Tyr diagnosis showed cross-kinase labels on substrates' own activation loops, for example AXL on Src-family loop sites. These are plausibly indirect. Treating them as mask-only was proposed and **rejected**, because it would change labels after seeing results. A hub-site reporting stratum was added instead. *Firmness:* not tested. Source: `docs/DESIGN.md` A26; `docs/features_tyr.md` §1.

**1.5 Resources that looked useful but leak** (measured reasons for the exclusions listed in the main text).
- *STRING (all channels, including "experiments"):* 34% of PSP kinase–substrate pairs have a score, against ≈4% of controls. STRING imports BioGRID "Biochemical Activity" records, which are kinase assays.
- *iKiP-DB "combined" GMT and activity tables:* they contain 61% of PSP pairs (A11).
- *COMPARTMENTS knowledge channel:* imports GO/UniProt curation; strongest single context feature (0.74), but only +0.003 M1 as an alarm run (0.776 vs 0.773). Ablation only.
- *PSP-derived site scores in SELPHI2.0:* dropping 4 such features raised M1 from 0.708 to 0.724 and lowered M3 from 0.316 to 0.265. The M3 gain was site popularity.
- *Flags for whether the centre site was observed:* 18% of human-native positives lack MS evidence, against 0.4% of unlabelled sites, so the flags would identify labels (A6).
- *Uncross-fitted label features:* knn_sim fitted on all labels reached AUROC 0.93–1.00, against 0.59–0.76 when cross-fitted.

*Firmness:* firm; each exclusion rests on a measured overlap or control. Sources: `docs/REPORT_internal.md` §3.2; `docs/RESEARCH_TREE.md` 1.4, 3.2; `docs/features_knn.md`.

**1.6 Sugiyama 2019 as labels; F8 as a feature.** The in vitro lysate screen was not used as labels, because a published benchmark reported that it hurt activity inference (Müller-Dott S, et al. Comprehensive evaluation of phosphoproteomic-based kinase activity inference. Nat Commun. 2025;16:4771. doi:10.1038/s41467-025-59779-y; main-text ref 2). The planned F8 feature from the same table was never built (empty group, `F8_unavailable`), although reviewers in round 2 had measured signal in it (about 6× chance) and no PSP curation. *Firmness:* not tested. Sources: `docs/LITERATURE.md` §2; `docs/model_r2.md` §2; `docs/reviews/r2_data_leakage.md`.

## 2. Features

**2.1 knn_sim (similarity to training substrates).** *Tried:* 16 cross-fitted columns measuring ESM and BLOSUM similarity of a site to the kinase's training substrates. In the sanity harness on fold 0 it looked complementary to the learned PSSM (rank average 0.789 vs PSSM 0.765). *Rule:* ΔM1 ≥ 0.005 with a passing canary. *Result in R2:* ΔM1 +0.003 averaged over 5 seeds (seed 0 +0.0054 [0.000, 0.012], seed 1 +0.0004). ΔM3 was +0.012 [0.002, 0.024], but ΔM2 was −0.018 [−0.030, −0.006]. *Why:* the per-kinase density and training-count columns are not comparable across kinases, so they damage the site-centric ranking. The sanity harness used random negatives and no paralog masks, so its numbers are not comparable to the official metric. *Firmness:* firm that it harms M2 in this form. Sources: `docs/model_r2.md` §6; `docs/features_knn.md`.

**2.2 F9: SH2 docking and other Tyr features.** *Tried:* Yaron-Barir SH2 matrices scored on other pY sites of the substrate (24 of 28 SH2 kinases), ITAM motifs and receptor × membrane terms. *Result in R2:* M1 0.7699 vs 0.7682 (+0.0017), below the 0.005 threshold. In the sanity harness, kinases that have an SH2 matrix gained +0.019 [−0.018, +0.057] (p = 0.33, 18 kinases). Kinases whose SH2 columns were all missing gained more (+0.044), so most of the apparent gain came from the missing-value pattern acting as a kinase-class indicator, not from docking. Those sanity CIs resampled kinases, not clusters, and are too narrow. *Firmness:* not shown to help. Sources: `docs/features_tyr.md` §3; `docs/model_r2.md` §6.

The same diagnosis rejected two explanations of the Tyr gap: paralog masking does not leave Y kinases short of negatives (median 779 E-hard negatives per Y kinase), and evidence class does not explain it (0.603 vs 0.596 within 17 Tyr kinases). Excluding hub sites left R2's Y M1 unchanged (0.627 vs 0.629). Sources: `docs/features_tyr.md` §1; `docs/strata_cv.md` §4.

**2.3 A32(a): priming-offset features.** *Motivation:* on primed positives, the priming-aware Cantley arm scored 0.945 against 0.870 for an R2 run. *Tried:* three label-free columns: an indicator that the kinase's preferred phospho-offset position is an observed phosphosite, a preference-weighted flank score, and the kinase's maximum pS/pT preference. *Pre-declared rule:* overall ΔM1 ≥ −0.002, ΔM2 ≥ −0.005, primed-positive ΔM1 ≥ +0.02 with CI > 0, and a passing canary. *Result:* ΔM1 +0.0037 [−0.001, 0.009], ΔM2 +0.010, canary 0.566 (narrow pass), primed ΔM1 −0.006 [−0.022, +0.008] (0.893 → 0.888; 8 kinases, 231 positives). *Why:* F1 already contains priming-aware Cantley percentiles, which ranked 11–53 by gain; the new columns ranked 86–167 of 175. The final R2 configuration already scored 0.893, so the gap was smaller than first thought. *Firmness:* firm against the rule. Source: `docs/strata_cv.md` §9.

**2.4 Phosphoproteomics features (F6).** The brief emphasised this modality (co-regulation, inhibitor-response deconvolution, CPTAC kinase state). *Result:* in the ablation curve F6a–c move M1 from 0.753 to 0.754, about +0.002 in total; in the final configuration, a variant without F6 scores 0.771 against 0.768 with F6 (seed 0; five-seed mean with F6 0.771). The features were kept because they passed the leakage checks. *Why (likely):* coverage is sparse and biased, since labelled sites are covered 2.1–2.5× more often than unlabelled ones, and missingness alone would give per-kinase AUROC of about 0.58–0.63. Inhibitor response is not direct phosphorylation, and some drug potencies were entered from memory (flagged as an assumption). Downstream, removing F6 widened the gap to PSP on unseen drugs (−0.08 at top100), so F6 may help there; that evidence is indirect. *Firmness:* not shown to help CV M1. Sources: `docs/REPORT_internal.md` §6.2, §8; `docs/features_phospho.md` §5.

**2.5 ESM-2 3B window embeddings in R3.** *Rule:* fold-0 ΔM1 ≥ +0.01. *Result:* 0.700 vs 0.695 alone (+0.005), and 0.725 vs 0.741 when combined with Cantley distillation (−0.015). *Firmness:* not shown to help. Source: `docs/model_r3.md` §2.

**2.6 Other weak feature ideas.** PPI random walk with restart scored AUROC 0.67 against random negatives but 0.54 against shuffled controls, so its signal is mostly how well studied a protein is. The count of localisation sources alone scored 0.656; such study-bias proxies were excluded (A11). A phospho channel did not help the learned PSSM (M1 0.654 vs 0.660). Sources: `docs/RESEARCH_TREE.md` 4.6; `docs/features_context.md` §6; `docs/REPORT_internal.md` §5.1.

## 3. Models and training regimes

The main text reports that one global model is best; the alternatives are below (`docs/regimes_labels.md` §3; seed 0).

**3.1 Per-kinase models.** *Tried:* separate LightGBM or L2-logistic models for kinases with at least 20 training positives, with the global model for the rest. *Result:* ΔM1 −0.0119 [−0.0168, −0.0085], ΔM2 −0.0188, ΔM3 −0.0390. The models lost most on the 78 kinases that could be fitted in every fold (global 0.8046 vs per-kinase 0.7788). *Why:* sharing across kinases is worth about +0.025 M1 even for data-rich kinases. *Firmness:* firm.

**3.2 Per-Manning-group models.** *Result:* ΔM1 −0.0074, ΔM2 −0.0223. Splitting the data by group lost more than group-specific capacity gained. *Firmness:* firm for one seed.

**3.3 Partial pooling.** *Tried:* a per-kinase ridge residual or residual boosting, both shrunk toward the Manning family. λ was chosen from {0, 3, …, 300, ∞} on fold 0. *Result:* λ = 300 was selected, which is effectively the global model (ΔM1 +0.0002 [−0.0001, +0.0005]). *Side finding:* a small λ (3–10) raised R2-alone M2 from 0.192 to 0.232 for an M1 cost of 0.001–0.003. That observation led to the post-TEST intercept experiment (§5.1). *Firmness:* firm that it gives no M1 gain.

**3.4 Ranking objectives for R2.** On fold 0, lambdarank grouped by kinase matched the binary loss on M1 (0.754) but collapsed M2 (0.086 vs 0.223), because its scores are not comparable across kinases. Site-grouped lambdarank scored M1 0.730 and M2 0.189. *Firmness:* not shown to help; fold 0 only. Source: `docs/model_r2.md` §4.

**3.5 R5: fine-tuned ESM-2 with a late-interaction head.** *Tried:* the brief asked for a cross-encoder, but in a timing pilot, inference alone took 3.6–14.3 GPU-h per fold (budget about 1.5 GPU-h). A ColBERT-style late-interaction model on a fine-tuned ESM-2 35M or 150M trunk was built instead. *Pre-registered gate:* fold-0 M1 ≥ 0.751 or M2 ≥ 0.241 (R3 + 0.01 / + 0.02). *Result:* the best configuration scored M1 0.7250 and M2 0.1668, against 0.741 and 0.221 for R3. Removing the late-interaction term improved it (M1 0.7270, M2 0.1985), and the late-interaction term alone collapsed to 0.678. *Why:* early stopping chose step 100–600 of 3,000 in every seed (about 7.4k training positives), and the 150M trunk had the best inner-validation AUROC (0.787) but a worse M1 (0.712): the trunk overfits training proteins. Token-level matching added no signal at this label volume. *Firmness:* firm for this training regime; freezing most of the trunk (LoRA) was proposed but not tried. Cantley distillation replicated on this architecture (+0.034 M1; positive canary 0.895 → 0.958). Source: `docs/model_r5.md`.

**3.6 A32(b): prior-debiased site score.** *Motivation:* M2 for kinases with 5–19 labels was 0.142 (final fusion), below Cantley's 0.162. *Tried:* subtract β·log(training count) from the site score, β ∈ {0, 0.1, 0.25, 0.5, 0.75, 1}. *Pre-declared rule:* adopt β > 0 only if pooled M2 falls by no more than 0.003 and the 5–19 bin improves. *Result:* β = 0 was best pooled and was chosen in all 5 nested folds. At β = 0.1, pooled M2 fell by 0.027 (0.263 → 0.236). The 5–19 bin moved only from 0.142 to 0.146, and the ≥100 bin fell from 0.621 to 0.526. *Why:* a global count prior is too blunt. *Firmness:* firm. Source: `docs/strata_cv.md` §9.

**3.7 Simple learned baselines (B2, B7).** These baselines were built to check that the learned models beat simple alternatives, not as candidates.
- *B2, learned PSSM with a family prior:* M1 0.660, below Cantley (paired Δ −0.056 [−0.071, −0.038]).
- *B7, logistic regression on Cantley, the learned PSSM, a PPI edge, localisation and per-kinase intercepts:* M1 0.726. Against Cantley on matched pairs, ΔM1 was +0.0047 [−0.0052, +0.0154], so it was no better on M1. It was better on M2 (+0.042 [0.023, 0.057]).
- *Canary check:* when canaries shared one group prior, B2's positive canary fell to 0.849 (fail). This led to A25, which gives each canary its own group.

*Firmness:* firm. On this benchmark, learning a motif matrix from PSP labels alone does not match the label-free atlas. Source: `docs/REPORT_internal.md` §5.1.

**3.8 Planned but not built** (*not tested*). Deferred in DESIGN §0: KolossuS 15B, an AF-Multimer pilot, unobserved S/T/Y sites (U_all), the L3 label variant and docking-motif features. Not implemented in R3: nnPU, hard-negative mining and a third loss term. Also not built: F6d (kinase-activity profiles) and the promised `scripts/run_all.sh`. A two-stage reranker for R5 was rejected because the frozen harness ranks unscored pairs last. Sources: `docs/DESIGN.md` §0; `docs/model_r3.md` §5; `docs/model_r5.md` §1; `docs/REPORT_internal.md` §10.

## 4. Evaluation, controls and calibration

**4.1 The pooled-OOF metric defect (A24).** Pooling out-of-fold scores across folds lets fold-varying offsets, such as per-kinase intercepts and count features, distort per-kinase AUROC. The label-free B0-prior control scored M1 0.40 instead of 0.50. An audit found the defect at 17:01 (T+3.7 h), after 47 of the eventual 135 leaderboard runs had been scored. *Fix:* rank-normalise within fold and kinase for M1, M3, M4 and M5 (M2 unchanged), and re-evaluate every run. *Effect:* B1 did not move (ΔM1 0.000), learned models gained +0.010 to +0.014 M1, and the R2-full margin over Cantley widened from +0.046 to +0.056. The artefact had also inflated the raw 1–4-label bin (R2 0.821, 0.764 after the fix). This is a post-freeze metric change, applied to all methods; TEST, with one model, is unaffected. Sources: `docs/DESIGN.md` A24; `docs/audit/audit_cv_r1.md` M-1, M-4; `docs/retro/final_retro.md` §1.3.

**4.2 Negative canaries: v1 failed as a null; v2 passes only pooled.** Negative canaries are fake kinases with randomly chosen "substrates", which any model should score near 0.5. Canary v1 drew random U_ms sites, which tend to lie on proteins with few phosphosites, so a label-free score (the negative count of U_ms sites on the protein) reached 0.76. Canary v2 draws real labelled sites, clustered by protein and uniform over sites; the same probe gives 0.489. Under v2 the final fusion passes pooled (0.539 ≤ 0.57). Per fold, however, 2 of 5 folds exceed 0.57 (0.575 and 0.676), and every model sits above 0.5, including a label-free R2 control (0.539 on fold 0). A sampling difference between canaries and E-hard negatives is the likely cause, but this is not proven. The proposed site-uniform v3 canary was not built. *Firmness:* defect fixed for v1; the v2 gate is weaker than intended. Sources: `docs/core_infra.md` §6; `docs/audit/audit_cv_r1.md` M-3; `docs/REPORT_internal.md` §6.7.

**4.3 Coverage metric as first designed.** Counting sites in some kinase's top 1% was vacuous: random scores already cover 98.5% of S/T and 64% of Y sites. It was replaced by a site-level assignment probability, reported next to a random baseline (A17). Source: `docs/reviews/r2_statistics.md`; `docs/DESIGN.md` A17.

**4.4 The ≥0.8 assignment tier (A31(b)).** The first assignment probability was a softmax over per-kinase calibrated scores. On labelled DEV sites, its [0.8, 1] tier had a claimed probability of 0.87 and an observed precision of 0.54 (n = 167). Y sites in that tier were 0 of 41 correct, and the Y Platt slope was negative in 4 of 5 folds. After isotonic recalibration, done separately for S/T and Y, only 2 DEV sites reach 0.8 and no Y site reaches 0.5. Earlier coverage figures (14,235 sites at 0.5 and 1,674 at 0.8) were withdrawn. Cantley-only coverage shows a ≥0.8 tier (1,828 sites), but it comes from the tail of an isotonic map fitted on 152 sites. *Firmness:* firm that the model has no reliable ≥0.8 tier. Sources: `docs/REPORT_internal.md` §4.5, §9; `docs/DESIGN.md` A31.

**4.5 Calibrated assignment probability as the M2 score (A31(a)).** The design planned M2 on the calibrated assignment probability. That score gave M2 0.205, against 0.263 for the raw fused score (nested selection 0.259). The raw score was adopted at 22:12, after both values had been seen. The choice was made on cross-fitted CV only, but it is a late change in a favourable direction on a headline endpoint (see §6). Sources: `docs/DESIGN.md` A31; `docs/retro/final_retro.md` §3.2.

**4.6 Cantley on single-residue kinases (A30).** Ten dual-specificity kinases have an atlas matrix for one residue type only. On CV their other residue was ranked last; on TEST they were left unscored in both Cantley arms (386 of 503 kinases scored), which can only help the baseline. The issue was found during TEST preparation, not when B1 was built. Source: `docs/DESIGN.md` A30.

**4.7 Stale cross-fitted stores.** R2's label-derived stores for folds 0 and 3 were built from an intermediate canary table (331/329 training canary pairs instead of 348/350), so the canary gate was not end-to-end for those folds; real-kinase features were unaffected. The stores were rebuilt at 17:30, before the final R2 seed runs. Source: `docs/audit/audit_cv_r1.md` M-2.

## 5. Post-TEST exploration (A34; CV only, nested; released model unchanged)

The A34 bar for "promising" was a nested-CV gain on the target stratum with CI > 0 and an overall M1 loss of at most 0.002. An independent auditor reproduced every number and overruled two "promising" verdicts by the agents who ran the experiments. Sources: `docs/posttest_cv.md`; `docs/posttest_tyr.md`; `docs/audit/audit_posttest.md`; `docs/REPORT_internal.md` §11.1.

**5.1 Per-kinase intercepts on the site score (M2 for understudied kinases).** The pre-declared unweighted variant made its target worse: 5–19-label bin −0.013 [−0.021, −0.004], pooled −0.006. A kinase-balanced variant was added after the first unweighted grid point had been seen. It raised the 5–19 bin by +0.085 [+0.048, +0.116] and pooled macro M2 from 0.263 to 0.299. It also lowered the ≥100-label bin by 0.124 (28 of 31 kinases worse), lowered per-site Hit@5 from 0.535 to 0.465, and cut the share of top-5 slots held by unlabelled kinases from 14% to 2.3%. *Why:* the intercepts move top-5 slots between kinases; they do not add information. *Firmness:* firm that the pre-declared variant fails; the balanced variant is a trade-off and would need a pre-registered test on new labels.

**5.2 Priming-gated blend with priming-aware Cantley.** On primed positives the blend gained +0.018 [+0.013, +0.024], but unprimed positives of the same kinases lost 0.006 [−0.010, −0.001]. All positives of the priming kinases moved by +0.001 [−0.003, +0.005], and M3 fell slightly. A label-free flag ("the priming position is an observed phosphosite") with no Cantley term reproduced the primed gain (0.909–0.917). Against negatives carrying the same flag, the gain was +0.002. *Why:* the primed stratum is defined by the same property that the blend rewards. *Firmness:* firm that this is a rank shift.

**5.3 Y-specific fusion weight.** Y M1 over 38 kinases changed by −0.0039 [−0.0097, +0.0010]. The weight was re-chosen at the global value in 3 folds and gave Δ < 0 in the other 2. *Firmness:* not shown to help.

**5.4 Tyr motif repair.** *Tried:* replace a kinase's Cantley Tyr row with its own, nearest or cross-family atlas row, chosen in inner folds. The motivation was atlas matrices that are anti-correlated with their own labels. *Result:* +0.0042 [−0.0103, +0.0167]. Inner-fold gains of +0.012 to +0.027 shrank to +0.0056 on held-out folds (3 of 5 positive). Of 38 Y kinases, 18 improved and 19 worsened; two predicted beneficiaries lost (TAM/Axl −0.053, Tec −0.031). The repair parameter γ hit the top of its grid (0.5) in every fold. *Why:* 38 kinases and about 1,150 Y positives cannot support selecting even one parameter. *Firmness:* not shown to help; Tyr specificity remains unsolved.

## 6. Process failures of the agent team

The points below summarise the fresh-agent retrospective (`docs/retro/final_retro.md`), with additions from other committed files where marked. They concern scheduling and discipline; the retrospective did not re-examine the science.

1. **The sealed TEST set was spent early, and no spare shard existed.** TEST finished at T+11.6 h of 24 (gate T+16 h); FREEZE was committed 5.3 h before its gate and 56 min before the PI's own deadline for CV exploration. The three post-TEST experiments on known weaknesses could therefore only be exploratory. No second sealed shard was reserved, although 11,609 substrate clusters would have allowed one.
2. **Most of the budget did not change the model.** The delivered configuration was the 19:19 fusion run (T+6.0 h). Of 135 timed CV runs, 62 (46%) started at or after 20:00, and the best of them reached M1 0.774 against 0.786. Compute was not the constraint: all 135 runs took 11.8 h of serial wall time.
3. **A required experiment was scheduled late.** The brief asked for a training-regime comparison. It was delegated at T+6.6 h, used 55 of 135 runs (38% of serial CV compute), produced no adoption, and closed one minute before FREEZE.
4. **Controls were not run through the harness first.** The A24 defect (§4.1) was found by an audit after 47 runs. Canary v1 (§4.2) had to be replaced.
5. **Searches continued past the noise floor.** Seed and knn_sim sweeps account for 74 of 404 run directories; knn_sim gained +0.003 and was not adopted.
6. **Code on the irreversible path was reviewed late.** The TEST-harness review started at 20:25, and its fixes were applied at 23:30, 34 min before FREEZE. It found a bypassable single-batch guard and a crash path that left no trace.
7. **A known stronger comparator was left out of FREEZE.** Priming-aware Cantley beat the model on primed CV positives (0.945 vs 0.893). This was recorded at 21:11 and in A32, but the arm was not frozen. As a result, the TEST primed-site claim holds only against the non-priming-aware arm (`docs/audit/audit_test.md` verdict #11).
8. **The protocol was still being amended while results arrived.** A23–A34, 12 of 34 amendments, were made in the last 8.6 h before and around the freeze, and six of them (A23–A28) were made within 40 min. Each followed a result. Most were conservative, for example A30, A23 and A24. A31(a) changed M2's score column late, in the favourable direction (§4.5). M2's definition changed four times (A18, A20, A27, A31(a)).
9. **A rule was broken after it was written.** The balanced intercept variant in §5.1 was added after the first result had been seen, eight hours after the lesson that forbids this. LESSON 14 ("no amendments after FREEZE") was itself written one minute after FREEZE, and A34 followed it.
10. **Memory was budgeted against the wrong limit.** The retrospective counts two OOM kills (17:50 and about 23:00) from budgeting against the node's 375 GB instead of the 128 GiB job cgroup; the second came after the memory rule was written. Agent documents record more: about 20:14, 20:40 and 21:05 (`docs/regimes_labels.md` §6), 20:50 (`docs/strata_cv.md`), two GPU OOMs in R5 (`docs/model_r5.md` §4), and the downstream harness on the 131 M-row score file (`docs/model_r2.md`).
11. **A mandatory baseline nearly missed FREEZE.** Pretraining of B3 (Phosformer; 4.3 GPU-h) slipped from 19:15 to 22:30 and then to 23:30, 34 min before FREEZE. The PI also fixed its negative ratio at the paper's 1:16, replacing the pre-registered choice between 1:16 and 1:50 (`docs/REPORT_internal.md` §5.2).
12. **Status records went stale.** `docs/DELEGATION.md` was stale at every retrospective; `docs/RESEARCH_TREE.md` carried 11 "in progress" markers, 8 for finished work.
13. **There was a long idle gap.** No file was written for about 6 h 15 min (03:17 to 09:32).
14. **The single batch omitted analyses that cannot now be added:** hub-excluded Y, external and temporal calibration, TEST coverage baselines, a paired T3 test, stratum CIs and the priming-aware arm (`docs/REPORT_internal.md` §7.9).

The retrospective's general finding: rules built into code (the cross-fit driver that raises `LeakageError`, the canary gate, SHA-pinned FREEZE inputs, a standing auditor) changed behaviour; rules kept as prose did not. It also records what worked: failed pre-registered gates (R5, knn_sim, F9, A32(a/b), A34) were accepted without re-arguing them, and four independent audits caught the A24 artefact, the stale stores, the overconfident ≥0.8 tier and the bypassable guard. A final review of the report made 7 factual fixes, none changing a scientific claim (`docs/review_report_v3.md`).

## Notes on sources

- The retrospective counts two OOM kills; agent documents record more (item 10). Some may be the same events.
- `docs/LESSONS.md` calls the TEST priming claim "retracted"; `docs/audit/audit_test.md` rates it "trusted with caveats" (comparator was the non-priming-aware arm only). This text follows the audit.
- `docs/REPORT_internal.md` §3.2 gives cross-fitted knn_sim AUROC as 0.59–0.76; the fold-0 table in `docs/features_knn.md` shows 0.653–0.763 for the columns it lists. The REPORT range is used.
