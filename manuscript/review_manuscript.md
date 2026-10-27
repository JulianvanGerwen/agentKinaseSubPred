# Independent review of `manuscript/manuscript.md` and `manuscript/supplement.md`

Reviewer: independent agent; I did not write the manuscript. Read-only on both files; no git; nothing under `sealed/`, `data/` or `external/` was opened. Reference formatting was not reviewed (another agent owns it).
Locators are `manuscript.md:<line>` and `supplement.md:<line>` as of 2026-10-27 12:08 CET.

## Verification summary

- **Recomputed from TEST batch files** (`results/test/721dc1977366/`):
  - All 9 primary deltas, CIs, p_boot and Holm p (`superiority.tsv`), and all 17 secondary tests (`secondary.tsv`).
  - Every standalone M1–M4 value and CI in Table 3 (`*_ci.tsv`); Table 3 runtimes (`results/test/LEDGER.md`).
  - All 55 Table 4 cells (`*_strata.tsv`, `*_ext_tiers.tsv`).
  - Tyr kinases above B0-pop: 10 of 23. M3 against B4: 52 / 27 / 40 ties (`*_per_kinase.tsv`).
  - T3 M1 for final, R2, R3, B7 and B2 (0.752 / 0.744 / 0.703 / 0.708 / 0.616), and Tyr T3 0.531 (`*_t3.json`, `final_t3_per_kinase.tsv`).
  - Reliability: n = 1,458, ECE 0.026, top-1 0.353. The ≥0.5 tier has n = 236 + 37 + 72 = 345, claimed 0.616 and observed 0.603 (`final_reliability.tsv`).
  - Coverage: 29,886 / 26,762 / 0 / 30,258 / 1,828 / 6, and 8.52% of 350,736 (`results/release/coverage.tsv`).
- **Recomputed from CV files:**
  - CV deltas (`results/models/r2/compare__*.tsv`) and the fusion sweep (`fusion_r3_sweep*.tsv`).
  - Table 1, Table 2 and Table 5 rows (`results/cv/LEADERBOARD.md`, `results/models/r2/cv_table.tsv`, `docs/regimes_labels.md`).
- **Supplement:** I checked 30+ supplement numbers (§1.1–1.3, §3.1–3.6, §5.1–5.4, §6.1–6.10) against `docs/regimes_labels.md`, `docs/strata_cv.md`, `docs/model_r5.md`, `docs/audit/audit_posttest.md` and `docs/retro/final_retro.md`. All match except item 1.5 below.
- **Figures 1–5:** I read them all. Fig 2 matches Table 3 and `superiority.tsv`, and its D panel shows the correct 0.156. Fig 3 matches Table 4. Fig 4 matches `final_reliability.tsv` and `coverage.tsv`. For Figs 1 and 5, see 2.E and 5.
- **Method mismatches:** a code-reading sub-agent checked the Methods. I spot-checked 4 of its 6 mismatches in the code (3.1, 3.2, 3.4, 3.5), and all 4 hold.
- **The flagged discrepancies are resolved in the Appendix.**

---

## 1. Must-fix errors

**1.1 Author summary misstates the Tyr result for the other methods** (`manuscript.md:15`). On Tyr kinases, Cantley scored 0.496, B3 0.534 and B4 0.575. All three are below the B0-pop floor of 0.585; only the final model (0.595) is above it (`*_strata.tsv`, residue_class=Y).
> Replace "For tyrosine kinases, our method and the others did little better than a simple rule that favours well-studied sites." with: "For tyrosine kinases, our method did barely better than a simple rule that favours well-studied sites, and the other methods did worse than that rule."

**1.2 The Abstract says "Accuracy did not drop from cross-validation to test"** (`manuscript.md:11`). This is false for Tyr kinases: 0.667 → 0.595, paired Δ −0.072 ± 0.033 (`docs/audit/audit_test.md` §4). It also invites a CV→TEST reading of M2, which the audit forbids (§0 #7, §7.1). A full replacement Abstract (297 words, within the PLOS limit of 300) is below. It also carries items 2.1, 2.2 and 2.11.
> "Most human phosphosites have no known kinase, and published kinase–substrate predictors are often evaluated on random splits over kinase–substrate pairs, which leak information through shared peptides, homologous proteins and other sites on the same protein. We built a benchmark that holds out substrate-protein clusters, cross-fits every label-derived feature by cluster, uses paralog-masked hard negatives and weights kinases equally. On it we trained a predictor that scores 503 human protein kinases against every mass-spectrometry-observed S/T/Y site of matching residue class, treating PhosphoSitePlus annotations as positive–unlabelled data. The final model is a late fusion of a global LightGBM pair model, which uses motif, structure, language-model, interaction, localisation, co-expression and phosphoproteomics features, and a neural kinase × site model pretrained by distilling label-free Cantley motif-atlas scores. All choices were made on 5-fold cross-validation (203 kinases; macro AUROC 0.786). The model was then scored once on a sealed test set of 2,353 substrate clusters. There, macro AUROC on hard negatives was 0.794 [95% CI 0.774, 0.811], site-centric Hit@5 0.363 (0.348 for the same kinases in cross-validation) and top-1% recall 0.363. The model beat the Cantley atlas and our retrained Phosformer and SELPHI2.0 on all nine pre-registered, Holm-corrected tests; the margin over SELPHI2.0 on top-1% recall was modest (+0.059 [0.010, 0.120]). Macro AUROC did not drop from cross-validation to test, and on unseen kinase families with unseen proteins it was 0.752. Tyrosine kinases scored 0.595, barely above a site-popularity floor of 0.585, and fell from cross-validation. No site reached an assignment probability of 0.8. Part of the margin over Cantley reflects a weak test draw for Cantley. Phosphoproteomics features added little, and for kinase-activity inference the predicted library was at parity with PhosphoSitePlus overall and worse on unseen drugs. AI agents ran the study; no human has yet checked the code or analyses."

**1.3 The Tyr explanation contradicts the project's own diagnosis and S1** (`manuscript.md:142`, `:199`).
- The main text attributes the Tyr gap to hub sites and "probable indirect in-vivo events".
- But for the models, excluding hub sites left R2's Tyr M1 unchanged on CV: 0.627 against 0.629 (`docs/strata_cv.md:136,211,225`: "Hubs are not the cause"). The same result is in `supplement.md:87`.
- The in-vivo label-noise explanation was tested and rejected: within 17 Tyr kinases, 0.603 vs 0.596 (`docs/features_tyr.md` §1.3; `supplement.md:87`).
- The "quarter of the gap" decomposition was done on the Cantley atlas (B1 arm), not on the model (`docs/features_tyr.md` §1, header).
- `REPORT.md:527` and `:959` carry the same wording; please tell the PI.

> Line 142: replace "On CV, about a quarter of the Tyr gap came from macro dilution … (FGFR1, BTK, Eph, VEGFR, TAM)." with: "A CV diagnosis of the Cantley atlas's Tyr gap attributed about a quarter of it to macro dilution, because 42% of Tyr kinases have three or fewer positives. It linked the rest to E-hard pools rich in promiscuous pTyr hub sites (one GLO1 site is labelled to 15 Tyr kinases) and to atlas matrices that anti-correlate with their own labels for several families (FGFR1, BTK, Eph, VEGFR, TAM). Within kinases, in-vivo-only labels were not harder than in-vitro labels (0.596 against 0.603). For the R2 model, however, removing hub sites left Tyr M1 unchanged on CV (0.627 against 0.629), so hubs do not explain the model's Tyr gap. We have no tested explanation for it."
>
> Line 199: replace "Hub sites and indirect in-vivo labels dominate the Tyr hard negatives, so part of the problem may lie in the endpoint." with: "Removing hub sites did not change the model's Tyr M1 on CV, so the cause is unknown, and the hub-excluded TEST stratum was not computed."

**1.4 The in-vivo-only regime explanation contradicts S1** (`manuscript.md:152`). `supplement.md:63` and `docs/regimes_labels.md:125` say the comparison is size-confounded: 2,808 pairs against 8,987.
> Replace "Training on in-vivo-only evidence was worse on every endpoint, which fits the idea that these labels are enriched for indirect events." with: "Training on in-vivo-only evidence was worse on every endpoint. However, that label set is about three times smaller than the in-vitro set (2,808 against 8,987 pairs), so the comparison mainly measures label count."

**1.5 The no-phospho comparison is not like-for-like** (`supplement.md:91`, and implicitly `manuscript.md:59`).
- The no-phospho variant (0.771) uses the final v3 configuration at seed 0.
- Its matched reference is v3 seed 0 = 0.768, or the five-seed mean of 0.7711 (`docs/model_r2.md:284-290,457-458`).
- It is not the ablation-configuration run `r2_full` (0.773).

> `supplement.md:91`: replace "a variant without F6 scores 0.771 against 0.773" with "in the final configuration, a variant without F6 scores 0.771 against 0.768 with F6 (seed 0; five-seed mean with F6 0.771)".
>
> `manuscript.md:59`: replace "The phosphoproteomics families together added about +0.002 M1, and a variant without them scored 0.771." with "Along the ablation curve the phosphoproteomics families together added about +0.002 M1. In the final configuration, removing them gave M1 0.771 against 0.768 with them (seed 0), inside the noise."

**1.6 Numeric errors** (small, but they are errors):
- (a) `manuscript.md:86` gives "+0.041 [0.007, 0.060]". The source value is 0.040451 (`secondary.tsv`), so this is a double-rounding error. Replace with "+0.040 [0.007, 0.060]". The same drift is in `docs/audit/audit_test.md` §0 #15.
- (b) `manuscript.md:71` (F6b M2 0.154), `:74` (full M2 0.208) and `:59` ("from 0.160 to 0.208"). The cited source `cv_table.tsv` and the A24 leaderboard (rows 40 and 15) give **0.153** and **0.207**. `REPORT.md:477,480` and `docs/model_r2.md:184,188` give 0.154 and 0.208. The caption says "A24" and cites `cv_table.tsv`, so use 0.153 / 0.207 and replace "from 0.160 to 0.208" with "from 0.160 to 0.207". Tell the PI about the REPORT drift.
  - Alternative: keep the REPORT values and cite only `REPORT.md` §6.2.

---

## 2. Overstatements and audit consistency (with replacement text)

**2.1** `manuscript.md:11` says "AI agents ran the study with limited human verification", and `:329` says "Human verification was limited to supplying inputs". There was no verification (`:310`). Use the last sentence of the 1.2 Abstract. For `:329`, write: "There was no human verification of code, analyses or text; the human only supplied inputs. The publishing researcher is responsible for checking the results before use."

**2.2** The Abstract and `:86` quote TEST M2 0.363 without the CV value for the same kinases (0.348). The audit requires it next to 0.363 (`audit_test.md` §0 #7). This is carried in the 1.2 Abstract. Results `:110` already has it.

**2.3** `manuscript.md:15` says "ranks the correct kinase better than three widely used methods". Two of the three are our retrains, and "widely used" has no source in the project files.
> Replace with: "ranks the correct kinase better than a widely used motif atlas and our own retrained versions of two published predictors".

**2.4** Heading `manuscript.md:84`, "TEST: the frozen model beats all three named baselines".
> Replace with: "TEST: the frozen model beats Cantley and our retrains of Phosformer and SELPHI2.0".

**2.5** `manuscript.md:193` says "The design choices that mattered most were not model choices. They were the cluster split, cross-fitting …". No experiment ablated these controls, so the claim is an unsupported comparative causal statement.
> Replace with: "We did not ablate the benchmark controls themselves, but several observations point to them as the main reason why controlled numbers fall below published ones. Label features that were not cross-fitted reached AUROC 0.93–1.00, against 0.59–0.76 when cross-fitted. 51% of Phosformer's test-positive peptides occur as training positives. SELPHI2.0's random-pair metric stayed at 0.947 on our split. Part of the baseline drop may also reflect our limited retraining budgets (Limitation 5)."

**2.6** `manuscript.md:195` says "Late fusion of the tree and neural models helped on every endpoint, on CV and on TEST." It did not help on CV T2: fusion 0.749 against R2 0.751 (`docs/model_r2.md:17`).
> Replace with: "Late fusion helped on every endpoint for seen kinases on CV and TEST, but not for held-out families on CV (T2: 0.749 against 0.751 for R2 alone)."

**2.7** `manuscript.md:195` says "… which agrees with the idea that understudied kinases borrow signal through shared descriptors." The supporting evidence is about data-rich kinases.
> Replace with: "One global model with kinase descriptors beat per-kinase and per-group models, even for the 78 kinases with at least 20 labels in every fold (0.805 against 0.779). Sharing across kinases therefore helps data-rich kinases too; we did not test directly whether understudied kinases gain through the descriptors."

**2.8** The post-TEST paragraph (`manuscript.md:189`) omits wording that the audit requires (`audit_posttest.md:101-108,168-172`).
> Replace "Per-kinase intercepts raised M2 … from well-studied kinases. A priming-gated blend … is an observed phosphosite." with: "The pre-declared per-kinase intercepts lowered M2 for 5–19-label kinases (−0.013 [−0.021, −0.004]). A balanced variant, added after the first result was seen, raised it. It did so only by taking top-5 slots from well-studied kinases and from kinases with no labels, whose share fell from 14% to 2.3%, and per-site Hit@5 fell from 0.535 to 0.465. A priming-gated blend with priming-aware Cantley raised primed positives by +0.018. However, a label-free indicator of whether the priming position is an observed phosphosite reproduced the gain, and against negatives with the same flag the gain was +0.002."

**2.9** Audit re-cuts are quoted as TEST results without a label. `audit_test.md:12-14` says such rows "must not be quoted as TEST results". The affected quotes are:
- `:86`: "97 of 116", and 52 vs 27;
- `:171`: S/T 0.615 / 0.603 (n = 345); Tyr 0.41 / 0.34; the 0.445 cap.
> After each one, add "(audit re-cut of the batch per-kinase or per-site tables)", or move them into a sentence that begins "An audit re-cut of the batch tables showed …".

**2.10** Spot checks (`manuscript.md:177`). TP53 S15→ATM and AKT1 S473→MTOR are textbook pairs. They are very probably training labels of the deployment model, which was refit on DEV+TEST, so they test recall rather than prediction. I could not confirm this without reading PSP; the PI should check.
> Replace "Spot checks showed both hits and misses. TP53 S15 was assigned to ATM at 0.70, which is correct." with: "Spot checks on well-known sites showed both hits and misses. These sites are probably training labels of the deployment model, so the checks test recall rather than prediction. TP53 S15 was assigned to ATM at 0.70."

**2.11** Coverage (`manuscript.md:177`). After "(30,258)", add: "so coverage at 0.5 does not separate the model from the label-free atlas."

**2.12** Runtime (`manuscript.md:185`). After "(2,043 s)", add: "These times exclude feature construction, which B3 does not need: cross-fitting alone costs about 21 CPU-minutes per fold, on top of the ESM-2 and structure features."

**2.13** `manuscript.md:82` says "The released Phosformer weights reached 0.796". This was on only 2% of the required pairs (`LEADERBOARD.md` row 1: coverage 0.020).
> Replace with: "The released Phosformer weights reached 0.796 on the 2% of required pairs they could score".

**2.14** `manuscript.md:144` says "Kinase-equal M1 held up across label bins: 0.759 at 5–19 labels". The 13 TEST kinases in this bin are its best-labelled members, which the same paragraph concedes for M2.
> Replace with: "Kinase-equal M1 held up across label bins on TEST (0.759 for the 13 kinases with 5–19 labels, which are the best-labelled of that bin, and 0.795 at ≥100 labels)".

**2.15** `manuscript.md:59` says "because its per-kinase scales are not comparable across kinases; the kinase descriptors repaired this". This is asserted in `docs/model_r2.md:193` but was not tested.
> Replace with: "probably because its per-kinase scales are not comparable across kinases; M2 recovered once the kinase descriptors were added".

### 2.E Agent-setup statements not recorded in project files

- **E1.** "the 'PI', Claude Opus 5.5" appears in `:5`, `:327`, Fig 5 and `supplement.md:3`. No project file records the lead agent's model. `PROMPT.md:1` names only the models for the team. The source is the publication brief, so the human must confirm it.
- **E2.** The flag says sub-agents were "(Claude Opus 5.5 and Sonnet 5.5)", which implies both models were used. The files do not record whether Sonnet was used at all, and `:306` concedes "which model each used was not recorded". The flag wording is mandated, so the human/PI decides. I suggest "(Claude Opus 5.5 or Sonnet 5.5; the model of each sub-agent was not recorded)".
- **E3.** `:310` says the human "did not … select models … and did not check the code or analyses before this manuscript was drafted". This is not recorded anywhere; the human must confirm it.
- **E4.** Fig 5 says "Fresh-agent retros at M1, M2, M3". `docs/retro/M1_retro.md`, `M2_retro.md` and `M3_retro.md` have no author line. Only `final_retro.md:1-3` and `LESSONS.md:2` record a fresh agent.
  > Replace with "Milestone retros (author not recorded); final retro by a fresh agent".
  - The milestone labels M1–M3 also collide with the endpoint names M1–M3.
- **E5.** Fig 5 says "Lead agent … recomputes headline numbers". This is not recorded; the recomputation was done by audit agents (`audit_test.md`, `review_report_v3.md`).
  > Replace with "commissions independent recomputation".
- **E6.** `:306` says agents "received a written brief with GPU, core and memory budgets". `docs/DELEGATION.md` lists only GPU and core budgets. The ≤60 GB memory cap was added to `docs/AGENT_BRIEF.md` only at 17:55, after the first OOM. Individual task briefs are not archived.
  > Replace with: "received a written task brief and a shared agent brief (`docs/AGENT_BRIEF.md`) with GPU and core budgets; a per-agent memory cap of 60 GB was added at 17:55. Task briefs were not archived."
- **E7.** `:298` says seeds were not kept "in a central list, which deviates from the brief". In fact `docs/AGENT_BRIEF.md` fixes the list [0, 1, 2, 3, 4].
  > Replace with: "The shared agent brief fixed the seed list [0–4]. Downstream analyses used the benchmarKIN and decoupler seeds 42 and 1, and the seeds actually used are recorded per run manifest."
- **E8.** `:298` gives "375 GB RAM". Jobs ran under a 128 GiB cgroup (`final_retro.md:179`). Add "(jobs were limited to a 128 GiB cgroup)".

---

## 3. Method mismatches (Methods against code)

- **3.1** E-hard (`:225`): the code requires a *human-native PSP* non-autophosphorylation label (`src/ksp/eval.py:17`, `src/ksp/labels.py:602`, `src/ksp/models/r2/sample.py:9-10`).
  > Replace "a non-autophosphorylation label" with "a human-native, non-autophosphorylation PSP label".
- **3.2** F1 (`:240`): the final models use only the score and the U_ms percentile (unprimed, soft-primed and hard-primed) plus the imputation distance (`src/ksp/models/r2/assemble.py:91-97`). R3 also uses universe percentiles.
  > Replace "scores and percentiles against the Ochoa reference (as published) and against U_ms" with "scores and percentiles against U_ms (the published Ochoa-reference percentile is used only as baseline arm B1(i))".
- **3.3** F4 (`:243`): the final R2 also uses a COMPARTMENTS *experiments* cosine, a combined cosine, an HCM NMF cosine and a Hein profile correlation (`assemble.py:118-119`).
  > Add "and the COMPARTMENTS experiments channel (its knowledge channel was excluded)".
- **3.4** R3 loss (`:262`): the code uses a balanced masked BCE over all mask-filtered kinase × site entries in a batch. Each batch holds 1,024 sites, 25% of them carrying a training positive, and positives are weighted by n_k^−0.5. No hard-negative selection is done (`src/ksp/models/r3/model.py:115-125`, `train.py:10`, `scripts/models/r3_cv.py:46`).
  > Replace "binary cross-entropy with masked hard negatives" with "a balanced binary cross-entropy over all mask-filtered kinase × site entries in each batch of 1,024 sites (25% carrying a training positive; positives weighted by n_k^−0.5)".
- **3.5** R3 Cantley descriptor (`:262`): the code uses 32 PCs per residue type, 64 in total, from unprimed scores on 3,000 S/T and 1,000 Y reference sites (`src/ksp/models/r3/data.py:112,125`).
  > Replace "a 32-dimensional PCA of the kinase's Cantley score profile" with "32 principal components per residue type (S/T and Y) of the kinase's unprimed Cantley scores on 3,000 S/T and 1,000 Y reference sites".
- **3.6** Canary gate (`:286`): the positive canary is gated on both E-hard and E-all (`eval.py:610,617`).
  > Replace "the gate was ≤0.57 and ≥0.90 on E-hard" with "the gate was ≤0.57 on E-hard for the negative canaries and ≥0.90 on both E-hard and E-all for the positive canary".
- **3.7** CPTAC residualisation (`:245`) also regresses out stromal, immune and purity scores (`docs/features_phospho.md:30-31`).
  > Replace "residualised on protein abundance" with "residualised on protein abundance and on stromal, immune and purity scores".
- **3.8** Fusion grid (`:268`): the ensemble sweep (`fusion_r3_sweep__r2_final_ens5.tsv`) covers w = 0.3–0.8. The endpoints at 0 and 1 are the standalone R3 and the R2 ensemble.
  > Replace "over a grid from 0 to 1" with "over w ∈ {0, 0.3, 0.4, …, 0.8, 1}, where 0 and 1 are R3 and the R2 ensemble alone".
- **3.9** R2 search (`:258`): "noise floor about ±0.005" conflicts with `docs/model_r2.md:110`, which gives fold-0 noise of about ±0.01 (0.005 appears at lines 140-141 for a single perturbation). The early-stopping AUROC uses sampled hard- and E-all-type negatives ("mixed").
  > Replace with: "early stopping used macro per-kinase AUROC of positives against the sampled hard and E-all-type negatives of an inner validation split …; fold-0 noise about ±0.01".
- **3.10** R3 is not refit: the best checkpoint is restored, so 20% of training clusters never train R3. R2, by contrast, is refit. Add one sentence to `:264`.
- **3.11** Omitted hyperparameters, from the code check; please add them to Methods or S1:
  - R2: min_data_in_leaf 100, feature_fraction 0.8, bagging 0.8, λ_L2 1.0, max_bin 255, ≤3,000 rounds, patience 100. Refit rounds per seed on TEST: 477 / 200 / 180 / 392 / 303. 623,373 training rows (12,223 positives). F7 = ESM PCA-16, Cantley PCA 8 (S/T) + 8 (Y), and residue class.
  - R3: d_site 192, h 128, d_kin 256, 4 heads. AdamW lr 3e-4, wd 0.05, warm-up 200, ≤6,000 steps, patience 10. Dropout 0.2, ESM-input dropout 0.3. Pretraining lr 1e-3, batch 1,024.
  - B3: focal γ = 2, lr 2e-5, batch 64; the 1:16 negatives are 50% hard.

---

## 4. Gaps against the human's request

- **4.1 Credit to the human brief.**
  - The brief asked for grouping by protein and homologs, splits on sites or proteins rather than pairs, leave-kinase-out tiers, kinase-balanced metrics, a training-regime comparison, species pooling, in vitro/in vivo analysis and priming kinases (`PROMPT.md`; `docs/AGENT_BRIEF.md`). `:27` presents the benchmark as the agents' own design.
    > Add: "The human brief asked for splits that group sites by protein, for leave-kinase-out tiers, for kinase-balanced evaluation and for a comparison of training regimes; the agents turned these requests into the benchmark below."
  - In Methods "Splits", credit the prior split designs of Phosformer-ST (CD-HIT 70% + Hamming-5) and KolossuS (phylogeny-aware) (`LITERATURE.md:15,17`).
- **4.2 Literature credit for method components.** LightGBM, the InfoNCE/contrastive loss, positive–unlabelled learning, stacking/late fusion, Platt and isotonic calibration, and leakage in pairwise-input prediction are not credited. None is in `LITERATURE.md`, so they must be added there and verified ("[verify]") before citing. The credits that are present check out against project files: the reviewer proposal (`docs/reviews/r1_external_PI.md:137`), Phosformer-ST, KolossuS, DeepKinZero, CoPheeMap, SELPHI2.0 and structuremap.
- **4.3 Definitions of the evaluable sets are missing** (`:35`).
  > Add: "K_eval: ≥3 TEST positives on ≥2 clusters; K_cv: ≥5 DEV positives on ≥2 clusters; K_any: ≥1 TEST positive; pseudokinases need ≥3 labels" (`results/splits/split_report.md`).
- **4.4 Important decisions not explained:**
  - (a) Why TEST was opened at T+11.6 h rather than at its gate. No rationale is recorded, and `final_retro.md` §1.2 calls it contrary to the brief. Say so explicitly.
  - (b) Why the priming-aware arm was left out of FREEZE, although it was known at 21:11.
  - (c) Why B3's negative ratio was fixed at 1:16. `docs/baseline_phosformer.md:108` records the PI instruction "no hyperparameter changes after seeing fold metrics".
  - (d) Why design review round 3 was skipped (`RESEARCH_TREE.md:43`).
  - (e) The fold-0 result that chose 1:50 negatives with 50% hard.
- **4.5 Conceptual explanation.** Results and Methods lack a short conceptual account of the final method.
  > Suggest at the start of the R2/R3 Methods: "R2 integrates many weak, heterogeneous context features for each kinase–site pair, but its scores are only partly comparable across kinases. R3 learns one embedding space for site windows and kinase descriptors, so its scores are comparable across kinases and it can score kinases with few labels through their descriptors. Fusing z-scored outputs combines the two."
- **4.6 Missing agent-setup details.**
  - How sub-agents were launched and how much context they shared.
  - The total number of agents: about 20 owners in `DELEGATION.md`, plus 11 design reviewers, plus auditors and retro agents.
  - The number of human interactions (`DATAINEED.md` requests).
  - That the PI read sub-agent summaries rather than logs.
  - Token or compute cost, if it is recorded anywhere.
- **4.7 PLOS Comp Biol elements.** Present: Title, Authors placeholder, Abstract (291 words), Author summary (187), Introduction, Results, Discussion, Materials and methods, Data and code availability, Acknowledgments, References, and Supporting information. Missing:
  - Funding statement;
  - Competing interests;
  - Author contributions (CRediT);
  - a working repository URL (`:323` "[repository URL]");
  - figure files in TIFF/EPS for submission.
- **4.8** The AI-exercise flag is present at `:5`, `:327` and `supplement.md:3` ✓. The agent setup is in Methods (`:300-319`) ✓. The tone is generally critical ✓, except for the items in section 2.

---

## 5. Minor

- `:33`: 18,638 L2 pairs is the total before restriction to U_ms. 14,747 L2 pairs were used (`split_report.md`, "Dropped from training"). Add "(14,747 after restriction to U_ms)".
- `:35`: define M2 precisely as "for each kinase, the fraction of its labelled sites at which it ranks in the top five filtered candidate kinases, averaged over kinases". The audit shows that macro and micro Hit@5 differ.
- `:55`: Table 1 gives B0-pop M4 as "—", but `LEADERBOARD.md` row 81 gives 0.562.
- `:59`: "doubled M3, from 0.125 to 0.236" → "nearly doubled".
- `:90`: the column header "Score 25.9 M pairs" does not fit B1, which scores 19.2 M pairs (LEDGER pre-flight).
- `:116`: T2 R2 is 0.751 in `docs/model_r2.md:17,491` but 0.750 in `LEADERBOARD.md` row 45 (A24 re-eval). Keep REPORT's value and tell the PI.
- `:86`: "the fusion's weakest margin was M2 over R3". It is the least significant one; the lowest CI bound is R2 M4 (0.0068). → "Its least significant margin was …".
- `:314`: "about 5 h before its own gate". FREEZE was 5.3 h before the gate and the batch completed 4.4 h before 05:20. → "FREEZE was committed 5.3 h before its gate".
- Fig 1 says "Ledger refuses a second batch". Add "(unless a logged override is used)". The TEST-harness review found the guard bypassable before the fix (A29).
- Fig 2D: the title says "CI lower bound below 0.02", but the caption says "the one modest test". Align them.
- `supplement.md:79`: drop "[verify title]". The title is verified in `docs/lit/01_methods_benchmarks.md:249`.
- `REPORT.md` items for the PI, not the manuscript:
  - `:40` gives 0.157, which should be 0.156.
  - `:477` and `:480` give 0.154 / 0.208 where the A24 values are 0.153 / 0.207.
  - `:527` and `:959` carry the Tyr hub/in-vivo wording (see 1.3).

---

## Appendix: resolution of the flagged discrepancies

1. **Upper CI of M2 against B4, 0.157 vs 0.15649.** `superiority.tsv` gives hi95 = 0.1564936, so the correct value is **0.156**. `REPORT.md:40` gives 0.157, a double-rounding error (0.15649 → 0.1565 → 0.157). The manuscript text does not quote this bound, and Fig 2D shows 0.156, which is correct. Fix REPORT only. The audit's "+0.041" for final vs R3 M2 (item 1.6a) is the same kind of error.
2. **B3 CV M1, 0.677 vs 0.679.** 0.677 is the A24 re-evaluation of the saved predictions by the frozen harness (`results/cv/LEADERBOARD.md:79`, generated 01:23; fold mean 0.689). 0.679 is the earlier run-time evaluation in `docs/baseline_phosformer.md:121` (22:27). The same document gives B3-pub 0.793 against 0.796 on the leaderboard, and M4 0.657 against 0.655, which is a consistent offset between harness versions. The manuscript uses 0.677 (`:53`, `:82`), which matches REPORT v3 (`REPORT.md:406`) and the audit (`audit_test.md` §4). This is correct and needs no change.
3. **Müller-Dott reference [2].** Every element is traceable to project files:
   - title: `docs/lit/01_methods_benchmarks.md:249`;
   - Nat Commun 16:4771: `docs/downstream.md:7` and `docs/lit/02_ks_annotation_sources.md:120`;
   - DOI 10.1038/s41467-025-59779-y: `LITERATURE.md:33`.
   
   `DATA.md:40` spells the name "Muller-Dott", which is cosmetic. Reference [2] as written is consistent with the files, and the "[verify title]" at `supplement.md:79` can be removed.
