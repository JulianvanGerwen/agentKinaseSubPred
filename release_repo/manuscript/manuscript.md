# Predicting direct kinase–substrate sites across the human proteome: a leakage-controlled benchmark and a late-fusion model, built and evaluated by AI agents

[authors to be added by the publishing researcher]

> **This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.**

---

## Abstract

Most human phosphosites have no known kinase, and published kinase–substrate predictors are often evaluated on random splits over kinase–substrate pairs, which leak information through shared peptides, homologous proteins and other sites on the same protein. We built a benchmark that holds out substrate-protein clusters, cross-fits every label-derived feature by cluster, uses paralog-masked hard negatives and weights kinases equally. On it we trained a predictor that scores 503 human kinases against every mass-spectrometry-observed S/T/Y site of matching residue class, treating PhosphoSitePlus annotations as positive–unlabelled data. The final model is a late fusion of a neural kinase × site model, pretrained by distillation of the label-free Cantley motif atlas, and a global LightGBM model of kinase–site pairs with motif, structure, language-model, interaction, localisation, co-expression and phosphoproteomics features. All choices were made on 5-fold cross-validation (203 kinases; macro AUROC 0.786). The model was then scored once on a sealed test set of 2,353 substrate clusters. There, macro AUROC on hard negatives was 0.794 [95% CI 0.774, 0.811], site-centric Hit@5 0.363 (0.348 for the same kinases in cross-validation) and top-1% recall 0.363. The model beat the Cantley atlas and our retrained Phosformer and SELPHI2.0 on all nine pre-registered, Holm-corrected tests; the margin over SELPHI2.0 on top-1% recall was modest (+0.059 [0.010, 0.120]). Macro AUROC did not drop from cross-validation to test, and on unseen kinase families with unseen proteins it was 0.752. Tyrosine kinases scored 0.595, barely above a site-popularity floor of 0.585, and fell from cross-validation. No site reached an assignment probability of 0.8. Part of the margin over Cantley reflects its weak test draw. Phosphoproteomics features added little, and for kinase-activity inference the predicted library was at parity with PhosphoSitePlus overall and worse on unseen drugs. AI agents ran the study; no human check of the code or analyses is recorded.

## Author summary

Protein kinases switch other proteins on and off by adding phosphate groups to specific serine, threonine or tyrosine residues. Mass spectrometry finds hundreds of thousands of such phosphorylated sites in human cells, but for most of them we do not know which of the roughly 500 kinases is responsible. Computer predictions can help, but they are often tested in ways that make them look better than they are. One reason is that very similar proteins can end up in both the training and the test data. We built a stricter test. Whole families of similar proteins were kept out of training, and the final test set was opened only once. On this test, our predictor ranks the correct kinase better than a widely used motif atlas and our own retrained versions of two published predictors. That conclusion holds for serine/threonine kinases. For tyrosine kinases, our method did barely better than a simple rule that favours well-studied sites, and the other methods did worse than that rule. The predictor also cannot yet say with high confidence which kinase acts on a given site. AI agents carried out this study under a human-written brief; every result should be checked before use.

---

## Introduction

About 95% of known human phosphosites have no assigned kinase, and annotated sites are concentrated on a few kinases. In PhosphoSitePlus (PSP) [1], the top 20 kinases account for about half of human kinase–substrate annotations, and 146 of 421 annotated kinases have fewer than five (repository `docs/lit/02_ks_annotation_sources.md`). Predictions of direct kinase–substrate (KS) relations are therefore used to interpret phosphoproteomics experiments, to choose candidate substrates for experiments, and to build the kinase–substrate libraries on which kinase-activity inference depends [2].

Two kinds of prior method dominate. Motif methods score a site against each kinase's intrinsic sequence preference. The Cantley lab measured these preferences with positional scanning peptide arrays for 303 Ser/Thr kinases [3] and 93 Tyr kinases [4], so its atlases are not trained on curated KS labels. Learned methods train on curated pairs. Examples are the transformer Phosformer [5] and the integrative random forest SELPHI2.0 from the Petsalaki lab [6], which combines motif, co-expression and co-regulation features. Papers on learned methods report high accuracy, but often under random splits over KS pairs. On Phosformer's released splits, 51% of test-positive peptides also occur as training positives under another kinase (repository `docs/LITERATURE.md` §1). The SELPHI2.0 paper reports a cross-validated AUC of 0.97, which falls to 0.75–0.83 on external data [6]. LinkPhinder drops from 0.97 on random splits to 0.53 on independent data [7]. Phylogeny-aware splits give more realistic values, such as AUROC 0.755 for S/T sites for the sequence-only KolossuS model [8].

Three biases make evaluation hard. First, 68% of KS pairs have another labelled site on the same protein, so a leave-site-out split leaks protein-level information (repository `docs/RESEARCH_TREE.md`). Second, sites labelled to other kinases are an attractive source of hard negatives, but many of them are substrates of the kinase's own paralogs. Third, metrics pooled over pairs are dominated by a few heavily studied kinases. PSP labels are also biased toward well-studied proteins and sites, so a model can score well by learning annotation popularity instead of kinase specificity.

The human-written brief asked for splits that group sites by protein rather than by KS pair, for benchmarks that leave out substrates or whole kinases, for evaluation that is not dominated by well-studied kinases, for a comparison of training regimes, and for attention to unorthodox kinases, such as those that need a priming phosphosite. We describe the benchmark that the agents built from these requests, the resulting predictor and its limits, and, because the study was also an experiment in whether AI agents can run such a project, the agent setup and what went wrong.

## Results

### Benchmark design

Positives are PSP pairs and there are no true negatives, so the task is positive–unlabelled (PU) learning [9] (Fig 1; details in Methods). Each of 503 kinases is scored against every site of matching residue class that was observed by mass spectrometry (U_ms). Of 11,609 substrate homology clusters, 2,353 (20.3%; 3,957 proteins, 2,205 test positives) were sealed as TEST, and DEV was split into five CV folds.

There are three tiers: T1, a seen kinase on unseen proteins (CV and TEST); T2, an unseen kinase family on seen proteins (CV); and T3, an unseen family on unseen proteins (TEST). The endpoints are macro averages over kinases. M1 (primary) is per-kinase AUROC against hard negatives (E-hard: sites labelled to kinases outside the kinase and its paralogs). M2 is site-centric Hit@5: for each kinase, the fraction of its labelled sites at which it ranks in the top five of the filtered candidate kinases, averaged over kinases. M3 is recall of positives within the kinase's top 1% of all candidate sites (E-all). M4 is M1 over all kinases with ≥1 test positive. The evaluable sets were frozen with the split: K_eval (≥3 TEST positives on ≥2 clusters; 119 kinases), K_any (≥1 TEST positive; 256) and K_cv (≥5 DEV positives on ≥2 clusters; 203). Pseudokinases need ≥3 labels.

**Fig 1. Method and evaluation design.** Top: inputs, the two models, late fusion, calibration and outputs. Bottom: splits, endpoints and the single sealed TEST batch. Labels, negatives and features feed both models. R3 uses the ESM-2 site windows and the kinase descriptors (ESM-2 embedding of the kinase domain, profile of Cantley scores, Manning group and family) rather than the full feature set. Asterisks mark label-derived features, which are cross-fitted by substrate cluster. The ledger blocks a second batch unless a logged override is used. File: `figures/Fig1_method.png` (PDF alongside).

### Cross-validation: the fused model leads on every endpoint

On pooled 5-fold DEV CV over 203 kinases, the late fusion of R2 (LightGBM) and R3 (neural model), each with five seeds, was best on every endpoint, and each component was weaker alone (Table 1; per-fold M1 0.767–0.794). Against Cantley arm (ii), on the pairs that both methods score, the paired differences were +0.063 [0.048, 0.074] for M1, +0.105 [0.082, 0.125] for M2 and +0.227 [0.189, 0.251] for M3. Fusion added ΔM1 +0.010 [0.005, 0.016] and ΔM2 +0.063 [0.051, 0.078] over the R2 ensemble.

**Table 1. Pooled 5-fold DEV cross-validation (K_cv = 203 kinases; A24 metrics).**

| Method | M1 | M2 Hit@5 | M3 | M4-CV |
|---|---|---|---|---|
| **Final: R2+R3 fusion** | **0.786** | **0.263** | **0.361** | **0.774** |
| R2, 5-seed ensemble | 0.776 | 0.200 | 0.312 | 0.767 |
| R3 (Cantley-distilled), 5 seeds | 0.754 | 0.223 | 0.294 | 0.731 |
| B7 simple-strong logistic regression | 0.726 | 0.196 | 0.179 | 0.708 |
| B1 Cantley, arm (i) / arm (ii) | 0.717 | 0.157 / 0.171 | 0.133 | 0.690 |
| B4 SELPHI2.0, our retrain (its own 34% of pairs) | 0.708 | 0.176 | 0.316 | 0.714 |
| B3 Phosformer, our retrain | 0.677 | 0.132 | 0.104 | 0.655 |
| B2 learned PSSM | 0.660 | 0.118 | 0.109 | 0.638 |
| B0-pop (site popularity floor) | 0.552 | 0.020 | 0.142 | 0.562 |

M4-CV is over the 339 kinases with ≥1 DEV positive. B1 arm (i) uses the published percentile and arm (ii) a percentile recalibrated on our site universe; they tie on M1 and M3. Sources: repository `docs/REPORT_internal.md` §5.1, §6.1; `results_summary/cv/LEADERBOARD.md`.

In the R2 ablation (Table 2), site-intrinsic structure and ESM-2 features nearly doubled M3, and the kinase descriptors added +0.014 M1 and raised M2. The learned PSSM lowered M2 until the kinase descriptors were added, probably because its per-kinase scales are not comparable across kinases. The phosphoproteomics families together added about +0.002 M1, inside the noise (S1 Text §2.4), which is a negative result for the modality that the brief emphasised. Two alarm runs showed no sign of circularity: removing the 18 compendium studies cited as PSP label sources left M1 unchanged (0.773), and a model trained on missingness patterns alone reached 0.522.

**Table 2. R2 ablation curve (seed 0, pooled CV, A24).**

| Step | M1 [95% CI] | M2 | M3 |
|---|---|---|---|
| F1 Cantley atlas | 0.723 [0.707, 0.737] | 0.136 | 0.125 |
| + F0 site structure and ESM-2 | 0.733 | 0.169 | 0.236 |
| + F2a learned PSSM | 0.738 | 0.094 | 0.239 |
| + F3 interaction screens | 0.748 | 0.143 | 0.260 |
| + F4 localisation | 0.753 | 0.141 | 0.269 |
| + F6a co-regulation | 0.755 | 0.149 | 0.274 |
| + F6b inhibitor response | 0.755 | 0.153 | 0.281 |
| + F6c CPTAC kinase state | 0.754 | 0.146 | 0.272 |
| + F5 co-expression | 0.759 | 0.160 | 0.283 |
| + F7 kinase descriptors (full) | 0.773 [0.755, 0.785] | 0.207 | 0.312 |

Source: repository `docs/model_r2.md` §5 and `results_summary/cv/LEADERBOARD.md` (A24 re-evaluation).

Six candidate improvements failed their pre-declared adoption rules and were not adopted (S1 Text §2–3): a fine-tuned ESM cross-encoder, a feature for similarity to training substrates, Tyr/SH2 docking features, priming-offset features, a prior-debiased site score and ESM-2 3B windows in R3.

### Retrained baselines are well below their published numbers

The brief asked for retrained baselines, so that all methods learn from the same splits, labels and negatives (Methods). As a fidelity check, our Phosformer code (B3), trained on Phosformer's own published training split (a 60% subsample, one epoch), gave AUROC 0.888 on its test split; the paper reports 0.885. This rules out a broken pipeline but not under-training on our folds. On our benchmark, B3 scored below the label-free Cantley atlas on every CV endpoint (M1 0.677 against 0.717). The released Phosformer weights reached 0.796 on the 2% of required pairs that they could score, but they were trained on PSP pairs that include our positives, so we report them only as a contaminated contrast. Our SELPHI2.0 port (B4) scored M1 0.708 on the 34% of pairs that its site features cover, against 0.728 for Cantley on the same pairs. We traced SELPHI2.0's published AUC of 0.97 mainly to random-pair negatives: on our cluster split, the same "positives versus random pairs" metric gave 0.947.

### TEST: the frozen model beats Cantley and our retrains of Phosformer and SELPHI2.0

The single TEST batch scored 12 frozen configurations over 25.9 million kinase–site pairs (Table 3, Fig 2). All nine pre-registered primary tests met the "beats" rule, which requires a Holm-adjusted p below 0.05 and a paired CI lower bound above 0 (Fig 2D). Eight had bootstrap p < 0.0002 (the resolution floor at B = 10,000) and Holm p ≤ 0.0018. Against Cantley, the M1 difference was +0.075 [0.046, 0.099], and an audit re-cut of the per-kinase tables showed that the final model was better on 97 of 116 matched kinases. The ninth test, M3 against B4, was modest: +0.059 [0.010, 0.120], Holm p 0.018, and in the same re-cut the final model was better on 52 kinases and B4 on 27 (40 ties). In the secondary family, the fusion beat both of its components on every endpoint (least significant: M2 over R3, +0.040 [0.007, 0.060], Holm p 0.013), so the fusion choice made on CV replicated.

**Table 3. TEST results (single batch, FREEZE 721dc197; K_eval = 119).** For standalone values, unscored pairs are ranked last; brackets are 95% cluster-bootstrap CIs (B = 1,000). Runtimes are wall-clock on a shared node and exclude feature building.

| Config | M1 | M2 Hit@5 | M3 | M4 | Train (s) | Score TEST pairs (s) |
|---|---|---|---|---|---|---|
| **Final (R2+R3)** | **0.794** [0.774, 0.811] | **0.363** [0.324, 0.397] | **0.363** [0.309, 0.401] | **0.799** | 1,136.6 | 316.9 |
| R2 alone | 0.778 [0.758, 0.795] | 0.264 [0.224, 0.293] | 0.312 [0.261, 0.348] | 0.781 | 417.9 | 298.9 |
| R3 alone | 0.775 [0.752, 0.793] | 0.322 [0.290, 0.362] | 0.292 [0.242, 0.330] | 0.769 | 718.7 | 3.1 |
| B7 simple-strong LR | 0.725 [0.700, 0.751] | 0.243 | 0.200 | 0.717 | 89.1 | 287.0 |
| B1 Cantley (i) / (ii) | 0.714 [0.685, 0.743] | 0.163 / 0.189 [0.155, 0.214] | 0.161 [0.129, 0.194] | 0.699 | 0 | 4.4 / 1.0 |
| B3 Phosformer, our retrain | 0.708 [0.686, 0.733] | 0.163 [0.145, 0.190] | 0.104 [0.077, 0.129] | 0.688 | 3,152.4* | 6,478.8 |
| B4 SELPHI2.0, our retrain | 0.705 [0.673, 0.733] | 0.233 [0.204, 0.266] | 0.304 [0.233, 0.351] | 0.701 | 737.0 | 2,043.4 |
| B2 learned PSSM | 0.685 [0.647, 0.717] | 0.155 | 0.138 | 0.664 | 3.3 | 39.7 |
| B0-pop (popularity floor) | 0.578 [0.543, 0.608] | 0.020 | 0.127 | 0.588 | 0 | 0.4 |
| B0-prior (training count) | 0.500 | 0.109 | 0.010 | 0.500 | 0 | 0.8 |

*B3 training time excludes about 4.3 h of masked-language-model pretraining. Scoring covers 25.9 M pairs, except for B1, which scores 386 of 503 kinases (19.2 M pairs). Sources: repository `results_summary/TEST_LEDGER.md`; `results_summary/test_721dc1977366/<config>_ci.tsv`.

**Fig 2. TEST endpoints and the nine pre-registered tests.** (A–C) Standalone M1, M2 and M3 with 95% cluster-bootstrap CIs (B = 1,000). The Cantley row uses arm (i) for M1 and M3 and arm (ii) for M2, the stronger arm for each endpoint. (D) Paired differences, final minus baseline, on the matched universe (B = 10,000) with Holm-adjusted p. A p-value of 0.0018 corresponds to the bootstrap floor of 0.0002 before adjustment. Orange marks the one modest test, the only one whose CI lower bound is below 0.02. Sources: repository `results_summary/test_721dc1977366/*_ci.tsv`, `superiority.tsv`. File: `figures/Fig2_test_endpoints.png`.

An independent audit agent reproduced every TEST number exactly (repository `docs/audit/audit_test.md`). Four of its caveats limit how the numbers can be read:

- *No CV→TEST optimism in M1.* On the 117 kinases in both K_cv and K_eval, CV M1 was 0.792 and TEST M1 0.793 (paired Δ +0.001 ± 0.012), consistent with the audited optimism bound of about 0.005 over the roughly 90 CV configurations.
- *The M2 rise is composition, not improvement.* Hit@5 rises steeply with a kinase's label count, and K_eval holds relatively more well-labelled kinases than K_cv. On the same 117 kinases, CV M2 was 0.348 against 0.363 on TEST, and the label-free B0-prior control moved the same way (0.064 to 0.109).
- *The margin over Cantley is partly a weak-baseline draw.* The margin grew from +0.063 on CV to +0.075 on TEST, but on the 115 kinases that Cantley scores in both, Cantley fell from 0.741 to 0.717 (−0.025 ± 0.013), almost entirely on Tyr kinases, while the final model was flat (−0.003).
- *B3 and B4 are our reimplementations* under this benchmark's split, labels and negatives; the claims apply to these retrains, not to the published systems as distributed. For B4, 66.6% of TEST pairs lacked SELPHI site features and used default values; a CV check put the cost at about 0.013 M1, which does not explain the gap.

### Generalisation to unseen kinase families

With a kinase's whole family held out of training (T2, on CV), R2 reached M1 0.751 and the fusion 0.749. On TEST (T3), with each kinase scored by the model trained without its family, the fusion reached 0.752, R2 0.744, B7 0.708, R3 0.703 and the learned PSSM 0.616. The fusion's T3 value is above Cantley's standalone T1 value of 0.714, but no paired T3 test was pre-registered, so this comparison is descriptive. T2 is not truly kinase-naive, because label-free Cantley percentiles of held-out kinases enter R3 pretraining. On CV T2, sites not labelled to any other kinase (55% of held-out positives) scored 0.703, against 0.757 in the matched CV control, whereas sites labelled to another kinase scored 0.798, against 0.799.

### Where the model works and where it does not

Fig 3 and Table 4 show TEST strata. The strata carry no CIs in the batch, so all statements here are descriptive.

**Fig 3. Per-kinase and stratified TEST performance.** (A) Per-kinase TEST AUROC on E-hard negatives against the number of DEV training labels (log scale), for the 119 K_eval kinases. Only 10 of 23 Tyr kinases score above their own site-popularity (B0-pop) AUROC. (B) TEST M1 within strata for the final model, Cantley arm (i) and B0-pop. Asterisks: primed and unprimed positives of the curated priming kinases. The comparison with Cantley here uses its non-priming-aware arm. Strata have no CIs. Sources: repository `results_summary/test_721dc1977366/{final,B0_pop}_per_kinase.tsv`, `{final,B1_i,B0_pop}_strata.tsv`. File: `figures/Fig3_kinases_strata.png`.

**Table 4. TEST M1 by stratum (macro over K_eval kinases in the stratum; descriptive).**

| Stratum | Kinases | Final | Cantley (i) | B3 | B4 | B0-pop |
|---|---|---|---|---|---|---|
| S/T kinases | 94 | 0.841 | 0.772 | 0.751 | 0.736 | 0.576 |
| **Tyr kinases** | 23 | **0.595** | 0.496 | 0.534 | 0.575 | **0.585** |
| 5–19 training labels | 13 | 0.759 | 0.646 | 0.610 | 0.694 | 0.639 |
| ≥100 training labels | 31 | 0.795 | 0.748 | 0.757 | 0.703 | 0.542 |
| Primed positives | 6 | 0.897 | 0.704 | 0.788 | 0.665 | 0.517 |
| Unprimed positives, priming kinases | 8 | 0.740 | 0.725 | 0.729 | 0.644 | 0.445 |
| Hub sites (≥5 non-paralogous kinases) | 85 | 0.863 | 0.712 | 0.731 | 0.802 | 0.788 |
| In vivo evidence only | 88 | 0.788 | 0.713 | 0.713 | 0.706 | 0.589 |
| In vitro evidence only | 96 | 0.827 | 0.730 | 0.715 | 0.751 | 0.615 |
| New since 2021 | 61 | 0.766 | 0.710 | 0.660 | 0.696 | 0.556 |
| External curated, tier A | 41 | 0.788 | 0.684 | 0.721 | 0.679 | 0.531 |

Sources: repository `results_summary/test_721dc1977366/<config>_strata.tsv`, `<config>_ext_tiers.tsv`; `docs/REPORT_internal.md` §7.5–7.6.

**Tyrosine kinases.** We found no demonstrated Tyr-kinase specificity. On TEST, Tyr kinases reached M1 0.595, barely above the site-popularity floor of 0.585; the model beat that floor on only 10 of 23 Tyr kinases, and Cantley was at chance (0.496). On unseen families (T3), Tyr kinases scored 0.531. Tyr-kinase M1 also fell from CV to TEST on the same kinases (0.667 to 0.595; paired Δ −0.072 ± 0.033), whereas S/T kinases rose (0.827 to 0.841). A CV diagnosis linked the Tyr gap of the Cantley atlas to macro dilution, promiscuous pTyr hub sites among the E-hard negatives and atlas matrices that anti-correlate with their own labels (repository `docs/features_tyr.md` §1). Neither hub sites nor evidence class explained the model's own Tyr gap on CV (S1 Text §2.2), and we have no tested explanation for it. The audit rated the Tyr result as not trustworthy as a skill claim.

**Understudied kinases.** On TEST, kinase-equal M1 held up across label bins (Table 4), and M4 was 0.807 for 16 kinases with no DEV labels. Site-centric ranking (M2) did not hold up. On CV, Hit@5 was 0.142 for kinases with 5–19 labels, below Cantley's 0.162, against 0.621 at ≥100 labels. The TEST value for that bin (0.208 against Cantley's 0.165) does not refute the CV finding, because TEST contains only 13 such kinases and they are the best-labelled of the bin.

**Priming kinases.** For the curated priming kinases (GSK3, CK1, CK2α), primed positives scored 0.897 on TEST, against 0.704 for Cantley arm (i), which ignores priming; on unprimed positives the two nearly tied (0.740 against 0.725). On CV, the priming-aware Cantley arm beat the model on primed positives (0.945 against 0.893), but that arm was not in FREEZE, so no TEST claim can be made against it.

**Popularity.** On hub sites, the final model scored 0.863, but the popularity floor alone reached 0.788. Gains on hub sites and on external tier B, where popularity is high (B0-pop 0.724), therefore come mostly from popularity. On pairs added to PSP since 2021, the margin over Cantley was smaller (+0.056) but still positive.

### Training regime and label variants

The brief asked whether one model for all kinases, one model per kinase, or something in between works best. Within the regimes tested, one global model was best on CV (Table 5), even on the 78 kinases with at least 20 labels in every fold (0.805 against 0.779), and partial pooling selected a shrinkage that is effectively the global model. Pooling mouse and rat labels with human labels raised M2 at no M1 cost. In-vivo-only labels were worse on every endpoint, but that label set is about three times smaller than the in-vitro set, so the comparison mainly measures label count. External curated pairs were not adopted: they gained +0.006 M1 on seed 0 but +0.001 over three seeds, and they had been pre-registered as ablation-only, to keep the external test stratum independent (S1 Text §1.3).

**Table 5. Regimes and label variants (R2, seed 0, pooled CV; Δ against the same-seed L2 global run).**

| Configuration | M1 | ΔM1 | M2 | ΔM2 | M3 |
|---|---|---|---|---|---|
| **L2 (human + mouse/rat), global (default)** | 0.768 | — | 0.192 | — | 0.298 |
| L1 human only | 0.771 | +0.003 | 0.163 | −0.029 | 0.295 |
| L2, in vitro evidence only | 0.765 | −0.003 | 0.145 | −0.047 | 0.276 |
| L2, in vivo evidence only | 0.739 | −0.029 | 0.103 | −0.090 | 0.225 |
| Per-kinase LightGBM (N ≥ 20, else global) | 0.756 | −0.012 | 0.174 | −0.019 | 0.259 |
| Per-kinase logistic regression | 0.755 | −0.013 | 0.171 | −0.021 | 0.261 |
| Partial pooling, ridge residual (λ = 300) | 0.768 | +0.000 | 0.199 | +0.007 | 0.298 |
| Per-Manning-group LightGBM | 0.761 | −0.007 | 0.170 | −0.022 | 0.298 |

Source: repository `docs/regimes_labels.md`; `docs/REPORT_internal.md` §4.4.

### Calibration and site assignment

The assignment probability is the calibrated probability that the top-ranked kinase is a true kinase of the site, given that the site is a substrate of a modelled kinase (Methods). On 1,458 labelled TEST sites, expected calibration error was 0.026 and top-1 accuracy 0.353 (Fig 4A), and no site reached 0.8. In an audit re-cut, the ≥0.5 tier (S/T sites only; n = 345) had a claimed probability of 0.615 and an observed precision of 0.603. No top-ranked Tyr kinase could exceed 0.445 (see Methods for an exception among lower-ranked kinases), and Tyr sites in the 0.2–0.5 tier were mildly overconfident (0.41 claimed, 0.34 observed). These are precisions on labelled sites and do not transfer to unlabelled sites. The ≥0.8 tier disappeared when an isotonic map, fitted separately for S/T and Tyr sites, fixed an overconfidence problem found on CV (S1 Text §4.4).

**Fig 4. Assignment reliability and coverage.** (A) Reliability of assignment_prob for the frozen DEV-trained model on labelled TEST sites, in bins of claimed probability (marker area ∝ sites per bin). The shaded region marks the tier that no site reaches. (B) Coverage counts of the deployment model (retrained on DEV+TEST) over the 350,736-site universe, at assignment probability ≥0.5 and ≥0.8, next to a Cantley-only and a random-score assignment. Coverage counts are model-based and are not precision. Cantley's ≥0.8 tier rests on an isotonic map fitted on very few labelled sites. The deployment model's calibrators were fitted on DEV and were not re-tested. Sources: repository `results_summary/test_721dc1977366/final_reliability.tsv`, `results_summary/coverage.tsv`. File: `figures/Fig4_assignment.png`.

### Coverage of the deployment model

For release, we refitted the frozen configuration on DEV and TEST labels (15,202 training pairs, including canaries) and scored all 350,736 observed sites against the 503 kinases of matching class (131,144,410 pairs). With the DEV-fitted calibrators, 29,886 sites (8.5%) had a top kinase at assignment probability ≥0.5, 26,762 of them without any PSP kinase, and none reached 0.8 (Fig 4B). A Cantley-only assignment covered a similar number of sites (30,258) and random scores covered 6, so coverage at 0.5 does not separate the model from the label-free atlas. Observed sites are about 19% of all S/T/Y residues in the proteome (repository `docs/features_structure_esm.md`); the others were not scored. The outputs are best used as ranked shortlists. Spot checks on three well-known sites showed hits and misses (TP53 S15: ATM at 0.70; AKT1 S473: MTOR third; EIF2S1 S52: PKC-family kinases first, although its training labels are EIF2AK1–4 and other kinases). All three lie on proteins with labels in the deployment training set, so these checks test recall of training labels, not prediction.

### Downstream use: kinase-activity inference

We tested the library predicted by a "clean" R2 variant (trained on DEV only; CV M1 0.771; Methods) in benchmarKIN [2], which evaluates kinase-activity inference on perturbation experiments [10] and CPTAC tumours; our harness reproduced the published PSP values (perturbation AUROC 0.782 against 0.79; tumours 0.664 against 0.662). On matched kinases, the library was at parity with PSP, not better (top-50: −0.013 [−0.043, 0.012] on perturbations, +0.011 [−0.027, 0.046] on tumours), and on drugs that the model had not seen it trailed PSP (top-50: −0.050 [−0.096, −0.001]). It was much better than a Cantley top-100 library (+0.176 [0.091, 0.263] on perturbations), made more kinases evaluable, and added slightly to PSP on tumours (+0.034 [0.005, 0.068]; repository `docs/downstream.md`). These comparisons were not in the pre-registered primary family, and their CIs are unadjusted.

### Runtime next to accuracy

The final model was the most accurate configuration, and it scored the TEST pairs about 20 times faster than B3 and 6 times faster than B4 (Table 3). These times exclude feature construction, which B3 does not need; cross-fitting alone costs about 21 CPU-minutes per fold, on top of the ESM-2 and structure features. The deployment run scored 131 million pairs in about 1,000 s, with a peak memory of 18.7 GB.

### Post-TEST exploration did not change the model

After the TEST batch, three exploratory experiments on DEV CV targeted M2 for understudied kinases, primed sites and Tyr specificity; the released model was not changed (S1 Text §5). None was a clean improvement. A Tyr-specific fusion weight and a nested motif repair failed. Pre-declared per-kinase intercepts lowered M2 for kinases with 5–19 labels; a balanced variant, added after the first result was seen, raised it only by taking top-5 slots from well-studied and unlabelled kinases. A priming-gated blend raised primed positives by +0.018, but a label-free flag for an observed priming phosphosite reproduced the gain, and against negatives with the same flag the gain was +0.002.

## Discussion

On a leakage-controlled benchmark, our late-fusion predictor beat the Cantley atlas and our retrained Phosformer and SELPHI2.0 on all nine pre-registered tests. The evidence is strongest for S/T kinases (TEST M1 0.841, with no CV→TEST optimism in M1). Under the benchmark's controls, both retrained learned baselines scored below the label-free atlas on M1. We did not ablate the controls themselves, but three observations point to them as the main reason why controlled numbers fall below published ones: label features that were not cross-fitted reached AUROC 0.93–1.00, against 0.59–0.76 when cross-fitted; 51% of the test-positive peptides on Phosformer's released splits occur as training positives (Introduction); and SELPHI2.0's random-pair metric stayed high on our split. Limited retraining budgets may also contribute (Limitation 5). We think this gap between published and controlled numbers is the most transferable result of the study.

The modelling results are modest. One global model with kinase descriptors beat per-kinase and per-group models, even for data-rich kinases; whether understudied kinases gain through the descriptors was not tested directly. Cantley distillation and late fusion both helped, although fusion did not help held-out families on CV. The phosphoproteomics modality, which the brief emphasised, added almost nothing to CV M1, possibly because co-regulation is confounded by protein abundance and study structure, because the inhibitor data cover only part of the kinome, or because effects at the kinase level do not resolve which of several sites a kinase phosphorylates directly. We did not test these explanations.

### Limitations

1. **Tyrosine kinases are not trusted.** Tyr-kinase M1 on TEST (0.595) was within 0.01 of a site-popularity floor, fell from CV, and was 0.531 on unseen families. Neither this model nor Cantley shows Tyr-kinase specificity here. Removing hub sites did not change the model's Tyr M1 on CV, so the cause is unknown, and the hub-excluded TEST stratum was not computed.
2. **There is no high-confidence assignment tier.** No site reached 0.8, and the ≥0.5 tier means about 60% top-1 precision on *labelled* S/T sites. The "which kinase?" output is a ranked shortlist.
3. **Site-centric ranking favours well-studied kinases.** The TEST rise in M2 reflects which kinases are evaluable, and on CV the model was below Cantley for kinases with 5–19 labels.
4. **Part of the margin over Cantley is a weak-baseline draw.** The CV margin (+0.063) is the more conservative estimate.
5. **The baselines are our retrains.** B3 used one seed, a fixed 1:16 negative ratio and at most three fine-tuning epochs; B4 used 10 of 100 bags, two substituted columns and default site features for 66.6% of TEST pairs. The original authors might do better on this benchmark. The priming-aware Cantley arm, which led on primed sites in CV, was not a TEST comparator.
6. **The labels are biased.** PSP reflects what has been studied, some labels, especially in-vivo-only ones, may reflect indirect relations (the label-count confound in the regime comparison prevents testing this), and there are no true negatives. All probabilities and precisions are lower bounds, and gains on hub sites are mostly popularity.
7. **One split, one batch, used early.** TEST was evaluated once, at about the halfway point of the 24 h budget, on one split. Strata carry no CIs, and there is no paired T3 test. Hub-excluded Tyr strata, calibration on the temporal and external sets and TEST coverage baselines were not computed and cannot be added without new data. Everything after the batch is exploratory.
8. **Some controls were approximate.** The negative-canary gate passed only when pooled (two of five folds exceeded it); a sampling difference between canaries and E-hard negatives is the likely but unproven cause.
9. **An agent-run study.** AI agents produced the design, code, analyses and text. Agent audits reproduced the headline numbers from saved predictions. No human re-derivation of the pipeline, code review or in-depth check of the biological plausibility of the predictions is recorded.

### What would establish real accuracy

The estimates are internally consistent, but they rest on one curated source and one split. Three steps would test them more strongly. The first is a prospective evaluation on KS pairs added to PSP or other curated resources after the 2025-05-19 release that we used: the model would be frozen now and scored later, with a pre-registered analysis that includes the Tyr, hub-excluded and calibration analyses missing here. The second is an evaluation on direct-substrate data that do not come from literature curation, such as kinase-specific in vitro or engineered-kinase screens, scored for S/T and Tyr kinases separately. The third is a re-run of the baselines by their authors, or with their full training budgets, on the same split. For the "which kinase?" use case, precision on unlabelled sites needs new experiments that test the top-ranked kinase for a random sample of unannotated sites.

## Materials and methods

### Data and site universe

Labels came from the PSP Kinase_Substrate_Dataset, release 2025-05-19 [1], with kinases keyed by UniProt accession [11]. We included 503 kinases (391 ST, 92 Y, 20 STY; 38 pseudokinases flagged) and dropped regulatory subunits. Residue class followed the Manning group [12] (Tyr kinase group → Y, else ST), with curated overrides; for example, DYRK, CLK, HIPK and GSK3 are ST for trans substrates, and MAP2Ks and TTK are STY. Labels outside a kinase's class were dropped, and sites were mapped to canonical sequences by 15-mer. The observed-site universe U_obs (350,736 sites) is the union of PSP sites, the provided compendium of perturbation phosphoproteomics, CPTAC [13] and Ochoa 2020 [14]; U_ms is the subset with ≥1 mass-spectrometry observation. Positives and negatives were both restricted to U_ms, and flags for whether the centre site was observed were not allowed as features, because 18% of human-native positives lacked mass-spectrometry evidence, against 0.4% of unlabelled sites, so observation status would otherwise act as a proxy for the label. Versions, checksums and licences of all datasets are in `DATA_SOURCES.md` and `docs/DATA_internal.md`.

### Labels

L2 labels are human-native PSP pairs plus mouse and rat pairs transferred to the human site with an identical 15-mer (18,638 pairs for 401 kinases; 14,747 after restriction to observed sites). We chose L2 over human-only labels because it raised M2 by 0.029 at no M1 cost (Table 5). Autophosphorylation (substrate accession = kinase accession) was kept in training and excluded from all primary endpoints. TEST positives are human-native PSP pairs only. The known-positive set KP(k) holds every pair of kinase k from any source and species, plus sites with an identical 15-mer; KP(k) sites are never negatives, and family-level UniProt notes and complex-subunit rows are used only as masks. A secondary external test used curated OmniPath [15], SIGNOR [16] and UniProt pairs on TEST proteins that are absent from PSP (PSP near-copies dropped; stricter tier A); these pairs were never used in training. A temporal stratum holds pairs absent from the 2021 PSP file bundled with kinase-library 1.8.0 [17].

### Negatives and the PU framing

For kinase k, candidate negatives are the U_ms sites of k's residue class, minus KP(k), minus a paralog mask and minus 15-mers within Hamming distance 2 of a KP(k) site. The paralog mask Par(k) covers kinases with ≥60% kinase-domain identity or the same Manning subfamily (same-family fallback: ≥50% identity); we added it because a design reviewer found that more than half of the SRC-family hard negatives were substrates of SRC paralogs. Two negative sets were pre-registered. E-hard, the primary set, contains sites with a human-native, non-autophosphorylation PSP label to a kinase outside k ∪ Par(k); these negatives test specificity among real substrates and are less sensitive to popularity than random sites. E-all contains every masked U_ms site in the evaluated split, scored exhaustively. Training negatives used the same masks, built from training-fold labels only.

### Splits and sealing

Substrate clusters are connected components of an mmseqs2 [18] all-against-all search (≥40% identity, ≥50% coverage) over proteins that carry a U_obs site, unioned with identical-15-mer edges (the largest of the 11,609 clusters holds 2.35% of labels). Grouping by protein and holding out whole kinases were requests of the human brief. We grouped by homology cluster rather than by single protein because protein-level features would otherwise carry homologs' labels across the split. We drew on the split designs of Phosformer-ST, which clustered kinases at 70% identity and separated peptides by Hamming distance [19], and of KolossuS, which uses a phylogeny-aware split [8]. TEST was chosen by simulated annealing (seed 0) to give each kinase with ≥5 labels about 20% of its positives (187 of 228 such kinases fell within 15–25%); DEV was split into five folds by cluster under the same constraints. For leave-kinase-out tiers, kinase clusters are components at ≥50% domain identity (KinCore alignments [20]) merged with Manning family (120 clusters, five folds). TEST labels were sealed with a recorded SHA-256; a test ensures that only the TEST harness reads them, and DEV and TEST share no proteins or clusters.

### Cross-fitting of label-derived features

Because most KS pairs have co-labelled sites on the same protein, every feature computed from labels was cross-fitted by substrate cluster. Inside each outer CV fold, training rows of such a feature came from inner 5-fold cross-fitting over training clusters, so a site's own protein and cluster never contributed to its own feature; eval rows used a randomly assigned inner-fold fit, so that their label coverage matched. One driver (`ksp.crossfit`) implements this rule for all such families and raises an error if a fitted label lies on a query protein or cluster. In a toy test with planted labels, 0 of 9,721 training positives were flagged with cross-fitting, against all of them without it. The same rule applies when the final model is trained on DEV and scored on TEST.

### Features

Features are kinase × site or kinase × protein matrices, loaded through a store that checks index hashes and gathered per kinase on demand.

- **F0, site structure and sequence:** pLDDT, exposure, disorder, solvent accessibility and secondary structure from AlphaFold DB v6 models [21,22] (structuremap exposure measures [23], re-implemented without the PAE adjustment); ESM-2 650M per-residue embeddings in full-protein context [24] (PCA-64 for R2).
- **F1, Cantley atlas [3,4]:** scores and percentiles against U_ms (the published Ochoa-reference percentile only in baseline arm B1(i)), in unprimed, soft-primed and hard-primed versions; the primed versions score flanking U_ms phosphosites as pS/pT/pY (soft versions weight them by evidence). 396 of 503 kinases have their own matrix; the others receive the nearest atlas matrix plus a distance feature.
- **F2a, learned PSSM (label-derived):** per-kinase ±7 log-odds with an empirical-Bayes family and group prior.
- **F3, interactions:** systematic screens only (BioPlex 3.0 [25], HuRI without its literature set [26], OpenCell [27], the Buljan kinase interactome [28] and the Human Cell Map [29]; 279,677 edges) as direct-edge flags and random-walk proximity, plus a label-derived proximity to the kinase's training substrates (F3-label). STRING [30] was excluded in all channels: 34% of PSP KS pairs have a STRING score, against about 4% of controls, and it imports kinase-assay records.
- **F4, localisation:** cosine similarity of compartment vectors from Human Protein Atlas v25 immunofluorescence [31], OpenCell, the Human Cell Map, Hein et al. [32] and the COMPARTMENTS experiments channel (not its knowledge channel); a combined cosine, an NMF cosine from the Human Cell Map and a profile correlation from Hein et al.; and UniProt topology.
- **F5, co-expression:** tissue RNA from HPA [33] and DepMap 26Q1 co-essentiality [34–36].
- **F6a, co-regulation (label-derived):** correlation of a site with the kinase's training substrates per study block of the compendium and per CPTAC cohort (top-3 and median Fisher z), excluding the same protein and cluster. CPTAC phosphosites were residualised on protein abundance and on stromal, immune and purity scores, following CoPheeMap [37], which showed that co-phosphorylation across tumours is dominated by abundance. Co-regulation with known substrates is also a SELPHI2.0 feature [6]; ours is cross-fitted.
- **F6b, inhibitor response:** ridge deconvolution of each site's response across 692 inhibitor conditions, mainly Hijazi et al. [38] and decryptM [39], on a drug × kinase occupancy matrix from Kinobeads [40,41], CurveCurator re-fits [42,43] and KINOMEscan affinities [44].
- **F6c, CPTAC kinase state:** correlation with the kinase's abundance and activation-loop phosphosites in 10 cohorts.
- **F7, kinase descriptors:** residue class, Manning group and family [12], a 16-component PCA of the ESM-2 embedding of the kinase domain, and 8 + 8 components (S/T and Y) of a PCA of the kinase's profile of Cantley scores.

Excluded because of measured or likely leakage (S1 Text §1.5): the iKiP-DB combined tables [45] (61% of PSP pairs), the Ochoa functional score [14] (trained on PSP sites), PSP regulatory-site and literature-count fields, study-bias proxies (HPA reliability, number of localisation sources), CoPheeMap/CoPheeKSA and other predictors' outputs, and the COMPARTMENTS knowledge channel (ablation only). A test fails any source that overlaps PSP KS pairs at more than three times chance. The in vitro lysate screen of Sugiyama et al. [46] was planned as a feature but never built, and was not used as labels because it hurts activity inference [2].

### Priming-dependent and atypical kinases

The priming kinases were fixed before any metric was seen: GSK3A/B (priming at +4), the CK1 family (−3), CK2α (pS in place of acidic residues) and PLK1. There is no separate model for them; priming enters through the primed F1 scores, the phospho channel of R3 and the soft-primed distillation targets. Explicit priming-offset features (A32a) failed their adoption rule (S1 Text §2.3). PIKKs (S/T-Q motif) and FAM20C (S-x-E, secretory; 87.5% of its sites are lumenal) rely on motif, localisation and topology features; they are reported as strata, but they are small (4 PIKKs and 1 FAM20C kinase on TEST).

### Why two models

R2 integrates many weak, heterogeneous context features for each kinase–site pair, but its scores are only partly comparable across kinases. R3 learns one embedding space for site windows and kinase descriptors, so its scores are comparable across kinases by construction, and it can score kinases with few labels through their descriptors. Fusing the z-scored outputs combines the two. This was the design rationale; we did not test each part of it separately.

### R2: global LightGBM pair model

R2 is one LightGBM model [47] over (kinase, site) rows with feature groups F0–F7, so understudied kinases share trees with well-studied ones through the kinase descriptors. Settings: binary objective, 127 leaves, learning rate 0.03, min_data_in_leaf 100, feature_fraction 0.8, bagging fraction 0.8, λ_L2 1.0, max_bin 255, at most 3,000 rounds, early-stopping patience 100. Early stopping used the macro per-kinase AUROC of positives against sampled hard and E-all-type negatives in an inner validation split (20% of training clusters, real kinases only); the model was then refitted on all training rows with the early-stopped number of rounds. Sample weights were n_k^−0.5, normalised (n_k, the kinase's number of positives), which tempers the dominance of well-studied kinases without fully equalising them. Each kinase had 50 negatives per positive, half of them E-hard type, with E-all-type sites filling the remainder where the hard pool was too small. A staged search on fold 0 (19 configurations over tree size, weight exponent α ∈ {0, 0.5, 1}, negative ratio, hard fraction and objective; fold-0 noise about ±0.01) kept the base setting of 1:50 with 50% hard negatives (fold-0 M1 0.754), because no alternative beat it by more than the noise (best alternative 0.751); ranking objectives did not beat the binary loss (S1 Text §3.4). The final R2 is a five-seed ensemble (seeds 0–4) averaged in logit space, trained on 623,373 DEV rows (12,223 positives) with 172 features; the refitted round counts for seeds 0–4 were 477, 200, 180, 392 and 303. OOF and TEST predictions are exhaustive over all candidate pairs (about 20 million per CV fold).

### R3: neural kinase × site model with Cantley distillation

The site encoder is a two-layer transformer over frozen ESM-2 650M embeddings of the ±7 window in full-protein context, an amino-acid embedding, a phospho-channel bit that marks flanking U_ms phosphosites (for priming; the centre bit is always 0) and a protein-mean token. The kinase encoder combines four blocks, each with descriptor dropout 0.2: the ESM-2 embedding of the kinase domain; 32 principal components per residue type (S/T, Y) of the kinase's unprimed Cantley scores on 3,000 S/T and 1,000 Y reference sites; Manning group; and Manning family. Groups and families without a training kinase in a fold are zeroed. A bilinear-plus-MLP head scores each pair; the bilinear compatibility between site and kinase embeddings is similar to zero-shot designs such as DeepKinZero [48]. The loss is a balanced binary cross-entropy over all mask-filtered kinase × site entries in each batch of 1,024 sites (25% carrying a training positive; positives weighted by n_k^−0.5), plus 0.5 × an InfoNCE term [49] over kinases for each positive site, with other true kinases and paralogs filtered out. Settings: site dimension 192, hidden dimension 128, kinase dimension 256, 4 attention heads, dropout 0.2, ESM-input dropout 0.3, AdamW (learning rate 3 × 10⁻⁴, weight decay 0.05), 200 warm-up steps, at most 6,000 steps, early stopping with patience 10 evaluations (every 100 steps) on the macro hard-negative AUROC of an inner validation split. The best checkpoint is restored and R3 is not refitted, so the 20% of training clusters held out for validation are never used to train R3, whereas R2 is refitted on them.

R3 was first pretrained for 3,000 steps (learning rate 10⁻³, batch 1,024) to predict label-free Cantley percentiles against the site universe (unprimed and soft-primed) on training-fold sites, for the 396 kinases with their own matrix. During fine-tuning, an auxiliary distillation loss of weight 0.2 was kept. The design reviewer in the external-PI role proposed the idea of pretraining on atlas percentiles. It follows Phosformer-ST [19], which is trained on Cantley percentiles, and KolossuS [8], which pretrains on the atlases. Distillation raised R3's M1 by +0.019 [0.010, 0.033] (pre-A24 harness) and lifted the positive canary from 0.83 to 0.93. Five seeds were trained.

### Late fusion

R3 was combined with R2 by late fusion rather than stacked into R2 as a feature: nested stacking would have needed 36 neural trainings with inner folds, against 6 for late fusion, and late fusion keeps R3's out-of-fold predictions out of R2's training rows. Each model's logits were z-scored within the prediction set (per CV fold, or over the TEST candidates). This puts both models on one scale and so makes the fused score comparable across kinases, which the site-centric endpoint needs. The fused score is 0.6 · z(R2) + 0.4 · z(R3), with the weight chosen on pooled CV M1 over w ∈ {0, 0.3, 0.4, …, 0.8, 1} (0 and 1 are R3 and the R2 ensemble alone). M1 was 0.754 at 0, 0.780 at 0.3, 0.786 at 0.5 and 0.6, 0.785 at 0.7 and 0.776 at 1. Nested leave-one-fold-out selection chose 0.6 in every fold, so the M1 selection optimism was 0. The site-centric score (score_site) uses a weight of 0.5, chosen on CV M2.

### Calibration and the meaning of assignment_prob

Calibrators were fitted on DEV out-of-fold predictions, cross-fitted by fold: per-kinase Platt scaling [50] for kinases with ≥30 OOF positives (61), a family slope with a per-kinase intercept (261), or a global slope (17). Out of sample, expected and observed positives agreed within bins of training-label count (for example 5,739 against 5,400 at ≥100 labels). Under PU labelling these probabilities estimate the chance that a pair is a curated PSP positive, which is a lower bound on the chance that it is a true substrate.

The assignment probability of site s is a softmax over the calibrated logits of all kinases of matching class (temperature 0.93, fitted by OOF log-likelihood), followed by an isotonic map [51] of the top kinase's probability, fitted separately for S/T and Tyr sites (A31b, added after an audit found the ≥0.8 tier overconfident). It is conditional on s being a substrate of a modelled kinase, estimates top-1 precision on PSP-like labelled sites, and gives no information on whether an unlabelled site is a substrate. M2 uses the raw fused score_site (CV M2 0.263; nested 0.259) rather than the calibrated assignment probability (0.205); we flag that this choice was made late, after both values were seen (A31a). The remaining probability mass is shared among the other kinases in proportion to their softmax values. In the released site table, the kinase with the highest `assignment_prob` is therefore not always the kinase with the highest score. This happens at 25,532 of 350,736 sites (7.3%). At 19 Tyr sites (INSR/IGF1R paralog pairs), a non-top kinase has `assignment_prob` 0.55. The 0.445 ceiling applies to the top kinase only. We found this while preparing the release and did not change the released files; the `prob` and `score_site` columns give the model's own ordering.

### Baselines

B1 is the Cantley atlas, not retrained because it is not fitted on KS labels, and computed with our reimplementation of kinase-library 1.8.0 [17] (agreement within 10⁻³) in three arms: (i) the published percentile against the Ochoa reference, (ii) a percentile recalibrated on U_ms, and (iii) a priming-aware score. Arms (i) and (ii) tie on M1, and arm (ii) was better on M2 (+0.0135 [0.007, 0.020]), so we compared against arm (i) on M1 and M3 and arm (ii) on M2 (A23), which strengthens the baseline. Ten dual-specificity kinases have a matrix for only one residue type (MAP2K1/2/4/5/6/7, TTK, TESK1, NEK10, PKMYT1); on TEST they were left unscored in B1 and dropped from matched comparisons, which can only help the baseline (A30). B2, the retrained motif analogue, is the learned per-kinase PSSM alone; B7 is a logistic regression on Cantley, the learned PSSM and minimal context.

B3 reimplements Phosformer [5] (6-layer, 12-head RoBERTa [52], d = 768; kinase domain + 11-mer; focal loss [53]). It was pretrained by masked language modelling for 50 epochs on a corpus from which 642 TEST substrates, 458 TEST homologs and 152 near-identical kinase domains had been removed. It was then fine-tuned on our folds with focal γ = 2, learning rate 2 × 10⁻⁵ and batch 64, at the paper's 1:16 negative ratio with 50% hard negatives (seed 0). The pre-registered plan was to choose between 1:16 and 1:50 on inner validation; the PI instead fixed 1:16 with the instruction "no hyperparameter changes after seeing fold metrics", after the fold-0 1:16 run had finished and before any 1:50 result existed. B4 ports SELPHI2.0 [6] with its random forest and full 72-configuration grid, but with 10 instead of 100 bags. Of its 44 model columns, 39 were regenerated, 3 were regenerated with a method substitution and 2 were substituted; label-derived columns were cross-fitted. The floors are the kinase's training count (B0-prior), site popularity (B0-pop) and LightGBM on feature missingness (B0-miss).

### Evaluation, metrics and statistics

All metrics come from one frozen harness (`ksp.eval`, commit 974d520). Kinases are weighted equally because annotations are concentrated on a few kinases (Introduction), so pooled metrics would mostly measure those kinases. Because pooled out-of-fold AUROC is distorted when score offsets differ between folds (B0-prior scored M1 0.40 instead of 0.50; S1 Text §4.1), CV M1, M3 and M4 are computed on scores rank-normalised within each fold and kinase, for all methods (A24); TEST has a single model and is unaffected. M2 uses filtered ranks, which remove the site's other true kinases and unlabelled paralogs of the true kinase. Head-to-heads use only pairs that both methods score, and for M2 only kinases that both score (A27).

CIs come from a bootstrap over substrate clusters (B = 1,000; B = 10,000 for TEST primary tests), with paired draws and a kinase-resampling sensitivity analysis. Exactly nine primary tests were pre-registered (final against B1, B3 and B4 on M1–M3), with Holm correction [54] over the nine; "beats" requires Holm p < 0.05 and a paired CI lower bound above 0. A secondary family (against R2, R3, B7, B2, and M4 against B1) was corrected separately. Twelve negative canary kinases (random real DEV-labelled sites, clustered by protein) and one positive canary (planted motif and interaction edge) ran through the whole pipeline, with gates of ≤0.57 on E-hard for negative canaries and ≥0.90 on E-hard and E-all for the positive canary; the final model passed when pooled (0.539, 1.000). Selection used the pre-declared rule (best CV M1, M2 as tie-break); about 90 configurations were evaluated on CV, about 55 of them selectable.

### The single TEST batch and the FREEZE/ledger mechanism

All TEST deliverables (12 configurations, 17 prediction files including T3) were pinned by SHA-256 in FREEZE.json, which also lists the nine primary tests and was committed before the batch. The harness wrote a STARTED block to an append-only ledger before reading the sealed labels and a COMPLETED block afterwards. It required an explicit `--real` flag, B = 10,000 and the commit of the frozen harness, and it blocks further batches once a COMPLETED block exists unless a logged override is used. All strata, kinase-resampling CIs and reliability tables had to be computed inside that batch (A29). The batch ran on 2026-10-27 from 00:10 to 00:57 (2,838 s); an audit agent then confirmed the pins and the DEV/TEST disjointness. TEST labels were afterwards released only for the deployment refit, and further experiments were restricted to DEV CV and marked exploratory (A34).

### Downstream evaluation

We used benchmarKIN [2] (Zenodo 12566560), with the perturbation compilation of Hernandez-Armenta et al. [10], Hijazi et al. [38], a tyrosine set and seven tumour cohorts from CPTAC. Libraries of the top N predicted substrates per kinase (N fixed in advance as the PSP size, 50 or 100) were compared with PSP, Cantley and shuffled controls, using z-score inference and matched kinases. The clean R2 variant excluded the 22 benchmark studies and all CPTAC features. Membership of U_ms still depends on the compendium and CPTAC, which we disclose as a residual channel.

### Seeds, runtime and reproducibility

Seeds came from the list [0–4] fixed in the shared agent brief: seed 0 for the split, annealing, ablations, B3, B4, regimes (external-label variant also seeds 1–2) and the bootstrap, and seeds 0–4 for the final R2 and R3. Downstream analyses used the benchmarKIN and decoupler seeds 42 and 1. Every run manifest records the seeds used, git hash, configuration, input checksums, runtime and fallbacks. Runtimes are wall-clock on a shared node (64 cores, 375 GB RAM, four RTX 4090 GPUs; jobs limited to a 128 GiB cgroup), where contention inflated GPU times about two-fold. The end-to-end script promised in the design was not written; the commands are listed in `docs/REPORT_internal.md` §12.

### Agent and sub-agent setup

The study ran in a 24 h budget (2026-10-26 13:20 to 2026-10-27 13:20 CET) under a human-written brief (repository `PROMPT.md`). The brief asked a lead agent to act as PI, delegate to sub-agents, keep a research record, have the design reviewed in up to three rounds, and rely on no result until an agent that did not produce it had checked it (Fig 5).

**Fig 5. Agent setup.** The human researcher supplied the brief, data and PhosphoSitePlus access. The lead agent ('PI') delegated to sub-agents with assigned paths and compute budgets; the model run by each agent is not recorded in the project files. Independent agents audited results that they did not produce, and code-level guards enforced cross-fitting, index consistency and the single TEST read. The bottom panel lists the main process failures found by the outside-view retrospective. Sources: repository `docs/DELEGATION.md`, `docs/AGENT_SETUP.md`, `docs/retro/final_retro.md`, `docs/LESSONS.md`. File: `figures/Fig5_agents.png`.

**Roles.** The lead agent wrote the design and its 34 amendments, kept the research tree, lessons file and delegation register, decided adoptions by pre-declared rules, and was the only agent that ran git. The brief made Claude Opus 5.5 and Sonnet 5.5 available; which model each agent ran is not recorded in the project files, and the model names in the statement at the top come from the publication brief. Four agents wrote the literature reviews. Round 1 of design review had seven reviewers (an external PI and six specialists), and about 80 of their 84 requested changes were accepted; round 2 had four reviewers. Round 3 was skipped, because round 2 found no critical issues at the design level and the PI wanted to avoid rework while six implementers were running. Each implementation package had one owner agent that wrote only in its assigned paths (`docs/DELEGATION.md`). Agents received a task brief (not archived) and a shared agent brief (`docs/AGENT_BRIEF.md`) with GPU and core budgets; a memory cap of 60 GB per agent was added at 17:55, after the first out-of-memory kill. Agents returned summaries of at most 600 words, which the lead agent was asked to read instead of raw logs. The total number of agent sessions, how they were launched, how much context they shared and their token cost are not recorded.

**Review by agents that did not produce the result.** Separate agents audited CV results twice, the code of the TEST harness, the TEST batch, the post-TEST work and the report's numbers, recomputing headline numbers from saved predictions with read-only data access. They found the pooled-metric defect (A24), stale cross-fitted stores, the overconfident ≥0.8 tier (A31b), a single-batch guard that could be bypassed (A29) and a post-TEST variant added after its first result was seen. Milestone retrospectives were written at three points (author not recorded); a fresh agent with no prior involvement wrote the final outside-view retrospective, from which the lessons were pruned to six rules (`docs/LESSONS.md`).

**The human's role.** The human wrote the brief, supplied the data, the two baseline papers and PhosphoSitePlus access, and resolved two of the data requests that agents logged (the Kinobeads tables and further PhosphoSitePlus files); ten other requests stayed open. The brief also supplied the main ideas for the split and the study design: grouping by protein, leave-substrate-out against leave-kinase-out benchmarks, the comparison of training regimes, species pooling, in vitro against in vivo labels, and attention to priming kinases. Whether the human checked any code or analysis is not recorded.

**What went wrong** (repository `docs/retro/final_retro.md`; full list in S1 Text §6):
- The final configuration existed at T+6.0 h; no later CV run improved M1, although 62 of 135 timed runs started after 20:00 (best 0.774) and compute was not the constraint (11.8 h of serial wall time for all CV runs).
- TEST was spent early: FREEZE was committed at 00:04, 5.3 h before its gate, and the batch completed at 00:57 (T+11.6 h of 24), with no recorded rationale. The retrospective calls this contrary to the brief, which asked for the best outcome rather than the fastest completion. No second sealed shard had been reserved, so post-TEST experiments could not be adopted.
- The protocol was amended while results arrived: 12 of 34 amendments (A23–A34) were made in the last 8.6 h before and around the freeze, six of them within 40 minutes, and A34 after TEST. Most strengthened baselines or fixed errors for all methods, but A31a changed the M2 score column late, in the favourable direction, after both values were seen.
- The priming-aware Cantley arm was known to beat the model on primed sites about three hours before the freeze but was not frozen; no reason is recorded.
- Also: fixes from the review of the TEST harness were applied 34 minutes before the freeze; the required regime comparison was delegated at T+6.6 h and used 38% of CV compute; there were two out-of-memory kills; the delegation status table was never kept current; and about 6 h passed with no file written after the 03:17 report draft.

The same retrospective lists as strengths the standing auditor, the sealed set with SHA pins and ledger, pre-registered gates that were honoured when they failed (six negative results accepted) and a report that leads with its caveats. Its main finding is that rules built into code (the cross-fitting error, the ledger, the pins) changed behaviour, whereas rules kept as prose were often ignored.

## Data and code availability

Code, the frozen evaluation harness, run manifests, CV and TEST metric files, figure scripts and the derived prediction files will be available at [repository URL]. The figures in this manuscript are generated by `manuscript/figures/make_figures.py` from committed summary files only. PhosphoSitePlus data are not redistributed; they are available from PhosphoSitePlus under its licence. Released predictions are derived from PhosphoSitePlus and are therefore distributed under CC BY-NC-SA (non-commercial). They include `predictions/kinase_top200.tsv.gz` (top 200 sites per kinase) and `predictions/site_top5.tsv.gz` (top five kinases per site with assignment_prob), with SHA-256 pins and limitations in `predictions/RELEASE_MANIFEST.json`. Metric summaries of CV and the TEST batch are in `results_summary/`, and the project's design, reports, audits and retrospectives are in `docs/`. Sources, versions, checksums and licences of all inputs are listed in `DATA_SOURCES.md`.

## Acknowledgments and statement on the use of AI

**This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.**

AI agents performed the literature review, design, design review, implementation, analysis, auditing and writing, including this manuscript. No human verification of code, analyses or text is recorded; the human supplied inputs. The publishing researcher is responsible for checking the results before use. We thank the developers of the public resources listed in the references and in `DATA_SOURCES.md`.

## Funding

[Funding statement to be added by the publishing researcher.]

## Competing interests

[Competing-interests statement to be added by the publishing researcher.]

## Author contributions

[CRediT roles of the human author(s) to be added by the publishing researcher.] Conceptualization of the task (the brief), provision of data and PhosphoSitePlus access: the human researcher. The brief also proposed protein- and sequence-clustered splits, leave-substrate-out against leave-kinase-out benchmarks, the comparison of training regimes, species pooling, in vitro against in vivo labels and the treatment of priming kinases. Methodology in detail, software, formal analysis, investigation, validation (independent agent audits), visualisation and writing: AI agents (a lead agent and sub-agents; see Methods).

## References

1. Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E. PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015;43(D1):D512–D520. doi:10.1093/nar/gku1267
2. Müller-Dott S, Jaehnig EJ, Munchic KP, Jiang W, Yaron-Barir TM, Savage SR, et al. Comprehensive evaluation of phosphoproteomic-based kinase activity inference. Nat Commun. 2025;16:4771. doi:10.1038/s41467-025-59779-y
3. Johnson JL, Yaron TM, Huntsman EM, Kerelsky A, Song J, Regev A, et al. An atlas of substrate specificities for the human serine/threonine kinome. Nature. 2023;613(7945):759–766. doi:10.1038/s41586-022-05575-3
4. Yaron-Barir TM, Joughin BA, Huntsman EM, Kerelsky A, Cizin DM, Cohen BM, et al. The intrinsic substrate specificity of the human tyrosine kinome. Nature. 2024;629(8014):1174–1181. doi:10.1038/s41586-024-07407-y
5. Zhou Z, Yeung W, Gravel N, Salcedo M, Soleymani S, Li S, et al. Phosformer: an explainable transformer model for protein kinase-specific phosphorylation predictions. Bioinformatics. 2023;39(2):btad046. doi:10.1093/bioinformatics/btad046
6. Maier BD, Petursson B, Lussana A, Petsalaki E. Data-driven extraction of human kinase-substrate relationships from omics datasets. Mol Cell Proteomics. 2025;24(8):100994. doi:10.1016/j.mcpro.2025.100994
7. Nováček V, McGauran G, Matallanas D, Vallejo Blanco A, Conca P, Muñoz E, et al. Accurate prediction of kinase-substrate networks using knowledge graphs. PLoS Comput Biol. 2020;16(12):e1007578. doi:10.1371/journal.pcbi.1007578
8. Jha K, Shonai D, Parekh A, Uezu A, Fujiyama T, Yamamoto H, et al. Deep learning-coupled proximity proteomics to deconvolve kinase signaling in vivo. bioRxiv [Preprint]. 2025. doi:10.1101/2025.04.27.650849
9. Elkan C, Noto K. Learning classifiers from only positive and unlabeled data. In: Proceedings of the 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD 2008). New York: ACM; 2008. p. 213–220. doi:10.1145/1401890.1401920
10. Hernandez-Armenta C, Ochoa D, Gonçalves E, Saez-Rodriguez J, Beltrao P. Benchmarking substrate-based kinase activity inference using phosphoproteomic data. Bioinformatics. 2017;33(12):1845–1851. doi:10.1093/bioinformatics/btx082
11. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025;53(D1):D609–D617. doi:10.1093/nar/gkae1010 (data: release 2026_03, human reference proteome UP000005640)
12. Manning G, Whyte DB, Martinez R, Hunter T, Sudarsanam S. The protein kinase complement of the human genome. Science. 2002;298(5600):1912–1934. doi:10.1126/science.1075762
13. Li Y, Dou Y, Da Veiga Leprevost F, Geffen Y, Calinawan AP, Aguet F, et al. Proteogenomic data and resources for pan-cancer analysis. Cancer Cell. 2023;41(8):1397–1406. doi:10.1016/j.ccell.2023.06.009
14. Ochoa D, Jarnuczak AF, Viéitez C, Gehre M, Soucheray M, Mateus A, et al. The functional landscape of the human phosphoproteome. Nat Biotechnol. 2020;38(3):365–373. doi:10.1038/s41587-019-0344-3
15. Türei D, Korcsmáros T, Saez-Rodriguez J. OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat Methods. 2016;13(12):966–967. doi:10.1038/nmeth.4077
16. Lo Surdo P, Iannuccelli M, Karis K, Meo E, Omidi P, Tosoni M, et al. SIGNOR 4.0: the 2025 update with focus on phosphorylation data. Nucleic Acids Res. 2026;54(D1):D682–D690. doi:10.1093/nar/gkaf1237
17. The Kinase Library. kinase-library: Python package, version 1.8.0 [software]. Available from: https://github.com/TheKinaseLibrary/kinase-library (documentation: https://kinase-library.phosphosite.org)
18. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol. 2017;35(11):1026–1028. doi:10.1038/nbt.3988
19. Zhou Z, Yeung W, Soleymani S, Gravel N, Salcedo M, Li S, et al. Using explainable machine learning to uncover the kinase–substrate interaction landscape. Bioinformatics. 2024;40(2):btae033. doi:10.1093/bioinformatics/btae033
20. Modi V, Dunbrack RL. Kincore: a web resource for structural classification of protein kinases and their inhibitors. Nucleic Acids Res. 2022;50(D1):D654–D664. doi:10.1093/nar/gkab920 (resource: https://dunbrack.fccc.edu/kincore/)
21. Varadi M, Bertoni D, Magana P, Paramval U, Pidruchna I, Radhakrishnan M, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Res. 2024;52(D1):D368–D375. doi:10.1093/nar/gkad1011
22. Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021;596(7873):583–589. doi:10.1038/s41586-021-03819-2
23. Bludau I, Willems S, Zeng WF, Strauss MT, Hansen FM, Tanzer MC, et al. The structural context of posttranslational modifications at a proteome-wide scale. PLoS Biol. 2022;20(5):e3001636. doi:10.1371/journal.pbio.3001636
24. Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023;379(6637):1123–1130. doi:10.1126/science.ade2574
25. Huttlin EL, Bruckner RJ, Navarrete-Perea J, Cannon JR, Baltier K, Gebreab F, et al. Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. Cell. 2021;184(11):3022–3040.e28. doi:10.1016/j.cell.2021.04.011
26. Luck K, Kim DK, Lambourne L, Spirohn K, Begg BE, Bian W, et al. A reference map of the human binary protein interactome. Nature. 2020;580(7803):402–408. doi:10.1038/s41586-020-2188-x
27. Cho NH, Cheveralls KC, Brunner AD, Kim K, Michaelis AC, Raghavan P, et al. OpenCell: endogenous tagging for the cartography of human cellular organization. Science. 2022;375(6585):eabi6983. doi:10.1126/science.abi6983
28. Buljan M, Ciuffa R, van Drogen A, Vichalkovski A, Mehnert M, Rosenberger G, et al. Kinase interaction network expands functional and disease roles of human kinases. Mol Cell. 2020;79(3):504–520.e9. doi:10.1016/j.molcel.2020.07.001
29. Go CD, Knight JDR, Rajasekharan A, Rathod B, Hesketh GG, Abe KT, et al. A proximity-dependent biotinylation map of a human cell. Nature. 2021;595(7865):120–124. doi:10.1038/s41586-021-03592-2
30. Szklarczyk D, Nastou K, Koutrouli M, Kirsch R, Mehryary F, Hachilif R, et al. The STRING database in 2025: protein networks with directionality of regulation. Nucleic Acids Res. 2025;53(D1):D730–D737. doi:10.1093/nar/gkae1113
31. Thul PJ, Åkesson L, Wiking M, Mahdessian D, Geladaki A, Ait Blal H, et al. A subcellular map of the human proteome. Science. 2017;356(6340):eaal3321. doi:10.1126/science.aal3321 (data: Human Protein Atlas v25.1, https://www.proteinatlas.org)
32. Hein MY, Peng D, Todorova V, McCarthy F, Kim K, Liu C, et al. Global organelle profiling reveals subcellular localization and remodeling at proteome scale. Cell. 2025;188(4):1137–1155.e20. doi:10.1016/j.cell.2024.11.028
33. Uhlén M, Fagerberg L, Hallström BM, Lindskog C, Oksvold P, Mardinoglu A, et al. Tissue-based map of the human proteome. Science. 2015;347(6220):1260419. doi:10.1126/science.1260419
34. Yun Y. Chronos parameters (Public 26Q1) [dataset]. DepMap, Broad Institute. figshare; 2026. doi:10.6084/m9.figshare.31660582.v2
35. Tsherniak A, Vazquez F, Montgomery PG, Weir BA, Kryukov G, Cowley GS, et al. Defining a cancer dependency map. Cell. 2017;170(3):564–576.e16. doi:10.1016/j.cell.2017.06.010
36. Dempster JM, Boyle I, Vazquez F, Root DE, Boehm JS, Hahn WC, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol. 2021;22:343. doi:10.1186/s13059-021-02540-7
37. Jiang W, Jaehnig EJ, Liao Y, Shi Z, Yaron-Barir TM, Johnson JL, et al. Deciphering the dark cancer phosphoproteome using machine-learned co-regulation of phosphosites. Nat Commun. 2025;16:2766. doi:10.1038/s41467-025-57993-2
38. Hijazi M, Smith R, Rajeeve V, Bessant C, Cutillas PR. Reconstructing kinase network topologies from phosphoproteomics data reveals cancer-associated rewiring. Nat Biotechnol. 2020;38(4):493–502. doi:10.1038/s41587-019-0391-9
39. Zecha J, Bayer FP, Wiechmann S, Woortman J, Berner N, Müller J, et al. Decrypting drug actions and protein modifications by dose- and time-resolved proteomics. Science. 2023;380(6640):93–101. doi:10.1126/science.ade3925
40. Klaeger S, Heinzlmeir S, Wilhelm M, Polzer H, Vick B, Koenig PA, et al. The target landscape of clinical kinase drugs. Science. 2017;358(6367):eaan4368. doi:10.1126/science.aan4368
41. Reinecke M, Brear P, Vornholz L, Berger BT, Seefried F, Wilhelm S, et al. Chemical proteomics reveals the target landscape of 1,000 kinase inhibitors. Nat Chem Biol. 2024;20(5):577–585. doi:10.1038/s41589-023-01459-3
42. Bayer FP, Gander M, Kuster B, The M. CurveCurator: a recalibrated F-statistic to assess, classify, and explore significance of dose–response curves. Nat Commun. 2023;14:7902. doi:10.1038/s41467-023-43696-z
43. Bayer FP, The M. CurveCurator: a recalibrated F-statistic to assess, classify, and explore significance of dose-response curves – example datasets [dataset]. Zenodo; 2023. doi:10.5281/zenodo.8399823
44. Davis MI, Hunt JP, Herrgard S, Ciceri P, Wodicka LM, Pallares G, et al. Comprehensive analysis of kinase inhibitor selectivity. Nat Biotechnol. 2011;29(11):1046–1051. doi:10.1038/nbt.1990
45. Mari T, Mösbauer K, Wyler E, Landthaler M, Drosten C, Selbach M. In vitro kinase-to-phosphosite database (iKiP-DB) predicts kinase activity in phosphoproteomic datasets. J Proteome Res. 2022;21(6):1575–1587. doi:10.1021/acs.jproteome.2c00198
46. Sugiyama N, Imamura H, Ishihama Y. Large-scale discovery of substrates of the human kinome. Sci Rep. 2019;9:10503. doi:10.1038/s41598-019-46385-4
47. Ke G, Meng Q, Finley T, Wang T, Chen W, Ma W, et al. LightGBM: a highly efficient gradient boosting decision tree. In: Advances in Neural Information Processing Systems 30 (NIPS 2017). Red Hook (NY): Curran Associates; 2017.
48. Deznabi I, Arabaci B, Koyutürk M, Tastan O. DeepKinZero: zero-shot learning for predicting kinase–phosphosite associations involving understudied kinases. Bioinformatics. 2020;36(12):3652–3661. doi:10.1093/bioinformatics/btaa013
49. van den Oord A, Li Y, Vinyals O. Representation learning with contrastive predictive coding. arXiv:1807.03748 [Preprint]. 2018. doi:10.48550/arXiv.1807.03748
50. Platt JC. Probabilistic outputs for support vector machines and comparisons to regularized likelihood methods. In: Smola AJ, Bartlett P, Schölkopf B, Schuurmans D, editors. Advances in Large Margin Classifiers. Cambridge (MA): MIT Press; 1999 [verify year: 1999 or 2000]. p. 61–74.
51. Zadrozny B, Elkan C. Transforming classifier scores into accurate multiclass probability estimates. In: Proceedings of the Eighth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD 2002). New York: ACM; 2002. p. 694–699. doi:10.1145/775047.775151
52. Liu Y, Ott M, Goyal N, Du J, Joshi M, Chen D, et al. RoBERTa: a robustly optimized BERT pretraining approach. arXiv:1907.11692 [Preprint]. 2019. doi:10.48550/arXiv.1907.11692
53. Lin TY, Goyal P, Girshick R, He K, Dollár P. Focal loss for dense object detection. In: 2017 IEEE International Conference on Computer Vision (ICCV). IEEE; 2017. p. 2999–3007. doi:10.1109/ICCV.2017.324
54. Holm S. A simple sequentially rejective multiple test procedure. Scand J Stat. 1979;6(2):65–70.

## Supporting information

**S1 Text. Approaches that failed or were not adopted, and post-TEST exploratory analyses.** This covers the cross-encoder R5, the similarity-to-training-substrates feature, Tyr/SH2 docking features, priming-offset features, the prior-debiased site score, ESM-2 3B windows, external curated labels, and the three post-TEST experiments with their audit. Provided separately as `manuscript/supplement.md`.
