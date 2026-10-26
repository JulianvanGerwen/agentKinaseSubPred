# Review r2: external PI (computational signaling)

**Reviewer:** fresh external PI. I did not write round 1.
**Reviewed:** DESIGN.md v1 (14:30), with its §15 change map. I also read PROMPT.md and the r1 reviews (all critical and major items), and checked the current code and data state (`src/ksp/splits.py`, `cache/labels/ks_pairs.parquet`, `cache/sites/universe.parquet`, `external/phosformer/`, docs/DELEGATION.md).
**Time of review:** about 14:45. G1 is at 16:50.

## Verdict
v1 is a large improvement and is scientifically coherent. Five r1 fixes are in place:
- cross-fitting;
- paralog-masked negatives;
- autophosphorylation removed from the endpoints;
- a sealed TEST set with freeze-then-batch scoring;
- M1–M5 endpoints with pre-declared superiority claims.

As written, I would trust the evaluation more than any published kinase–substrate benchmark.

The remaining problems are about delivery, not concept:
1. **No feature contract.** The cross-fitting rule is the single most important correctness device. It is stated as a principle, but there is no interface for it, and four agents are building label-derived features in parallel before the folds exist.
2. **G3 is overloaded.** Nested R3 stacking alone means about 36 R3 trainings.
3. **Three loose ends could still bias results:**
   - the canaries do not exercise the leak path they are meant to detect;
   - one input channel would hand the model the site's own MS-observation status;
   - the stratified hard-negative mix breaks the single log ρ_k prior correction that the coverage numbers rely on.

I propose 8 changes. Two items are cuts.

## Quick checks I ran (read-only, ≤4 cores)
Inputs were L2 pairs with an included kinase, in class and in the universe: 17,681 pairs over 398 kinases.

| Check | Result | Used in |
|---|---|---|
| Label-count bins (1–4 / 5–19 / 20–99 / ≥100) | 110 / 150 / 101 / 37 kinases. Of the 503 included kinases, 105 have zero L2 labels. | context |
| Positives with a same-kinase sibling on the same protein (kinases with 10–100 labels) | **65.9%**; mean 1.87 sites per (kinase, protein) | R2PI-2 |
| L2 positive pairs on sites **not** in U_ms | **16.6%** (2,935 pairs). U_obs \ U_ms has only 3,382 sites, and 67% of them are labelled. | R2PI-4 |
| Share of sites with no compendium or CPTAC value (F6 all NA): positives vs unlabelled U_obs sites | 36.7% vs 36.0%. The flag has AUROC 0.525 overall and 0.455 inside U_ms. | R2PI-4 (good news) |
| Annotated (L2-labelled) sites as a fraction of U_ms | 2.7%. A 50% hard-negative mix therefore oversamples that stratum about **18×** and undersamples the rest about 0.5×. | R2PI-5 |
| LightGBM, synthetic 400k×160 with 40% NA, 4 threads, 100 trees | 6.8 s at 63 leaves and 11.6 s at 255 leaves, i.e. about 70–115 s per 1,000 trees | feasibility |
| Phosformer MLM corpus on disk | `external/phosformer/hf_dataset_curated/corpus/{kinases,substrates}.fa` is present, so B3 pretraining is feasible. | R2PI-8 |

The LightGBM measurement suggests §10's "15 min at 64 threads" is about 10× pessimistic for a 63–255-leaf model. Model fitting is not the bottleneck. Feature construction and integration are.

## (1) Were round-1 critical and major points incorporated correctly? (spot-check of 8)

| r1 item | v1 location | Assessment |
|---|---|---|
| PI-2 / DL-1: cross-fitting for training rows; same protein and cluster excluded | §7 header | **Correct** as a rule. The interface is missing (R2PI-1), and the canary does not test it (R2PI-2). |
| KB-1 / STAT-2 / ALG-2: paralog-masked E-hard and training negatives | §5 | **Correct.** It takes the union of the identity rule (≥60%) and the subfamily rule, plus a Hamming-≤2 mask. It is conservative, which is right. Minor: the unmasked E-hard sensitivity arm (DL-2, STAT-2) was dropped. Keep it as a one-line secondary; it costs nothing. |
| DL-3: autophosphorylation by accession, excluded from endpoints | §3 | **Correct** for positives. The wording should also exclude *k*'s own unlabelled sites as negatives for *k* (DL-3 item 2). The R3/B3 ablation that drops in-domain self pairs from training was dropped. That is acceptable, since evaluation already excludes self pairs. |
| ALG-1: case–control prior correction | §6.6, §8.5 | **Partially correct.** Item 4 (non-uniform negative sampling invalidates log ρ_k) was not carried over, and §5's hard-negative mix is exactly that case (R2PI-5). |
| KB-2: non-circular downstream benchmark | §6.7 | **Partial.** Two parts landed: dropping benchmark studies by PMID/PXD, and dropping CPTAC features. Three parts were lost: leave-drug/target-out, the "ours minus PSP" and "PSP ∪ ours" library arms, and separate reporting of upstream MAP2K/RAF kinases (R2PI-7). |
| PRIOR-3 / CMP-2: Phosformer fair retrain | §8.1 B3 | **Mostly correct.** v1 reproduces the MLM pretraining and re-tunes the negative ratio. The paper's augmentations and the shared negative generator are not stated (R2PI-8). |
| STAT-1: endpoints for understudied kinases | §6.3 M4/M5 | **Correct.** STAT-1(b), selection-aware nested CV for M5, was not adopted. Acceptable if M5 is labelled "optimistic (selection on same CV)" in the report. |
| PI-9: site-evidence strata for the coverage claim | §6.4, §6.6 | **Not landed.** The §15 map points to §6.5, but neither §6.4 nor §6.6 has LT-vs-MS-only or compendium-coverage strata. Folded into R2PI-5. |

## (2) Coherence and feasibility against the §0 gates

- **G1 (16:50): feasible**, with two notes for core-infra. These are not scope changes.
  - (a) The per-kinase 15–25% TEST band cannot always be met. Clusters are indivisible, and some cases have no solution at all: a kinase with 7 positives can only get 14% or 29%. The optimiser must use a **soft penalty**, and the acceptance criterion should be stated. Example: ≥90% of kinases with ≥10 positives inside the band, and every K_eval kinase between 10% and 35%.
  - (b) Canary labels must be in the label table **before** the split freeze, so every feature agent computes features for them (R2PI-2).
- **G2 (20:50): feasible only if the label-derived features (F2a, F3-label, F6a) arrive through one shared cross-fitting driver.** If each of three agents writes its own outer×inner fold loop against folds that are frozen two hours from now, the likely outcomes are a late G2 or a silent inconsistency. The second is worse (R2PI-1).
- **G3 (02:20): overloaded.** In 5.5 h the plan has R3, nested stacking, B3, B4, the regime grid, four label variants and LKO T2.
  - Nested stacking alone means 5 outer × 5 inner + 5 outer + 5 inner + 1 final ≈ 36 R3 trainings. At the budgeted ≤1 GPU-h each, that exceeds the GPU-hours available on GPUs 1–3 overnight.
  - B4 has no owner in DELEGATION.md.
  - See R2PI-3.
- **G4–G6:** feasible as written.

## (3) The single biggest risk
**Integration of the cross-fitted, label-derived features.** These features (learned PSSM, co-regulation with training substrates, RWR proximity to training substrates, and R3 stacked scores) carry most of the expected gain over Cantley for kinases with ≥20 labels. They are also where leakage would inflate CV.

Three facts combine:
- they are being built by three agents without a common interface;
- the folds are not frozen yet;
- the end-to-end alarm (canaries) uses independent random sites, which have essentially no same-protein siblings. Real positives have one 66% of the time.

So the most likely failure is not a bad model. It is a CV number that looks good because of a leak no alarm can see, followed by a CV-vs-TEST gap after the freeze. R2PI-1 and R2PI-2 together close this.

## Proposed changes (8)

### R2PI-1 (critical): §7 header and §11.1, add a fit/transform contract and one cross-fitting driver
**Change.** Add the following to §7 after the cross-fitting rule:

> Every label-derived feature family exposes `fit(train_pairs: DataFrame[kinase_acc, site_id], ctx) -> Model` and `Model.transform(kinase_accs, site_ids) -> float32 array`. These functions contain no fold logic and read no split files.
>
> A single driver, `ksp.crossfit` (owner: core-infra), loops over:
> - the outer folds;
> - the inner folds from `results/splits`;
> - DEV→TEST.
>
> It writes `cache/features/<family>/outer{f}.{train|val}.npy`, where train rows come from inner-fold fits and val rows from the outer-train fit, plus manifests. Feature agents deliver `fit`/`transform` and a unit test. They never write fold-specific outputs themselves.
>
> The DL-1 toy test is run through the driver for **each** family. In that test, two same-protein sites with identical profiles are planted, and the test fails if one sees the other.

Make this the first item in each feature agent's brief now: F2a, F3-label and F6a are due by G1 + 1 h.

**Why.** It turns the most important correctness rule into one piece of code that can be tested. It also removes the dependency of three agents on split files that do not exist yet, since they can develop against random folds. Without it, G2 depends on three independent fold loops agreeing.

### R2PI-2 (major): §11.3, make the canaries exercise the same-protein leak path
**Change.**
- **Sampling.** Sample each negative canary's labels **protein-clustered**:
  - draw proteins from U_ms proteins, weighted like real substrate proteins;
  - then draw 1–4 sites per protein from the empirical per-(kinase, protein) site-count distribution (mean 1.87), until 10–100 labels are reached.
- **Timing.** Inject the canaries into the label table at split freeze (G1).
- **Evaluation.** Score the canaries under two negative sets:
  - E-all-type negatives (U_ms, same masks), where the ≤0.57 threshold applies;
  - E-hard, reported for information only. Random-site positives against annotated-site negatives mixes popularity into the result.

**Why.** 65.9% of real positives have a same-kinase sibling on the same protein. Independent random sites almost never do. So the current canaries cannot detect the leak that motivated DL-1 (sibling co-regulation and the PSSM including a sibling's window). With clustered canaries, a broken cross-fit would push canary CV AUROC well above 0.57.

### R2PI-3 (major): §8.4 and §12 G3, use late fusion instead of nested R3 stacking, and de-load G3
**Change.**
1. **Late fusion by default.** Train R3 on the 5 outer folds plus one DEV-final run: 6 trainings.
   - The fused score is `σ(a·logit_R2 + b·logit_R3 + c)`, with (a, b, c) fitted on outer-fold out-of-fold predictions.
   - Optionally, use separate (a, b) per label-count bin. That is about 3–8 parameters, so the optimism is negligible.
   - Nested stacking (R3 as an R2 input) is kept as an option only if one R3 fold measures ≤15 GPU-min **and** late fusion is already in place by G3 − 2 h.
2. **Assign owners for the G3 items now.** B4 SELPHI2.0-lite gets a CPU owner who starts at G2. LKO T2 and the label variants are run with **R2 only**, seed 0, and without re-fusing R3.
3. **Regime search is a coordinate search on fold 0.** The order is: negative ratio × mix (6 runs), then α (3), then objective (3). After that, run 5-fold on the top 3, logged in REGISTRY.tsv. The full 54-cell factorial is not run.
4. **Cuts:**
   - B6 KolossuS: it is a contaminated comparator and none of the three pre-declared claims uses it;
   - L2_noauto: evaluation already excludes self pairs; deferred.

**Why.** Nested stacking means about 36 R3 trainings, which does not fit next to B3 on GPU 0 and the embeddings on GPUs 2–3. Late fusion gives essentially the same accuracy benefit when there is one base learner and one neural model.

### R2PI-4 (major): §8.4, §7 F1 and §9, never encode the site's own observation status
**Change.** Add to §7:

> No feature may encode whether the scored site itself is observed by MS (U_ms membership, PSP MS/LT counts, Ochoa presence).

In the R3 phospho channel and the F1 priming-aware score, the **central position is always encoded as unphosphorylated**. Only flanking positions carry U_ms status. Add an assert in the R3 dataloader, and a toy test in which one centre bit is flipped and the output must not change.

**Why.** 16.6% of L2 positive pairs are on sites outside U_ms. E-all-type training negatives are U_ms by construction. A centre bit would therefore be a perfect positive flag for one in six positives. On deployment it would also push LT-only sites up for every kinase.

The good news is that F6 missingness is balanced: 36.7% of positives vs 36.0% of unlabelled sites have no compendium or CPTAC value. So the general NA handling is not a comparable shortcut. Only this channel needs fixing.

### R2PI-5 (major): §5 training negatives, §6.4, §6.6 and §8.5, correct the prior properly under the hard-negative mix and stratify coverage by evidence
**Change.**
1. **Weights.** With an h% hard mix, weight each negative by π_stratum / s_stratum (normalised). Here π is the stratum's share of k's masked universe and s is its share of the sample. The weighted negatives then represent uniform sampling, and `logit + log ρ_k` is exact.
2. **Grid.** Add "weighted vs unweighted" to the {25%, 50% hard} grid.
   - If the unweighted mix wins on M1-CV, it may be used for rankings.
   - Released **probabilities** (§6.6) must then come from the weighted model or from a calibrator fit on **exhaustive** out-of-fold U_ms predictions. The calibrator must not be fit on the sampled negatives.
3. **Strata.** Add these strata to §6.4 and §6.6 (PI-9 had not landed):
   - site evidence: PSP LT ≥1 vs MS-only;
   - compendium coverage: absent / <10 / ≥10 conditions.
   Coverage counts are reported per stratum, with the low-evidence stratum separate.

**Why.** Annotated sites are about 2.7% of U_ms, so a 50% hard mix oversamples them about 18× and undersamples the rest about 0.5×. A single log ρ_k is then wrong by about +2.9 nats in the annotated stratum and about −0.7 nats in the rest, as far as the model can tell the strata apart through popularity proxies.

The direction matters. The Priority-2 headline ("sites with no PSP kinase reaching p ≥ 0.8") is in the stratum whose probabilities would be **over**-stated.

### R2PI-6 (major): §6.2 and §6.8, pre-declare which score each endpoint uses, and a selection rule that can pick the M2 model
**Change.** Add to §6.8:

> (a) The **within-kinase ranking** (substrate libraries, M1, M3, M4, M5) is chosen by M1-CV. When candidates differ by less than the MDD (0.015), the tie is broken by M3-CV.
>
> (b) The **cross-kinase score** (M2 and the §6.6 probabilities) is chosen by **M2-CV** among three candidates:
> - calibrated R2-binary;
> - site-grouped R2;
> - late-fused R2+R3.
>
> (c) Each baseline is scored on M2 with its native cross-kinase output: the Cantley percentile and the Phosformer probability.
>
> (d) Both choices are written to FREEZE.json.

**Why.**
- With M1-CV as the only selection criterion, the site-grouped objective (PI-5), which was added specifically for M2, can never be selected.
- The prior-corrected probability adds about log n_k per kinase. That helps a micro site-centric ranking but hurts the within-kinase-averaged Hit@5 that M2 actually measures.
- So the score fed to M2 is a real design choice, and it must be fixed before TEST.

### R2PI-7 (major): §6.7, complete KB-2 (leave-drug-out, and arms that measure added value)
**Change.**
1. **Leave-drug/target-out.** For the perturbation benchmark, the clean variant also drops every compendium condition whose drug's target set (from the F6b map) includes a benchmark target kinase. 17 Hijazi inhibitors recur in other compendium studies.
2. **Library arms.** At each fixed N, add "ours minus PSP pairs" and "PSP ∪ ours" beside the existing PSP-only, Cantley and shuffled arms. Make the shuffle degree-preserving, within kinase.
3. **Upstream kinases.** Report the upstream, direct-substrate-poor targets (MAP2K1/2, BRAF, RAF1, MAP3K*) separately. State in advance that the libraries are not tuned to them.

**Why.**
- Without (1), F6b transfers the benchmark signal through drug identity.
- Without (2), an "ours" library that contains PSP mostly measures PSP.
- This is the only evidence for the user's third use case (kinase activity), so it should not be over-claimed.

### R2PI-8 (minor): §8.1 B3, write down the fairness contract and a fallback
**Change.**
- **Shared setup.** B3 shares with R3 the labels (L2), the splits, the negative generator (KP and paralog masks), the seed and the early-stopping criterion (M1-CV).
- **Augmentations.** B3 includes the paper's augmentations as one block:
  - kinase up-sampling (to 50 and 1,118 positives);
  - independent ±5 shifts of the kinase-domain boundaries;
  - 5% kinase-token masking.
- **Fallback.** If MLM pretraining has not converged by G2 (20:50), fine-tune from the checkpoint at that point and document the loss curve. Do not let pretraining push past G3.
- **Corpus.** The pretraining corpus is already in `external/phosformer/hf_dataset_curated/corpus/kinases.fa`. Record its checksum in DATA.md as the pretraining input.

**Why.** Leaving out the paper's augmentations, or using a different negative generator, would make the "beats Phosformer" claim, one of the three pre-declared superiority claims, contestable. The fallback keeps GPU 0 from becoming a G3 blocker.

## Cut or add, in one line each
- **Cut:**
  - nested R3 stacking (default);
  - B6 KolossuS;
  - L2_noauto;
  - the full regime factorial;
  - F5 and F6d, already deferred: keep them deferred.
- **Add:**
  - the cross-fit driver and contract;
  - clustered canaries;
  - the centre-position mask;
  - stratum-weighted negatives for probabilities;
  - an M2-specific selection rule;
  - leave-drug-out downstream.
- **None needed:** labels and universe (§2–3), split construction (§4.1, §4.3–4.6), M1/M3/M4 definitions, the superiority and Holm rule (§6.8), the PPI exclusion list (§7 F3), and the process rules (§11.4–11.8).

## Risks I could not resolve
1. PSP in-vivo labels may come from compendium studies, since there are no PMIDs per KS row. Only the in-vitro stratum is a clean check of F6. Same as r1.
2. E-hard's label source for TEST negatives (human-native only, or all L2/KP sources) is not stated explicitly in §5. I assume "any L2 label to a kinase outside {k} ∪ Par(k), after KP masking". Please state it.
3. For the 105 included kinases with zero labels, accuracy rests on Cantley plus kinase descriptors. LKO T2 measures this on CV, but T3 (unseen family and protein) will have few kinases. Expect wide CIs.
4. The NFS-shared tree and convention-only sealing remain.
