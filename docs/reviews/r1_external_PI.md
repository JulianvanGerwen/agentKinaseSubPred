# Review r1 — external PI (computational signaling biology)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Reviewer role: external PI. I judged the overall strategy: whether this design is likely to produce the most accurate, broadly useful human kinase–substrate (KS) predictor with ≈21 h left on 64 cores and 4×RTX 4090. I also checked whether the endpoints match the three user needs (signaling exploration, phosphoproteomics interpretation, kinase activity estimation).
Documents read: PROMPT.md, DESIGN.md v0, LITERATURE.md, docs/data_profile.md, docs/data_acquisition.md, docs/lit/01 (skimmed the others).

## Overall verdict
The scientific framing is strong. It uses protein-cluster splits, a positive–unlabelled (PU) framing with two declared negative sets, kinase-equal macro metrics, fold-aware label-derived features, a leak-audited PPI and a locked test set. This is a better evaluation design than anything in the published KS literature.

The main risk is not scientific; it is **execution within 21 h**. The design lists:
- 6 baselines,
- 9 feature families (≈20 sub-features),
- 5 training regimes,
- 6 label variants,
- 3 LKO tiers,
- 2 universes,
- and a downstream benchmark,

with no time boxes. In practice, feature engineering (F3–F6) and baseline re-implementations (B3/B4) will eat the budget, and the global model will be under-tuned at freeze time.

Two technical issues would silently inflate CV numbers and need fixing before any model is trained:
1. self-inclusion in label-derived features (no cross-fitting for training rows);
2. same-protein co-regulation in F6a/F3c.

The endpoints are close to right, but the two most user-relevant ones (site-centric kinase assignment and kinase-activity inference) are secondary, and the E-all co-primary is the one most exposed to ascertainment bias.

## Quick checks run (read-only, ≤4 cores)
Scripts are in the scratchpad: `sim.py` and `temporal.py`. Results:
- **Kinases evaluable on TEST.** With a 20% protein split of HMR substrate proteins (5 random seeds; this ignores clustering, so the real numbers are slightly lower):
  - 162 ± 8 kinases have ≥3 human-native test positives, 115 have ≥5 and 67 have ≥10;
  - about 2,900 test positive pairs fall on about 2,580 annotated test sites;
  - 58 of the ≥3 kinases have <20 training labels;
  - the top 20 kinases hold 54% of test positives;
  - in one seed, **40 of 151 kinases with ≥3 test positives have them on ≤2 proteins**, and 13 on a single protein. For these, the "per-kinase AUROC" is effectively 1–2 independent observations (often autophosphorylation clusters).
- **Paralog contamination of E-hard.** I estimated the share of a kinase's E-hard negatives that are labelled to a paralog (same gene-symbol stem, a crude proxy):
  - CDK4/6/7/9 ≈13%, CDK5 12%, CDK2 10%, MAPK8–14 9–12%;
  - the median kinase is 0.15%, and 22 kinases exceed 5%.
  - So for CMGC, E-hard penalises genuine shared specificity.
- **Free temporal hold-out exists.** `external/kinase_library/.../Kinase_Substrate_Dataset_count_07_2021.txt` is the full PSP KS file from July 2021 (21,403 rows).
  - Against the 2025 release it gives **1,717 new human×human pairs** (704 in vivo only, 601 both, 412 in vitro only) across 262 kinases (147 with ≥3, 96 with ≥5).
  - 448 of these pairs are on proteins with no 2021 annotation; 507 are on sites already annotated to another kinase in 2021.
  - Only 26 pairs were removed.
- **Covariate shift to the coverage target.** 93% of KS-labelled human sites have LT_LIT ≥1, and only 2.8% are "≤1 MS record and no LT". Among all PSP human sites, 52.5% are in that low-evidence class. Every test positive comes from the well-studied tail, while Priority 2 (coverage) is about the other half.

## Proposed changes (prioritised)

### PI-1 (critical) — §12, add §12a "Time-boxed plan and cut list"
**What to change.** Replace the ordered list with a schedule that has gates. Hours are wall-clock from now (≈T+3 h).

| Window | Deliverable (parallel agents) | Gate |
|---|---|---|
| +0–2.5 h | Labels L2, kinase list, universe, splits and tests; eval harness and metric tests; B1 Cantley and B2 PSSM on CV | Harness reproduces hand-computed toy metrics; B1 macro AUROC E-hard reported |
| +0–4 h (in parallel) | Feature store v1: F0 (one-hot, ESM-2 650M, pLDDT/RSA), F1 priming-aware, F2a cross-fitted, F6a, F6b, F3 screen-PPI, F4. **GPU0:** ESM embeddings. **GPU1:** R3 | Each family passes the planted-leak toy test |
| +4–7 h | **MVP:** R2 LightGBM on F0–F2 + F6 + F3/F4, α-weighting, 1 seed | MVP beats B1 on E-hard macro AUROC on CV (paired bootstrap) |
| +4–9 h (in parallel) | B3 Phosformer (ESM-2-initialised, no MLM pretraining) on GPU2; B4 SELPHI2.0-lite on CPU; R3 → stacked into R2 | Hard stop at +9 h for baselines |
| +7–13 h | Ablation curve (leave-one-family-out from the MVP, not add-one); R1 vs R2; label variants L1/L2/L2+ext; negative ratio | — |
| +13–16 h | LKO T2; unorthodox-kinase analysis; calibration | **Freeze at +16 h** (about T+19 h) |
| +16–19 h | Single TEST run, temporal hold-out, benchmarKIN, U_obs predictions, independent result review | — |
| +19–21 h | REPORT.md | — |

**Explicitly deferred unless the schedule is ahead:**
- F5 GTEx co-expression, F6c, F6d, docking-motif ELM regexes;
- U_all scoring;
- L3;
- B6 KolossuS-15B (use the 650M variant if `pip install` works, otherwise skip);
- the R4 hierarchical model, unless R1 beats R2 for kinases with ≥100 labels;
- 5 seeds for ablations (use 1 seed plus a paired cluster bootstrap; 5 seeds only for the final ≤3 candidates).

**Why.** The current plan has about 40 model trainings × 5 seeds × 5 folds before LKO and finals, and no point at which a usable model must exist. An MVP gate guarantees a reportable result. Leave-one-family-out ablations from the full model answer the same question as an add-one curve, with fewer runs per decision.

### PI-2 (critical) — §7 header, F2a, F3 (fold-aware proximity), F6a, F6d; §11
**What to change.** Add a rule: "Every label-derived feature for a *training* row is computed by **cross-fitting**: an inner 5-fold split over DEV-train substrate clusters, so a row's feature never uses labels from its own substrate cluster. Validation and test rows use all training labels of the outer fold. Additionally, label-derived features for site *s* **exclude all training substrates on the same protein or substrate cluster as *s***, for all rows."

Extend the §11 planted-leak test so that it also fails if the train-row and test-row feature distributions differ for a null kinase, using a KS-test on a random-label kinase.

**Why.**
- "Leave-site-out" is not enough, for three reasons:
  - Sites on the same protein co-vary in perturbation data (protein abundance, multi-site regulation), so for training rows F6a picks up the kinase's other labelled sites on the same protein.
  - Test proteins never have such neighbours.
  - The learned PSSM (F2a) and the "proximity to training substrates" feature (F3) include the row's own label.
- LightGBM then learns to trust these features more at training time than they deserve at test time. The CV-vs-test drop alarm would only reveal this after the freeze.
- This is the same mechanism by which SELPHI2.0's co-regulation features leak; the design must not reproduce it.

### PI-3 (major) — §6, endpoint set and kinase inclusion
**What to change.**
1. **Kinase inclusion for the primary macro:** ≥3 human-native test positives on **≥2 distinct substrate clusters**, excluding autophosphorylation from the positive count.
   - The kinase list is frozen and committed when the split is made, before any model is trained.
   - About 130–140 kinases qualify, versus about 160 if the ≥2-cluster rule is dropped.
   - Report the ≥3-any set as secondary.
2. **Utility endpoints:** promote two endpoints to co-primary, as the two user-facing endpoints:
   - (a) **kinase-balanced site-centric Hit@5**: for each true kinase, average Hit@5 over its test sites, then macro-average over kinases. This is the phosphoproteomics-interpretation endpoint, and the kinase-frequency prior cannot game it.
   - (b) **benchmarKIN perturbation AUROC** of "PSP + top-N predicted" vs "PSP only".
3. **Demote E-all AUPRC** to secondary, and add **E-matched**: E-all negatives reweighted to match the positives' distribution over PSP site-evidence bins (LT_LIT 0/1/2–4/≥5 × MS 0–1/2–9/≥10) and compendium coverage (absent/<10/≥10 conditions). Report macro AUPRC-lift under E-matched.
4. **Pre-register "beats baseline"** as: paired cluster-bootstrap Δ (ours − baseline) on the primary endpoint, with the 95% CI excluding 0. Report the same Δ within the label-count bins 0–4 and 5–19, so that "harmonious" is tested rather than asserted.

**Why.**
- My check shows a quarter of the ≥3 kinases have positives on ≤2 proteins, which is pure noise at kinase-equal weighting.
- E-all rewards studied-ness: KS sites are 93% LT ≥1, and LightGBM's NA handling of compendium/PPI missingness is a popularity proxy. A diagnostic popularity model does not stop a co-primary from being won by popularity.
- Users mostly ask "which kinase did this?" (site-centric) and "is kinase X active?" (benchmarKIN); the endpoints should reflect that.

### PI-4 (major) — §3 test labels and §6, add a temporal hold-out
**What to change.** Add **T-temporal**: the PSP 07/2021 KS file is the label snapshot. Positives are human-native pairs that are new in 05/2025 (1,717 pairs, 147 kinases with ≥3).
- **T-temporal-dev:** new pairs on DEV proteins. The model is retrained on DEV labels present in 2021 only. This is the "new kinase for a known protein/site" setting.
- **T-temporal-test:** new pairs on TEST proteins, evaluated at the final run with the standard final model.
- The same E-hard and E-all definitions apply. Report it alongside TEST, and use it to check the calibration claims behind the coverage numbers (estimated precision ≥0.5 or 0.8).

**Why.** It is prospective, costs nothing to obtain, needs only the existing loader, and contains 507 pairs on sites already annotated to other kinases. That is exactly the multi-kinase disambiguation case.

**Caveat.** About 41% of the new pairs are in vivo only and may come from perturbation studies in the compendium. Report results by evidence class.

### PI-5 (major) — §8 R2/R3 and calibration
**What to change.** Compare two objectives for R2:
- the existing kinase-grouped binary/lambdarank;
- a **site-grouped listwise objective**: lambdarank grouped by site over the kinases of matching residue class, positives = labelled kinases, on annotated sites only. Use it either as a second model or as an auxiliary loss in R3 (softmax over kinases).

Then define the released "cross-kinase comparable" score as the calibrated site-grouped probability. Keep the within-kinase percentile for substrate lists.

**Why.**
- Per-kinase isotonic calibration with a "prevalence correction" needs the unknown true substrate count per kinase, which is unidentifiable under PU. Cross-kinase comparability will therefore be poor, and the site-centric use case is where users feel that.
- A site-grouped head is cheap: about 13k annotated sites × ≤430 kinases. It directly optimises the interpretation task, and the E-hard construction already tests it.

### PI-6 (major) — §7 F6b, prioritise the perturbation modality and broaden its kinase map
**What to change.** Before building F6, create `results/labels/condition_kinase_map.tsv` for all 1,937 compendium conditions:
- drug → targets from the Hijazi matrix (60 drugs);
- for the other compendium drugs, especially the decryptM series: targets from the decryptM paper's own target annotations or ChEMBL/DTC pChEMBL ≥7. Record the source per row; if a drug cannot be mapped in 1 h, leave it unmapped;
- genetic perturbations (85) → the perturbed gene;
- ligand/GF → no mapping.

F6b then uses "k inhibited/perturbed vs not". Rank feature-building effort as: F6b > F6a > F1-priming > F3 > F4 > F5 > F6c/F6d.

**Why.**
- Inhibitor and knockout response is the only modality that is causal and kinase-specific. It is what lets the model beat motif-only methods for kinases with shared motifs (CDK/MAPK, AGC).
- With the Hijazi matrix alone it covers 183 of 1,937 conditions. The decryptM dose series (790 conditions) are the richest data and are currently unmapped, because the Kinobeads route is blocked.

### PI-7 (major) — §7 F8 and data acquisition, acquire Sugiyama 2019 now
**What to change.** Sugiyama et al. 2019 (Sci Rep, open access; in vitro lysate screen of about 385 kinases, about 175k pairs) is listed as a feature, but it is not in data_acquisition wave 1. Download it in the next 30 min as a site × kinase binary or score feature, with an ablation.

Also add one R3 experiment: **pretrain on Cantley-atlas percentiles** (Phosformer-ST label rule: >90 → positive, <50 → hard negative) **before fine-tuning on L2**. Timebox it to 2 GPU-h.

**Why.**
- Both are kinase-specific signals that exist for under-studied kinases independently of PSP. This is the most direct lever for the "harmonious across kinases" goal and for the label-count bins 0–4 and 5–19.
- KolossuS, the cleanest published comparator, is trained on exactly these two sources.

### PI-8 (major) — §8 R0 baselines, scope them to what makes the comparison fair
**What to change.**
- **B3:** reimplement Phosformer as a 6-layer encoder **initialised from ESM-2 (35M or 150M)**, with no MLM pretraining, trained on our folds. Add the released weights for inference only. Budget: ≤3 GPU-h total.
- **B4 "SELPHI2.0-lite":** their RF hyperparameters and their feature *definitions*, computed from our fold-aware feature store, with the deviations documented. Budget: ≤2 h. Do not port their pipeline.
- **Add B7 "motif-only ensemble":** logistic regression on B1 and B2 (priming-aware). This shows how much the non-sequence modalities add over the best sequence-only model we can build, which is the scientific claim.
- **B5 NetworKIN** stays as published.
- **B6** is dropped, or run as the 650M variant only.

**Why.** PROMPT asks for fair retraining of these three benchmarks (Cantley, Phosformer, SELPHI2.0), not faithful reproduction of their engineering. Without a time cap, B3/B4 are the most likely budget sink.

### PI-9 (major) — §6 and §5, stratify by site evidence (coverage validity)
**What to change.** Report every primary and utility endpoint, and the E-matched metric, stratified by:
- (i) PSP site-evidence bin of the positive (LT ≥1 vs MS-only);
- (ii) compendium coverage (absent / <10 / ≥10 conditions);
- (iii) PPI or localisation data available vs missing.

Add a **missingness-only control** on E-hard: R2 trained on NA-indicator features only.
- If its macro AUROC is ≥0.60, the missingness pattern is acting as a label proxy. Mitigate by imputing with fold-mean and adding explicit indicators only for families where the ablation shows they help on E-matched.

The Priority-2 coverage claim (§6 Coverage) is made per stratum, and the low-evidence stratum is reported separately.

**Why.** Labelled sites are 93% LT ≥1, while 52% of the target universe has ≤1 MS record and no LT. Without this stratification the coverage numbers extrapolate precision from a population that does not resemble the sites being covered.

### PI-10 (minor) — §5 E-hard and §6, add a paralog-aware variant and family-level metrics
**What to change.**
- Add **E-hard-fam**: E-hard minus sites labelled to any kinase with ≥60% kinase-domain identity to *k*, or in the same Manning family. Report it as a sensitivity analysis.
- Also report **family-level** site-centric Hit@k, where kinases are collapsed to Manning family by max score.

**Why.**
- 9–13% of E-hard negatives for CDK and MAPK members are paralog-labelled, and PSP often assigns family-ambiguous kinases (CDK1 vs CDK2, PKC isoforms).
- Family-level accuracy is what many users can act on, and it shows how far isoform-level claims are supported.

### PI-11 (minor) — §2 U_all and §8 inference
**What to change.**
- Defer U_all (1.87M × ~430 ≈ 800M pairs) to an optional final step with a sequence-only model; the released core is U_obs.
- Restrict U_obs to sites with PSP evidence ≥1 LT or ≥2 MS, or ≥1 compendium/CPTAC observation, so that single unconfirmed MS localisations are flagged rather than mixed in.
- **Pre-register top-N for activity libraries:** choose N (50/100/200 or within-kinase percentile ≤1%) on the benchmarKIN Hernandez benchmark, then report CPTAC and Hijazi at that N. Note that Hijazi is contaminated, because F6b uses it.

**Why.** This saves several CPU-hours and avoids selecting N on the evaluation data.

### PI-12 (minor) — §4.2 and §11, safeguard the locked test
**What to change.** The TEST kinase list (PI-3), the frozen config hash, and the code commit hash are written into `results/splits/test_lock.json` before the freeze. `eval_test.py --final` refuses to run if the config hash differs from the frozen one. For LKO, T3 uses the same lock.

**Why.** This is a cheap guard against the most common agent failure mode (quiet re-runs after seeing test numbers), and it makes the "one evaluation" claim auditable.

## Answers to the brief's questions
- **Most accurate, broadly useful within budget?** Yes, if PI-1 (time boxes and the MVP gate), PI-2 (cross-fitting) and PI-6 (a broad perturbation→kinase map) are adopted. Accuracy gains over Cantley will come mainly from:
  - priming-aware motif plus learned PSSM;
  - inhibitor/knockout response;
  - co-regulation with known substrates.
  They will not come from PPI, localisation or co-expression; expect ≤0.01–0.02 macro AUROC from these, given 0.3–1.5% PPI coverage.
- **Over-engineered:**
  - 6 label variants;
  - 5 regimes including R4;
  - F5/F6c/F6d/docking motifs;
  - U_all;
  - B3 with pretraining, B4 full port, B6-15B;
  - 5 seeds everywhere.
- **Missing:**
  - cross-fitting and same-protein exclusion;
  - the temporal hold-out;
  - a site-grouped objective and cross-kinase calibration;
  - the Sugiyama feature;
  - a drug/genetic → kinase map for the full compendium;
  - evidence-stratified metrics for the coverage claim;
  - a pre-registered "beats baseline" rule.
- **Primary endpoints:** keep macro AUROC under E-hard (with the tightened kinase set). Add kinase-balanced site-centric Hit@5 and benchmarKIN as co-primary utility endpoints. Move E-all AUPRC to secondary and add E-matched.

## Risks I could not resolve
1. **PSP in vivo labels may be derived from the same perturbation studies** used as features. PSP KS rows carry no PMIDs, so the overlap cannot be audited directly. The in-vitro-supported stratum is the only clean check for F6.
2. **Gene-symbol orthology mapping** of rodent kinases is not sequence-verified. This is a low risk, but it matters for small kinase families.
3. **Family-ambiguous PSP annotations** (generic "PKC" or "CDK" assignments curated to one isoform) cap the attainable isoform-level accuracy. We cannot estimate the noise rate without manual curation.
4. **The benchmarKIN gold standard partially overlaps PSP training labels.** The "PSP + predicted" vs "PSP only" contrast mitigates this but does not remove it.
5. **Agent integration risk.** With about 6 parallel agents writing features into one store, schema drift is likely. One owner of the feature-store schema and a per-family contract test are advisable, but I cannot judge staffing from the documents.
