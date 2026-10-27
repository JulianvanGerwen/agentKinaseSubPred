# Clarity edit and trim: change log (2026-10-27, 12:28–12:52 CET)

Files edited in place: `manuscript/manuscript.md`, `manuscript/supplement.md`. I did not change tables, headings, code spans that were already there, figure/file names, citation numbers, the AI-exercise flag, the AI-use statement or the Author summary. The six-pattern rules come from the edit-for-clarity skill.

## Word counts

| Part | Before | After |
|---|---|---|
| Main text (Introduction to the end of Methods; table rows excluded, captions included) | 9,346 | 8,486 (−860, −9%) |
| Abstract | 299 | 299 (limit 300) |
| Author summary | 199 | 199 (unchanged) |
| Supplement (whole file) | 4,917 | 4,944 (clarity edits only, no cuts) |

**The target of about 7,000 was not reached.** Under the "never cut" list (limitations, captions, tables, Methods parameters, literature credits, key decisions, agent-setup content, every caveat), I could not find more than about 860 words of repetition to remove. Clarity fixes also add words, since inserting "that" or unpacking stacked nouns makes phrases longer. The PI would need to approve any of these further cuts (rough savings in brackets):
- (a) Move the R2/R3 hyperparameter lists and the feature-source lists in Methods into an S2 table (~450).
- (b) Replace the "What went wrong" bullets with one sentence and a pointer to S1 Text §6, which already lists every item (~180).
- (c) Move the Downstream and Coverage results to S1 and keep one sentence each (~280).
- (d) Remove the numbers from the TEST audit-caveat bullets that Limitations 3–5 repeat (~120).
- (e) Move the late-fusion sweep values, the B3 corpus counts and the canary-gate details to S1 (~120).

## Clarity changes, by pattern

### manuscript.md

**Pattern 1, rare word combinations (5):** "amendments landed" → "were made"; "fixes … landed" → "were applied"; "rules compiled into code" → "rules built into code"; "popularity-heavy external tier B" → "external tier B, where popularity is high"; "Its caveats govern how the numbers can be read" → "Four of its caveats limit how the numbers can be read".

**Pattern 2, inanimate subject with a human-only verb (13):**
- "The ledger refuses a second batch" → "blocks" (Fig 1 caption), and "refuses further batches" → "blocks" (Methods).
- "Two alarm runs found no sign" → "showed no sign".
- "97 of 116 matched kinases favoured the final model" → "the final model was better on 97 of 116 matched kinases", and the same fix for 52/27 kinases.
- "The Cantley-lab atlases measured these preferences" → "The Cantley lab measured these preferences".
- "Learned methods report" → "Papers on learned methods report"; "SELPHI2.0 reports" → "The SELPHI2.0 paper reports".
- "the ≥0.5 tier … claimed 0.615 and observed 0.603" → "had a claimed probability of 0.615 and an observed precision of 0.603".
- "says nothing about whether" → "gives no information on whether".
- "The same retrospective credits" → "lists as strengths".
- "clusters … never train R3" → "are never used to train R3".
- "Standalone values rank unscored pairs last" → "For standalone values, unscored pairs are ranked last" (Table 3 caption).

**Pattern 3, nouns used as verbs (0):** I found none that are non-standard. "port", "pin" and "key by" are normal usage for this readership.

**Pattern 4, missing "that"/"which" (6):** "proteins carrying" → "proteins that carry"; "pairs they could score" → "pairs that they could score"; "libraries that kinase-activity inference depends on" → "libraries on which kinase-activity inference depends"; "a phospho-channel bit marking" → "a bit that marks"; "filtered ranks, removing" → "filtered ranks, which remove"; "proteins whose labels the deployment model was trained on" → "proteins with labels in the deployment training set".

**Pattern 5, stacked nouns (22):**
- Abstract: "global LightGBM pair model" → "global LightGBM model of kinase–site pairs"; "label-free Cantley motif-atlas scores" → "distillation of the label-free Cantley motif atlas".
- "ESM-2 kinase-domain embedding" / "ESM-2 domain embedding" → "ESM-2 embedding of the kinase domain" (Fig 1 caption, F7, R3).
- "Cantley score profile" → "profile of Cantley scores" (Fig 1 caption, F7).
- "perturbation phosphoproteomics compendium" → "compendium of perturbation phosphoproteomics".
- "centre-site observation flags" → "flags for whether the centre site was observed".
- "planted-label toy test" → "toy test with planted labels".
- "Human Cell Map NMF cosine" → "NMF cosine from the Human Cell Map"; "Hein profile correlation" → "profile correlation from Hein et al.".
- "compendium study block" → "study block of the compendium"; "HPA tissue RNA" → "tissue RNA from HPA".
- "Cantley universe percentiles" → "Cantley percentiles against the site universe"; "external-PI design reviewer" → "design reviewer in the external-PI role"; "R3 phospho channel" → "phospho channel of R3".
- "frozen harness commit" → "commit of the frozen harness"; "Hernandez-Armenta perturbation compilation" → "perturbation compilation of Hernandez-Armenta et al."; "seven CPTAC tumour cohorts" → "seven tumour cohorts from CPTAC".
- "batch per-kinase/per-site tables" → "per-kinase tables" / "audit re-cut"; "TEST-harness code (review)" → "(review of the) code of the TEST harness", in 2 places.
- "design-level critical issues" → "critical issues at the design level"; "similarity-to-training-substrates feature" → "feature for similarity to training substrates"; "per-agent memory cap" → "memory cap of 60 GB per agent".

**Pattern 6, long subject–verb or verb–object gaps (7):**
- Abstract: the late-fusion sentence was restructured so the long feature list comes last.
- R3: the site-encoder sentence ("takes [list] … and passes them through") was restructured.
- B3: the corpus/fine-tuning sentence was split.
- B4: the column-port sentence was split.
- Results, calibration: the three-clause audit sentence was split.
- Late fusion: the z-scoring rationale was moved next to the z-scoring step.
- Data: the U_ms restriction sentence was reorganised so the reason follows the action.

### supplement.md

**Pattern 1 (5):** "does not starve the Y negative pool" → "does not leave Y kinases short of negatives"; "unseen-drug gap" → "gap to PSP on unseen drugs"; "fixes landed" → "were applied"; "landed" (twice) → "were made"; "sampler-population difference" → "sampling difference between canaries and E-hard negatives" (this matches the wording of main-text Limitation 8).

**Pattern 2 (5):** "The Tyr diagnosis found" → "showed"; "the round-2 review had measured" → "reviewers in round 2 had measured"; "a timing pilot measured" → "in a timing pilot, inference alone took"; "tier claimed 0.87 and achieved 0.54" → "had a claimed probability of 0.87 and an observed precision of 0.54"; "An independent audit … overruled" → "An independent auditor … overruled".

**Pattern 5 (3):** "Centre-site observation flags" → "Flags for whether the centre site was observed"; "fold-0 sanity harness" → "sanity harness on fold 0"; "B3 Phosformer pretraining" → "Pretraining of B3 (Phosformer; …)".

**Patterns 3, 4 and 6:** no instances found.

## What was cut or moved (manuscript.md)

**Results**
- **Benchmark design:** merged the Fig 1 lead-in into a pointer, and moved the L2 label counts (18,638 / 401 / 14,747) to Methods › Labels.
- **Cross-validation:**
  - Removed the restatement of Table 1 (M1/M2/M3 of the fusion, R2 0.776, R3 0.754).
  - The fusion-weight and nested-M2 sentences now live only in Methods (Late fusion, Calibration).
  - Removed the restatement of the Table 2 ablation steps; the F6 final-configuration values now sit behind a pointer to S1 §2.4.
- **Failed improvements:** per-item numbers moved behind a pointer to S1 §2–3.
- **Retrained baselines:**
  - Removed the B3 architecture recap and the B4 column counts, which Methods already gives.
  - Moved the B1/B2 definitions and the 72-configuration grid to Methods › Baselines.
- **TEST:** removed the final-model row restated from Table 3.
- **Tyr:**
  - The Cantley-gap diagnosis (about a quarter from macro dilution; 42%; the GLO1 site labelled to 15 Tyr kinases; FGFR1/BTK/Eph/VEGFR/TAM) now sits behind a pointer to `docs/features_tyr.md` §1.
  - Hub and evidence-class numbers moved behind a pointer to S1 §2.2.
- **Understudied:** removed the Table 4 bin values and merged a sentence that appeared twice.
- **Regimes:**
  - Removed the per-kinase −0.012, λ = 300 and +0.029, which Table 5 gives.
  - The 2,808 / 8,987 pair counts now sit behind a pointer to S1.
- **Calibration:** the first-softmax values 0.87/0.54 and 0/41 now sit behind a pointer to S1 §4.4.
- **Coverage:**
  - Removed the release-file description, which Data availability gives.
  - Removed the TP53 refit detail; the caveat that the checks test recall of training labels is kept.
- **Downstream:** the 0.04–0.06 range and the evaluable-kinase counts now sit behind a pointer to `docs/downstream.md`.
- **Runtime:** removed the per-method seconds (Table 3 has them) and the node description (Methods › Seeds).
- **Post-TEST:** per-experiment CIs, the 14% → 2.3% shift and the 0.535 → 0.465 change now sit behind a pointer to S1 §5.

**Discussion**
- Removed repeats of Results values (B4 0.708, AUC 0.97, 0.947, 78-kinase numbers, T2 0.749/0.751, the canary).

**Methods**
- **R2:** the ranking-objective M2 0.086 now sits behind a pointer to S1 §3.4; the settings are in compact list form.
- **Priming:** the A32a numbers now sit behind a pointer to S1 §2.3.
- **Baselines:** removed "B1 thus scores 386 of 503", which the Table 3 footnote gives.
- **Evaluation:** replaced "top 20 hold half" with a pointer to the Introduction.
- **TEST batch:** removed "recomputed every endpoint", which Results states.

**Agent setup (shortened only)**
- The reviewer-specialty list and the package list now sit behind a pointer to `docs/DELEGATION.md`.
- Removed the register counts (about 20 owner roles, 4 literature agents, 11 reviewers).
- Removed a sentence that repeated the AI statement ("Per the statement above…").
- Compressed the "What went wrong" bullets and added a pointer to S1 §6. All five points are kept.

## Number check (all digit strings, set comparison)

- **supplement.md:** identical; 333 distinct numbers before and after.
- **manuscript.md:** 30 numbers disappeared. Each one is still in the supplement, a table or a repository file:
  - **Rounded copies of Table 3 runtimes:** 1,137 / 317 / 6,479 / 2,043 / 287 (Table 3 gives 1,136.6 / 316.9 / 6,478.8 / 2,043.4 / 287.0).
  - **In S1 Text:** 0.54, 0.87, 0.004, 0.008, 0.030, 0.086, 0.465, 0.535, 0.596, 0.627, 0.629, 2,808, 8,987, 0.0010, 0.0017, 0.0039, 0.0042, 0.0097, 0.0103, 0.0167.
  - **In `docs/downstream.md` and REPORT.md:** 73, 149, 274, 0.04, 0.06.
- **New digit strings (6):** 1.3, 1.5, 2.2, 2.4, 3.4, 4.1. These are only S1 section numbers in the new pointers.
- **Other checks:**
  - All 58 table rows are byte-identical.
  - All 46 headings are unchanged.
  - No existing code span was changed. Two new pointers were added: `docs/downstream.md` and `docs/features_tyr.md`.
  - The flag appears twice, as before.

## Citation check

- No reference lost a mention.
- The set of references cited only once is unchanged.
- No numbers were changed.
- **Already out of order before this edit (for the PI):** [49] (Phosformer-ST) is first cited in Methods › Splits and sealing, before [19]–[48]. Fixing it needs either renumbering (Phosformer-ST → 19, with 19–48 shifted) or an earlier first mention of [49]. I did neither, because renumbering was not allowed.
- All other references appear in order of first citation.
