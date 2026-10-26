# r1 review — prior work & baselines (role: prior_work_baselines)

Reviewer: independent agent, round r1, 2026-10-26.
Scope: are B1–B6 the right set; are the retraining plans fair and faithful; are we comparing on the same
pairs/coverage; what we must be able to show to claim we "beat" each.
Inputs read: `PROMPT.md`, `DESIGN.md` (v0), `LITERATURE.md`, `docs/lit/01_methods_benchmarks.md`,
`docs/data_profile.md`, `docs/data_acquisition.md`, `data/paper1.pdf` (SELPHI2.0 — note: the task brief has
paper1/paper2 swapped; paper1 = SELPHI2.0, paper2 = Phosformer), `data/paper2.pdf` (Phosformer),
and the already-downloaded `external/{kinase_library,selphi2,phosformer,kolossus,benchmarkin,ochoa2020}`.

All numbers below were computed in this review from the repo's own cached data
(`cache/psp_human_mapped_rows.parquet`, `external/*`), read-only, on CPU.

---

## 0. Measurements that drive the recommendations

| Quantity | Value | Source of number |
|---|---|---|
| Human-native PSP pairs (S/T/Y, canonical-mapped) | 14,644 over 420 kinases | computed |
| Distinct annotated human-native sites | 10,166 (S 6,299 / T 2,212 / Y 1,655) | computed |
| Kinases with ≥15 human-native pairs (⇒ ≥3 expected TEST positives at 20%) | **149** | computed |
| …of which 3–5 / 5–10 / 10–20 / ≥20 expected test positives | 43 / 44 / 34 / 28 | computed |
| Cantley atlas (kinase-library 1.8.0): matrices | 311 S/T + 93 Tyr (309 ser_thr, 74 tyr, 21 dual) | `kinome_information.tsv` |
| Pairs scoreable by Cantley, **residue-aware** | **97.3%** of all human-native pairs; **98.6%** within the 149 evaluable kinases | computed |
| Evaluable kinases with no atlas matrix at all | **2** (BCR/ABL, LRRK1) | computed |
| benchmarKIN NetworKIN 3.1 file | 207 kinase ids, 44,253 sites, 7,211 proteins | `external/benchmarkin/networkin_human_predictions_3.1.tsv` |
| Our positives whose **site** even appears in that NetworKIN file | **52.0%** (within evaluable kinases) | computed |
| Ochoa 2020 feature table (source of most SELPHI features) | 116,258 sites / 10,588 proteins | `external/ochoa2020/MOESM4,5` |
| Our positives covered by the Ochoa table | **44.6%** (S 49.4 / T 38.5 / **Y 33.7**) | computed |
| SELPHI2.0 feature columns | 43 (+4 id cols) in `data/psp/pos_set_features_20240910.tsv` | inspected |
| …regenerable now from already-downloaded data | **37/43** (24 from Ochoa tables, 4 derived from Ochoa cols, 2 Cantley PSSM, 4 Hijazi co-reg, 2 HPA RNA, 2 UniProt/COMPARTMENTS) | computed |
| …needing substitution / acquisition / identification | `coreg293` (Mertins 2016, not downloaded), `gtex` (GTEx v8, not downloaded), `score`/`sel.1`/`sel.2` (undocumented) | computed |

**The single most consequential finding:** the three baselines have wildly different coverage of our
evaluation universe (Cantley ~99%, Phosformer/KolossuS ~100% by construction, SELPHI-retrained ~45%,
NetworKIN ~52% of sites × 207/420 kinases), and DESIGN §8 has no rule for how to compare across that.
Without such a rule, every headline comparison is confounded by coverage, which is exactly the flaw the
SELPHI2.0 paper has (it compared against NetPhos on 17 kinases) and that `01_methods_benchmarks.md` §5.1(8)
identifies. This is PRIOR-1.

**The second:** DESIGN §7 F1 and `01_methods_benchmarks.md` §5.4 both assume coverage is "the main axis on
which we can beat" Cantley (the lit note predicts "~100+ of 421 test kinases with no atlas matrix"). On the
kinases that actually enter the primary endpoint this is false — 147/149 are in the atlas and 98.6% of their
pairs are residue-aware scoreable. We will have to beat Cantley on **accuracy**, and the coverage win only
exists in the long tail of kinases that have too few labels to be evaluated. The design should say so
plainly so that nobody builds a narrative on a coverage advantage that the metric cannot show.

---

## 1. Proposed changes

### PRIOR-1 — Add a baseline coverage protocol and a matched-universe comparison rule
**Severity: critical.** **Section: §6 (new subsection "Baseline coverage protocol"), referenced from §8 R0.**

**Change.** Add to §6:

> For every baseline *B* and every evaluation table, we report three things:
> (i) **Coverage** `cov(B)` = fraction of evaluated (k,s) pairs that *B* can score, broken down by residue
> and by kinase-label stratum, with the reason for each gap (kinase absent / site absent / feature absent).
> (ii) **Matched-universe metrics (primary for any head-to-head claim):** both our model and *B* are scored
> on the intersection universe `U_B` = pairs scoreable by *B*, with kinases re-filtered to ≥3 test positives
> *inside* `U_B`. This is the only comparison from which a "we beat *B*" claim may be made.
> (iii) **Full-universe metrics (secondary):** *B* evaluated on the full universe with unscored pairs
> assigned the worst rank within kinase (ties broken at random, seed-averaged). This is a deployment-realistic
> view, not a fairness view, and is always reported next to `cov(B)`.
>
> Following KinPred (doi:10.1371/journal.pcbi.1008681) all baselines are mapped to our reference proteome and
> our ID space before scoring; mapping losses count towards `cov(B)` and are reported, not silently dropped.

**Why.** Measured coverage (table above) spans 45%–99% across B1–B6; NetworKIN covers only 52% of our
positive sites and 207 of 420 kinases, and a faithful SELPHI retrain covers 44.6% of positives (33.7% of
Tyr). Comparing macro-AUROC across such different universes measures coverage, not accuracy, in either
direction. Explicitly declaring (ii) as the claim-bearing comparison prevents both over- and under-claiming,
and (i)+(iii) preserve the Priority-2 coverage story.

---

### PRIOR-2 — Split B1 into three declared arms; move the nearest-kinase imputation out of "as-published"
**Severity: critical.** **Section: §8 R0 B1 (and §7 F1).**

**Change.** Replace the B1 bullet with:

> **B1 Cantley Kinase Library (no label training).** Four arms, all reported:
> - **B1-raw**: log2 PSSM score (within-site comparator).
> - **B1-pub**: percentile against the *published* reference distributions (Ochoa 82,735 S/T sites;
>   5,431 Tyr sites with localisation prob > 0.75), i.e. exactly as the Kinase Library ships it.
> - **B1-recal**: percentile obtained by ranking the raw score **within kinase over our own U_obs**. This uses
>   **no labels** and is therefore still "as-published" in information terms; it removes the mismatch between
>   their reference phosphoproteome and our evaluation universe. **B1-recal is the primary Cantley arm** for
>   cross-kinase (macro) metrics; B1-pub is reported to show the size of the calibration effect.
> - **B1-cal**: per-kinase Platt/isotonic calibration of B1-raw fitted **on training-fold labels only**. This
>   *is* training and is labelled as such; it is the apples-to-apples comparator for our calibrated output and
>   must never be presented as "the published method".
>
> Each arm is run priming-naive and priming-aware (`pp=True`, lowercase s/t/y for neighbouring known
> phosphosites), with `st_fav` on, and with/without the 15 non-canonical Tyr matrices — these are
> pre-declared sub-arms, not post-hoc choices.
>
> **B1-ext** (nearest-atlas-kinase imputation by domain identity, currently written into §7 F1 and §8 B1) is a
> *modification* of the baseline and is reported separately from B1. It affects 2 of the 149 evaluable
> kinases (BCR/ABL, LRRK1) and ~1.4% of pairs, so it must not be used to make B1 look better-covered than it is.

**Why.** The published percentile is calibrated on a different site universe than ours; conflating that
mismatch with predictive accuracy would let us "beat" Cantley by a calibration artefact (or lose to it by
one). The three-way split makes the information content of each arm explicit. Separating B1-ext matters
because PROMPT.md's Priority-2 is coverage and we must not launder a coverage claim through the baseline.
The nearest-kinase imputation is also, strictly, label-free transfer learning — it belongs to *our* method
family, not to Johnson/Yaron-Barir's.

---

### PRIOR-3 — Pin the Phosformer reimplementation to the paper's actual recipe, and fix the two things that would make the retrain unfair
**Severity: critical.** **Section: §8 R0 B3.**

**Change.** Replace B3 with an explicit contract (all values from `data/paper2.pdf` §2.1–2.4):

> **B3 Phosformer (reimplemented, CC BY-NC-ND avoided by reimplementation).**
> Architecture: RoBERTa-style encoder, **6 attention layers, 12 heads, 768 dims**; input
> `<s>`11-mer`</s><s>`unaligned kinase-domain sequence`</s>`; `<pad>` where the site is within 5 residues of a
> terminus; classification head reads the embedding **at the central peptide residue**; softmax binary output.
> Fine-tuning: **focal loss γ=2**, AdamW **lr 2e-5**, early stopping patience 5.
> Augmentations (all kept, and ablated as one block): resample under-represented kinases (their rule:
> kinases >50 positives upsampled to 1,118, <50 positives upsampled to 50); random ±5-residue shifts of both
> kinase-domain boundaries, independently; 5% kinase-token `<mask>`; multilevel negative sampling (n easy +
> m hard per positive).
> Initialisation, three arms: (a) **released MLM checkpoint** (`external/phosformer/hf_phosformer`,
> pre-trained 50 epochs, lr 5e-5, on 295,320 NCBI-nr kinase domains + 5,338 substrate proteins) — the faithful
> arm; (b) ESM-2 (35M/150M) encoder; (c) random init. We do **not** redo their 8×A5000 pretraining.
> Two fairness fixes:
> - **(F1) Re-tune the negative budget on our CV.** Their 1:16 total pos:neg with their n:m easy:hard split
>   was tuned on *their* data and *their* 1:1 test set. Freezing it on our folds would handicap the baseline
>   relative to our models, whose negative ratio §5 says we tune. n, m and the total ratio are tuned on DEV CV
>   with the same budget we give our own models, and the chosen values are reported.
> - **(F2) Match the window.** Phosformer is native ±5 (11-mer); our F0/R3 use ±7/±15. B3 is run at ±5
>   (faithful) **and** at our window, so that a window-length advantage is not reported as an architecture win.
> Negatives (easy and hard) are regenerated **within training folds only**, from our site universe.
> **B3-pub**: the released *fine-tuned* checkpoint run zero-shot, labelled "trained on PSP/Phospho.ELM/
> PhosphoNetworks across 13 organisms → overlaps our test labels; optimistic ceiling, not a comparator."

**Why.** The design currently says "we re-run their recipe … or initialise from ESM-2, documented", which
leaves the two decisions that determine the result undetermined. (F1) is the main fairness risk: the paper
chose 1:16 by benchmarking on a *balanced* 1:1 pos:hard-neg test set (paper2 §2.1, §3.1), a regime that does
not exist in our evaluation; carrying their number over unchanged would be an unfair retrain. (F2) matters
because their own Discussion names input scope as their main limitation. Also keep their kinase-upsampling:
it is the closest published analogue of our §8 R2 `α` weighting and directly serves PROMPT.md's "harmonious
across kinases" requirement — so it should be a shared knob, ablated the same way in B3 and R2.

---

### PRIOR-4 — Make the SELPHI2.0 retrain concrete: 37/43 features are already regenerable; declare the 45% universe; keep its leaks but audit them
**Severity: major.** **Section: §8 R0 B4 (and §7 "Excluded as features").**

**Change.** Replace B4 with:

> **B4 SELPHI2.0 (retrained; feature-engineering baseline).** Model class kept exactly: scikit-learn
> RandomForest, grid search over `max_depth ∈ {10,20,50,70,80,100}`, `min_samples_split ∈ {8,10,12}`,
> `n_estimators ∈ {150,500,1000,1500}`, bootstrap=True, **100 bagged runs with resampled negatives, final
> probability = mean**; negative ratio arms 10×/20×/50× random pairs (they report 50×) **plus** one arm with
> our residue-matched hard negatives.
> Feature regeneration status (43 columns of `data/psp/pos_set_features_20240910.tsv`):
> - **24 join directly** from data already in `external/ochoa2020/`: 23 columns of MOESM4
>   (`is_disopred, disopred_score, exp3d_ala_ddG_effect, exp3d_acid_ddG_effect, log10_hotspot_pval_min,
>   isHotspot, isInterface, adj_ptms_w21, netpho_max_all, netpho_max_KIN, paxdb_abundance_log10, w0_mya,
>   w3_mya, quant_top1, quant_top5, PWM_max_mss, ACCpro, SSpro, SSpro8, sift_{min,mean,ala,acid}_score`)
>   plus `functional_score` from MOESM5.
> - **4 derived from Ochoa columns**: `isProteinDomain`←`domain_uniprot`, `isUniprotRegion`←`region_uniprot`,
>   `isELMLinearMotif`←`elm_*`, `isEV_ala_prediction_epistatic5`←`EV_ala_prediction_epistatic`.
> - **2 from the Cantley matrices we hold**: `pssm_ser_thr`, `pssm_tyr` (min–max normalised per kinase, as in
>   `alussana/selphi2.0-pssms`).
> - **4 from the provided compendium**: `coreg.ntera`, `coreg.hl`, `coreg.mcf7`, `isKinaseCoreg` — Hijazi 2020
>   is PMID 31959955, **183 conditions present in the compendium**, with `data/kinhib_inhib_comb__1.csv` for
>   inhibitor→kinase mapping.
> - **2 from HPA v23** (`rna.tiss`, `rna.cell`), **2 from UniProt/COMPARTMENTS** (`isProteinKinaseDomain`,
>   `isCytoplasmic`). → **37/43 regenerable with data in hand.**
> - **`coreg293`**: Mertins 2016 not downloaded → substitute CPTAC BRCA phospho×proteome (we hold it) and
>   record the deviation. **`gtex`**: GTEx v8 not downloaded → acquire, else substitute HPA consensus tissue
>   RNA and record. **`score`, `sel.1`, `sel.2`**: not defined in the paper → resolve against the Zenodo
>   feature table before the run; if unresolved, drop them and report the 40-feature variant.
> Mandatory deviations (documented as such, not as "their method"):
> 1. `coreg.*` are computed **per training fold** with the target site *and all other sites of the target
>    protein* removed from the kinase's KSEA regulon (CoPheeKSA's wording, doi:10.1038/s41467-025-57993-2).
> 2. Our substrate-cluster-grouped splits replace their `StratifiedKFold` over pairs.
> 3. Positives/negatives follow our label policy (L2 default).
> Universe: features derived from the Ochoa table exist for **44.6% of our human-native positives (Y: 33.7%)**.
> B4 is therefore compared under PRIOR-1(ii) on `U_B4`, and a **B4-lite** arm using only features computable
> proteome-wide (PSSM + co-reg + RNA + localisation + UniProt flags) is reported on the full universe.
> Leakage: `functional_score`, `PWM_max_mss` and `netpho_max_KIN` are themselves derived from known
> kinase–substrate annotations or from another predictor. They are **kept** in the faithful arm (it is their
> method) with one ablation arm removing them, and the delta is reported as the leakage estimate.

Also extend the §7 "Excluded as features" list from `Ochoa 2020 functional score` to explicitly include
`PWM_max_mss`, `PWM_nkinTop*`, `netpho_max_KIN`, `ptmdb_coreg_kinases*` and `pubmed_counts` from the same
table — these are label-derived or ascertainment proxies and are easy to import by accident when we join the
Ochoa table for the legitimate columns.

**Why.** The design says "its feature set regenerated fold-aware … adapted to our features where its raw
inputs are unavailable", which is unverifiable and leaves the fairness of the comparison to improvisation.
The audit above shows the retrain is ~86% regenerable today, which makes a *faithful* retrain achievable —
and it also exposes the 45% universe limit, which must be declared in advance rather than discovered during
evaluation. The in-advance declaration also protects against the temptation to evaluate SELPHI only where it
happens to look weak (Tyr, where its coverage is 33.7%).

---

### PRIOR-5 — KolossuS: establish the released checkpoint's training set before calling it a clean comparator
**Severity: major.** **Section: §8 R0 B6.**

**Change.** Replace B6 with:

> **B6 KolossuS (as-published, 650M variant).** Run via `kolossus-cli --model-small` with
> `kolossus-extract` ESM-2-650M embeddings; the 15B model is not required and is not in `external/kolossus/
> hf_weights/`. **Before quoting any number, establish from the preprint methods and the repo what the
> released checkpoint was fine-tuned on.** The weights shipped are
> `saved_model_cs_ws_ft_humansty.epoch_{55,59}.pth`; the `ft_humansty` suffix is consistent with fine-tuning
> on **human S/T/Y literature labels**, which would make this checkpoint contaminated with respect to our test
> set in exactly the way GPS 6.0 and the Phosformer release are. If contaminated, it is labelled
> "contaminated upper bound" and the published 0.755 (S/T) / 0.838 (Tyr) PSP numbers are **not** quoted as a
> bar we must clear; if it is confirmed to be atlas+Sugiyama only, it is our cleanest published comparator and
> those numbers stand as the reference point. Either way, report coverage (sequence-only ⇒ ~100%) and add an
> arm retraining the architecture on our folds if the as-published arm turns out contaminated.

**Why.** `docs/lit/01_methods_benchmarks.md` §4.5 calls KolossuS "the cleanest published comparator for our
setting" and names 0.755/0.838 as "the number to beat". That rests entirely on the claim that it never saw
PSP labels. The filename of the checkpoint we have downloaded contradicts that reading strongly enough that
we must check before building a headline claim on it. This is a cheap check (read the preprint methods) with
a large downside if skipped.

---

### PRIOR-6 — Fix the baseline set: add B0 (prior-only) and B7 (motif+minimal context); keep GPS 6.0; demote NetworKIN; drop LinkPhinder; move KSTAR downstream
**Severity: major.** **Section: §8 R0 (whole block) and §12 step 2.**

**Change.**
> - **B0 Prior-only baselines, reported as a first-class row in every results table** (not only as a §11
>   diagnostic): (a) kinase-frequency prior (rank kinases by training-label count, identical for every site);
>   (b) site-popularity prior (rank sites by number of distinct kinase annotations / observation count,
>   identical for every kinase). Any method whose margin over B0 on E-hard is small is reporting study bias.
> - **B7 Cantley ⊕ minimal context**: per-kinase logistic regression on {B1-recal percentile, co-expression,
>   localisation similarity, systematic-PPI flag}, fitted fold-aware. Cheap (minutes) and it is the honest bar
>   that SELPHI and NetworKIN actually represent: if the full F3–F6 stack cannot beat B7, the stack is not
>   earning its complexity or its runtime.
> - **GPS 6.0 — keep as a baseline** (offline package at gps.biocuckoo.cn; not yet in `external/`). It was the
>   strongest comparator in the SELPHI2.0 external benchmark (AUC 0.81 vs 0.83). Run as-published; label
>   "trained on PSP/GPS-gold ⇒ our test labels are in its training set; optimistic ceiling". Do not retrain
>   (577 transfer-learned models, no training pipeline released).
> - **NetworKIN — keep but demote** to the coverage-restricted table and the downstream benchmarKIN
>   evaluation only. Two independent disqualifications from the headline table: the benchmarKIN file covers
>   **207 kinase ids / 44,253 sites / 52% of our positive sites**; and its network term is STRING, which our
>   own leakage audit (`LITERATURE.md` §4) shows is contaminated with curated KS edges (34% of PSP pairs have
>   a STRING score vs ~4% of controls). Report it as "coverage-limited and label-contaminated".
> - **LinkPhinder — exclude**, with a one-line justification in §8: predictions download is dead, and its own
>   independent evaluation collapsed to AUROC 0.532 (with GPS 0.509, NetworKIN 0.511, NetPhorest 0.510,
>   Scansite 0.502). Our B0 reproduces its failure mode at zero cost, so nothing is lost. Cite it in the
>   report as the cautionary case for pair-random splits.
> - **KSTAR — not a KS-prediction baseline.** Move to §6 "Downstream": use it (and its NetworKIN-pruning
>   libraries) as a comparator *library* inside benchmarKIN, and adopt its paralog-redundancy finding
>   (homologous kinases share up to 80% of site evidence, pruned to 18%) as the target for the §4.5 audit.

**Why.** The current B1–B6 list omits the two baselines that most constrain our claims: a prior-only control
(the LinkPhinder lesson — on PSP random splits, degree and study bias reach 0.97) and a motif+context
baseline (the actual competitor to our integrative stack). It simultaneously spends effort on NetworKIN,
whose coverage we have now measured at 52% of sites and whose network input we have already banned as a
feature for leakage, and names KolossuS/NetworKIN as peers of Cantley when their evaluation status differs.

---

### PRIOR-7 — State that the primary endpoint is a macro over ~149 kinases, and add a minimum-pool rule
**Severity: major.** **Section: §6 (Kinase-centric / Primary endpoint).**

**Change.** Add to §6:

> With the ≥3-test-positive rule, a 20% cluster split and human-native test labels, the primary endpoint is a
> macro average over an expected **~149 kinases** (kinases with ≥15 human-native pairs; measured on the
> provided PSP release), of which **43 fall in the 3–5 test-positive bin** where per-kinase AUROC has a
> standard error of roughly 0.19–0.25. Consequences, pre-registered:
> - the evaluable-kinase count and its composition are reported in every table;
> - the 3–5 bin is always shown as its own stratum, and no per-kinase claim is made for it;
> - a kinase enters the macro only if its E-hard negative pool has **≥50 sites** of the matching residue class
>   in the evaluated split (otherwise it is excluded and counted; for reference, the expected TEST pool is
>   ~2,033 annotated sites, of which ~331 are Tyr);
> - all method comparisons are made on the **paired cluster bootstrap of the macro metric**, never by counting
>   per-kinase wins.

**Why.** §6 defines the primary endpoint without saying how many kinases it averages over; the answer (~149,
a third of the 420 with any label) materially changes how the result reads and how much resolution we have.
The ~331 Tyr test sites also mean the Tyr stratum and the Tyr-specific baselines (Yaron-Barir atlas,
KolossuS-Tyr) will be the noisiest part of the whole evaluation — worth knowing before, not after.

---

### PRIOR-8 — Pre-register explicit superiority criteria for "we beat X"
**Severity: major.** **Section: §6 (new subsection "Superiority criteria").**

**Change.** Add:

> We claim to beat baseline *B* only if **all** of the following hold, computed once at final evaluation:
> 1. On the matched universe `U_B` (PRIOR-1 ii), the paired cluster-bootstrap 95% CI for Δ(macro per-kinase
>    AUROC, E-hard) excludes 0, in our favour.
> 2. No deterioration on the co-primary (macro per-kinase AUPRC, E-all) beyond its CI.
> 3. `cov(ours) ≥ cov(B)` on the full universe, stated numerically.
> 4. The win is not carried by well-studied kinases: Δ is reported per label-count stratum (0 / 1–4 / 5–19 /
>    20–99 / ≥100) and must be ≥0 in the 5–19 and 20–99 strata.
> 5. Claims under leave-kinase-out (T2/T3) are made separately and never generalised from T1.
> 6. For each contaminated baseline (Phosformer-released, GPS 6.0, possibly KolossuS, NetworKIN-via-STRING)
>    the claim is phrased as "we beat a baseline that has seen our test labels", which is a stronger statement
>    and must be worded that way rather than hidden.
> Additionally, prior numbers are never cited as the bar: Phosformer's 0.885 is on a 1:1 pos:hard-negative
> test set in which (our own audit, `scripts/phosformer_split_leakage_check.py`) 51.1% of test-positive
> peptides are training positives under another kinase; SELPHI2.0's 0.97 is pooled pair-level CV with
> label-derived co-regulation features; Phosformer-ST's 0.992 is distillation of the Cantley PSSMs. Only
> numbers we recompute on our splits enter comparison tables.

**Why.** PROMPT.md requires us to "beat established benchmarks" and forbids changing metrics after seeing
results. Without a written criterion, "beat" will be decided after the numbers arrive. Criterion 4 is the one
that enforces the brief's "harmonious across kinases" requirement at the level of the claim, not just the
reporting.

---

### PRIOR-9 — Strengthen B2 so it is not a straw man
**Severity: minor.** **Section: §8 R0 B2 (and §7 F2a).**

**Change.** B2 becomes: per-kinase log-odds PSSM from training positives with a **Dirichlet prior from the
kinase's Manning family profile**, prior weight `w` tuned on DEV CV; background from residue-matched
unlabelled U_obs sites; minimum 5 training positives, else fall back to the family-level PSSM and flag.
Report B2 both with `w` tuned and with `w=0` (naive counts), and report how many kinases fall back.

**Why.** §8 B2 currently says only "per-kinase PSSM learned from training labels". An unregularised
count PSSM on a kinase with 5 positives is close to noise; beating it would be meaningless, and the bulk
of our evaluable kinases (87 of 149 have <10 expected test positives) sit exactly in that regime. §7 F2a
already specifies a family Dirichlet prior — B2 should inherit it so that the baseline and the feature are
the same object.

---

### PRIOR-10 — Record the licence and contamination status of each baseline in the design, not only in lit notes
**Severity: minor.** **Section: §8 R0 (preamble) and a new row set in `DATA.md`.**

**Change.** Add a table to §8 with one row per baseline: {code licence, weights licence, trained-on, does its
training set intersect our TEST labels (yes/no/unknown), redistribute in our repo (yes/no)}. Known values:
Phosformer code CC BY-NC-ND (⇒ reimplement, do not redistribute modified weights; the *article* is CC BY);
Phosformer-ST CC BY-NC-ND; `alussana/selphi2.0-pssms` and `rohitsinghlab/KolossuS` have **no licence file**
(run locally, never redistribute); SELPHI2.0 main repo Apache-2.0; kinase-library CC BY-NC-SA 3.0; CoPheeMap
MIT; KinPred MIT. Contamination: Phosformer-released **yes** (PSP + Phospho.ELM + PhosphoNetworks, 13
organisms); GPS 6.0 **yes**; NetworKIN **yes, indirectly** (STRING); KolossuS **unknown — PRIOR-5**;
Cantley **no** (PSPA arrays only; PSP enters only their evaluation, and the percentile reference is the
Ochoa MS phosphoproteome, which overlaps our universe but carries no K–S labels).

**Why.** PROMPT.md requires DATA.md to carry source/version/licence per dataset; baselines are datasets and
code with redistribution constraints, and the contamination column is the thing a reader needs to interpret
every number in the comparison table. Keeping it in `docs/lit/` only means it will not reach REPORT.md.

---

### PRIOR-11 — Hold library size fixed in the benchmarKIN downstream comparison
**Severity: minor.** **Section: §6 (Downstream).**

**Change.** Add: the downstream comparison uses **matched substrate-set sizes** — top-N predicted substrates
per kinase for N ∈ {10, 25, 50, 100} — against PSP-only, Cantley-top-N and NetworKIN-top-N libraries, and is
reported both on all kinases each library covers **and** on the common-kinase subset. Report the
kinase-coverage of each library alongside its AUROC.

**Why.** Müller-Dott 2025 found that libraries matter more than inference methods, and that on a common
17–31-kinase subset NetworKIN (0.77) and OmniPath (0.76) beat the curated combination (0.75) — i.e. the
headline ordering flips with coverage. An unmatched comparison here would let us "win" by shipping a larger
library, which is not evidence of accuracy. §6 currently says "top-N predicted substrates per kinase" without
fixing N or the kinase set.

---

### PRIOR-12 — Add one "their-labels" arm per retrained baseline to separate label quality from model quality
**Severity: minor.** **Section: §8 R0 (B3, B4) and §12 steps 5–6.**

**Change.** Run B3 and B4 twice: once with our L2 labels, and once with each method's **native** label policy
mapped onto our splits (Phosformer: PSP+Phospho.ELM+PhosphoNetworks, 13 organisms, no evidence filter;
SELPHI2.0: human×human PSP only, no evidence filter, 14,542 pairs). Report both. (Phospho.ELM is behind an
academic click-through per `docs/data_acquisition.md`; if it is not obtained, use the PSP+PhosphoNetworks
subset and say so.)

**Why.** If we win, the reader cannot otherwise tell whether it is our model or our label curation. PROMPT.md
asks for a fair retrain of these benchmarks and separately asks us to explore label quality; this arm answers
both at a cost of two extra training runs (Phosformer ≈ a few GPU-h; SELPHI RF ≈ CPU-hours).

---

## 2. Risks I could not resolve

1. **SELPHI2.0's stated universe is internally inconsistent with the data we hold.** The paper says
   predictions span 238,374 phosphosites "listed in PhosphoSitePlus … and had a functional score assigned to
   them", but the Ochoa 2020 functional-score table has 116,258 rows (10,588 proteins). Either they used a
   different/extended Ochoa product or the 238,374 figure means something else. This changes B4's universe by
   a factor of two. Resolve from their Zenodo `feature_table.tsv` header / `positive-set-20240910.tsv` before
   fixing `U_B4`; do not guess.
2. **`score`, `sel.1`, `sel.2`** in their feature file are undocumented in the paper; I could not determine
   what they are. If they turn out to be label-derived, B4's leakage profile is worse than estimated.
3. **KolossuS checkpoint provenance** (PRIOR-5) — unresolved from the files in `external/kolossus/`; the
   README and HF card say only "MIT".
4. **GPS 6.0 has not been downloaded**, so its kinase set (349 claimed) and its HGNC mapping are unverified;
   its coverage of our evaluable kinases is an estimate, not a measurement.
5. **Cantley percentile self-reference**: sites in our test set that are also in the Ochoa reference
   distribution contribute to their own kinase's background. The effect is small but I did not quantify it;
   B1-recal (PRIOR-2) sidesteps it for the primary arm.
6. **Phospho.ELM licence** blocks an exact reproduction of Phosformer's label set (PRIOR-12).
7. **Whether 149 evaluable kinases is enough** for the stratified claims in PRIOR-8 criterion 4: the 5–19 and
   20–99 strata will hold roughly 40–50 kinases each, so a stratum-level Δ CI will be wide. I did not simulate
   the power; a statistics reviewer should.
