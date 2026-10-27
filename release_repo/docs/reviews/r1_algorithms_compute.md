> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Review r1 — algorithms & compute

Reviewer role: ML engineer (PLMs, GBDTs, multitask/metric learning, PU learning, protein–peptide models).
Target: `DESIGN.md` v0 (2026-10-26). Inputs read: `PROMPT.md`, `DESIGN.md`, `LITERATURE.md`,
`docs/lit/01_methods_benchmarks.md`, `docs/lit/04_other_modalities.md`, `docs/data_profile.md`.

Verdict: the design's *evaluation* scaffolding is strong; its *modelling* sections (§5 negatives, §8 regimes,
calibration) and its *compute* section (§10) are the weak parts. §10 in particular is off in both directions —
it over-budgets the things that are free and under-budgets the only thing that is actually expensive (pair-feature
materialisation). Everything below is backed by measurements I ran on this machine.

---

## 0. Measurements I ran (read-only, this box: 4× RTX 4090 24 GB, 64 cores, 375 GB RAM)

Scripts in `<agent-scratchpad>/`
(`bench_esm.py`, `bench_coreg.py`, `bench_lgbm2.py`).

| What | Measured | Implication |
|---|---|---|
| ESM-2 35M, fp16, 31-mer windows, 1 GPU | **19,700 seq/s** (B=2048–4096); 43 TFLOP/s effective ≈ 52% MFU | — |
| ESM-2 35M, fp16, L=300 | 2,005 seq/s | — |
| → extrapolated ESM-2 650M (×19.6 non-embed params, same MFU) | ≈ **1,000–1,300 windows/s**; full 20.4k-protein proteome at **full-length context ≈ 8–15 GPU-min** | PLM embeddings are *not* a budget item |
| → extrapolated ESM-2 3B | ≈ 250–300 windows/s; full proteome ≈ **40–70 GPU-min** | 3B is cheap too — the argument against it is not compute (CMP-3) |
| Pairwise-complete Pearson, real `cache/compendium_wide.parquet` (189,962 × 1,937, 95.4% missing), masked-matmul on 1 GPU, fp32, min-n=8 | **0.71 G pairs/s**; full **190k sites × 13.2k labelled sites in 3.6 s**, 10.4 GB peak VRAM | §10's "co-regulation is the heaviest, vectorised on GPU if needed" is wrong by ~3 orders of magnitude. It is free (CMP-4) |
| LightGBM train, 300k rows × 120 feats, 600 trees/63 leaves, 16 threads | **7.5 s** | a small-config R2 fold is seconds |
| LightGBM train, 400k rows × 160 feats, 1000 trees/127 leaves, 64 threads (cores shared with other agents at the time) | **888 s** | training cost is very config- and contention-sensitive: budget 10 s – 15 min per fold; run ≤16 threads/job, 4 jobs in parallel |
| LightGBM predict, 64 threads, 120 feats, 600 trees / 160 feats, 1000 trees | **1.19 / 0.94 M rows/s** → **129M pairs in 1.8–2.3 min** (16 thr: 7.9 min; 4 thr: 31 min) | tree scoring is cheap; §10's "~1 h CPU" is accounting for the wrong thing |
| PSP HMR label structure (`cache/psp_human_mapped_rows.parquet`) | 19,477 HMR pairs / 14,644 human-native; 442 / 420 kinases ≥1 | — |
| Kinases with ≥3 **human-native** test positives at a 20% protein-level holdout, 5 seeds | **mean 163 (range 149–184)** of 442; ≥5 pos: 115; ≥10 pos: 68 | the §6 primary endpoint is computed on ~37% of the kinome (ALG-8) |
| Paralog contamination of the §5 hard-negative pool | MAPK1/MAPK3 share **398** sites (Jaccard **0.43**); at 1:20 + 50% hard, CDK1 draws **~241** CDK2-only sites and MAPK1 **~103** MAPK3-only sites as "negatives" | ALG-2 |

---

## Proposed changes

### ALG-1 (critical) — §5 "Training negatives" + §8 "Calibration": decouple the training ratio from the deployment prevalence, and correct the prior explicitly

**Problem.** §5 fixes 1:20 pos:neg. The deployment prevalence is ~4 orders of magnitude lower: a median S/T
kinase has ~10 HMR positives against a U_obs S/T universe of ~255k sites (π ≈ 4e-5); a Y kinase ~10 against
~45k Y sites. Training at 1:20 shifts every score's log-odds by `log(0.05/4e-5) ≈ +7.1 nats`. §6's co-primary
(macro AUPRC under E-all) and §6's coverage endpoint ("precision ≥ 0.5 / ≥ 0.8") are both *absolute*-probability
claims, so they are governed entirely by this constant, not by the model. §8 currently repairs it with "per-kinase
isotonic … with prevalence correction", which is under-specified and, per ALG-7, not estimable.

**Change.** In §5 and §8:
1. State the sampling scheme as **case–control sampling with known rates** and record, per kinase *k* and per fold,
   the negative sampling rate `ρ_k = n_sampled_neg_k / n_unlabelled_k`. Apply the closed-form correction
   (Elkan 2001 / King & Zeng 2001) to the training logit before any calibrator:
   `logit_corrected = logit_train + log ρ_k`. This is exact under random sampling and must be applied *before*
   the kinase-level calibrator, so the calibrator only has to fix curvature, not a 7-nat offset.
2. Raise the default ratio to **1:50** and tune on {1:20, 1:50, 1:200}. 1:200 (≈3.3M rows) costs about 1.5 min
   per fold with a small tree config but could take tens of minutes with large trees (see §0). Use early stopping
   for that arm.
3. **Compute all CV ranking metrics on the full unlabelled universe, not on the sampled negatives.** At 1.19 M
   rows/s this is 1.8 min for the whole 129M-pair grid, so there is no reason to approximate E-all by a sample.
   Add one sentence to §6 saying E-all is evaluated exhaustively.
4. For the mined-negative arm (ALG-5), the sampling is no longer random, so `log ρ_k` no longer applies —
   in that arm the correction must come from a held-out *uniformly* sampled negative set. State this.

**Why.** Without (1) the "precision ≥ 0.8" coverage deliverable is not interpretable, and R2 vs R3 comparisons
under AUPRC are confounded by whatever ratio each happened to use.

---

### ALG-2 (critical) — §5 "Training negatives": exclude paralogue and same-cluster sites from hard negatives

**Problem.** §5 defines hard negatives as "sites labelled to other kinases", with no paralogue exclusion. §4.5
already recognises paralogue redundancy, but only for *evaluation*. Measured on our own label table:
MAPK1 and MAPK3 share 398 of their 769/563 annotated sites (Jaccard 0.43); CDK1/CDK2 share 124. At the §5 default
(1:20, 50% hard), **CDK1's negative set would contain ~241 CDK2-annotated sites and MAPK1's ~103 MAPK3-annotated
sites** — for a pair with 43% label overlap, the large majority of those are true positives mislabelled as the
hardest negatives. This systematically teaches the model to *separate* paralogues that are biochemically
inseparable, which is the opposite of what §4.5's audit asks for and will hurt exactly the high-n kinases that
dominate E-hard.

**Change.** In §5, add to the negative sampler, in this order:
- exclude any site annotated (in the training fold, any label variant, any species) to a kinase with **≥70%
  kinase-domain identity** to *k*, or in the same mmseqs 50%-identity kinase cluster used by §4.4;
- exclude any site in the **same substrate cluster** as a training positive of *k* (paralogous substrate sites);
- exclude any site whose **15-mer window is within Hamming distance 2** of a training positive window of *k*
  (the Phosformer-ST dedup rule, applied to negatives rather than to splits);
- report the number and fraction excluded per kinase in `results/leakage/`.
Mirror the first bullet in the **E-hard evaluation negatives**, otherwise E-hard penalises correct paralogue
predictions; keep the unfiltered version as a declared sensitivity arm.

---

### ALG-3 (major) — §8 R1/R2/R4: replace three regimes with one partial-pooling family indexed by a shrinkage parameter

**Problem.** As written, R1 ("LightGBM/logistic per kinase with ≥20 training labels") is a straw man. Only 126
kinases have ≥20 human-native pairs (147 HMR); a LightGBM on 20 positives + 400 negatives is pure variance, so
"R1 < R2" is a foregone conclusion and the experiment in §12.7 buys no information. R4 ("blended with weight
depending on n_k") leaves the blending rule unspecified, which is the only part that matters. PROMPT.md
explicitly asks for "something in between" to be *tested*; this needs a parameterisation, not three labels.

**Change.** Rewrite §8 R1/R4 as a single hierarchical model:
- R2 (global pair model over all kinases) stays as the pooled component, unchanged.
- **R1′ = per-kinase ridge logistic residual** on a small, fixed, ~8-feature vector (Cantley score/percentile,
  learned-PSSM log-odds, F6a co-regulation, F6b inhibitor contrast, F3 RWR proximity, F4 localisation cosine,
  F0 disorder/pLDDT, R3 score), fitted on the R2 logit as an offset, with the coefficient vector shrunk toward
  the **Manning-family mean coefficient vector**: penalty `λ‖β_k − β_family‖²`. Logistic, not trees — with ~10–100
  positives a tree model cannot estimate interactions, and the shrinkage target is only definable for a linear
  parameterisation.
- **R4 = the blend with an empirical-Bayes weight** `w_k = n_k / (n_k + λ)`, so that `λ → ∞` recovers R2 exactly
  and `λ → 0` recovers fully per-kinase fitting. Tune the single scalar λ (and the family-shrinkage λ) on DEV CV.
- Report the fitted λ and the implied `w_k` per label-count bin. That number ("per-kinase information starts to
  pay off above n_k ≈ X") is a genuinely publishable result and directly answers PROMPT.md §18.
This is cheaper than the current plan (one λ sweep instead of three model families). The λ grid needs only refits of
the per-kinase ridge models on top of a fixed R2 offset, which takes seconds, so the grid costs almost nothing.

---

### ALG-4 (major) — §8 R3: the neural model is one sentence; replace it with a concrete, zero-shot-robust spec

**Problem.** §8 R3 currently says "a site encoder (ESM-2 embedding of the ±15 window, or a fine-tuned small ESM)
and a kinase encoder (ESM-2 kinase-domain embedding + Cantley vector), combined by a bilinear / MLP head".
Three decisions in there are likely wrong, and one is unspecified but decisive:
- **±15 window-only embeddings.** ESM-2 is trained on whole proteins; a free-floating 31-mer is out of
  distribution, and it discards exactly the signal that motif PSSMs cannot capture (docking motifs sit 10–40 aa
  away for Cy/RxL, anywhere in the protein for MAPK D-motifs, per `docs/lit/04` §4.3). Full-length context costs
  **8–15 GPU-min for the whole proteome** (measured extrapolation), so there is no reason to use window-only.
- **"or a fine-tuned small ESM"** left as an option. With ~15.8k DEV positives, full fine-tuning of a 650M encoder
  is ~41,000 parameters per positive example. KolossuS's own ablation (atlas-only 0.754 → +Sugiyama fine-tune
  0.762 AUROC) says the marginal return on scaling/fine-tuning the PLM here is ~0.01 AUROC.
- **Kinase encoder = mean-pooled ESM-2 domain embedding + Cantley vector.** With 442 kinases this is a 442-point
  lookup table and the model will memorise kinase identity; nothing in the architecture forces the descriptors to
  be used, so T2/T3 (unseen kinase family) performance will collapse. This is the single biggest risk to the
  zero-shot part of the brief.

**Change.** Specify R3 in §8 as:
- **Site representation = concatenation of three blocks**: (i) ±7 one-hot (20 aa + gap) **plus an explicit
  phospho-channel** marking neighbouring observed phosphosites (the §9 priming requirement — PLMs have no
  phospho token, so this must be an added channel, not a sequence hack); (ii) frozen ESM-2 650M per-residue
  embedding at the site, taken **from full-length protein context** (sliding 1,022-aa windows, site most central);
  (iii) frozen ESM-2 mean-pool over ±25 and over the whole protein (the docking-motif / protein-context channel).
  Blocks (ii)–(iii) are precomputed once and cached as fp16 (~4.6 GB for S/T/Y residues).
- **Kinase representation = a block-structured descriptor vector**: [Cantley PSSM flattened → PCA 32d |
  specificity-determining residues one-hot from the KinCore alignment | Manning group/family one-hot |
  frozen ESM-2 kinase-domain mean-pool → PCA 64d]. The Cantley block is the only non-label-derived, directly
  specificity-informative descriptor and covers ~396/442 kinases.
- **Descriptor dropout (the key trick): during training, zero out whole descriptor blocks independently with
  p = 0.2, and additionally zero the family/group one-hot with p = 0.2.** This is standard modality dropout and
  costs nothing; without it, nothing prevents identity memorisation and the T2/T3 numbers are not interpretable.
  Pre-register a "kinase-identity-free" ablation (train with the domain-embedding block always dropped) as the
  honest zero-shot ceiling.
- **Head**: low-rank bilinear (site_proj ∈ R^128, kinase_proj ∈ R^128, score = ⟨W_s x, W_k z⟩ + MLP([x; z; x⊙z]))
  rather than a pure cosine head (KolossuS's cosine head is what limits it to 0.755; it cannot express
  "kinase k dislikes residue R at +1").
- **Default = frozen encoders.** LoRA fine-tuning (rank 8, top 6 layers only) of ESM-2 650M is a *pre-registered
  ablation* with a decision rule, not the default. Budget: see CMP-1.

---

### ALG-5 (major) — §8 R3 loss: add a listwise objective over kinases, PU-aware negatives, and mined hard negatives

**Problem.** §6 declares a site-centric endpoint (Hit@1/@5/@10, MRR over the ~430 candidate kinases for a site),
but **no model in §8 optimises it**. R3 is described only as "trained multi-task across all kinases with hard
negatives", i.e. pointwise binary loss, which optimises neither of §6's primary endpoints (within-kinase ranking)
nor the site-centric one. §5's static 50/50 easy/hard split is also the weakest of the negative schemes surveyed
in `docs/lit/01` §5.5.

**Change.** Specify the R3 objective as a weighted sum of three terms, with the weights tuned on CV:
1. `L_pair` — binary cross-entropy with `pos_weight` (not focal loss alone; focal γ=2 at 1:50 imbalance
   down-weights the easy negatives that carry the prevalence information needed by ALG-1).
2. `L_site` — **InfoNCE / softmax cross-entropy over the kinase axis**: for each training positive site *s* with
   true kinase *k*, softmax over {k} ∪ {sampled 63 kinases of matching residue class, excluding paralogues of *k*
   per ALG-2}. This is the loss whose population optimum is the §6 site-centric ranking.
3. `L_kinase` — sampled listwise loss over the site axis within a kinase (e.g. `rank_xendcg`-style or a sampled
   softmax over 1 positive + 255 negatives), which matches the §6 primary (per-kinase AUROC/AUPRC).
Plus:
- **Periodic hard-negative mining**: at the end of each epoch, score a 1% random slice of the full unlabelled
  universe with the current model and refresh 25% of each kinase's negative pool from the top-scoring unlabelled
  sites, after applying every ALG-2 exclusion. Measured inference throughput makes a full re-mining pass cost
  seconds; this is the cheapest remaining source of AUPRC-at-the-top.
- **PU-awareness**: mined hard negatives are enriched for unlabelled *true* positives by construction. Add
  **nnPU** (Kiryo et al. 2017) with class prior π estimated per kinase group, or at minimum down-weight mined
  negatives by `(1 − π̂)`, and run the choice {plain BCE | nnPU | label-smoothed negatives ε=0.05} as a declared
  3-arm ablation. Given 95% of PSP's supported human phosphosites have no kinase annotation (`data_profile` §2),
  the "unlabelled ≠ negative" assumption is not a technicality here.
- Apply the same three-term objective to **R2** where possible: LightGBM `lambdarank` grouped by kinase gives
  term 3; grouping by *site* gives term 2. §8 currently mentions "binary loss vs lambdarank grouped by kinase"
  only — add the site grouping as a third arm. Caveat to record: LightGBM's lambdarank needs group-sorted data
  and degrades with very large groups, so cap group size by subsampling negatives within group.

---

### ALG-6 (major) — §8 R3→R2 stacking: "out-of-fold within DEV" is not sufficient; specify nested OOF

**Problem.** §8 says R3's score "feeds R2 (stacked; out-of-fold within DEV to avoid leakage)". The naive reading
— train R3 on DEV\j, predict fold j, use that column everywhere — leaks: when R2 is trained on DEV\j, the R3
feature values for those same rows were produced by a model that saw their labels. The stacked model then learns
"trust R3", and the fold-j R3 column (produced without fold-j labels) is distributed differently from the training
column. This is the classic stacking leak and it produces a CV>TEST drop that §11's leakage alarm would flag only
*after* the compute is spent.

**Change.** Write the exact procedure into §8:
> For each outer DEV fold *j*: (a) partition DEV\j into 4 inner folds; train R3 on each 3-of-4 inner subset and
> predict the held-out inner fold, giving the R3 feature column for all rows of DEV\j; (b) train R3 once on all
> of DEV\j and predict fold *j*, giving the R3 column for fold *j*; (c) train R2 on DEV\j using the inner-OOF
> column and evaluate on fold *j*. Total: 5 × (4 + 1) = **25 R3 trainings**. Apply the identical scheme to any
> other stacked component (the learned-PSSM F2a score, the SELPHI-style RF, Phosformer).

Add the cheaper alternative as a declared fallback: **skip stacking and blend** — fit a 2–3 parameter logistic
blend of (R2 logit, R3 logit, Cantley percentile) on the inner folds only; costs 5+5 trainings and cannot leak
through the feature matrix. Given the measured R3 cost (CMP-1), the nested version is affordable; say which one
is the default so the comparison is not decided ad hoc.

---

### ALG-7 (major) — §8 "Calibration": per-kinase isotonic is not estimable for most kinases

**Problem.** §8 specifies "per-kinase isotonic / Platt on out-of-fold CV predictions". Median kinase has ~10 HMR
positives *in total*; within DEV-CV the number of OOF positives per kinase is that times 0.8. Isotonic regression
fitted on 8 positives is a step function with ≤8 jumps and will produce wildly overconfident calibrated
probabilities; the §6 coverage deliverable ("number of sites with ≥1 prediction at estimated precision ≥0.5 and
≥0.8") then rests on an unestimable quantity. Only 62 kinases have ≥50 human-native pairs, which is roughly where
isotonic starts to beat a 2-parameter fit.

**Change.** Replace the §8 calibration sentence with a tiered rule:
- **n_k^OOF ≥ 50** (≈62–76 kinases): per-kinase isotonic, as now.
- **5 ≤ n_k^OOF < 50**: **beta calibration** (3 parameters, Kull et al. 2017) with its parameters shrunk toward a
  pooled estimate fitted across all kinases in the same label-count bin and kinase group — i.e. a hierarchical
  calibrator, not 400 independent ones.
- **n_k^OOF < 5**: use the pooled bin-level calibrator only; do not fit anything per-kinase.
- In all tiers, apply the ALG-1 `log ρ_k` prior correction *before* the calibrator.
- Report **ECE and a reliability curve per label-count bin** (not per kinase), and restrict the §6 "precision ≥0.8"
  coverage claim to the ≥50 tier, reporting the other tiers as percentile-threshold coverage only (§6 already has
  "top 1% per kinase" — make that the headline for low-n kinases and say so).

---

### ALG-8 (major) — §6/§8: tie the kinase weighting α to the metric it is supposed to optimise, and state the effective per-kinase weight

**Problem.** §8 R2 sets "sample weights equalising kinases (∝ 1/n_k^α, α ∈ {0, 0.5, 1} tuned)". Two issues.
(a) The weight is specified on positives, but §5 already equalises kinases through the negative sampler
("1:20 per kinase, min 50 negatives"), so the *effective* per-kinase weight in the loss is
`w_k · n_k · (1 + max(20, 50/n_k))` — not `n_k^{1−α}`. For a kinase with 3 positives the min-50 floor alone gives
it 53 rows vs 20 for an unweighted 1:20 draw, a 2.6× implicit up-weight that interacts with α in an undocumented
way. (b) Tuning α against the §6 primary is circular in a specific way: I measured that **only ~163 kinases
(range 149–184 over 5 seeds) have ≥3 human-native test positives**, i.e. the "kinase-equal" primary endpoint is
in fact an average over the best-studied 37% of the kinome, so α=1 optimises a population the metric does not
cover.

**Change.**
- In §5/§8, define and log one quantity, `W_k` = the total loss weight kinase *k* receives (positives + negatives,
  after sampling and weighting), and tune **that** on a declared target (e.g. `W_k ∝ n_k^{1−α}`), rather than
  tuning α on top of an uncontrolled sampler.
- In §6, state explicitly that the primary endpoint is computed over ~150–185 kinases and report that count with
  every macro number; add as a **co-primary for coverage** a metric in which every test positive contributes —
  the §6 site-centric MRR/Hit@k, macro-averaged over *label-count bins* rather than restricted to kinases with
  ≥3 positives. Select α on the stratified table (the 1–4 and 5–19 bins), not on the overall macro number.
- (Flagging for the statistics reviewer: with 3–5 positives, a per-kinase AUROC has an SE of roughly 0.1–0.15,
  so the 1–4 bin is close to uninformative per kinase even if the macro average over ~50 such kinases is not.)

---

### CMP-1 (major) — §10 "Compute plan": rewrite with measured numbers; the real constraint is pair-feature materialisation, not GPU

**Problem.** §10 is wrong in both directions. It says co-regulation is "the heaviest" (measured: **3.6 s**), that
ESM-2 650M embeddings are "<1 GPU-h" (measured extrapolation: **8–15 GPU-min**, with full-protein context), that
LightGBM is "minutes per fold" (measured: 7.5 s for a small config, up to ~15 min for a large one under core
contention — so this one §10 estimate is roughly right), and that 130M-pair inference is "about 1 h CPU" (tree
scoring: **1.8–2.3 min**). What it does not mention is the only real cost: a 129M-row × ~150-float32 feature matrix is
**78 GB**, and several of those columns (F1 PSSM, F3 RWR, F4 cosine, F5 correlations, F6a/b/c/d) are *pair-level*
and must be produced for every one of the 129M rows, not looked up.

**Change.** Replace §10 with:
1. **Inference architecture (normative, not advisory): stream by kinase.** Loop over the ~430 kinases; for each,
   build the ≤300k-row × ~150-column block on the fly (180 MB), score, and write only the top-N per kinase plus
   the within-kinase percentile. Peak RAM stays under 2 GB; no 78 GB intermediate ever exists. All pair-feature
   kernels (PSSM scoring, co-regulation aggregation, RWR, cosines) are expressible as dense matmuls against a
   single kinase's descriptor column, so one GPU is enough and the four can run four kinases in parallel.
2. **Consequence to exploit:** because full-universe scoring costs minutes, (a) E-all can be evaluated
   exhaustively rather than sampled (ALG-1), (b) hard-negative mining over the *whole* universe is affordable
   (ALG-5), and (c) **U_all (1.87M S/T/Y × 430 = 800M pairs) is ~11 min of tree scoring**, so §2's "U_all scored
   with sequence/context features only" can instead be scored with the full model, with kinase-specific features
   set to their missing-value encoding. That directly serves PROMPT.md's Priority 2.
3. Measured/estimated budget for the whole project:

| Step | Device | Estimate |
|---|---|---|
| ESM-2 650M per-residue, full proteome, full-length context | 1 GPU | 8–15 min (fp16; cache 4.6 GB fp16 for S/T/Y) |
| ESM-2 650M kinase domains (442 × ~300 aa) | 1 GPU | < 10 s |
| Cantley PSSM scores, 396 kinases × 300k sites (as one-hot @ matrix matmul) | 1 GPU | < 1 min |
| F6a co-regulation, full 190k × 13.2k pairwise-complete, **per study block**, **per CV fold**, leave-protein-out | 1 GPU | ~4 s per pass; < 10 min for 5 folds × ~180 blocks |
| R2 LightGBM, one CV fold (330k–3.3M rows × ~150 feats) | 16 CPU | 10 s (63 leaves/600 trees) – 15 min (127 leaves/1000 trees, measured 888 s under contention); scales ~linearly with rows, so 1:200 ratios need early stopping / `max_bin` 63 |
| R2 full hyperparameter + α + λ sweep (5 folds × ~100 configs, 4 parallel 16-thread jobs) | 64 CPU | 2–12 h depending on tree size; tune tree size first on one fold |
| R3 frozen-encoder head, one fold, 10 epochs over ~800k pairs | 1 GPU | 3–8 min |
| R3 nested-OOF stack (25 trainings, ALG-6) × 5 seeds | 4 GPUs | 1–3 h wall |
| R3 LoRA-fine-tuned-650M variant (ablation only), one fold | 1 GPU | 1.5–2.5 h |
| B3 Phosformer MLM pretraining reproduction (CMP-2) | 4 GPUs | 1.5–2 h wall (≈6 GPU-h) |
| B3 Phosformer fine-tune, 5 folds × 5 seeds | 1–2 GPUs | 3–6 GPU-h |
| B4 SELPHI2.0 feature regeneration + 100-bag RF | 64 CPU | 1–2 days (feature pipeline, not model) |
| Final inference, U_obs × K ≈ 100–130M pairs (residue-class matched) | 64 CPU + 1 GPU | 15–40 min incl. feature construction |
| Final inference, U_all × K ≈ 800M pairs | 64 CPU + 1 GPU | 2–4 h |

4. **GPU allocation** (replacing §10's one-line version): GPU0 = R3 development (interactive); GPU1 = embeddings
   then R3 seeds; GPU2 = B3 Phosformer (pretrain then fine-tune, runs unattended for a day); GPU3 = feature
   kernels (co-regulation, PSSM, RWR) and full-universe inference. Keep 16 of 64 cores free for the CPU-side
   feature pipeline so LightGBM sweeps do not starve it.

**Why it matters for accuracy, not just bookkeeping:** §10 as written would lead agents to subsample the
co-regulation features and to approximate E-all, both of which cost accuracy, to save compute that does not exist.

---

### CMP-2 (major) — §8 B3: reproduce Phosformer's MLM pretraining rather than substituting ESM-2 init; define the fairness contract

**Problem.** §8 B3 says "we re-run their recipe and pretraining on kinase domains, **or** initialise from ESM-2,
documented". Those are not interchangeable: initialising a 6-layer/768 RoBERTa from ESM-2 benchmarks ESM-2's
pretraining corpus, not Phosformer. PROMPT.md requires a *fair* retrain of this benchmark, so the pretraining must
be reproduced, and the only reason given for not doing so would be cost.

**Cost estimate (so the decision is not made by default).** Phosformer pretrained on 295,320 kinase domains
(86.9 M aa) + 5,338 substrate proteins (3.85 M aa) ≈ 90.8 M tokens/epoch × 50 epochs = **4.5e9 tokens**, over a
6-layer × 768-dim encoder (≈42.5 M non-embedding params). At 6ND that is **1.16e18 FLOPs**; at the ~43 TFLOP/s
effective throughput I measured on this box for a short-sequence encoder, that is **≈ 6.4 GPU-hours**, i.e.
**1.5–2 h wall on 4×4090 with DDP**. This is less than one hyperparameter sweep. There is no cost argument for
skipping it.

**Change.** In §8, replace B3 with three named arms and a fairness contract:
- **B3a (headline)**: our clean-room reimplementation of the architecture, with their MLM pretraining reproduced
  on this hardware (~6 GPU-h, one-off, checkpoint cached and checksummed in DATA.md), fine-tuned on **our** labels,
  **our** folds, **our** negative generator.
- **B3b (ablation)**: identical architecture and fine-tuning, initialised from ESM-2 35M/650M instead — this
  isolates "does the kinase-domain-specific MLM corpus matter?", which is a real and currently open question.
- **B3c (contaminated reference)**: the released checkpoint, inference only, labelled as having seen our test
  labels. Report but never compare on equal footing.
- **Fairness contract to write into §8**: B3a/B3b and R3 must share the identical label set, splits, residue-class
  policy, negative generator and sampling ratio, seed list, and early-stopping criterion; the *only* differences
  permitted are architecture, inputs and pretraining corpus. Keep their reported augmentations (kinase-boundary
  ±5 shift, 5% kinase-token masking, under-represented-kinase upsampling) in B3a since they are part of the method,
  and test them separately in R3.
- Licence note already in `docs/lit/01` §6 (CC BY-NC-ND) stands: reimplement, do not redistribute modified weights.

---

### CMP-3 (minor) — §7 F0 / §10: fix ESM-2 650M as the default, test ESM-C 600M as the single alternative, drop 3B and 15B

**Problem.** §7 F0 says "ESM-2 650M … (later 3B or ESM-C 600M)" and §10 says "ESM-2 3B: optional" — an open
option with no decision rule, which in practice means an agent will spend a day on it.

**Evidence against scaling here.** (i) KolossuS used ESM-2 **15B** and reached AUROC 0.755 (S/T) on PSP labels,
and its own fine-tuning only moved 0.754 → 0.762 — the ceiling is the label/biology, not the encoder.
(ii) Raw ESM-2 zero-shot on this task is 0.503 (KolossuS), i.e. the PLM prior contributes ~nothing without
supervision, so a better prior is not the lever. (iii) Vieira et al. 2025 report small transfer gains 650M→3B→15B.
(iv) A practical argument the design should record: **3B gives 2,560-d embeddings, which cannot enter a GBDT over
129M pairs** (that is 1.3 TB of feature matrix) — they can only enter after PCA to ~64 dims or via R3, at which
point the extra width is discarded anyway. Compute is *not* the objection (3B is 40–70 GPU-min for the proteome);
expected value is.

**Change.** In §7 F0 and §10:
- Default: **ESM-2 650M**, frozen, full-length context.
- One pre-registered alternative: **ESM-C 600M** (1,152-d, MIT, same cost class, better per-parameter quality;
  note it needs `<plm-env>` transformers ≥5 for the native port — do not run the repo's `trust_remote_code`
  file).
- **Decision rule, written before running:** promote the alternative only if it improves DEV-CV macro AUROC
  (E-hard) by **> 0.005** with paired cluster-bootstrap p < 0.05 on the same folds and seeds.
- **Skip ESM-2 3B and 15B and ESM-3** — state the reason (ceiling evidence above + the dimensionality constraint),
  so the decision is documented rather than re-litigated.

---

### CMP-4 (minor) — §7 F6a / §10: co-regulation is free; upgrade the feature rather than approximating it

**Problem.** §7 F6a specifies "the correlation of *s*'s profile with the profile of each training substrate of *k*,
aggregated (mean of top-m and median Fisher-z)… computed per data block… combined with study-size-aware
weighting", and §10 hedges it as the heaviest step. The aggregation ("mean of top-m") is an approximation chosen
under an assumed compute constraint that does not exist.

**Change.** With the measured 3.6 s for a full 190k × 13.2k pairwise-complete correlation at min-n=8, specify in
§7 F6a that the correlations are computed **exactly, per study block, per CV fold, with full leave-protein-out**
(exclude every site on the substrate protein of *s* from *k*'s reference set, not just *s* itself — this is the
leak CoPheeKSA fixed only at site level and SELPHI2.0 did not fix at all; see `docs/lit/01` §1.4), and that the
per-kinase aggregation keeps a richer summary than top-m: `{mean, median, max, n_eff, 90th percentile, Fisher-z
weighted mean, and the same restricted to each of the 3–5 largest study blocks}`. Also state the implementation
(masked matmul: `n = M_a M_bᵀ`, `Σx = A M_bᵀ`, `Σxy = A Bᵀ`, …, 6 GEMMs, `r` masked where `n < 8`) so the agent
does not reach for a pandas/scipy loop, which would be ~10⁴× slower and is the realistic failure mode here.
Same kernel serves F5 (CPTAC co-abundance) and F6c/F6d unchanged.

---

## Which changes give the most accuracy per hour (my ranking)

1. **ALG-2** (paralogue exclusion from negatives) — hours of work, removes a systematic label error affecting the
   highest-n kinases, which are the ones the §6 primary actually averages over.
2. **CMP-4 + ALG-1** (exact leave-protein-out co-regulation; exhaustive E-all) — the literature's clearest
   complementary-signal modality (CoPheeKSA), now shown to cost seconds. Highest expected ΔAUPRC per hour.
3. **ALG-5** (listwise loss + mined hard negatives) — directly optimises two of the three declared endpoints;
   ~a day of work, minutes of GPU.
4. **ALG-4 descriptor dropout** — one line of code; without it the T2/T3 zero-shot numbers are not meaningful,
   and zero-shot coverage is half the brief.
5. **ALG-3** (λ-parameterised partial pooling) — converts §12.7 from a foregone conclusion into a result.
6. **CMP-1** (stream-by-kinase inference) — unlocks 2, 3 and the U_all coverage deliverable.
7. **ALG-7 / ALG-1 calibration** — does not change ranking accuracy at all, but the §6 coverage deliverable is
   uninterpretable without it.
8. **CMP-2** (Phosformer pretraining, 6 GPU-h) — buys comparability, not accuracy; run it unattended on GPU2.
9. **Not worth it**: ESM-2 3B/15B, AF-Multimer per-pair, KolossuS 15B inference. (B6 in §8 should be dropped to
   the 650M KolossuS variant only, which `docs/lit/01` §5.4 B6 already notes fits one 4090.)

---

## Risks I could not resolve

- **The 0.755 ceiling.** KolossuS's phylogeny-aware PSP number (0.755 S/T / 0.838 Y) may be close to the ceiling
  imposed by PSP label noise rather than by models. If so, most of the modelling choices above will move the
  metric by <0.02 and the honest differentiators will be *coverage* (understudied kinases, T2/T3) and the
  downstream benchmarKIN result. I cannot estimate the label-noise floor from here; a Johnson-style
  "report-count vs motif-match" stratification on our own labels would bound it and is cheap.
- **Whether R3 adds anything over R2 at all.** If F1 (Cantley) + F2a (learned PSSM) already saturate the sequence
  channel, R3's marginal value is only in T2/T3. The design cannot know this before running, but §12 should order
  the experiments so that R2-with-F1/F2a is measured before the R3 compute is committed — §12 step 3 before step
  5 already does this; keep it.
- **nnPU prior estimation.** π (the fraction of unlabelled sites that are true positives for *k*) is not
  identifiable from our data without an assumption; the Sugiyama/iKiP in vitro screens could bound it per kinase
  but they are themselves biased. I would run nnPU as an ablation with a swept π rather than a fitted one, and
  would not make it the default.
- **LightGBM lambdarank at this scale.** I did not benchmark `lambdarank` with kinase-sized groups
  (up to 300k items); it may be impractical and need group subsampling. Measure before committing §8 to it.
- **Full-length-context ESM throughput** is extrapolated from the 35M model, not measured on 650M (the 650M
  weights are not cached and I did not download them). The extrapolation is by non-embedding parameter count at
  constant MFU; attention cost at L≈550 adds perhaps 15%. I would treat the 8–15 GPU-min figure as ±2×, which does
  not change any conclusion.
- **Compendium density.** 95.4% of the compendium matrix is missing and 46% of KS sites are absent entirely
  (`data_profile` §3). The co-regulation features will be NA for half the training positives, and §7's
  "missingness is not given as a numeric count feature" rule is right but means the model sees a systematically
  different feature distribution for annotated vs unannotated sites. I flagged the mechanism; quantifying the
  residual bias is a data/statistics question, not one I can close here.
