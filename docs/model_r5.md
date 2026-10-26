# R5: fine-tuned ESM-2 late-interaction kinase–site model

Owner: R5 agent, 2026-10-26. Harness frozen at commit 974d520 (DESIGN A28); every metric on this page comes from
`ksp.eval`. Labels only through `ksp.labels.load_dev()`.

Paths: code `src/ksp/models/r5/{data,model,train}.py`, drivers `scripts/models/r5_{pilot,kinase_embed,cv}.py` and
`r5_queue.sh`, cache `cache/r5/`, results `results/models/r5/`, runs `results/runs/r5_*`.

## 1. What was asked, and what the compute allowed

The brief asked for a **cross-encoder** that jointly attends over the kinase domain and the substrate window, with a
fallback to **late interaction** if full cross-encoding could not score the ~21.3M candidate pairs of a CV fold
within about 1.5 GPU-h.

**Time pilot first** (`scripts/models/r5_pilot.py` → `results/models/r5/pilot_timing.json`, one RTX 4090, bf16,
inference only, joint sequence = 300 kinase-domain + 31 window + 3 special tokens = 334):

| Model | batch | pairs/s | GPU-h for 21.3M pairs (1 fold) |
|---|---|---|---|
| esm2_t12_35M | 64 | 1,637 | **3.6** |
| esm2_t12_35M | 256 | 1,383 | **4.3** |
| esm2_t30_150M | 64 | 447 | **13.2** |
| esm2_t30_150M | 256 | 414 | **14.3** |

A full cross-encoder is therefore **2.4×–9.5× over budget for inference alone**, before any training, and ×5 folds.

The two permitted fallbacks:
- **Two-stage (top-K rerank).** Rejected. `ksp.eval.evaluate` scores the *exhaustive* E-all universe and ranks
  unscored pairs last (`missing="worst"`); M3 (top-1% recall over E-all) and the E-all secondary metrics would be
  computed against a universe the model never saw, and the leaderboard already blanks M3 for runs with partial
  E-all coverage. The harness is frozen, so it cannot be adapted. A two-stage R5 would not be comparable to R2/R3.
- **Late interaction (ColBERT-style).** Adopted. It keeps a *token-level* kinase↔substrate interaction but
  factorises the two encoders, so each side is encoded once per fold and only the MaxSim head runs per pair.
  Measured cost of scoring a whole fold: **4.4 s** for 21.28M pairs (0.0012 GPU-h), i.e. three orders of magnitude
  inside budget.

## 2. Model

**Site tower (fine-tuned ESM-2).** ESM-2 (35M or 150M) over the ±15 window (`window31` of
`cache/index/sites.parquet`; 31 residues + CLS/EOS, terminus overhang padded and attention-masked). Three
**zero-initialised** additive inputs are injected at the embedding layer, so training starts from the unmodified
pretrained trunk:
1. the frozen ESM-2 650M per-residue embedding of the same residue **in full-protein context** (`cache/r3`, exactly
   R3's site input) — the fine-tuned trunk refines R3's representation rather than replacing it;
2. a phospho-channel embedding (a flanking residue is itself a U_ms phosphosite). **The centre site's own status is
   never used** (forced 0; DESIGN A6/DL2-7);
3. the protein-context vector (ESM-2 650M protein mean), broadcast over the window.

ESM-2's `token_dropout` scaling (×0.88) and its `position_ids` construction are reproduced explicitly, because the
`inputs_embeds` path in `transformers` skips both.

Outputs: per-token late-interaction vectors `S` (B, 33, 64) and a pooled site vector `z` (centre ‖ masked mean).

**Kinase tower.** Per-residue ESM-2 650M embeddings of the catalytic domain, cached once
(`scripts/models/r5_kinase_embed.py` → `cache/r5/kin_tok.npy`, 503 × 320 × 1280 fp16, built in 4 s). On top of them
a trainable Linear → 2-layer transformer (d=192) → `m=8` learned query tokens attend over the domain → `Q` (K, 8, 64).
The R3 descriptor blocks (kinase-domain mean embedding, Cantley score-profile PCA, Manning group/family, descriptor
dropout p=0.2) give a global kinase vector `u` (K, 128). Caching the kinase *tokens* is what makes the kinase tower
cost ~1% of a step instead of dominating it.

**Head.**

```
score(s,k) = gamma * sum_{j<8} max_t <Q[k,j], S[s,t]> / sqrt(64)     <- late interaction (new)
           + <W z_s, u_k>/sqrt(h) + MLP(z_s * u_k) + b(u_k)          <- R3's bilinear pair head
```

The second line makes R5's head a strict superset of R3's, so the late-interaction term is an addition that can be
ablated (`use_li` / `use_pair`).

**Losses (unchanged from R3, re-used from `ksp.models.r3.model`).** Balanced masked BCE (positives weighted
∝ n_k^-0.5) + 0.5 × filtered listwise InfoNCE over kinases per positive (denominator excludes the site's other true
kinases, masked entries and Par(k)) + optional Cantley-percentile distillation. Targets, masks (KP, Hamming ≤ 2,
Par(k)) and the training/validation protocol come from `ksp.models.r3.data.build_targets` and
`ksp.splits.train_valid_split_inner` verbatim, so R5 and R3 are directly comparable.

**Canaries.** Canaries have no kinase domain, so they get random-but-fixed pseudo-domain token blocks (16 tokens,
seeded, scaled to the real token s.d.) — the same role `ksp.canary.canary_descriptors` plays for R3. They carry no
sequence information.

**Efficiency note (identical metric, cheaper).** The inner-validation early-stopping metric only ever reads
label-carrying rows (positives, and negatives restricted to rows carrying some positive), so R5 scores only those
rows for it. The value is identical to R3's; the cost is ~50× lower.

## 3. Results — R5 does NOT clear the adoption gate; it is not carried to 5 folds

Pre-registered rule and config list: `results/models/r5/SELECTION.md` (written before any fold-0 R5 metric was
read). Gate from the PI brief: fold-0 **M1 ≥ 0.751** (= R3 fold-0 0.741 + 0.01) **or M2 ≥ 0.241** (= 0.221 + 0.02).

**Fold 0, `ksp.eval.prepare('cv', fold=0)`, 180 evaluated kinases, 21,283,695 pairs scored, 0 missing:**

| Config (3 seeds) | M1 | M2 Hit@5 | M3 top-1% | M4 | inner-val AUROC | canary neg / pos (E-hard) | gate |
|---|---|---|---|---|---|---|---|
| **R3 `reg1x5:w7_nce05_cpt_aux` (reference, 5 seeds)** | **0.741** | **0.221** | 0.267 | – | – | – | – |
| `r5li:li35_base` | 0.6915 | 0.1506 | 0.1451 | 0.6957 | 0.7695 | 0.503 / 0.895 | – |
| **`r5li:li35_cpt_aux`** (best R5) | **0.7250** | 0.1668 | 0.2106 | 0.7261 | 0.7853 | 0.519 / **0.958** | – |
| `r5li:li150_cpt_aux` | 0.7118 | 0.1749 | 0.1852 | 0.7155 | 0.7869 | 0.516 / 0.935 | – |
| *(time pilot: li35_base, 200 steps, 1 seed)* | 0.6861 | 0.1440 | 0.1483 | 0.6859 | 0.7528 | 0.458 / 0.863 | pilot |

**Best R5 = M1 0.7250, i.e. 0.016 BELOW R3 and 0.026 below the gate.** (The best *diagnostic* variant, `noli`,
reaches 0.7270 and is still 0.024 below the gate; it is not selectable in any case.) The best M2 over any R5 config is 0.1749
(150M), 0.066 below the 0.241 M2 gate. Neither arm of the gate is met, so per the pre-registered rule the 5-fold
run was not performed and nothing is carried forward. All runs are in `results/cv/REGISTRY.tsv` (`model=R5`).

**Diagnostic ablations** (fold 0, 3 seeds, *not eligible for selection*, exactly as R3's `noctx` was not):

| Ablation | M1 | M2 | M3 | M4 | inner-val |
|---|---|---|---|---|---|
| `li35_cpt_aux` (full head) | 0.7250 | 0.1668 | 0.2106 | 0.7261 | 0.7853 |
| `li35_cpt_aux_lionly` (late interaction only, R3 pair head removed) | 0.6778 | 0.1426 | 0.1885 | 0.6724 | 0.7561 |
| `li35_cpt_aux_noli` (late-interaction term removed) | **0.7270** | **0.1985** | **0.2509** | **0.7292** | 0.7992 |

### Why it falls short — three measured observations

1. **The fine-tuned trunk overfits almost immediately.** Early stopping selected step **100–600** of a 3,000-step
   budget in *every* seed of *every* config (`best_step` in each run's `summary.json`). There are only ~7.4k
   training positives on the fold-0 training proteins, against 35M (or 150M) trunk parameters. The
   inner-validation AUROC of the best R5 (0.785) is far above its fold-0 M1 (0.725): the model generalises to
   held-out *sites* of training proteins much better than to held-out *proteins*, which is what M1 measures.
2. **More trunk capacity makes it worse, not better.** The 150M trunk has the highest inner-validation score
   (0.787, the best of all configs) and the *second-worst* M1 (0.712). That ordering — inner-validation up, M1
   down — is the signature of the same overfitting, and it rules out "the trunk was too small" as the explanation.
3. **The late-interaction term contributes nothing — it is slightly harmful.** Removing it (`noli`) *improves*
   every endpoint: M1 0.7250 → 0.7270, M2 0.1668 → 0.1985, M3 0.2106 → 0.2509. Late interaction *alone*
   (`lionly`, R3's pair head removed) collapses to M1 0.678. So the model's whole score comes from the R3-style
   bilinear pair head, and the token-level kinase↔substrate matching that was the point of the brief adds no
   signal at this label volume. Even the best ablation (0.7270) is still 0.024 short of the gate, so this is not
   a near miss that a tuned late-interaction head would close.

The Cantley-distillation result from R3 **replicates cleanly** here: +0.034 M1 (0.6915 → 0.7250), +0.066 M3, and it
fixes the positive canary (0.895 → 0.958), exactly the behaviour reported in docs/model_r3.md §2.1. This is
independent confirmation of R3's PI-6/7 finding on a different architecture.

### Canaries

Canaries were injected in every run (`load_dev(with_canaries=True)`), and `ksp.eval.canary_check` was computed on
every fold-0 run. For the best R5 (`li35_cpt_aux`): negative-canary mean AUROC on E-hard **0.519 ≤ 0.57 (PASS)**,
positive canary **0.958 ≥ 0.90 (PASS)**; E-all values 0.720 / 0.987. The un-distilled base fails the positive
canary (0.895), as R3's un-distilled base did (0.830). **The gate is defined on pooled CV**, and no pooled R5 run
exists, so these are per-fold (fold-0) values and are reported as such.

## 4. Runtime (RTX 4090, one GPU per job, 4 threads per job, ≤ 3 jobs)

| Step | Time |
|---|---|
| Kinase-domain token cache, 503 domains, one-off | **4 s** |
| Substrate ±15 token/phospho cache, 347,354 sites, one-off | ~40 s (CPU) |
| Frozen-650M window store → GPU, per job | 47 s |
| Training, one CV fold, 3 seeds, 35M (`li35_base` / `li35_cpt_aux`) | 623 s / 779 s |
| Training, one CV fold, 3 seeds, 150M | 1,191 s |
| — of which Cantley distillation pretraining | ~66 s per seed |
| **Scoring all 21.28M candidate pairs of a fold (3-seed ensemble)** | **12.3 s (35M) / 36.9 s (150M)** |
| `ksp.eval` on one fold | ~36 s |
| **Total per fold, end to end** | **11 min (35M) / 21 min (150M)** |
| Peak GPU memory | 18.2–22.3 GB (one job per card) |
| Peak host RAM per job | ~13 GB (9.4 GiB frozen-650M window store + targets); 3 jobs ≈ 40 GB, inside the 60 GB agent budget |

The headline efficiency claim holds: **0.0034 GPU-h to score a whole fold** against the ~1.5 GPU-h budget and the
3.6–14.3 GPU-h a full cross-encoder would have needed — a factor of ~10³. A full 5-fold R5 run would have cost
about 70 GPU-min. The architecture is cheap; it is the accuracy that does not justify it.

Three runs are recorded as `status: failed` manifests, none of them hidden: one aborted before training because a
manifest input (the substrate token cache) did not exist yet, and two OOMed holding three seed models plus the
9.4 GiB frozen window store on one 24 GB card. The `noli` ablation was then rerun successfully with
`PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`. GPUs 1–3 only; GPU 0 was left to Phosformer throughout.

## 5. Deviations and assumptions

1. **Late interaction instead of a full cross-encoder**, for the measured reason in §1. The kinase and the substrate
   still interact at token level, but not through joint self-attention.
2. **Kinase-domain embeddings are computed on the isolated domain sequence** (the `emb_domain_iso` convention of
   `src/ksp/features/esm_kinases.py`), not in full-protein context.
3. **Long domains** (35 of 503, e.g. SRPK1/2, MASTL, CDC7) keep the N- and C-terminal 160 residues and drop the
   middle insert. **5 of 503 kinases have no annotated catalytic domain** and fall back to the first 320 residues;
   both are flagged in `cache/r5/kin_meta.json`.
4. **GPUs 1–3** as assigned (GPU 0 is Phosformer's), ≤ 12 cores (3 jobs × 4 threads).
5. **3 seeds per run**, not R3's 5, for time. This is conservative: it can only understate R5 against R3's 5-seed
   numbers.
6. **No git** (agent rule), so manifests carry the hash the manifest module queries.

## 6. Conclusion and what would be needed to revisit

**R5 is a negative result and is not carried forward.** The engineering goal was met — a token-level
kinase↔substrate interaction that scores a 21.3M-pair fold in 12 s — but the accuracy goal was not: the best
configuration is M1 0.725 against R3's 0.741 on the same fold, and the late-interaction term itself contributes
nothing (removing it gives 0.727). Nothing from R5 should enter the late fusion or the final model.

Two findings are worth keeping regardless:

1. **Full cross-encoding is not affordable in this harness.** Measured: 3.6–4.3 GPU-h/fold (ESM-2 35M) and
   13.2–14.3 GPU-h/fold (150M) for inference alone, against a ~1.5 GPU-h budget. Any future cross-encoder work
   needs either a two-stage design *and* a harness change to evaluate partial coverage correctly (the frozen
   `ksp.eval` ranks unscored pairs last), or a factorised encoder as used here.
2. **Cantley distillation replicates on a second, independent architecture**: +0.034 M1, +0.066 M3, and it fixes
   the positive canary (0.895 → 0.958). R3 saw +0.019 M1 and the same canary repair. This strengthens the PI-6/7
   finding well beyond a single model.

If this direction is revisited with more time, the measured failure mode points at the training regime, not the
architecture: every seed early-stopped by step 100–600 of 3,000, and the larger trunk raised inner-validation while
lowering M1. The things to try first would be freezing most of the trunk (LoRA or last-k layers only), a much lower
trunk learning rate with a longer schedule, and an inner-validation split grouped by protein rather than by
substrate cluster so that early stopping measures the same generalisation M1 does. None of that was run here, and
none of it is claimed to work.

## 7. Reproduce

```
CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_pilot.py            # cross-encoder vs late-interaction timing
CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_kinase_embed.py     # cache/r5/kin_tok.npy (4 s)
CUDA_VISIBLE_DEVICES=1 python -I scripts/models/r5_cv.py train --fold 0 --config li35_cpt_aux --n-seeds 3
python -I scripts/models/r5_report.py folds                            # the fold-0 table above
# 5-fold run (NOT performed: gate not met) would be:
bash scripts/models/r5_queue.sh li35_cpt_aux "1 2 3 4"
python -I scripts/models/r5_cv.py pool --config li35_cpt_aux --runs <5 fold run_ids>
```
