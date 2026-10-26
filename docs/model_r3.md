# R3: neural kinase–site model (DESIGN §8.4, ALG-4/5, A6/A13/A18/A20)

Owner: R3-neural agent, 2026-10-26. Status: 5-fold CV and LKO-T2 done. Every metric on this page comes from `ksp.eval`.

## 1. Headline results (DEV CV, pooled out-of-fold predictions, `ksp.eval.prepare('cv')`)

Selected configuration: `reg1x5:w7_nce05`, meaning a ±7 window, InfoNCE weight 0.5 and a 5-seed ensemble.

| Endpoint | R3 pooled 5-fold CV | 95% CI (cluster bootstrap, B=1000) |
|---|---|---|
| **M1-CV = M5**: macro AUROC on E-hard over K_cv (203 kinases) | **0.729** | 0.710–0.747 |
| **M2**: filtered Hit@5, site-centric | **0.176** | 0.162–0.193 |
| **M3**: recall in the top 1% of E-all | **0.247** | 0.213–0.272 |
| **M4**: macro AUROC on E-hard over K_any (339 kinases) | **0.728** | 0.707–0.748 |

Secondary endpoints: AUROC on E-all 0.832, MRR 0.150, Hit@1 0.093, Hit@10 0.262. No (kinase, site) pairs were missing.

- Run: `results/runs/r3_w7_nce05_cv_pooled__20261026-151359__s0/`.
- Copies: `results/models/r3/r3_w7_nce05_cv_pooled_{summary.json,per_kinase.tsv}`.
- Strata: `results/runs/<that run>/strata.tsv`.

**Per fold** (M1 / M2 / M3):

| fold | M1 | M2 | M3 |
|---|---|---|---|
| 0 | 0.704 | 0.136 | 0.182 |
| 1 | 0.755 | 0.204 | 0.294 |
| 2 | 0.723 | 0.213 | 0.237 |
| 3 | 0.738 | 0.183 | 0.250 |
| 4 | 0.737 | 0.176 | 0.219 |

**LKO tier T2** (unseen kinase family, seen proteins; `results/splits/kinase_folds.tsv`). Each held-out kinase was scored by the model that never saw its fold's labels, and the results were pooled over the 5 kinase folds:

| Endpoint | Value |
|---|---|
| M1 (K_cv) | **0.664** |
| M3 | **0.203** |
| M4 (K_any) | **0.663** |

- Run: `results/runs/r3_w7_nce05_lko_pooled__20261026-151849__s0/`.
- Per LKO fold, M1 was 0.689, 0.687, 0.642, 0.682 and 0.613.
- **M2 is not reported for LKO.** Each LKO model scores only its held-out kinases, so a cross-kinase ranking would mix models (with `missing='drop'` the rank is among held-out kinases only).
- Losing about 0.065 of M1 against T1 shows the kinase descriptors (ESM domain, Cantley profile, group) carry real zero-shot signal.

**Context.** These are other agents' runs, read from `results/runs/*/metrics.json` and not reviewed by me:

| Model | Scope | M1 | M2 | M3 |
|---|---|---|---|---|
| B1 Cantley, universe percentile | pooled CV | 0.719 | 0.171 | 0.133 |
| B1 Cantley, imputed arm | pooled CV | 0.724 | — | — |
| B2 learned PSSM | pooled CV | 0.665 | — | — |
| R2 LightGBM | fold 0 only | 0.753 | 0.186 | 0.257 |

On the same fold 0, R3 scores 0.704 against R2's 0.753. So R3 alone is weaker than R2 on M1. Its main gain over Cantley is M3: 0.247 against 0.133. Under A13, R3 enters R2 by late fusion.

## 2. Model

### Inputs (no label-derived features)

**Site.** For each position of the ±7 window (±10 in one config), the encoder takes:
- the frozen **ESM-2 650M per-residue embedding from full-protein context**;
- an amino-acid embedding;
- a **phospho-channel** bit: 1 if that *flanking* residue is itself a U_ms site. The centre's own bit is forced to 0, as A6 requires.

It also takes one protein-context token, the ESM-2 protein mean (`esm2_650m/protein_mean` via `ksp.store`).

- **Window embeddings had to be recomputed.** `cache/esm2_650m` stores S/T/Y residues only, so the other window positions were missing.
  - `scripts/models/r3_window_embed.py` re-runs `ksp.features.esm_embed.embed_sequences` with identical settings (HF `facebook/esm2_t33_650M_UR50D`, last layer, fp16, 1022-aa windows with stride 511) on the 18,683 proteins that carry U_ms sites.
  - It keeps 3,656,284 residues within ±10 of the 347,354 U_ms sites.
  - Cost: 5.7 GPU-min (2.9 min wall on GPUs 2 and 3). Output: `cache/r3/win_emb.npy`, 9.4 GB fp16.
  - Site windows are indexed by an explicit `site_id` array (`win_sites.site_ids.npy`) and row matrix (`win_rows.npy`), and matched to `cache/index/sites.parquet` by ID.
  - Checks: centre-residue cosine against the struct+ESM cache is ≥ 0.99999 (mean 1.0000, n=5000), and window residues agree with `window31` with 0 mismatches.

**Kinase.** The 503 real kinases plus 13 canaries (`ksp.canary.kinases_with_canaries()` order, which is the `ksp.eval` kinase_idx order). There are four descriptor blocks:
1. **ESM-2 kinase-domain embedding** (`esm2_650m/kinase_domain`, z-scored, plus a missing flag; 5 kinases have no domain).
2. **Cantley block.** `cache/features/motif` contains no PSSM *vector* store, only kinase × site score matrices. I therefore used each kinase's unprimed log2 Cantley score profile over a fixed random reference of 3,000 S/T and 1,000 Y U_ms sites (seed 0), taken from `motif/cantley_unprimed_score`, or from `cantley_imputed_unprimed_score` for kinases without their own matrix. Each profile is z-scored per kinase, then reduced by PCA to 32 dimensions per residue type, with missing, imputed and distance flags. This is label-free. Coverage: 415 kinases have an S/T block (104 imputed) and 101 have a Y block (16 imputed).
3. **Manning group** embedding.
4. **Manning family** embedding.

Group and family embeddings are zeroed for groups and families with no training kinase in a run.
- **Descriptor dropout p = 0.2**, applied per kinase and per block independently at every step.
- Canaries use `canary_descriptors()` random vectors in the ESM block and have no Cantley block.

### Architecture
- **Site encoder:** a 2-layer pre-norm transformer (d=192, 4 heads); the output is [centre ; masked mean] passed to an MLP giving h=128.
- **Kinase encoder:** the sum of the block projections, then LayerNorm and an MLP giving h=128.
- **Head:** a bilinear term ⟨W z_s, u_k⟩/√h, plus MLP(z_s ⊙ u_k), plus a kinase bias b(u_k).
- Each step scores a batch of sites against **all** kinases (B × K logits).

### Loss
The total loss is the BCE term plus 0.5 × the InfoNCE term.

- **BCE term:** balanced masked BCE. Positive entries are weighted w_k ∝ n_k^-0.5 (n_k = the kinase's fit-row positives, normalised to a mean of 1 over positives); the negative term is a plain mean.
- **Targets** are built with `ksp.labels`:
  - positives are L2 labels (`load_dev('L2', with_canaries=True)`, U_ms only);
  - negatives are candidates of matching residue class (`kinases.residue_class`);
  - **masked** entries (ignored) are KP(k), identical or Hamming ≤ 2 15-mers, and Par(k)-labelled sites, all via `labels.NegativeMasks`.
- **No eval-fold label enters the training targets.** The masks are built from `load_dev('all')` restricted to the training proteins (and, for LKO, the training kinases).
- **InfoNCE term:** a filtered listwise InfoNCE over kinases for each positive (s, k).
  - The denominator is {k} ∪ {j : candidate negative at s, j ∉ Par(k)}.
  - Other true kinases of s and masked entries are excluded.
- **Kinases with no training positive** are not trained (all entries ignored), but they are still scored zero-shot from their descriptors.

### Sampling, early stopping and ensembling
- **Batches:** B = 1024 sites. 25% are sites carrying at least one training positive (these supply the E-hard-type negatives); the rest are uniform over training rows. AdamW at lr 3e-4 with weight decay 0.05; dropout 0.2 and ESM-input dropout 0.3.
- **Early stopping** (patience 10 × 100 steps) uses **only the inner validation split of the training folds**: `ksp.splits.train_valid_split_inner`, 20% of substrate clusters.
  - The metric is macro AUROC (hard negatives) over non-canary kinases with ≥ 3 inner-validation positives.
  - The CV evaluation fold is never used for stopping.
  - The best step is typically 200–500.
- **Ensembling:** 5 seeds per training run (0..4, the project seed list), with logits averaged. This raised inner-validation AUROC from about 0.76 to 0.779 on fold 0.

Code:
- `src/ksp/models/r3/{data,model,train}.py`
- `scripts/models/r3_cv.py`
- `r3_window_embed.py`, `r3_dev_sweep.sh`, `r3_queue.sh`
- diagnostics: `r3_canary_diag.py`, `r3_canary_strat.py`

## 3. Protocol and configurations tried (`results/cv/REGISTRY.tsv`)

| Stage | Data used | Configs | Result |
|---|---|---|---|
| Time pilot (untuned base: lr 1e-3, p_lab 0.5, 1 seed) | fold-0 eval **seen** | 2 (`pilot0:*`) | M1 0.687 / 0.692. The best step was 200, meaning heavy overfitting, so the base was re-tuned. |
| Base hyper-parameter dev sweep | inner validation of fold-0 training rows only; eval fold untouched | 10 (`dev:*`, split=inner_dev) | All within 0.750–0.764. `reg1` (0.762) was chosen as the base, with a 5-seed ensemble on top. |
| **Config selection on fold 0 (M1-CV)** | fold-0 eval | 3 | **w7_nce05 0.7035** > w10_nce05 0.7007 > w7_nce0 0.6907 → w7_nce05 |
| All 5 folds | — | selected config | see §1 |
| Diagnostic ablation, not eligible for selection | 5 folds | 1 (`w7_nce05_noctx`) | no protein token and no phospho channel: M1 0.719, M2 0.182, M3 0.223 |

In total, **5 configurations touched the fold-0 eval metric**, plus 10 inner-only dev runs and 1 post-hoc ablation. The ablation shows the protein and phospho context adds +0.010 M1 and +0.024 M3. The w10 config was better on M2 and M3 (0.145 and 0.196 on fold 0 against 0.136 and 0.182), but M1 is the pre-declared selection criterion.

## 4. Canaries (`load_dev(with_canaries=True)`; reported separately, excluded from all endpoints)

| Pooled CV | neg mean AUROC E-hard | neg mean AUROC E-all | pos AUROC E-hard | pos AUROC E-all | `canary_check` |
|---|---|---|---|---|---|
| R3 (selected) | **0.668** | **0.582** | 0.825 | 0.803 | neg **FAIL** (> 0.57), pos **FAIL** (< 0.9) |
| R3 no protein/phospho context | 0.624 | 0.521 | 0.845 | 0.864 | neg FAIL (E-hard) |

Per-canary results: `results/models/r3/r3_w7_nce05_cv_pooled_per_kinase.tsv`. Negative canaries have 7–56 DEV positives each.

**Diagnosis: the negative canary is confounded by protein phosphosite density.** The canary sampler draws proteins uniformly and then 1–4 sites per protein. Its sites are therefore strongly enriched on proteins with few U_ms sites, a label-free property that any model can see.

Label-free, kinase-agnostic scores evaluated with `ksp.eval` (`scripts/models/r3_canary_diag.py` → `results/models/r3/canary_diagnostic.tsv`):

| Score | neg E-hard | neg E-all |
|---|---|---|
| −(number of U_ms sites on the protein) | **0.764** | **0.737** |
| −(flanking U_ms count) | 0.629 | 0.574 |
| −(protein length) | 0.601 | 0.684 |

All of these "fail" the 0.57 threshold without using a single label.

Stratifying R3's canary AUROC by the protein's U_ms count (`scripts/models/r3_canary_strat.py`; positives-weighted within 16 density strata) brings the negative-canary mean down to:
- **E-all 0.477** and **E-hard 0.561** for the selected model;
- 0.462 and 0.560 for the no-context ablation.

Both are ≤ 0.57.

R3's training targets structurally exclude eval-fold labels:
- target rows for non-training proteins are −1;
- masks are built from training-protein labels only;
- early stopping uses inner validation only.

So I read the raw canary failure as a **canary-design artefact, not leakage**. The density-stratified check is my own diagnostic, built on `ksp.eval` primitives, and should be confirmed by core-infra or an independent reviewer. The suggested fix for core-infra is to sample canary sites with probability matched to the real-label protein-density distribution, or to accept canaries against a density-matched null.

**Positive canary (0.80–0.83).** The planted motif alone (an R-x-x-[ST]-[FILMV] indicator) reaches 0.963 E-hard and 0.988 E-all. R3 has no PPI channel, so the planted-PPI half is invisible to it. It also learns this 2-position motif for a random-descriptor kinase only partially before the global early stop, which is reached at about 300 steps and is driven by the real kinases. This is a real limitation: per-kinase specificity is under-fit for kinases identified only by descriptors.

## 5. Runtime (RTX 4090, one GPU per job, ≤ 4 CPU threads)
- **One-time window embedding:** 5.7 GPU-min.
- **One CV fold** (5 seeds): about 87 s of training, 1.2 s to score about 21.3M (kinase, site) pairs, and about 37 s of `ksp.eval`. That is **about 2.4 min wall per fold**, well under the 30-min target. Peak GPU memory is 11.4 GB, because the full window-embedding table is held on the GPU.
- **One LKO fold:** about 3.3 min (scoring all DEV sites).
- **Pooled evaluations:** CV with 1000 bootstraps took 290 s on CPU; LKO took 234 s.

## 6. Deviations and assumptions (stated, not silent)
1. **GPUs.** DESIGN §8.4 says GPUs 0–1; I used GPUs 2–3 as assigned.
2. **Cantley descriptor.** It is a score-profile PCA (see §2) rather than a flattened PSSM, because only score matrices exist in `cache/features/motif`.
3. **Extra inputs.** ALG-4 suggests a one-hot window plus a protein-context mean. I included both (an aa embedding and a protein-mean token) on top of the per-residue ESM window. Specificity-determining-residue one-hots were not used (optional in F7).
4. **Untested ALG-5 extras.** nnPU, label smoothing, hard-negative mining and the third loss term (L_kinase) were not implemented. Hard-type negatives come from the 25% labelled-site sampling instead.
5. **Git hash.** `ksp.manifest.git_hash` is overridden in `r3_cv.py` to `"not_queried(agent no-git rule)"`, so the agent never runs git. Manifests therefore carry no git hash.
6. **Fold-0 pilots.** The first two pilot runs (untuned base) looked at fold-0 eval metrics before the base hyper-parameters were re-tuned on inner validation. They are counted in REGISTRY.
7. **Kinase coverage.** Kinases without training positives in a fold (about 165–183 of 516) are scored zero-shot. In LKO, held-out kinases get no targets at all, positive or negative.
8. **Predictions.** Written by `RunManifest.save_predictions` to `cache/runs/<run_id>/predictions.parquet` with columns kinase_idx, site_idx, kinase_acc and site_id (categoricals) and score (raw logit). There is one row for every eval-fold U_ms site × every kinase of matching residue class, about 21.3M per fold and 106M pooled. The logit is comparable across kinases, so it also serves the cross-kinase site ranking (M2). No separate `score_site` is written, because a softmax over kinases at a site gives the same ranking.

## 7. Reproduce
```
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_window_embed.py run --worker 0 --n-workers 2   # + worker 1 on GPU 3
python -I scripts/models/r3_window_embed.py merge
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --fold F --config w7_nce05          # F = 0..4
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --lko J --config w7_nce05           # J = 0..4
python -I scripts/models/r3_cv.py pool --config w7_nce05 --runs <5 fold run_ids>
python -I scripts/models/r3_cv.py pool-lko --config w7_nce05 --runs <5 lko run_ids>
```

## 8. Suggested next steps
- **Late fusion with R2** (A13) on the OOF scores above.
- **Fix the canary sampler** (core-infra); see §4.
- **Under-fitting of descriptor-only kinases** (positive canary 0.83 against a motif ceiling of 0.96). Candidates:
  - per-kinase early stopping or longer training with a lower site-encoder learning rate;
  - a low-rank ESM projection to curb site memorisation;
  - a mined-negative / nnPU arm (ALG-5).

  Any of these must be logged as new configurations.
- **Final model** (DEV → TEST): same configuration with 5 seeds. Training on all of DEV takes about 2 min.
