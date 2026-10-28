# R3: neural kinase–site model (DESIGN §8.4, ALG-4/5, A6/A13/A18/A20)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Owner: R3-neural agent, 2026-10-26.

Status: **round 2 is final**, using canary v2. The best variant is **`reg1x5:w7_nce05_cpt_aux`**. Its out-of-fold predictions (5 CV folds plus T2) for late fusion with R2 are listed in `results/models/r3/FINAL_OOF.json`.

Every metric on this page comes from `ksp.eval`.

## 1. Final results: canary v2, pooled out-of-fold over DEV (`ksp.eval.prepare('cv')`)

| Endpoint | **R3 final (cpt_aux)** [95% CI, cluster bootstrap B=1000] | R3 base (w7_nce05) | Paired Δ final − base [95% CI], p |
|---|---|---|---|
| **M1-CV = M5** (K_cv, 203 kinases) | **0.740** [0.725, 0.757] | 0.721 | **+0.019** [+0.010, +0.033], p=0.002 |
| **M2** Hit@5 (site-centric, filtered) | **0.223** [0.211, 0.238] | 0.183 | **+0.040** [+0.025, +0.058], p=0.002 |
| **M3** top-1% recall | **0.282** [0.253, 0.310] | 0.233 | **+0.049** [+0.028, +0.078], p=0.002 |
| **M4** (K_any, 339 kinases) | **0.733** [0.719, 0.753] | 0.722 | +0.011 [+0.002, +0.025], p=0.014 |
| **LKO T2: M1** (K_cv, held-out kinase families) | **0.686** | 0.650 | +0.036 (not bootstrapped) |
| **LKO T2: M3** | **0.226** | 0.203 | +0.023 |
| **LKO T2: M4** | **0.676** | 0.655 | +0.021 |

The paired test is `ksp.eval.paired_bootstrap` (`results/models/r3/paired_cptaux_vs_base_v2.tsv`).

**M5 by training-label bin** (K_cv; macro AUROC on E-hard):

| Bin (n kinases) | CV base | CV final | T2 base | T2 final |
|---|---|---|---|---|
| 5–19 (83) | 0.713 | **0.727** | 0.651 | **0.698** |
| 20–99 (89) | 0.716 | **0.736** | 0.654 | **0.674** |
| ≥100 (31) | 0.756 | **0.784** | 0.637 | **0.688** |

In the K_any bin with 1–4 labels (79 kinases), the base scores 0.704 and the final 0.696, so there is no gain there.

The PI's hypothesis is **supported for T2**: the largest gain is +0.047 in the 5–19 bin. It is **only partly supported for CV**: the 5–19 gain (+0.014) is smaller than the ≥100 gain (+0.028), and kinases with 1–4 labels do not improve.

**Per fold, final variant** (M1 / M2 / M3):

| fold | M1 | M2 | M3 |
|---|---|---|---|
| 0 | 0.741 | 0.221 | 0.267 |
| 1 | 0.756 | 0.260 | 0.280 |
| 2 | 0.750 | 0.228 | 0.275 |
| 3 | 0.748 | 0.243 | 0.303 |
| 4 | 0.752 | 0.236 | 0.300 |

**Canary v2 gate** (`eval.canary_check`; the gate is the negative-canary mean on E-hard ≤ 0.57, and the positive canary ≥ 0.90 on both E-hard and E-all):

| Pooled CV | neg E-hard (gated) | neg E-all (not gated) | pos E-hard | pos E-all | Result |
|---|---|---|---|---|---|
| **Final (cpt_aux)** | **0.501** | 0.687 | **0.929** | **0.974** | **neg PASS, pos PASS** |
| Base w7_nce05 | 0.486 | 0.647 | 0.830 | 0.836 | neg PASS, pos FAIL |

Per fold, the final variant's negative E-hard mean was 0.498, 0.522, 0.456, 0.530 and 0.667. Fold 4 is high, but each fold holds only about 2–11 positives per canary; the pooled value is the gated one.

Distillation **fixes the positive canary**, which went from 0.83 to 0.93/0.97. Pretraining teaches the site encoder motif structure, so the kinase-specific fine-tuning learns the planted R-x-x-S/T-Φ motif even though the kinase is identified only by a random descriptor.

**Pooled runs:**
- final CV: `results/runs/r3_w7_nce05_cpt_aux_cv_pooled__20261026-162205__s0`
- final LKO: `r3_w7_nce05_cpt_aux_lko_pooled__20261026-164851__s0`
- base CV under canary v2: `r3_w7_nce05_cv_pooled__20261026-160122__s0`
- base LKO, round 2: `r3_w7_nce05_round2_lko_pooled__20261026-164850__s0`

Copies of the summaries and per-kinase tables are in `results/models/r3/*_{summary.json,per_kinase.tsv}`.

## 2. Round-2 experiments (fold-0 selection by M1-CV; 4 new eval-fold configs, all logged in `results/cv/REGISTRY.tsv`)

| Config (fold 0, canary v2) | M1 | M2 | M3 | Decision |
|---|---|---|---|---|
| base `w7_nce05` (re-run with canary v2) | 0.695 | 0.146 | 0.180 | reference |
| `w7_nce05_cpt`: Cantley pretraining | 0.739 | 0.217 | 0.238 | — |
| **`w7_nce05_cpt_aux`**: pretraining + auxiliary distillation 0.2 during fine-tuning | **0.741** | 0.221 | 0.267 | **selected** (highest M1) |
| `w7_nce05_esm3b`: ESM-2 3B windows | 0.700 | 0.135 | 0.203 | not adopted (+0.005, below the +0.01 rule) |
| `w7_nce05_cptaux_esm3b`: 3B + cpt_aux | 0.725 | 0.206 | 0.258 | not adopted (−0.015 against cpt_aux) |

- **Base under canary v2.** The base fold-0 M1 moved from 0.7035 (round 1) to 0.695. The label tables changed with canary v2: `dev_pairs.parquet` was re-written at 15:34.
- **Selection margin.** The cpt versus cpt_aux margin on M1 is within noise (0.0017). cpt_aux was chosen by the pre-declared M1 rule; it is also better on M3.
- **Two crashed runs.** Two 3B attempts crashed: a protein-mean dimension bug, and a job I killed. They are recorded as `status: failed` manifests: `r3_reg1x5_w7_nce05_esm3b_fold0__20261026-155930__s0` and `..._cptaux_esm3b_fold0__20261026-160048__s0`. The bug was fixed and both configs were re-run.

**Configuration count for R3 on eval folds (REGISTRY `split=cv`):**
- round 1: 5, including 2 untuned pilots;
- round 1 diagnostic ablation: 1 (`noctx`, not eligible for selection);
- round 2: 4;
- total distinct configs: 10.

There were also 10 inner-validation-only dev runs (`split=inner_dev`).

### 2.1 Cantley-atlas distillation (PI-6/7)
- **Targets:** label-free Cantley **universe percentiles divided by 100**, two channels:
  - `motif/cantley_unprimed_pct_universe`;
  - `motif/cantley_primed_soft_pct_universe`.
- **Coverage:** own-matrix kinases only (396 kinases; imputed rows are not used). Targets are NaN out of residue class and for canaries.
- **Sites:** **training-fold U_ms sites only.** In CV these are the DEV sites of the 4 training folds. In LKO they are all DEV sites; see §5 item 3.
- **Distillation head:** a separate head on the shared site and kinase encoders, with a bilinear term, a pair MLP and a kinase bias, giving 2 outputs. Its loss is soft-label BCE on every (site, own-matrix kinase) entry of a batch, about 400k targets per step.
- **Stage 1 (pretraining):** 3000 steps, lr 1e-3 with linear decay, B=1024 uniformly sampled training sites, descriptor dropout active. This takes about 35 s per seed on a free GPU.
- **Stage 2 (fine-tuning):** exactly the round-1 recipe: masked BCE plus 0.5 × filtered InfoNCE, early stopping on inner validation. In `cpt_aux`, the distillation loss is kept with weight 0.2.
- **Predictions** always come from the main head. The distillation head is discarded.
- Code: `train.cantley_targets`, `train.distill_loss`, `model.R3Net.distill`.

### 2.2 ESM-2 3B windows
- Built by `scripts/models/r3_window_embed.py --variant 3b`, ±7 window, giving 2,962,692 residues × 2560 (15.2 GB fp16) in `cache/r3/esm3b/`.
- Took 10.8 min wall on GPUs 2 and 3, about 20 GPU-min, running concurrently with training.
- Centre-residue cosine against `cache/esm2_3b` is ≥ 0.999995, and residues mismatch 0 times.
- Only the window embeddings change. The protein mean and kinase-domain embeddings stay 650M.
- Peak GPU memory is 17.2 GB. The 3B variant was not adopted.

## 3. Model (unchanged from round 1 except §2.1)

### Site
- For each position of the ±7 window, the inputs are the frozen ESM-2 650M per-residue embedding **from full-protein context**, an amino-acid embedding, and a **phospho-channel** bit (a flanking residue is a U_ms site; the centre bit is forced to 0 per A6).
- There is one protein-context token (the ESM-2 protein mean).
- The encoder is a 2-layer pre-norm transformer with d=192 and h=128.
- **Window embeddings were recomputed.** `cache/esm2_650m` stores S/T/Y residues only, so `scripts/models/r3_window_embed.py` re-ran `ksp.features.esm_embed` with identical settings. It kept 3,656,284 residues within ±10 of the 347,354 U_ms sites (5.7 GPU-min, 9.4 GB fp16, `cache/r3/`). Centre cosine against the cache is ≥ 0.99999.
- Window data are indexed by an explicit `site_id` array and matched to `cache/index/sites.parquet` by ID.

### Kinase
The 503 real kinases and 13 canaries, in `kinases_with_canaries()` order. There are four blocks, each with descriptor dropout p=0.2:
1. **ESM-2 kinase-domain embedding** (`esm2_650m/kinase_domain`, z-scored, with a missing flag).
2. **Cantley score profile.** No PSSM vector store exists, so this is each kinase's unprimed log2 score over 3,000 S/T and 1,000 Y reference U_ms sites, then PCA to 32 dimensions per residue type, with missing, imputed and distance flags.
3. **Manning group.**
4. **Manning family.** Groups and families with no training kinase in a run are zeroed.

Canaries take random descriptors from `canary_descriptors()`.

### Head and loss
- **Head:** a bilinear term plus a pair MLP plus a kinase bias. Each step scores B sites against all K kinases.
- **BCE:** balanced masked BCE, with positives weighted ∝ n_k^-0.5.
- **InfoNCE:** filtered listwise InfoNCE over kinases. The denominator excludes other true kinases, masked entries and Par(k).
- **Targets:** positives are L2 labels on U_ms sites. Candidates are of matching residue class. Masks are KP, Hamming ≤ 2 and Par(k), via `labels.NegativeMasks` built from **training-protein labels only**.
- **Untrained kinases:** kinases with no training positive are not trained, but are scored zero-shot.

### Training
- B=1024 per step, 25% of them sites carrying a training positive.
- AdamW at lr 3e-4 with weight decay 0.05; dropout 0.2 and ESM-input dropout 0.3.
- Early stopping (patience 10 × 100 steps) uses **only** the inner validation split of the training folds (`splits.train_valid_split_inner`).
- 5-seed ensemble (seeds 0–4), with logits averaged.

### Predictions
- Stored with `RunManifest.save_predictions` at `cache/runs/<run_id>/predictions.parquet`, columns kinase_idx, site_idx, kinase_acc, site_id and score (raw ensemble logit).
- **CV:** every U_ms site of the fold × every kinase of matching residue class, about 21.3M rows per fold.
- **LKO:** held-out kinases × all DEV U_ms sites.
- The logit is comparable across kinases, so it also serves M2.

## 4. Runtime (RTX 4090, one GPU per job, ≤ 4 threads per job, 2 jobs at a time)

| Step | Time |
|---|---|
| Final variant, one CV fold, 5 seeds | 5.5–7.5 min wall |
| — of which pretraining (alone on the GPU) | about 35 s per seed |
| — of which fine-tuning | about 20 s per seed |
| — of which `ksp.eval` | about 40 s |
| Final variant, one LKO fold | about 10 min |
| Base variant, one CV fold | about 2.4 min |
| Pooled CV evaluation (B=1000) | about 5 min |
| Paired bootstrap | about 10 min |

- Peak GPU memory is 12.1 GB (650M) and 17.8 GB (3B).
- All round-2 runs took about 75 min wall on GPUs 2 and 3.

## 5. Deviations and assumptions
1. **GPUs.** I used GPUs 2–3 as assigned; DESIGN §8.4 says GPUs 0–1.
2. **Cantley descriptor** is a score-profile PCA (see §3).
3. **Distillation in LKO.** Held-out kinases' *label-free* Cantley percentiles over DEV sites are part of pretraining. This is the same information as their Cantley descriptor and is available at deployment for any atlas kinase, but T2 for the final variant is therefore "unseen labels" rather than "no kinase-specific signal of any kind". Overall, 396 of 503 kinases (79%) have their own atlas matrix.
4. **Git hash.** Per PI instruction, the round-1 override of `ksp.manifest.git_hash` was removed. Round-2 manifests carry the git hash (`581e25c…`). Round-1 manifests show `not_queried(agent no-git rule)`.
5. **Untested ALG-5 extras.** nnPU, hard-negative mining and the third loss term (L_kinase) are not implemented.
6. **Canary v1** (round 1, now superseded): it was confounded by protein phosphosite density. Evidence is in `results/models/r3/canary_diagnostic.tsv` and `canary_stratified_*.tsv`, and `scripts/models/r3_canary_{diag,strat}.py`. Core-infra replaced it with v2.

## 6. Round 1 (history, canary v1)
- **Fold-0 selection among w7_nce05, w10_nce05 and w7_nce0** gave M1 of 0.7035, 0.7007 and 0.6907; w7_nce05 was chosen.
- **Pooled base, canary v1:** M1 0.729, M2 0.176, M3 0.247. LKO T2 M1 was 0.664.
- **No-context ablation** (no protein token, no phospho channel): M1 0.719.
- **Summaries:** `results/models/r3/r3_w7_nce05_canaryv1_cv_pooled_summary.json` and `r3_w7_nce05_round1_lko_pooled_summary.json`.

## 7. Reproduce
```
# embeddings (once): 650M (+-10) and, optionally, 3B (+-7)
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_window_embed.py run [--variant 3b] --worker 0 --n-workers 2  # + worker 1 on GPU 3
python -I scripts/models/r3_window_embed.py merge [--variant 3b]
# final variant
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --fold F --config w7_nce05_cpt_aux   # F = 0..4
CUDA_VISIBLE_DEVICES=2 python -I scripts/models/r3_cv.py train --lko J --config w7_nce05_cpt_aux    # J = 0..4
python -I scripts/models/r3_cv.py pool     --config w7_nce05_cpt_aux --runs <5 fold run_ids>
python -I scripts/models/r3_cv.py pool-lko --config w7_nce05_cpt_aux --runs <5 lko run_ids>
python -I scripts/models/r3_paired.py --a <final 5 runs> --b <base 5 runs> --name cptaux_vs_base_v2
```

Code: `src/ksp/models/r3/{data,model,train}.py`, `scripts/models/r3_*.py`, `r3_queue.sh`, `r3_dev_sweep.sh`.

## 8. Next steps
- **Late fusion with R2** (A13) on `FINAL_OOF.json`.
- **Final model:** the same configuration trained on all of DEV, with pretraining on DEV sites only. About 7 min on one GPU.
- **Possible gains:**
  - longer or stronger distillation, e.g. adding the imputed-matrix kinases as soft targets;
  - per-kinase early stopping;
  - an nnPU arm.

  Any of these must be logged as new configurations.

## 9. G4 TEST-batch deliverables (2026-10-26, 21:25–22:05; docs/G4_RUNBOOK.md)

**Model.** The frozen final configuration `reg1x5:w7_nce05_cpt_aux` (PI name "r3_final_distill"):
- trained on **all of DEV**, with Cantley pretraining on DEV U_ms sites only;
- seeds 0–4, logits averaged;
- early stopping on the inner validation split of DEV (seed 0).

Only TEST-protein **sites** are used for scoring (`cv_fold == -1`); no TEST label was read. Mode: `train.run(test=True)`.

**Assembled run: `r3_final_distill_test__20261026-220424__s0`.** The manifest records git hash `13c47c2…`, `train_wall_s`, `inference_wall_s`, the component run ids and the pre-flight results.

| File | Content | Rows | `--check-predictions` | SHA-256 |
|---|---|---|---|---|
| `cache/runs/r3_final_distill_test__20261026-220424__s0/test_predictions.parquet` | final model, all 503 kinases × TEST candidates | 25,922,452 | ok: true (503 kinases full, 0 partial) | `6223c72c3a93d8082c957753c541dc0737a8e9b5921e90ce19e29af71828d384` |
| `cache/runs/r3_final_distill_test__20261026-220424__s0/test_t3_predictions.parquet` | T3: each kinase scored by the model trained without its LKO family fold | 25,922,452 | ok: true (503 kinases full) | `9c6c403ea27690d8f5c337b47525a4cb84452ad0cc9be1cf31de78a4b78b8874` |

**Columns** are kinase_idx, site_idx and score (the raw 5-seed mean logit; higher means more likely). There is **no `score_site`**: R3's logit is used for the cross-kinase ranking, and the fused `score_site` is produced by `scripts/models/r2_fuse_test.py`.

My candidate set matches `ksp.eval_test.test_candidate_pairs()` exactly (25,922,452 pairs, 0 left-only, 0 right-only).

**Component runs** (each in `cache/runs/<id>/component_predictions.parquet`):

| Component | Run id | Kinases | train_wall_s | inference_wall_s |
|---|---|---|---|---|
| final | `r3_final_distill_test_final__20261026-212607__s0` | 503 | 718.7 | 3.1 |
| LKO fold 0 | `r3_final_distill_test_t3lko0__20261026-213814__s0` | 102 | 719.7 | 2.4 |
| LKO fold 1 | `r3_final_distill_test_t3lko1__20261026-212607__s0` | 100 | 793.1 | 2.3 |
| LKO fold 2 | `r3_final_distill_test_t3lko2__20261026-215020__s0` | 97 | 808.5 | 2.3 |
| LKO fold 3 | `r3_final_distill_test_t3lko3__20261026-213926__s0` | 103 | 731.7 | 1.1 |
| LKO fold 4 | `r3_final_distill_test_t3lko4__20261026-215144__s0` | 101 | 393.2 | 1.3 |

The five T3 trainings sum to 3,446 s of training and 9.4 s of inference.

**Runtimes** are wall-clock times on GPUs 2 and 3, which were shared with other agents' jobs: up to 20.7 GB of 24.5 GB in use, at 100% utilisation. They are inflated about 2× against an idle GPU (≈ 6 min per training). `train_wall_s` includes data loading, target building, the GPU upload, pretraining and fine-tuning for all 5 seeds.

**FREEZE entry draft:** `results/models/r3/FREEZE_entry_r3_final_distill.json`. Its role is "secondary"; the PI sets the role. Code: `scripts/models/r3_cv.py test|assemble-test`, `scripts/models/r3_test_queue.sh`.

## 10. G5 deployment (prepared and smoke-tested 2026-10-27 00:15; the real run waits for the PI's "GO deploy")

**Command** (one GPU, about 15 GB RAM):

    CUDA_VISIBLE_DEVICES=<g> OMP_NUM_THREADS=4 .venv/bin/python -I scripts/models/r3_cv.py deploy

**What it does.**
- **Config:** the frozen `r3_final_distill` config (`reg1x5:w7_nce05_cpt_aux`, seeds 0–4). Overrides are refused.
- **Labels:** `ksp.deploy.load_deploy_labels("L2")` for positives and `("all")` for the KP/Par/Hamming masks. Canaries are off, the module default. `deploy.check_test_completed()` is called before any GPU work, and the command was verified to refuse while the ledger has no COMPLETED block.
- **Training rows:** DEV+TEST **U_ms** sites (A6). Distillation pretraining uses the same rows. Early stopping uses an inner-validation 20% of substrate clusters.
- **Scoring:** every **U_obs** site × every real kinase of matching residue class, 131,144,410 pairs.

**U_obs window embeddings.** The U_ms cache lacks the 3,382 non-U_ms U_obs sites, so the windows were rebuilt with identical ESM-2 650M settings in `cache/r3/uobs/`:
- 350,736 sites and 3,691,633 residues;
- 5.8 GPU-min;
- centre cosine against the S/T/Y cache ≥ 0.99999, and 0 residue mismatches.

The phospho channel still marks flanking **U_ms** sites only.

**Output:** `cache/release/r3_deploy_<ts>/uobs_predictions.parquet`.
- Columns: kinase_idx int32, site_idx int32, score float32 (raw 5-seed mean logit).
- Written in chunks of 20k sites (zstd), ordered by site_idx and then kinase_idx.
- A `README.json` sits next to it with the SHA-256, `train_wall_s` and `inference_wall_s`. The manifest goes to `results/runs/r3_deploy__<ts>__s0/`.

**Smoke test** (`--smoke`: DEV labels and DEV rows, 1 seed, 200+200 steps; output in the scratchpad; run `r3_deploy_smoke__20261027-001319__s0`):
- 131,144,410 rows over 503 kinases;
- unique keys, all finite, residue class correct, and the count equals the expected U_obs × class total;
- train 32 s, inference (score plus write) 6.5 s;
- peak RSS 14.0 GB, peak GPU memory 12.6 GB.

**Real deployment run** (GO from the PI at 00:58; TEST batch COMPLETED for FREEZE 721dc197):

| Item | Value |
|---|---|
| Run | `r3_deploy__20261027-005817__s0`, git hash `37a3e6a…` |
| Output | `cache/release/r3_deploy_20261027-005816/uobs_predictions.parquet` |
| SHA-256 | `2f62988cf798f15ace6cd9925dce77e45ca2c6cf297549203734d256e9d002ba` |
| Rows | 131,144,410: 350,736 U_obs sites × 503 kinases |
| Checks | unique keys, all finite, residue class correct |
| Training labels | DEV+TEST L2, 14,747 positives on 347,354 U_ms training rows |
| Inner-validation ensemble AUROC (hard negatives) | 0.795 |
| `train_wall_s` | 446.6 |
| `inference_wall_s` | 11.8 (score 4.7 + write 7.1) |
| Peak memory | 13.6 GB RSS, 13.3 GB GPU |

The U_obs input-cache rebuild is recorded in the manifest under `input_cache_note`.
