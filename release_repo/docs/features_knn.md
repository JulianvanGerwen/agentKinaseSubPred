> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# knn_sim: similarity to a kinase's training substrates (label-derived, cross-fitted)

Owner: feat-knn. Family `knn_sim` = `ksp.features.knn.KnnSim`, registered in `ksp.crossfit.REGISTRY`.
DESIGN §7 cross-fitting rule / §16 A1: labels enter only through `fit(train_pairs, ctx)`.

## Files

| Path | What |
|---|---|
| `src/ksp/features/knn.py` | The family, the site-representation loader (`SiteReps`, `get_reps`), `columns()`, `load_tau()`, and the GPU engine (`KnnSim.dense`) |
| `scripts/features/knn_build_reps.py` | Builds the label-free representations (97 s, plus `--cat` 205 s) |
| `scripts/features/knn_tune_tau.py` | Chooses τ once, on CV fold 0 inner data → `results/features/knn/tau.json`, `tau_grid.tsv` |
| `scripts/features/knn_sanity.py` | Non-official fold-0 sanity check → `results/features/knn/sanity_fold0_*` |
| `scripts/features/knn_xfit.py` | Dense cross-fitted stores for R2 → `cache/features/knn/xfit/cv/fold{f}/knn_sim.npy` |
| `scripts/features/knn_test.py` | Unit tests: `python -m pytest -q scripts/features/knn_test.py` (9 tests, ~9 s). Kept in my paths; core-infra may move it to `tests/` |
| `cache/features/knn/reps_*.npy` | `reps_site_ids.npy` (row IDs; matched to the index **by ID**), `reps_esm_centre.npy`, `reps_esm_win7.npy`, `reps_esm_cat7.npy`, `reps_meta.json` |

## Representation spaces (all label-free)

The source is R3's ESM-2 650M per-residue embeddings in full-protein context (`cache/r3/win_emb.npy`, ±10 window, 347,354 U_ms sites; docs/model_r3.md §3). Non-U_ms sites have no representation, so their features are NaN.

| space | definition | similarity |
|---|---|---|
| `esmc` | Centre-residue embedding (1280-d). Centred on the mean of its residue type (S/T or Y) over U_ms, then L2-normalised | cosine |
| `esmcat` (default window space) | ±7 window, concatenated. Each position's embedding is projected on a global residue PCA-128 (fit on 200k random window residues; 63.9% variance), L2-normalised and scaled by √(w_p/Σw). Padding → 0. 1920-d | Σ_p (w_p/Σw)·cos_p: a centre-weighted mean of per-position cosines |
| `esmw` (optional) | ±7 window, mean-pooled over in-protein positions; centred and L2-normalised | cosine |
| `blosum` | One-hot 15-mer through a block-diagonal BLOSUM62 metric | Σ_p w_p·B62(q_p, t_p) / Σw. Padding and non-standard residues give 0 |

Weights are w_p = 1/(1+|p|/4) for p = −7..+7, with the centre included.

**Space (4), R3's site encoder, is skipped.** R3 saved only predictions (`cache/runs/r3_*/predictions.parquet`), not the encoder output.

**Window space choice.** Mean-pooled (`esmw`) and concatenated (`esmcat`) were both tuned on the same fold-0 inner data. `esmcat` won on lr macro AUROC (0.748 vs 0.726), so the default set is `knn.SPACES = ("esmc", "esmcat", "blosum")`. `esmw` stays available through `KnnSim(spaces=...)`.

## Features

For kinase k and query site s (residue type r ∈ {S/T, Y}), the set T is k's training sites of type r, **minus every training site on s's own substrate cluster** (from `ctx.clusters`, or the protein if it has no cluster).

This exclusion is done inside `transform`, so it holds even outside the driver. The driver additionally raises `LeakageError` if any such labels are fitted at all.

Columns per space `sp`:

| column | definition |
|---|---|
| `knn_<sp>_max` | max similarity to T |
| `knn_<sp>_top3` | mean of the top 3 (fewer if \|T\| < 3) |
| `knn_<sp>_kde` | log( Σ_t exp(sim/τ) / (\|T\| + c) ), with c = 1 (shrinks kinases with few labels) |
| `knn_<sp>_kde_bg` | log( mean_b exp(sim/τ) ) over a fixed label-free background: 2,000 random U_ms sites of type r (seed 0), minus s's cluster. Site-only |
| `knn_<sp>_lr` | kde − kde_bg: log density ratio of k's substrates against random sites of the same residue type |

There is one shared column, `knn_log_ntrain` = log1p(#training sites of k of type r), counted before exclusion. With the default spaces that makes 16 columns.

NaN means one of two things: k has no remaining training site of type r, or s has no representation. kde_bg is set to NaN wherever kde is NaN.

Canary kinases work like real kinases, because the family is keyed by `kinase_acc`.

## τ (chosen once; `results/features/knn/tau.json`)

**Protocol.**
- Outer fold-0 training labels (`crossfit.outer_split(0)`).
- Inner folds are `splits.inner_folds` over outer-training U_ms sites plus label sites, seed 0, as the driver assigns them. Inner fold 0 is held out.
- Fit on the other inner folds: 7,088 labels.
- Score 110 kinases with ≥5 held-out positives, each against 300 random held-out in-class sites with no DEV label of that kinase.
- Select τ by the max macro per-kinase AUROC of `knn_<sp>_lr`.
- The outer fold-0 eval proteins were not used.

| space | grid | τ | lr macro AUROC |
|---|---|---|---|
| esmc | 0.01–0.5 | **0.1** | 0.743 |
| esmcat | 0.01–0.5 | **0.1** | 0.748 |
| esmw | 0.01–0.5 | **0.1** | 0.726 |
| blosum | 0.2–8 | **4.0** | 0.774 (plateau from τ = 2 to 8) |

The full grid is in `tau_grid.tsv`.

## Leakage tests

**Unit test** (`scripts/features/knn_test.py`). Planted toy: homologous protein pairs share near-identical site embeddings and labels, and there is no signal across clusters.

| variant | AUROC (max / kde / lr) |
|---|---|
| leaky: fit on all labels, no exclusion | **1.00 / 1.00 / 1.00** |
| fit on all labels, in-transform exclusion | 0.48 / 0.47 / 0.49 |
| through `crossfit()` | **0.47 / 0.47 / 0.49** |
| positive control (planted shared direction), through the driver | 0.99 / 1.00 / 0.99 |

The test file also checks four other things:
- the driver raises `LeakageError` on overlapping rows;
- the engine matches a brute-force numpy computation on CPU and GPU, with multi-chunk paths and exclusions;
- `transform` matches `dense`;
- the BLOSUM metric matches a direct sum.

**Real data**, fold-0 outer-training rows (training positives plus 200 random negatives per kinase; macro AUROC):

| variant | esmc_lr | esmcat_lr | blosum_lr | esmc_max | blosum_kde |
|---|---|---|---|---|---|
| leaky (fit on all, no exclusion) | 0.999 | 1.000 | 0.961 | 1.000 | 0.929 |
| in-transform exclusion, no driver | 0.763 | 0.762 | 0.736 | 0.654 | 0.711 |
| driver | 0.759 | 0.759 | 0.730 | 0.653 | 0.706 |

## Sanity check (NON-OFFICIAL; CV fold 0; `knn_sanity.py`)

This is not `ksp.eval` and uses no paralog or Hamming masks.

**Setup.**
- Positives: DEV L2 labels on fold-0 proteins.
- Negatives: 500 random fold-0 in-class U_ms sites per kinase, with no DEV label of that kinase.
- Kinases: 94 with ≥5 positives (2,066 positives, 49k rows).
- knn and the learned PSSM (`f2a_pssm`) go through `crossfit()` with the same rows, seed and inner folds. Eval rows use the A1 inner-fit rule.

Per-kinase AUROC, macro:

| feature | macro | median | <20 train | 20–99 | ≥100 | Y class | frac. kinases > PSSM |
|---|---|---|---|---|---|---|---|
| knn_esmc_lr | 0.744 | 0.778 | 0.644 | 0.745 | 0.835 | 0.721 | 0.46 |
| knn_esmcat_lr | 0.737 | 0.772 | 0.625 | 0.740 | 0.836 | 0.719 | 0.46 |
| knn_esmw_lr | 0.712 | 0.753 | 0.631 | 0.716 | 0.781 | 0.737 | 0.36 |
| knn_blosum_lr | 0.750 | 0.757 | 0.602 | 0.767 | 0.850 | 0.685 | 0.54 |
| knn_blosum_kde | 0.718 | 0.724 | 0.570 | 0.736 | 0.813 | 0.657 | 0.28 |
| knn_esmc_kde | 0.681 | 0.703 | 0.635 | 0.686 | 0.714 | 0.677 | 0.21 |
| knn_esmc_max | 0.658 | 0.667 | 0.627 | 0.663 | 0.674 | 0.680 | 0.19 |
| knn_esmc_kde_bg (site-only) | 0.503 | 0.514 | 0.542 | 0.494 | 0.486 | 0.597 | 0.07 |
| **f2a_pssm** | **0.765** | 0.789 | 0.685 | 0.773 | 0.827 | 0.781 | — |
| f2a_pssm_ph | 0.774 | 0.798 | 0.709 | 0.777 | 0.832 | 0.797 | 0.56 |
| Cantley pct (label-free reference) | 0.785 | 0.803 | 0.704 | 0.799 | 0.833 | 0.719 | 0.56 |

All columns are in `sanity_fold0_summary.tsv`; per-kinase values are in `sanity_fold0_per_kinase.tsv`.

**Reading.**
- **lr ≫ kde/max.** Dividing by the background density matters: it removes generic "typical-site" similarity.
- **Comparison with the PSSM.** Alone, knn lr is slightly below the learned PSSM (0.74–0.75 vs 0.765), and below it for kinases with fewer than 20 labels. It is better for kinases with ≥100 labels (blosum_lr 0.850 and esmc_lr 0.835, vs 0.827).
- **Complementary to the PSSM.** The mean within-kinase Spearman correlation with `f2a_pssm` is 0.24 for esmc_lr, 0.28 for esmcat_lr and 0.57 for blosum_lr. A naive within-kinase rank average of esmc_lr and the PSSM gives 0.789, against 0.765 for the PSSM alone; esmcat gives 0.785.
- **The ESM spaces add the most new information.** blosum mostly re-learns the motif.
- **kde_bg is a site-popularity-type score.** It is ≈0.50 overall but 0.60 on Y kinases, so it should not be used alone (see provenance).

## Speed (RTX 4090, GPU 2, fp16 matmul)

- **Dense pass:** one fit plus `dense()` over **503 kinases × 347,354 U_ms sites** with 4 spaces takes **28.8 s**, including the NaN-filled output. The target was ≤ 2 min.
- **Fit** is a pure grouping and takes 0.02 s.
- **One-time setup:** loading and centring the representations takes about 55 s per process (NFS read). It is cached at module level, so later fits don't repeat it.
- **Dense xfit for one outer fold** (5 inner blocks of about 55k training sites plus eval sites, × 516 kinases) takes 32–39 s; fold 0 took 90 s because it includes the one-time load. That is about 5 s of transform per block, plus the 5.8 GB float16 write.
- **Built:** CV folds 0–4 (28 GB in total) on 2026-10-26.
- **Spot check of the fold-0 store:** 0% of outer-training positives have `knn_esmc_max` ≥ 0.999, i.e. there are no self-matches.
- **Engine:** an fp16 query × unique-training-site matmul per space, with own-cluster masking. Kinases are bucketed by training-set size (padding ≤ 1.25×), then gathered and reduced (max / topk / logsumexp in fp32).
- **CPU:** fp32 also works, which is what the unit tests use. It is much slower at full scale.

## Use in R2

**Option A: dense stores** (recommended, because R2 has no GPU). Run `CUDA_VISIBLE_DEVICES=2 python -I scripts/features/knn_xfit.py --folds 0 1 2 3 4 [final]`.
- It writes `cache/features/knn/xfit/cv/fold{f}/knn_sim.npy`.
- The format is float16 `(16, 516, 350736)`, plus `.block.npy` and `.json`. The layout is the same as `ksp.models.r2.xfit`:
  - kinase order = `canary.kinases_with_canaries()`;
  - site order = `cache/index/sites.parquet`;
  - eval rows use the A1 inner fit.
- `final` uses all DEV labels, with eval rows = TEST-protein U_ms sites. Whether that run should use `eval_mode` inner or all is the PI's call, so it was not run.

**Option B:** add `"knn": (lambda: knn.KnnSim, "site", knn.columns())` to R2's `FAMILIES`. `transform` accepts DenseExpand's pair chunks: it restricts the matmul to the requested kinases' training sites and memoises the background per site set. This needs a GPU to be fast.

## Assumptions and caveats

1. **Residue types are compared separately.** A S/T query is compared only with k's S/T training sites, and a Y query only with Y sites. An ST-class kinase's few Y labels form a separate slot, which R2 never scores because the classes don't match.
2. **The kde normalisation uses |T| after exclusion.** As a result, a positive whose own cluster holds other labels of k has a slightly smaller |T| and a slightly lower kde. This is the toy kde of 0.47 rather than 0.50. The bias works against the positive, so it is conservative.
3. **The background is label-free.** It is drawn from all U_ms sites, including other folds, and has no KS information, so this is not leakage.
4. **Only DEV labels were used.** τ and the window-space choice used fold-0 inner data only. The sanity numbers are fold 0 only, are not project metrics, and do not touch the official evaluation sets.
