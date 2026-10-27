> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Motif features: F1 (Cantley atlas, label-free) and F2a (learned PSSM, label-derived)

Owner: feat-motif. Validation: `results/features/motif/validation.md`.

## Files

| Path | What |
|---|---|
| `src/ksp/features/motif_cantley.py` | Clean-room Cantley scoring engine. It reads the kinase-library 1.8.0 **data files only** and never imports or executes package code. It also does nearest-atlas imputation. |
| `src/ksp/features/motif_learned.py` | F2a learned PSSM: `fit_pssm`, `score`, `score_pairs`, `universe_percentile`, `top_positions`, and the `LearnedPSSM` class implementing the `ksp.crossfit` contract (registered name `f2a_pssm`). |
| `scripts/features/motif_cantley_build.py` | Builds all F1 arrays (151 s). |
| `scripts/features/motif_cantley_validate.py` | Reproduction check against the published outputs. |
| `scripts/features/motif_extract_supp.py` | Streams the strict-OOXML supplementary xlsx files into parquet. |
| `scripts/features/motif_validate.py` | Non-official sanity checks: AUROC, site-centric rank, learned-PSSM logos, crossfit smoke test. |
| `cache/features/motif/cantley_*.npy` + `.json` | F1 arrays saved via `ksp.store.save_matrix`. Load with `ksp.store.load_matrix("motif/<stem>")`. `cantley_meta.json` describes them all. |
| `results/features/motif/cantley_kinase_map.tsv` | Matrix name → UniProt → gene, and whether each matrix is used. |
| `results/features/motif/cantley_imputation_neighbours.tsv` | Nearest atlas kinase and identity for every imputed kinase. |

## F1: Cantley score (not label-derived)

### Matrix-to-kinase mapping

- Sources: Johnson 2023 (S/T, 311 matrices in 1.8.0) and Yaron-Barir 2024 (Tyr, 78 canonical + 15 non-canonical `*_TYR` matrices).
- The UniProt accession comes from the `UNIPROT_ID` column of `kinome_information.tsv`.
- The build asserts that this mapping agrees with every `cantley_matrix_name` row of `results/labels/kinase_synonyms.tsv` (377 names cross-checked).
- It also asserts `PDK1` → O15530 (PDPK1). `PDHK1` maps to Q15118 (PDK1).
- All 404 matrices map to index kinases, one accession per matrix type.
- Coverage: 396 index kinases have at least one own matrix.

### Score definition

This is identical to kinase-library 1.8.0.
- **Scored positions:** S/T −5..+4 (9 rows); Tyr −5..+5 (10 rows). The centre and the other window positions are not scored.
- **Score:** log2 score = Σ_p log2M[p, aa], using the shipped `*_all_log2_matrices.tsv`.
- **Out-of-range residues:** padding `_`, `X` and `U` contribute 0.
- **Unprimed scoring:** every flank is upper-cased, so S/T/Y use the unmodified columns.
- **Rounding:** round to 3 dp; for S/T matrices add log2 st_fav[centre S|T] and round to 3 dp again.
- **Percentile:** 100·#{ref ≤ s}/n_ref (searchsorted, side='right'), rounded to 2 dp.

### Arrays

- Shape: float32 (503 kinases × 350,736 sites), in index order.
- Name pattern: `cantley_{variant}_{metric}`.

| variant | meaning |
|---|---|
| `unprimed` | As published. |
| `primed_hard` | Flank S/T/Y that is itself a U_obs site is scored with the phospho column (s/t/y). This equals kinase-library `pp=True` with those residues lower-cased. |
| `primed_soft` | A flank contributes (1−w)·log2M[upper] + w·log2M[lower]. w is the fraction of the 4 independent MS evidence channels that observed the flank site: PSP MS_LIT+MS_CST ≥1, compendium, CPTAC, Ochoa. w is 0 if the flank is not a U_obs site. |

Notes on the priming variants:
- The query site's own status is never used.
- The soft weight uses no KS information and no low-throughput literature channel.

| metric | meaning |
|---|---|
| `score` | log2 score. |
| `pct_published` | Arm (i). Percentile against the kinase's Ochoa-2020 reference shipped with kinase-library (82,735 S/T / 5,431 Y sites, unprimed). For primed variants this matches what kinase-library returns with `pp=True`. |
| `pct_universe` | Arm (ii). Label-free percentile within U_ms sites of the same residue type (S+T, or Y), computed within the same variant: a primed score is ranked against primed U_ms scores. |

The three arms of B1 are therefore:
- (i) `cantley_unprimed_pct_published`
- (ii) `cantley_unprimed_pct_universe`
- (iii) `cantley_primed_soft_pct_universe`, or `primed_hard`

### NaN rule

A cell is NaN if either:
- the kinase has no matrix for the site's residue type (S/T → ser_thr, Y → tyrosine); or
- the matrix type lies outside the kinase's `residue_class`, taken from `cache/index/kinases.parquet` at build time.

The second case masks 8 non-canonical `*_TYR` matrices of ST-class kinases on Y sites: BMPR2, LIMK1/2, PDK1/3/4, PINK1, TNNI3K. STY-class kinases (MAP2Ks, TESK1, PKMYT1, NEK10, TTK, ...) use their S/T matrix on S/T sites and their Tyr matrix on Y sites. **If residue classes change, rerun the build (2.5 min).**

### Imputation (kept separate from the published arm)

Imputed arrays follow the pattern `cantley_imputed_{variant}_{metric}`.
- They contain **only** the rows of kinases that lack an own in-class matrix: 104 S/T and 16 Y kinase-type slots.
- Each row is copied from the nearest atlas kinase with an in-class matrix.
- "Nearest" means highest kinase-domain identity:
  - domain = UniProt domain coordinates from the index, or the full length for 5 kinases without coordinates;
  - global BLOSUM62 alignment with free end gaps;
  - identity = identical residues / length of the shorter domain;
  - candidates are pre-ranked by alignment score, and the top 8 are re-ranked by identity.
- `cantley_imputed_distance` (503 × 2, columns [ST, Y]) = 1 − identity. Median identity is 0.38 (range 0.24–0.99).
- Example neighbours:
  - CSNK2A3 → CSNK2A1 (0.99)
  - MAP3K13 → MAP3K12 (0.87)
  - TIE1 → TEK (0.84)
  - TRRAP → ATR (0.24)
- Imputed scores are weak (sanity AUROC about 0.70 vs 0.76 for own matrices). The distance feature lets models discount them.

### Reproduction (details in validation.md)

| Target | Agreement |
|---|---|
| kinase-library scored Ochoa reference parquets, 82,735 × 311 and 5,431 × 93 | 98% exact, **100% within 0.001** |
| Johnson 2023 Suppl. Table 3 percentiles, 89,784 × 303 | 99.999% within 0.1 percentile points for the 275 matrices not revised since publication |
| Yaron-Barir 2024 Suppl. Table 3 | 99.0% within 0.1 for 87 of 93 matrices |

The 2% of non-exact values are 3-dp rounding ties. They come from the 4-decimal matrix files, which the package also uses.

Matrices revised after publication:
- S/T (28): CK1D, CDK5, CDK3, SRPK2, GRKs, ... (listed in validation.md);
- non-canonical Tyr (6): BMPR2_TYR, PINK1_TYR, MKK7_TYR, PDHK1/4_TYR, LIMK2_TYR.

Release 1.8.0 is used as "as published".

## F2a: learned PSSM (label-derived; MUST be cross-fitted)

### Usage

```python
from ksp.features.motif_learned import fit_pssm, score
m = fit_pssm(train_pairs[["kinase_acc", "site_id"]], prior="family", phospho=False)  # ~2.6 s, all 503 kinases
X = score(m)            # 503 x 350,736 float32 log-odds, NaN outside the kinase's residue class (~6 s)
```

### Model

- **Window and alphabet:** ±7 window including the centre (so the S-vs-T preference is learned). Alphabet = 20 aa + `_` + `X`.
- **Phospho channel (optional):** adds `s/t/y` for flank S/T/Y that are U_obs sites. The query site's own status is never used.
- **Residue types:** one model per (kinase, residue type ∈ class). Labels outside the class are dropped and counted.
- **Background:** frequencies of U_ms windows of the same residue type.
- **Empirical-Bayes hierarchy:** Dirichlet pseudocounts, with leave-kinase-out priors:
  - q_group = shrink(group counts − k, bg, γ = 10)
  - q_family = shrink(family counts − k, q_group, β)
  - q_k = shrink(C_k, prior, α)
- **Prior options:** `prior='family'`, `'group'` or `'none'`. Kinases with 0 labels fall back to family, then group, then background. `model.level` reports which level was used.
- **Tuning:** α and β come from a quick internal 5-fold CV over training pairs, grouped by substrate cluster (from `ctx.clusters`) or by protein. The objective is the macro (per-kinase) held-out log-likelihood.
- **Tuned values on DEV:** α = 64, β = 256 (phospho channel: α = 128).

### Crossfit class

`ksp.crossfit.REGISTRY["f2a_pssm"]` (`LearnedPSSM`) returns these columns:

| Column | Meaning |
|---|---|
| `f2a_pssm`, `f2a_pssm_pct` | Log-odds score, and its percentile within U_ms of the same residue type. |
| `f2a_pssm_ph`, `f2a_pssm_ph_pct` | The same, with the phospho channel. |
| `f2a_pssm_level` | Prior level used: 0 kinase, 1 family, 2 group, 3 background. |
| `f2a_pssm_log_ntrain` | log1p of the kinase's training-label count. |

The class was smoke-tested through `crossfit.crossfit()` with real substrate clusters (5 inner folds + eval in 66 s, no NaN).

### Sanity (non-official, CV fold 0)

- Macro AUROC over 94 kinases: family prior 0.766, phospho channel 0.782, Cantley 0.785.
- Shrinkage helps over no prior at 10–49 labels: 0.756 vs 0.708.
- The learned PSSM and Cantley are complementary.
- The learned logos recover:
  - R at −3/−2 for PKA;
  - P at +1 for CDK1 and MAPK1;
  - acidic +1..+3 for CK2;
  - Q at +1 for ATM;
  - pS at +4 for GSK3B, only with the phospho channel.

## Assumptions and caveats

1. **Hard priming treats every U_obs flank site as phosphorylated.** This overstates co-occurrence. The soft version weights flanks by MS evidence breadth, which is a proxy, not stoichiometry.
2. **The published percentile reference contains Ochoa sites.** Many of them are in U_ms, and some carry labels. This is a known mild circularity of the published method.
3. **Percentiles are stored as float32, not float16.** float16 resolution near 100 is 0.06 percentile points, which would collapse the top 0.1%, the range that matters for ranking.
4. **Priming does not help everywhere.** It improves priming-dependent kinases (GSK3A/B +0.1 AUROC, CK1 family) but costs about 0.02 on basophilic kinases (AKT1, PKC, CHEK1). Feed it as a separate feature.
5. **The learned-PSSM sanity numbers are non-official.** Validation used `ksp.labels.load_dev` (DEV only) and CV fold 0. Nothing in `src/` reads labels: the caller passes training pairs.
