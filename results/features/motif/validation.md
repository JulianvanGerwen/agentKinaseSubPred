# Motif features (F1 Cantley, F2a learned PSSM): validation

Scripts:
- `scripts/features/motif_cantley_validate.py` → `reproduction_summary.json`, `reproduction_per_kinase.tsv`
- `scripts/features/motif_validate.py` → `validation_auroc_cantley.tsv`, `validation_site_centric_cantley.tsv`, `validation_learned.tsv`, `motif_top_positions.tsv`, `validation_runtime.json`

Supplementary tables were extracted with `scripts/features/motif_extract_supp.py` into `cache/features/motif/validation/`.

## 1. Reproduction of the published Cantley scores

| Target | Sites × kinases | Agreement |
|---|---|---|
| kinase-library 1.8.0 scored Ochoa reference, S/T log2 score | 82,735 × 311 (25.7 M values) | **98.0% exact** at 3 dp; **100% within 0.001**; Pearson 0.99999999 |
| same, Tyr | 5,431 × 93 (0.5 M) | **97.8% exact**; **100% within 0.001** |
| percentiles of those reference sites (ours vs. from the shipped scores) | same | 98.3% exact; max \|Δ\| 0.13 (S/T), 0.44 (Y) percentile points |
| Johnson 2023 Suppl. Table 3, S/T percentiles | 89,784 × 303 (27.2 M) | all: 43.5% exact, 90.9% within 0.1. **275 unrevised kinases: 48.0% exact, 99.999% within 0.1**, median per-kinase Spearman 0.99999997 |
| Yaron-Barir 2024 Suppl. Table 3, Tyr percentiles (non-canonical sheet) | 7,315 × 93 | all: 63.3% exact, 92.6% within 0.1. **87 unrevised kinases: 67.6% exact, 99.0% within 0.1**, Spearman ≈ 1.0 |
| Johnson 2023 Suppl. Table 4, ATM percentile | 37 sites | max \|Δ\| 0.012; 100% within 0.1 |

### Why agreement is not 100% exact

**Matrix precision.** The 1.8.0 log2 matrices are stored with 4 decimals. Summing 9–10 of them can move the unrounded score by up to ±0.0005, which flips the 3-dp rounding when the score lies near a half-step. The kinase-library package uses these same files, so on new sequences our scores equal what the package itself returns.
- Our unrounded scores are never more than 0.001 from the published reference (max 0.0009997).
- With unrounded scores and an unrounded reference, the agreement with the paper's percentiles rises to about 80–89% exact. That points to full-precision matrices used for the paper tables.
- We keep the package convention: 3-dp rounding plus st_fav.

**Matrices revised after publication.** Agreement is low (Spearman < 0.999) for:
- 28 S/T matrices: CK1D, GRK7, VRK2, IKKA, CLK3, PKN3, LATS1, GRK5, MST1, SKMLCK, SMG1, PASK, GRK1, MLK3, MEKK6, CAMKK2, WNK3, MNK2, PBK, GRK4, CHAK1, HUNK, CDK5, NIM1, CDK3, ICK, SRPK2, CDK13;
- 6 non-canonical Tyr matrices: BMPR2_TYR, PINK1_TYR, MKK7_TYR, PDHK1_TYR, LIMK2_TYR, PDHK4_TYR.

The 1.8.0 data files differ from the paper's matrices for these kinases. 1.8.0 is used as the "as published" source.

**Residual ≤0.1 differences in percentile.** These are one or a few reference sites (1/82,735 or 1/5,431 of the distribution) that sit at rounding ties.

Score definition confirmed empirically: with st_fav the reference matches 98.0% exactly, without it 60.2%.

## 2. Cantley sanity AUROC (NON-OFFICIAL)

Setup:
- Positives: DEV L2 labels (`load_dev('L2')`, U_ms only, kinases with ≥5 in-class positives).
- Negatives: 2,000 random U_ms sites of the kinase's residue type, excluding any DEV label of that kinase.
- No paralog masks.

These are sanity checks, not project metrics.

| Kinases | n | score / pct_published / pct_universe | primed_hard (pct_univ) | primed_soft (pct_univ) |
|---|---|---|---|---|
| own matrix | 240 | 0.758 / 0.758 / 0.758 | 0.752 | 0.752 |
| imputed (nearest atlas kinase) | 7 | 0.702 / 0.695 / 0.700 | 0.651 | 0.668 |
| own S/T + imputed Y (MAP2Ks, TTK) | 5 | 0.690 / 0.613 / 0.645 | 0.501 | 0.517 |

Within one kinase, AUROC cannot separate the published arm from the universe-percentile arm, because both are monotone in the score. The arms differ only when kinases are compared with each other, as in the site-centric check below.

Selected kinases (score AUROC → primed_hard / primed_soft):
- PRKACA 0.880 → 0.878 / 0.878
- CDK1 0.879 → 0.884 / 0.888
- MAPK1 0.896 → 0.890 / 0.891
- CSNK2A1 0.886 → 0.874 / 0.884
- ATM 0.929 → 0.929 / 0.931
- SRC 0.657 → 0.678 / 0.676

Priming helps the priming-dependent kinases clearly:
- **GSK3B 0.674 → 0.788 / 0.788**
- **GSK3A 0.769 → 0.860 / 0.874**
- **CSNK1E 0.614 → 0.758 / 0.725**
- CSNK1A1 0.710 → 0.728 / 0.741
- CSNK1D 0.701 → 0.724 / 0.719

Priming slightly hurts basophilic kinases: AKT1 0.879 → 0.861, PRKCA 0.769 → 0.745, CHEK1 0.803 → 0.771. Use it as a separate feature, not a replacement.

**Site-centric check (M2-like, sanity).** For each labelled site, we take the filtered rank of the true kinase among all own-matrix kinases of the site's type (about 300 for S/T). Macro over 334 kinases, 11,693 pairs:

| Ranking by | Hit@5 | Hit@10 | MRR |
|---|---|---|---|
| raw score | 0.158 | 0.230 | 0.114 |
| pct_published | 0.144 | 0.226 | 0.112 |
| pct_universe | 0.156 | 0.229 | 0.117 |
| primed_soft pct_universe | 0.144 | 0.216 | 0.109 |

## 3. Learned PSSM (F2a; NON-OFFICIAL held-out = DEV CV fold 0 proteins)

- Fitted on the other 4 folds' L2 DEV labels.
- Positives: fold-0 labels.
- Negatives: 1,000 random U_ms sites on fold-0 proteins.
- 94 kinases with ≥5 held-out positives.
- Macro AUROC:

| Training-label bin | n kin | PSSM family prior | group prior | no prior | family + phospho channel | Cantley pct (own or imputed) |
|---|---|---|---|---|---|---|
| 0 | 1 | 0.569 | 0.430 | 0.500 | 0.467 | 0.519 |
| 1–9 | 5 | 0.642 | 0.645 | 0.639 | 0.668 | 0.667 |
| 10–49 | 47 | 0.756 | 0.759 | 0.708 | 0.770 | 0.786 |
| 50+ | 41 | 0.798 | 0.797 | 0.777 | 0.817 | 0.804 |
| all | 94 | 0.766 | 0.766 | 0.732 | 0.782 | 0.785 |

The 0 and 1–9 bins hold 1 and 5 kinases, so their AUROCs are noisy. Per-kinase rows are in `validation_learned.tsv`.

- Shrinkage toward family or group helps clearly over no prior at 10–49 labels: 0.756 vs 0.708.
- The phospho channel helps on well-labelled kinases: AKT1 0.898 → 0.915, SRC 0.734 → 0.751, PRKCA 0.785 → 0.811.
- The learned PSSM is complementary to Cantley rather than better:
  - better on SRC (0.734 vs 0.652), AKT1 and CAMK2A;
  - worse on MAPK1/3, CDK2 and the CK1 family.

**Tuned shrinkage** (internal 5-fold CV by substrate protein, macro held-out log-likelihood):
- family prior: α = 64, β = 256 (γ = 10 fixed);
- with the phospho channel: α = 128.

### Top-weighted cells (log2-odds)

From `motif_top_positions.tsv`, learned PSSM on all DEV L2 vs the Cantley matrix:

| Kinase | Learned PSSM (all DEV L2) | Cantley matrix | Expected |
|---|---|---|---|
| PRKACA (n=571) | −2R 2.60, −3R 2.55, +1F 1.30, −2K 1.15 | −2R 3.21, −3R 3.12, −2K 1.95 | R at −3/−2 ✔ |
| CDK1 (582) | +1P 2.58, −2P 1.22, +3K 1.16, +3R 1.06 | +1P 3.64, +3K 2.43, −2P 1.64 | P at +1 ✔ (+3 K/R) |
| CSNK2A1 (451) | +1D 2.41, +3E 2.25, −1D 1.72, +5D 1.61, +3D 1.60, +2D 1.54 | +3E 3.19, +1D 3.18, +4s/t 2.69 | acidic +1..+3 ✔ |
| ATM (199) | +1Q 3.95 | +1Q 3.06, −1y 1.95 | Q at +1 ✔ |
| GSK3B (284) | unprimed: +1P 2.04, **+4S 1.76**; phospho channel: **+4s 2.09**, +1P 2.05, −4s 1.67, +4t 1.31 | **+4s/t 6.06** | pS/pT at +4 ✔ (only with phospho channel / pp) |
| MAPK1 (483) | +1P 2.58, −2P 2.27, −3P 1.19 | +1P 3.40, −2P 1.60 | P at +1 (PxSP) ✔ |
| CSNK1D (71) | weak (+2M, +2W, ...); phospho: −3s 1.47 | −3s/t 4.23 | pS at −3 only weakly learned (few labels) |
| PLK1 (226) | +1M 1.92, −2D 1.61, +1F 1.56 | −2D 1.99, +1F 1.66 | D/E at −2, hydrophobic +1 ✔ |
| SRC (428) | −3D/E, +3P, −1V/I (weak; Tyr motifs are weak) | −1I 2.12, +3F/Y 1.86 | — |

ERK docking: not applicable (motif-only features).

Caveat for the phospho channel: rare flank codes (s/t/y at ±6/±7) can get spurious log-odds of about 1–1.5 (e.g. PRKACA +7y). They are shrunk, but they are noisy.

## 4. Runtime (≤12 cores, no GPU)

| Step | Time |
|---|---|
| Cantley build, all 350,736 sites × 404 matrices × 3 variants × 3 metrics, plus nearest-atlas imputation and the 18 arrays written via `ksp.store` | 151 s |
| Reproduction validation (25.7 M + 27.2 M + 0.7 M values) | 25 s |
| `fit_pssm`, all 503 kinases incl. α/β grid CV: family prior | 2.6 s |
| `fit_pssm`: + phospho channel | 3.7 s |
| `fit_pssm`: group / none prior | 0.6 s |
| `score(model)`, 503 × 350,736 | 6.2 s |
| `crossfit.crossfit(LearnedPSSM)`, 5 inner + 1 eval fit/transform (two models each, with universe percentiles) | 66 s |
| Extracting Johnson Suppl. Table 3 (1.7 GB XML) | about 6 min, once |
