# Tyrosine-kinase features (F9) and the diagnosis of the Y E-hard gap

Owner: feat-tyr. DEV CV only; no TEST access. Labels read exclusively through `ksp.labels.load_dev()`;
every feature below is **non-label-derived** and saved through `ksp.store`, so no `ksp.crossfit` family was
needed.

| Path | What |
|---|---|
| `src/ksp/features/tyr.py` | SH2 matrix reader/scorer, UniProt adaptor-domain parser, ITAM/ITIM motif logic |
| `scripts/features/tyr_diag.py` | Diagnosis (§1). Outputs `results/features/tyr/diag_*.tsv/json` |
| `scripts/features/tyr_build.py` | Builds all F9 arrays into `cache/features/tyr/` (16 s) |
| `scripts/features/tyr_eval.py` | Sanity evaluation (§3): single-feature E-hard AUROC + TK-restricted LightGBM ablation |
| `results/features/tyr/` | `diag_per_kinase.tsv`, `diag_evidence.tsv`, `diag_pool_composition.tsv`, `diag_crosskinase_auroc.tsv`, `diag_summary.json`, `eval_single_feature.tsv`, `eval_lgbm_per_kinase.tsv`, `eval_lgbm.json`, `kinase_adaptor_domains.tsv`, `tyr_build_meta.json` |

All numbers below are pooled out-of-fold on the frozen 5-fold DEV CV, E-hard negatives, per-kinase AUROC
macro-averaged, using the official B1 arm `motif/cantley_unprimed_pct_published` (A23) unless stated.

---

## 1. Diagnosis: the Y E-hard problem is real, but ~1/4 of the headline gap is macro-mean dilution

Headline: macro E-hard AUROC 0.534 for the 65 residue-class-Y kinases vs 0.732 for the 265 ST kinases
(`results/baselines/simple/summary.md`). Decomposing that 0.198 gap:

### 1.1 Small-n dilution (artefact; explains ~0.095 of the gap)

| n_pos cut | n Y kinases | Y macro | n ST kinases | ST macro | gap |
|---|---|---|---|---|---|
| ≥1 (headline) | 65 | 0.534 | 265 | 0.732 | 0.198 |
| ≥5 | 39 | 0.593 | 169 | 0.749 | 0.156 |
| ≥20 | 10 | 0.629 | 72 | 0.781 | 0.152 |

27 of 65 Y kinases (42 %) have ≤3 evaluation positives but hold only **4.2 %** of all Y positives, so they
move the macro mean a lot and the biology not at all. Eight Eph receptors and three VEGFRs each have exactly
one positive, and they are *systematically* below chance (EPHA1 0.037, EPHA4 0.108, EPHB1 0.132, KDR 0.124,
TYRO3 0.198). That is not noise around 0.5 — see §1.3.

**At matched n_pos ≥ 20 the gap is still 0.152, so three quarters of the deficit is genuine.**

### 1.2 Paralog masking does NOT starve the Y negative pool — rejected

Median E-hard negatives per Y kinase = **779** (ST: 5,091). The floor is the 9-member Src family
(~280 negatives after Par(k) masking), which is still ample for a stable AUROC. The Y pool is small simply
because DEV contains only **810** labelled Y sites in total, not because masking removes them.

### 1.3 Label noise by evidence class does NOT explain it — rejected

Y labels *are* more in-vivo-only than S/T labels (36.5 % vs 22.9 % of pairs, `evidence_native`). But the
within-kinase paired comparison on the 17 TKs with ≥5 positives of each class shows essentially no
difference:

**macro AUROC 0.603 (in vitro / both) vs 0.596 (in vivo only), Δ = +0.007** (`diag_evidence.tsv`).

In-vivo-only Y labels are not harder for the atlas than in-vitro ones. The hypothesis is rejected as stated.

### 1.4 What is actually wrong (A): E-hard for Y is dominated by promiscuous docking hubs

| | Y pool | S/T pool |
|---|---|---|
| sites in the DEV E-hard pool | 810 | 5,007 |
| fraction whose substrate is itself a protein kinase | **20.5 %** | 9.6 % |
| fraction of pairs on kinase proteins | **23.8 %** | 9.8 % |
| max kinases labelled to one site | **15** | 19 |
| fraction of sites with ≥4 labelled kinases | 5.6 % | 4.6 % |
| fraction of sites sharing a 15-mer with another pool site | 3.6 % | 1.7 % |
| fraction of the pool on its 10 most-labelled proteins | 12.2 % | 5.2 % |

The most-shared Y sites are the canonical pTyr docking hubs:
GLO1 Y136 (15 kinases), STAT5A Y694 (9), PLCG1 Y783 (9), ITGB1 Y795 (8), SHC1 Y427 (7),
PTK2/FAK Y925/Y576/Y861/Y194 (5–6 each).

These sites are each kinase's **positives and simultaneously every other kinase's E-hard negatives**. Most
tyrosine kinases genuinely do phosphorylate them; no motif model can or should separate them. This, not
motif weakness, is the dominant structural cause.

Inspecting the sub-chance families makes the mechanism concrete
(`scripts/features/tyr_diag.py` prints these; see also the scratch listing in §1.1):
- All of EPHA2/A5/A6/A7/A8, FLT1 and FLT4 carry **GLO1 Y136** as a positive — one site assigned to seven
  RTKs, almost certainly one promiscuous screen.
- AXL is labelled on six Src-family **activation-loop** sites (`RLIEDNEYTARQGAK` and variants: SRC Y419,
  LCK Y394, FYN Y420, LYN Y397, HCK Y411, YES1 Y426), all in-vivo-only. These are SFK
  *autophosphorylation* sites; "AXL activates SFKs" is a pathway statement, not a direct KS event.
  12 such pairs over 6 sites and 6 kinases are in the Y evaluation set.
- MERTK/TYRO3/AXL share STAT5A Y694 / STAT5B Y699 and the adjacent-duplicate windows
  ERRFI1 Y394/Y395 and TNK2 Y859/Y860.

### 1.5 What is actually wrong (B): the Tyr atlas separates families only weakly

Cross-kinase AUROC matrix (`diag_crosskinase_auroc.tsv`): row = scoring kinase's atlas matrix,
column = another kinase's positive set (shared sites removed), 18 well-populated TKs.

- mean off-diagonal AUROC **0.593**; 48 % of cells ≥ 0.6
- **within Manning family 0.508** — SRC vs FYN vs LCK vs LYN are indistinguishable. This is correct biology
  and Par(k) masking already handles it.
- **between family 0.600** — real but modest separation.

The matrix is strongly asymmetric, i.e. it is mostly a *row* effect ("does this matrix like its own
substrates against a generic pY background") rather than true pairwise discrimination. Rows with distinctive
recognition work (INSR 0.43–0.85, IGF1R 0.56–0.82, SYK 0.59–0.88 — NPxY/YxxM and ITAM). Rows that are
**anti-correlated with their own labels** are FGFR1 (mean 0.42), BTK (0.39), and by §1.1 the Eph, VEGFR and
TAM/Axl families. For those kinases the atlas matrix is worse than useless and imputation would be better.

### 1.6 Conclusion and the PI decision this implies

The Y deficit is real. Roughly: **1/4 macro-mean dilution, and the rest split between (A) an E-hard pool made
of genuinely promiscuous docking hubs and indirect in-vivo events, and (B) genuine atlas weakness for a
specific set of receptor families.** It is *not* paralog masking and *not* in-vivo-vs-in-vitro label quality.

Because (A) is a property of the endpoint rather than of any model, **no feature family can close this gap**.
Two options for the PI, neither of which I may take unilaterally (they change a frozen endpoint):
1. report a Y-specific stratum that excludes sites labelled to ≥ N kinases (a "hub site" stratum), as a
   diagnostic next to the frozen M1/M5; and/or
2. treat `kinase ∉ Par(k)` labels on a substrate's own activation loop as mask-only (A7-style) rather than as
   E-hard negatives, which would remove the AXL→SFK-A-loop class of indirect events.
Both would be endpoint changes and must be pre-registered, not chosen after seeing results.

---

## 2. The F9 feature family

### 2.1 Sources

| source | what | license/provenance |
|---|---|---|
| `external/cantley_yaronbarir2024/41586_2024_7407_MOESM7_ESM.xlsx`, sheet `sh2_norm_scaled_matrices` | **76 SH2 domains** × 6 positions (−2, −1, +1…+4 relative to the pY) × 20 aa, position-normalised and scaled to random amino-acid frequency (1.0 = no preference) | Yaron-Barir 2024 Nature suppl. Table 5, recorded in DATA.md |
| `external/uniprot/human_tsv/UP000005640_annotations.tsv.gz` (`Domain [FT]`, `Keywords`) | SH2 / SH3 / PID(PTB) domain counts, Transmembrane and Receptor keywords | UniProt 2026_03 |
| `cache/features/context/loc/vec_{hpa,hcm}_protein` | substrate plasma-membrane / cell-junction score | existing F4 stores |
| `external/uniprot/.../UP000005640_reviewed_isoforms.fasta.gz` | canonical sequences for ITAM/ITIM context | UniProt 2026_03 |

**Coverage of the SH2 panel.** 24 of the 28 SH2-bearing human kinases have a matrix
(ABL1/2, BLK, BMX, BTK, CSK, FER, FES, FGR, FRK, FYN, HCK, ITK, LCK, LYN, MATK, PTK6, SRC, SRMS, SYK, TEC,
TXK, YES1, ZAP70; SYK and ZAP70 via their tandem N/C domains, taken as the max). **JAK1, JAK2, JAK3 and TYK2
are missing** — their SH2-like domains were not profiled. The other 51 matrices are adaptors
(GRB2/7/10/14, SHC1–4, CRK/CRKL, NCK1/2, PIK3R1–3, PLCG1/2, PTPN11, VAV1–3, SHIP1/2, …) and are used only as
a kinase-independent "does this neighbourhood recruit SH2 adaptors at all" background. **No PTB/PID matrices
exist in any source I have**, so the PTB arm of the brief (IRS/SHC recruitment for INSR, NTRK, EGFR) could not
be built; see DATAINEED.

### 2.2 Phospho-information rule (important)

The SH2 docking features read the U_ms (mass-spec observation) status of the **other** Y sites on the
substrate protein. The query site's own observation status is never used, so this is exactly the rule already
accepted for the F1 `primed_soft` variant (DESIGN §7 F1) and does not violate A6/DL2-7. No KS label of any
kind enters these features, which is why none is a `ksp.crossfit` family.

### 2.3 Arrays (`ksp.store.load_matrix("tyr/<stem>")`)

**kinase × site** (503 × 350,736 float32; NaN outside {Y site} × {Y/STY-class kinase}):

| stem | definition | non-NaN kinases |
|---|---|---|
| `sh2_dock_max` | max SH2 log2 score over the **other** observed (U_ms) pY sites of the same protein, excluding the query site | 24 |
| `sh2_dock_local` | log2 Σ_j 2^score_j · exp(−|Δpos_j| / 50 aa): a distance-weighted "effective local pY-docking affinity" | 24 |
| `sh2_dock_near` | max over other observed pY within ±50 aa (NaN if none) | 24 |
| `sh2_dock_delta` | `sh2_dock_max` − SH2 score of the query site's own window (does the kinase prefer a *neighbour* as its docking anchor?) | 24 |
| `sh2_dock_spec` | **z-score of `sh2_dock_max` across all 76 SH2 domains at the same site** — the kinase-specific component, with the generic "good docking nearby" / site-popularity component divided out | 24 |
| `rtk_pm` | substrate plasma-membrane + cell-junction score, exposed **only** on transmembrane-receptor kinase rows (NaN elsewhere): an explicit RTK × membrane interaction | 54 |
| `itam_gated` | site ITAM order, exposed only on non-receptor TK rows (Syk, Tec, Src, Csk, Abl, Fer, Ack, FAK families) | 26 |

**site vectors** (350,736 × 4; NaN on non-Y sites):

- `site_sh2` = `n_other_py_protein`, `n_other_py_50aa`, `adaptor_sh2_max`, `adaptor_sh2_local`
  — **`n_other_py_*` and `adaptor_*` are site-popularity proxies** (B0-pop alone scores 0.583 on Y kinases,
  above Cantley's 0.534). They are kept in a separate column block so any F9 gain can be shown not to be
  popularity, and they must be ablated.
- `site_itam` = `itam_half` (Yxx[LIV]), `itam_full` (tandem Yxx[LIV] 6–12 aa away), `itim` ([ILVS]xYxx[LV]),
  `itam_order` (0/1/2).

**kinase vector** `kinase_desc` (503 × 6): `n_sh2`, `n_sh3`, `n_ptb`, `is_tm`, `is_receptor`,
`has_sh2_matrix`. (No human protein kinase carries a UniProt PID/PTB domain, so `n_ptb` is identically 0 and
should be dropped.)

---

## 3. Sanity evaluation on the DEV CV (TK kinases, E-hard)

### 3.1 Single feature alone (`eval_single_feature.tsv`)

Macro per-kinase E-hard AUROC, Y/STY kinases, restricted to the kinases where the feature is defined.

| feature | n kinases | macro (all) | macro (n_pos ≥ 5) |
|---|---|---|---|
| `f1_cantley_pct` / `_score` (reference) | 68 | 0.544 | **0.593** |
| **`f9_sh2_dock_max`** | 23 | 0.533 | **0.571** |
| `f9_n_other_py_50aa` (popularity proxy) | 73 | 0.481 | 0.554 |
| `f9_sh2_dock_local` | 23 | 0.531 | 0.554 |
| `f3_ppi_dist` | 69 | 0.560 | 0.551 |
| `f9_adaptor_sh2_local` | 65 | 0.494 | 0.545 |
| `f9_sh2_dock_near` | 23 | 0.528 | 0.543 |
| `f9_itam_gated` | 26 | 0.515 | 0.532 |
| **`f9_sh2_dock_spec`** (specificity-controlled) | 23 | 0.491 | **0.525** |
| `f9_sh2_dock_delta` | 23 | 0.517 | 0.525 |
| `f9_itam_half` / `f9_itam_order` | 73 | 0.486 | 0.519 |
| `f9_rtk_pm` | 32 | 0.637 | 0.509 |
| `f9_itim`, `f3_ppi_edge_any`, `f4_loc_cos_*`, `f4_site_topo_code` | — | — | 0.48–0.51 |

**SH2 docking is the strongest non-motif feature available for tyrosine kinases** (0.571 against the Cantley
atlas's 0.593) — but the honest control kills most of the story: once the across-domain z-score removes the
generic "this Y sits in a pY-dense, SH2-friendly neighbourhood" component, `sh2_dock_spec` falls to 0.525.
Alone, the signal is therefore mostly *kinase-non-specific* pY-cluster density, not SH2 sequence specificity.

### 3.2 TK-restricted LightGBM (`eval_lgbm.json`, `eval_lgbm_per_kinase.tsv`)

Fold-aware 5-fold CV over the frozen substrate-cluster folds, 40 TKs with ≥5 positives. Positives are the
fold-training evaluation positives; negatives are sampled 1:50 (50 % E-hard-type) from fold-training U_ms Y
sites after the KP / Par / Hamming≤2 masks; out-of-fold scores rank-normalised within fold and within kinase
(A24). Baseline = Cantley (score + pct) + F3 (PPI edge/studies/RWR/dist) + F4 (loc cosines, topology).
Averaged over seeds 0, 1, 2 (mean per-kinase seed sd 0.033).

| variant | n feats | macro (40 TKs) | n_pos ≥ 10 | n_pos ≥ 20 |
|---|---|---|---|---|
| `base` (Cantley + F3 + F4) | 10 | 0.5682 | 0.5942 | 0.5912 |
| `base + f9_sh2only` (the 5 SH2 docking matrices only) | 15 | **0.6009** | **0.6126** | 0.6165 |
| `base + f9_nopop` (all F9 except popularity proxies) | 27 | 0.5843 | 0.6089 | **0.6217** |
| `base + f9_all` | 31 | 0.5991 | 0.5976 | 0.6166 |
| `base + f9_nosh2` (F9 without SH2 docking) | 26 | 0.5905 | 0.5977 | 0.6064 |
| `base + f9_poponly` (popularity proxies only) | 14 | 0.5830 | 0.5942 | 0.5961 |

Paired per-kinase deltas, bootstrap over kinases (B = 10,000; note this is **not** the official
`ksp.eval.bootstrap`, which resamples substrate clusters and would be wider):

| comparison | n | mean Δ | 95 % CI | p_boot | improved |
|---|---|---|---|---|---|
| f9_sh2only − base, all TKs | 40 | +0.033 | [+0.011, +0.054] | 0.002 | 29/40 |
| f9_sh2only − base, **kinases that actually have an SH2 matrix** | 18 | +0.019 | [−0.018, +0.057] | **0.33** | 11/18 |
| f9_nopop − base, SH2 kinases | 18 | +0.022 | [−0.016, +0.060] | 0.26 | 11/18 |
| f9_nopop − base, all TKs | 40 | +0.016 | [−0.007, +0.039] | 0.17 | 23/40 |
| f9_poponly − base, all TKs | 40 | +0.015 | [−0.011, +0.040] | 0.26 | 24/40 |

**Read this carefully.** Splitting the `f9_sh2only` gain by whether the kinase has an SH2 matrix at all:

| | n | base | +f9_sh2only | Δ | improved |
|---|---|---|---|---|---|
| has an SH2 matrix | 18 | 0.592 | 0.611 | +0.019 | 11/18 |
| SH2 columns are all-NaN | 22 | 0.548 | 0.592 | **+0.044** | 18/22 |

The *larger* apparent gain is on kinases whose SH2 columns are entirely NaN. They cannot be using SH2
docking; the NaN pattern instead acts as an "is this an SH2-bearing non-receptor TK" indicator and changes
the jointly-trained model for everyone. **That part of the gain is a model-capacity/NaN-indicator artefact,
not biology.** On the kinases that genuinely have the feature the gain is +0.019 and not significant.

Two results are robust and useful, though:
- `f9_poponly` adds nothing at n_pos ≥ 10 (+0.000), so the gain is **not** site popularity.
- Every F9 variant beats `base` at n_pos ≥ 20 (0.591 → 0.616–0.622), i.e. the direction is consistent on the
  best-measured kinases.

### 3.3 Per Manning family (seed-averaged; full table in `eval_lgbm_per_kinase.tsv`)

| family | n kin | n pos | base | +f9_sh2only | +f9_nopop |
|---|---|---|---|---|---|
| Src | 9 | 601 | 0.626 | 0.636 | 0.640 |
| Abl | 2 | 139 | 0.655 | 0.595 | 0.593 |
| Syk | 2 | 66 | 0.649 | 0.682 | 0.704 |
| EGFR | 1 | 61 | 0.538 | 0.563 | 0.556 |
| Jak | 3 | 43 | 0.472 | 0.494 | 0.525 |
| Tec | 3 | 37 | 0.533 | 0.573 | 0.562 |
| InsR | 2 | 34 | 0.655 | 0.685 | 0.647 |
| FGFR | 3 | 28 | 0.400 | 0.507 | 0.495 |
| FAK | 2 | 22 | 0.679 | 0.704 | 0.633 |
| Met | 2 | 21 | 0.614 | 0.651 | 0.667 |
| Axl | 2 | 17 | 0.485 | 0.476 | 0.480 |

Consistent with §1.5: the families whose atlas matrices are anti-correlated with their labels (FGFR 0.400,
Jak 0.472, Axl 0.485) are where any added feature helps most, simply because the motif prior is actively
wrong there. Syk and Tec — the ITAM/SH2-dependent kinases where the SH2 hypothesis predicts a gain — do
improve (+0.03 to +0.06). Abl loses, which I cannot explain and do not want to explain away.

---

## 4. Recommendation for R2

**Include (cheap, fold-independent, non-label-derived):**
- `tyr/sh2_dock_max`, `tyr/sh2_dock_local`, `tyr/sh2_dock_spec` — 3 columns, 24 kinases. Best non-motif Y
  feature; gain unproven on its own kinases (p = 0.33) but directionally consistent at n_pos ≥ 20.
- `tyr/kinase_desc` columns `is_tm`, `is_receptor`, `has_sh2_matrix` (drop `n_ptb`, identically 0). These are
  genuinely informative *kinase class* descriptors and are most of what §3.2 shows the model exploiting.
- `tyr/site_itam` `itam_order`, and `tyr/itam_gated` — tiny but positive for Syk/Tec.

**Include only behind an ablation flag:** `tyr/site_sh2` `n_other_py_*` / `adaptor_sh2_*` (popularity
proxies; they add nothing at n_pos ≥ 10 and will inflate the negative canaries).

**Do not expect these to close the Y gap.** §1.4 and §1.6 are the actionable findings: the residual Y deficit
is mostly a property of the E-hard definition on a 810-site pool full of promiscuous docking hubs and
indirect in-vivo events, plus atlas matrices that are anti-correlated with labels for FGFR/Eph/VEGFR/TAM/BTK.
The highest-value next steps are a PI decision on the hub-site stratum (§1.6) and nearest-atlas imputation
for the anti-correlated families, not more features.

**Caveats I want on the record.**
- The §3.2 bootstrap resamples kinases, not substrate clusters, so its CIs are narrower than the official
  `ksp.eval.bootstrap` would give; nothing here has been run through `ksp.eval`/`ksp.manifest` and no F9
  configuration has been scored against the canaries. Any R2 run that adopts F9 must do both.
- `tyr_eval.py` is a self-contained sanity harness, not the official metric path. Its baseline numbers
  (0.568 macro over 40 TKs) are not comparable to the headline B1 numbers, which use `ksp.eval`, the full
  K_cv kinase set and `missing='worst'`.
- `sh2_dock_*` uses other sites' U_ms membership. If the PI rules that out (a stricter reading than the F1
  `primed_soft` precedent), the whole family dies.
