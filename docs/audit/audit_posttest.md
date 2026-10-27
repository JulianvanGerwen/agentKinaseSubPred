# Audit: post-TEST exploratory CV (A34): M2 intercepts (a), priming gate (b), Tyr (c)

Auditor: independent reviewer. Run 2026-10-27 02:34 CET to 03:16 CET (wall time ≈ 41 min, most of it full re-runs).
Scope: `docs/posttest_cv.md`, `docs/posttest_tyr.md`, their scripts and their artefacts, and REPORT §11.1.
- Nothing under `sealed/` and no TEST labels were read. No git. No existing file was modified.
- All re-runs went to the scratchpad: each script was imported with its `OUT` redirected.
- Audit scripts are in the session scratchpad (`a_saved.py`, `a_full.py`, `a_uw_ci.py`, `b_extra.py`, `rerun.py`). They are not part of the repo.
- Threads: ≤8 in total. Peak RSS per process was ≤ 27.6 GB (Tyr re-run; m2 24.1 GB, priming 23.1 GB). These are above the docs' "~15–16 GB" claims but inside the budget.

## Verdicts

| experiment | claim | verdict |
|---|---|---|
| (a) per-kinase intercepts, M2 | "PROMISING under A34": 5–19 bin +0.085 [+0.048, +0.116], pooled M2 +0.037 | **CONFIRMED-WITH-CAVEATS**. The numbers and the nesting are correct. But the gain exists only for the post-hoc "balanced" variant: with the pre-declared plain variant, the target stratum gets *worse*. The macro gain is a reallocation: pair-level (micro) Hit@5 falls by 0.070. The 0-label "dark" kinases, which M2 cannot score, are pushed out of the top 5. |
| (b) priming-gated Cantley blend | "PROMISING under the literal A34 rule": primed +0.018 [+0.013, +0.024] | **CONFIRMED-WITH-CAVEATS**. The numbers and the nesting are correct, and the author's "reallocation" reading holds. It is stronger than stated: a label-free 0/1 indicator ("the residue at the priming offset is a U_ms phosphosite"), with no Cantley term at all, reproduces the gain. Against negatives that carry the same flag, the gain is +0.002. |
| (c) Tyr: Y-specific weight (A) and motif repair (B) | both fail the A34 bar; neither adopted | **CONFIRMED** (negative result). Wording fixes are needed for the per-family claim and for the "~+0.019 without nesting" sentence. |

## Common checks (all three)

**Baseline reproduces the frozen fusion.** Base CV numbers are M1 0.786408, M2 0.262579, M3 0.360863. These equal `results/cv/LEADERBOARD.md` row 2 (`r2_fusion_r3__20261026-191901__s0`: M1 0.786, M2 0.263, M3 0.361), and each script asserts equality with the run's `metrics.json`:
- (a): difference 0.0;
- (b): M1 and M2 differences 0.0;
- (c): the fusion is rebuilt from R2 ens5 + R3 pooled, and overall M1 = 0.7864081926 = the run value.

**Bootstrap unit.** All three experiments use the harness's substrate-cluster bootstrap (`ksp.eval._cluster_weights` / `_boot_matrix`, B = 1000, seed 0, paired on identical weights).
- (b) and (c) use the harness cluster universe `_universe_clusters`, i.e. the clusters of positives plus E-hard negatives.
- (a) uses only the clusters of the positive sites. This is the only universe that matters for M2, and per-cluster weights are ≈ Poisson(1) either way.
- I re-derived the (a) CIs with my own cluster weights, and also with a kinase-resampling bootstrap. Both agree; see below.

**No TEST or sealed access.** No script reads anything beyond `ksp.eval.prepare('cv')`, the DEV CV run artefacts and label-free Cantley matrices or runs.

**One shared second-order caveat (stacking).** The OOF scores of the training folds g ≠ f came from base models whose training included fold f. The second-stage parameters (b_k, γ, w_Y, matrix picks) therefore see fold-f labels *indirectly*, through the base model. This is standard cross-validated stacking, and the frozen fusion weights (0.6 / 0.5) carry the same caveat. It is not a leak of the kind A34 forbids, but it means "nested" covers the second stage only.

---

## (a) Per-kinase intercepts on score_site (M2)

### Nesting: real
- `SiteSet.fit` only touches pairs whose site fold is in `train_folds`. The candidate masks for those pairs use only those sites' own label sets.
- Inner selection fits on 3 folds and scores the 4th through a fold-restricted view (`views[g]`). The outer intercepts are refit on the 4 training folds.
- I refit the 5 selected outer configurations myself. The maximum |Δb| and |Δτ| against `intercepts_nested.npz` are both 0.0.
- **Leakage probe.** I randomly permuted the true-kinase labels of every held-out fold-0 pair and refit the fold-0 intercepts. The maximum |Δb| is 0.0, so held-out labels cannot enter the fit.
- Selection is by inner pooled M2, not by the target bin. Ties go to the larger λ. `nested_selection.tsv` reproduces the selections (balanced λ = 3 in folds 0 and 4, λ = 10 in folds 1–3).

### Recomputed headline (from `site_level_*.parquet`, and by refitting)

| stratum (K_cv) | claimed base → nested, Δ [CI] | recomputed Δ, cluster CI | kinase-resample CI |
|---|---|---|---|
| 5–19 (83) | 0.142 → 0.227, +0.085 [+0.048, +0.116] | 0.1420 → 0.2271, +0.0850 [+0.0481, +0.1155] | [+0.055, +0.118] |
| 20–99 (89) | +0.048 [+0.033, +0.065] | +0.0480 [+0.0330, +0.0647] | [+0.026, +0.070] |
| ≥100 (31) | 0.621 → 0.497, −0.124 [−0.137, −0.110] | −0.1240 [−0.1368, −0.1100] | [−0.155, −0.092] |
| pooled M2 (203) | 0.263 → 0.299, +0.037 [+0.020, +0.050] | +0.0369 [+0.0200, +0.0499] | [+0.018, +0.056] |
| M1 / M3 | unchanged | Δ = 0 exactly (M1 and M3 read `score`, not `score_site`) | — |

All claimed numbers reproduce. Full re-run of the script into the scratchpad: **byte-identical** `grid_cv5.tsv`, `nested_selection.tsv`, `nested_vs_base.tsv`, `sanity.tsv`, and `summary.json` identical except runtime (1265 s wall).

### Post-hoc "balanced" variant: the claim depends on it
The plain (unweighted) softmax was the declared method. The balanced variant was added after one plain grid point had been seen on the *pooled* CV. Joint nested selection does not remove that researcher degree of freedom: the decision to add a variant was made from all-fold results.

**With the pre-declared plain variant only.** I used the saved inner scores, refit, and applied. Selection is λ = ∞ in folds 0, 2 and 4 and λ = 100 in folds 1 and 3.

| stratum | Δ vs base [cluster 95% CI] |
|---|---|
| **5–19 (target)** | **−0.013 [−0.021, −0.004]** |
| 20–99 | −0.007 [−0.012, −0.001] |
| ≥100 | +0.019 [+0.013, +0.026] |
| pooled M2 | −0.006 [−0.009, −0.000] |

So the pre-declared method is significantly *worse* on its own target stratum, and fails A34. The positive finding belongs to the balanced variant alone.

That said, the balanced gain is not a fragile tuning artefact:
- every balanced λ beats base on pooled M2 (0.278–0.303 against 0.263);
- the inner-to-outer optimism is only 0.0055.

The variant itself is the post-hoc part. Its objective is the macro-over-kinases form of M2, i.e. it is fitted to the metric.

### Cost on high-label kinases, and does the gain survive other views?

**Per-kinase direction.**

| bin | improved | worse | unchanged | median Δ |
|---|---|---|---|---|
| 5–19 | 39 | 4 | 40 | 0 |
| 20–99 | 54 | 7 | 28 | +0.04 |
| ≥100 | 1 | 28 | 2 | −0.129 |

**Pair-weighted (micro) Hit@5** over the 8,416 K_cv pairs, i.e. "a random labelled site": **0.535 → 0.465 (−0.070)**. Micro Hit@1 goes 0.317 → 0.260. With all 8,756 real pairs: 0.520 → 0.454.

So the macro gain is a transfer of hits away from the ≥100 kinases, which hold 5,400 of the 8,416 pairs. It is not better site-level assignment overall. Family-level Hit@5 rises only +0.010 (0.480 → 0.490), against +0.037 at kinase level, so much of the gain is reshuffling within families.

**M2 filtering and the intercepts.** The intercept objective uses exactly the M2 (A18) filtered candidate set: class-matched real kinases, minus the site's other true kinases and Par(k). The labels come from training-fold sites only, so there is no leak. There are two interactions, though:
1. **Dark kinases are pushed down at no cost to M2.**
   - The 130 kinases with 0 DEV labels get τ·b ≈ **−0.72** (mean; range −1.64 to −0.36). The 5–19 bin gets +0.47 and ≥100 gets −0.43.
   - These 0-label kinases are candidates in every filtered set but are never scored by M2 (K_cv needs ≥5 labels).
   - At labelled sites, their share of the top-5 slots falls from **14.0% to 2.3%**. The 5–19 share goes 8.9% → 23.5%, ≥100 goes 51.7% → 34.4%, and 1–4 goes 5.7% → 7.8%.
   - Part of the M2 gain therefore comes from demoting the *most* understudied kinases, which the metric cannot see. This runs against the motivation "M2 for understudied kinases".
2. **The sign of check 7 is misreported by omission.** Spearman(τb, log count) = +0.39 over all 503 kinases, which is driven by the 0-label kinases. Within K_cv it is **−0.79**. Both should be reported.

Control: an intercept-only ranking (score_site × 1e-6, used only to break ties) gives M2 0.002. The gain needs the score and is not a pure prior. The placebo (−0.051) also reproduces.

### Required wording (REPORT §11 / docs/posttest_cv.md)
- Replace "PROMISING under A34" with this: "Promising under A34 only for a variant added after the first result. The pre-declared unweighted intercepts lower the 5–19 bin (−0.013 [−0.021, −0.004])."
- Report the micro Hit@5 drop (0.535 → 0.465) and the ≥100 drop (−0.124) next to the gain. Do not call it an improvement in kinase assignment. Call it a macro-M2 reallocation towards 5–99-label kinases.
- State that 0-label kinases are pushed down (top-5 share 14% → 2%) and that M2 cannot see this.
- Check 7: give the within-K_cv correlation (−0.79) as well as the overall +0.39.
- The caveat "handled by joint nested selection" overstates the fix. Nesting controls λ and variant *selection*, not the post-hoc *introduction* of the variant.

---

## (b) Priming-gated blend with priming-aware Cantley

### Nesting: real
- γ*_f is the argmax of the mean per-fold priming-kinase M1 (`auroc_hard_f{g}`) over g ≠ f.
- Under A24, a fold's AUROC depends only on that fold's scores and labels, so fold-f labels never enter γ*_f.
- The z-scores are computed within (fold, kinase) over all candidate rows and use no labels.
- The primed stratum is not used for selection.
- The selections (0.2, 0.3, 0.3, 0.2, 0.3) follow from `grid.tsv`.

### Recomputed numbers
I rebuilt the score independently and applied my own A24 ranks over the 9 K_cv priming kinases:

| row | claimed | recomputed |
|---|---|---|
| primed positives | 0.893 → 0.911 | 0.8926 → 0.9110 |
| unprimed positives | 0.826 → 0.821 | 0.8264 → 0.8207 |
| all positives, priming kinases | 0.847 → 0.848 | 0.8472 → 0.8483 |

All three match `strata_nested_vs_base.tsv`: +0.018 [+0.013, +0.024], −0.006 [−0.010, −0.001] and +0.001 [−0.003, +0.005].

Overall endpoints, from `summary.json` paired bootstrap:
- M1: +0.00005 [−0.00015, +0.00023];
- M3: −0.0017 [−0.0042, −0.0002];
- M2: unchanged.

Full re-run into the scratchpad: **identical** `grid.tsv`, `nested_selection.tsv`, `strata_nested_vs_base.tsv`, and `summary.json` identical except runtime (784 s wall).

### Is "+0.018 primed" just rank shifting? Yes, and the mechanism is built into the stratum's definition
- The stratum "primed" means that the residue at the kinase's priming offset is an observed (U_ms) S/T phosphosite (`ksp.strata.primed_mask`).
- The Cantley arm is run with `pp=True` on U_ms-lower-cased flanks. It therefore rewards exactly that property.

On average 11.5% of the E-hard negatives of the priming kinases carry the same flag. Split by that flag (macro over 8 kinases):

| primed positives compared against | base | nested | Δ |
|---|---|---|---|
| all E-hard negatives (the reported stratum) | 0.8926 | 0.9110 | +0.018 |
| **flag-matched negatives** (offset also a U_ms phosphosite) | 0.7674 | 0.7691 | **+0.002** |
| unflagged negatives | 0.9081 | 0.9288 | +0.021 |
| *(flagged negatives against unflagged negatives)* | 0.702 | 0.772 | +0.070 |

**Label-free control.** I used score = z(fusion) + c·1[priming-offset residue is a U_ms phosphosite], with no Cantley term:

| c | primed AUROC | all-positive AUROC |
|---|---|---|
| 0.2 | 0.909 | 0.849 |
| 0.3 | 0.917 | 0.849 |

This equals or beats the nested Cantley blend (0.911 / 0.848).

The primed gain is therefore almost entirely "move every site with a phosphorylated priming position up". Within a matched comparison it does not discriminate primed true substrates from primed non-substrates. The author's reading ("reallocation, not a net improvement") is confirmed. Also:
- the all-positive row is +0.001 (CI spans 0);
- the unprimed row is −0.006;
- M3 falls by 0.0017 with its CI below 0.

### Notes
- The earlier A32 primed rule was ≥ +0.02 with CI > 0. This result (+0.018) would fail it. A34's bar is looser, which is legitimate but worth one clause in the REPORT.
- "Primed-positive gains are positive in every K_cv priming kinase" is true for all 8, but CSNK2A2 moves only 0.9976 → 0.9977. Two of the 8 kinases have 1 (CSNK1G2) and 4 (CSNK2A2) primed positives.

### Required wording
- Say "meets the literal A34 rule, but the gain is a stratum-definition effect". Include the flag-matched Δ (+0.002) and the indicator control.
- Always print the primed, unprimed and all-positive rows together, plus the M3 Δ with its CI below 0. The doc does this already; the REPORT must too.
- Do not describe it as "better at primed substrates".

---

## (c) Tyr specificity (`docs/posttest_tyr.md`)

### Nesting: real
- w_Y and γ are chosen by the inner-fold mean of within-fold Y AUROC.
- The per-kinase own / nbr / nbr_xfam / none pick uses the Cantley-channel AUROC restricted to inner-fold positives and negatives (`np.isin(fold[rp], inner)`).
- The Cantley channel is a label-free published matrix, z-scored over all of a kinase's rows.
- Held-out labels enter only `outer_*` and the final evaluate.

### Recomputed numbers
These are from `explore.json` and `per_kinase_Y.tsv`. A full re-run into the scratchpad gave an `explore.json` **identical** to the committed one except runtime (735 s wall).

| | claimed | recomputed |
|---|---|---|
| baseline Y M1 (38 kinases) | 0.6387 [0.5925, 0.6674] | mean of the per-kinase table 0.6387 |
| A: Y Δ | −0.0039 [−0.0097, +0.0010] | per-kinase mean 0.6347 (Δ −0.0040) |
| B: Y Δ | +0.0042 [−0.0103, +0.0167] | per-kinase mean 0.6429 (Δ +0.0042) |
| overall M1 Δ, A / B | −0.0005 / +0.0009 | from the evaluate |
| ST drift | 0 | 0 |

The inner/outer transfer table reproduces (B: mean inner gain +0.0187, mean outer gain +0.0056, 3 of 5 folds positive).

### Wording issues
1. **The per-family claim is selective.**
   - B improves 18 of 38 Y kinases and **worsens 19** (1 unchanged). Only the first figure is stated.
   - "The gains land where §1.5 predicted" holds for FGFR (+0.060).
   - Two of the *predicted* anti-correlated families get worse: **TAM/Axl −0.053** (MERTK −0.098) and Tec, the BTK family, −0.031 (BTK itself 0.000).
   - Ret (−0.052), PDGFR (−0.038) and Trk (−0.035) also lose and are not shown.
   - Eph and VEGFR have no K_cv kinase, so that part of the prediction is untested.
   - REPORT §11.1's "the direction and the per-family pattern do match the prediction" should become: "mixed: FGFR/Abl/Jak gain; the predicted TAM and Tec families lose; 18 up / 19 down".
2. **"Reported without nesting … would have claimed ~+0.019"** (doc §3.2 and REPORT §11.1). This was not measured. It is the mean inner-fold gain on the 42-kinase Y+STY fold-mean metric, not a non-nested pooled Y M1 Δ. Say "the inner folds showed +0.019", or run the non-nested fit.
3. **"0 / 5 folds improved" for A.** In 3 folds w_Y = 0.6 was re-chosen (Δ = 0 by construction), and in the other 2 it was negative. State that.
4. **Selection and reporting sets differ.** Selection uses the 42 Y+STY kinases; the headline uses the 38 pure-Y kinases. This is disclosed, but worth one clause.
5. **γ boundary.** γ = 0.5 was selected in all 5 folds, which is the top of the grid. This does not change the verdict (the CI already straddles 0), but "fails" holds only for γ ≤ 0.5.
6. "R2 and R3 are close on Y (0.629 vs 0.624)". These are the single-seed R2 and R3 values from strata_cv §4, not the ens5 R2 used in the fusion. This is minor.

---

## Summary of required REPORT changes
1. (a): name the post-hoc variant in the headline. Give the pre-declared result (5–19 Δ −0.013 [−0.021, −0.004], pooled −0.006). Give the micro Hit@5 drop (−0.070), the ≥100 loss (−0.124) and the demotion of 0-label kinases. Frame it as a macro-M2 trade-off, not an M2 improvement.
2. (b): the primed gain is reproduced by a label-free "phospho-priming-site present" indicator, and is +0.002 against flag-matched negatives. Report it as stratum rank-shifting, alongside unprimed −0.006, all-positive +0.001 and M3 −0.0017.
3. (c): fix the per-family sentence (18 up / 19 down; TAM and Tec lose), and reword the "~+0.019 without nesting" sentence as an inner-fold figure.
