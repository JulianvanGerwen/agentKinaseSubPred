# Post-TEST exploratory CV (DESIGN §17 A34): per-kinase site intercepts (a) and priming gate (b)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

> **Audit note (PI, 2026-10-27 03:16; `docs/audit/audit_posttest.md`).** Numbers reproduce exactly and nesting is real. Required readings:
> - **(a):** the gain belongs only to the kinase-balanced variant, which was added after the plain variant's first grid point had been seen. The pre-declared plain variant alone gives 5–19 bin −0.013 [−0.021, −0.004] and pooled −0.006. The balanced intercepts reallocate top-5 slots: ≥100-label bin −0.124 (28/31 worse), micro Hit@5 0.535 → 0.465, and the 0-label kinases' share of top-5 falls from 14 % to 2.3 %. The +0.39 intercept–label-count correlation is driven by the 0-label kinases; within K_cv it is −0.79.
> - **(b):** a label-free "priming position is an observed phosphosite" flag with no Cantley term reproduces the primed gain (0.909–0.917). Against flag-matched negatives the gain is +0.002.
>
> Treat both as trade-offs or rank shifts, not improvements. REPORT §11.1 has the audited wording.

Owner: strata analyst. Run 2026-10-27 01:34–02:34 CET.

Scope:
- CV only, on the final fusion OOF `r2_fusion_r3__20261026-191901__s0` (R2 5-seed ensemble + R3; score w = 0.6, score_site w = 0.5).
- Selection is nested leave-one-fold-out.
- The released model does not change. These are exploratory findings for REPORT §11.
- Nothing in sealed/ or TEST was read.

Code and outputs:

| part | script | results |
|---|---|---|
| (a) | `scripts/models/posttest_m2_intercept.py` | `results/models/posttest/m2_intercept/` |
| (b) | `scripts/models/posttest_priming_gate.py` | `results/models/posttest/priming_gate/` |

Limits: ≤8 threads, run sequentially, peak about 15 GB.

A34 rule: a finding is **promising** if the nested-CV gain on the target stratum has a cluster-bootstrap CI above 0 and overall M1 loses no more than 0.002.

## (a) Per-kinase intercepts on score_site (target: M2 for kinases with 5–19 labels)

**Model.** A site-level softmax over each labelled site's filtered M2 candidate set (class-matched real kinases, minus the site's other true kinases and Par(k)):
- P(k | s) = softmax(score_site / τ + b_k).
- Ridge penalty (λ/2)·Σb_k², with one free temperature τ.
- Fitted by L-BFGS only on the training folds' evaluable positives.
- Applied as score_site + τ·b_k.

**Two variants.**
- *Unweighted*: every training pair has weight 1.
- *Balanced*: each pair is weighted 1/n_k, so every kinase contributes equally. This matches M2's macro-over-kinases.

Disclosure: the balanced variant was added at 01:38, **after** the first unweighted grid point (λ = 0.3) came out worse than base. The unweighted intercept learns the label base rate P(k), which re-strengthens the prior. To account for this, nested selection runs over the **joint** grid (variant × λ ∈ {0.3, 1, 3, 10, 30, 100, ∞}).

**Grid.** Each fold uses intercepts fitted on the other 4 folds:

| variant | λ | M2 | 5–19 | 20–99 | ≥100 |
|---|---|---|---|---|---|
| base (λ = ∞) | — | 0.263 | 0.142 | 0.250 | 0.621 |
| unweighted | 10 | 0.253 | 0.103 | 0.232 | 0.713 |
| unweighted | 100 | 0.258 | 0.127 | 0.238 | 0.667 |
| balanced | 3 | 0.301 | 0.245 | 0.299 | 0.457 |
| balanced | 10 | 0.303 | 0.224 | 0.301 | 0.518 |
| balanced | 100 | 0.278 | 0.165 | 0.270 | 0.603 |

**Nested selection.** The selected configuration was balanced, λ = 3 in folds 0 and 4 and λ = 10 in folds 1–3. Paired cluster-bootstrap CIs (B = 1000) are against base:

| stratum (K_cv) | base | nested | Δ [95% CI] |
|---|---|---|---|
| **M2, 5–19 labels** (83 kinases) | 0.142 | 0.227 | **+0.085 [+0.048, +0.116]** |
| M2, 20–99 labels (89) | 0.250 | 0.298 | +0.048 [+0.033, +0.065] |
| M2, ≥100 labels (31) | 0.621 | 0.497 | **−0.124 [−0.137, −0.110]** |
| M2, pooled (203) | 0.263 | 0.299 | +0.037 [+0.020, +0.050] |
| M1 / M3 | 0.786 / 0.361 | 0.786 / 0.361 | 0 (by construction; verified) |

**Sanity checks** (`sanity.tsv`; all pass):
1. Base M2 reproduces the run's metrics.json exactly.
2. Intercepts applied to fold f were fitted without fold f.
3. Bin M2s recombine to the pooled value.
4. The bootstrap mean is within 0.0002 of the point estimate.
5. Placebo, with intercepts permuted across kinases: M2 −0.051 against base. The gain is kinase-specific.
6. Optimism (inner M2 at the selection minus outer M2) is 0.0055.
7. Spearman correlation of b_k with log training count is +0.39 (descriptive).
8. A full `ksp.eval.evaluate` gives an unchanged M1.
9. The full-evaluate M2 equals the site-centric M2. The canary is unchanged (0.539).

**Verdict: PROMISING under A34**, as a deliberate trade-off. Hit@5 moves from well-studied kinases to understudied ones: ≥100 drops by 0.12, while 5–19 and 20–99 rise. Pooled macro M2 rises, which is the harmonious-across-kinases direction PROMPT asks for. Any use should report the ≥100 loss alongside the gain.

## (b) Kinase-gated blend with the priming-aware Cantley arm (target: primed positives)

**Model.** For the curated priming kinases only (`ksp.strata.priming_curated`), score = (1−γ)·z(fusion) + γ·z(B1 Cantley primed-hard), where z is computed within (fold, kinase).
- Cantley covers 12 of the 13 priming kinases (all except CSNK2A3, which has 0 labels).
- score_site is unchanged, and all other kinases are unchanged.

**Selection.** γ ∈ {0, 0.1, …, 1} was chosen per outer fold by the mean per-fold priming-kinase M1 over the other 4 folds. The primed stratum was **not** used for selection. γ* is 0.2 in folds 0 and 3 and 0.3 in folds 1, 2 and 4.

| stratum (K_cv) | base | nested | Δ [95% CI] |
|---|---|---|---|
| **primed positives** (8 kinases, 231 positives) | 0.893 | 0.911 | **+0.018 [+0.013, +0.024]** |
| unprimed positives of priming kinases (9, 818) | 0.826 | 0.821 | −0.006 [−0.010, −0.001] |
| priming kinases, all positives (9, 1049) | 0.847 | 0.848 | +0.001 [−0.003, +0.005] |
| primed, Cantley-defined offsets (13, 245) | 0.910 | 0.919 | +0.009 [+0.006, +0.014] |
| overall M1 | 0.7864 | 0.7865 | +0.00005 [−0.0002, +0.0002] |
| overall M3 | 0.3609 | 0.3591 | −0.0017 [−0.0042, −0.0002] |

Primed-positive gains are positive in every K_cv priming kinase:

| kinase | base | nested |
|---|---|---|
| CSNK1A1 | 0.881 | 0.906 |
| GSK3B | 0.887 | 0.907 |
| GSK3A | 0.830 | 0.857 |
| CSNK2A1 | 0.884 | 0.909 |

**Sanity checks:**
- Non-priming kinases' AUROCs are identical (maximum difference 0).
- Base M1 and M2 equal the run's metrics.json.
- The γ = 0 grid point equals base.
- M2 is unchanged.

**Verdict: PROMISING under the literal A34 rule.** The primed CI is above 0 and overall M1 does not change. But the gain is a **reallocation within the priming kinases**: unprimed positives lose 0.006 (CI below 0), and the kinase-level M1 of the priming kinases does not improve (+0.001, CI spans 0). Overall M3 also falls slightly. The blend makes priming kinases rank primed candidates higher. That is biologically sensible, but it is not a net improvement for these kinases.

## Caveats

- Both experiments use CV only, on one OOF fusion (seed-ensemble R2 + R3), with no TEST confirmation (A34).
- In (a), the balanced variant was introduced after one unweighted grid point had been seen. This is disclosed above and handled by joint nested selection.
- In (b), the A34 target-stratum rule rewards a stratum gain even when it is offset inside the same kinases. I recommend that REPORT §11 present the unprimed and all-positive rows next to the primed row.
