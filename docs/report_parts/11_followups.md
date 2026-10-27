## 11. Follow-ups and remaining ideas

The TEST batch, FREEZE and ledger are done (`results/test/LEDGER.md`, 2026-10-27 00:57), and so is the Phosformer fidelity check (`results/baselines/phosformer/fidelity_check.json`); both are removed from this list.

**Any further evaluation must use new data, a new split, or a logged `--override-completed` (A29).** The sealed TEST set has been spent.

| Priority | Idea | Rationale / evidence |
|---|---|---|
| — | ~~Deployment model (DEV+TEST) plus U_obs coverage~~ | Done (G5, §9, `results/release/`) |
| 1 | Tyr specificity: nearest-atlas imputation for anti-correlated matrices (FGFR, Eph, VEGFR, TAM, BTK); a pre-registered A-loop/hub sensitivity analysis; PTB-adaptor recruitment | TEST Y 0.595 vs popularity floor 0.585; T3 Y 0.531 (`audit_test.md` §4); `features_tyr.md` §1.5–1.6. **Post-TEST exploration done (A34, `docs/posttest_tyr.md`): both cheap score-level fixes fail the A34 bar on CV** — Y-specific fusion weight −0.0039 [−0.0097, +0.0010]; nested nearest-atlas motif repair +0.0042 [−0.0103, +0.0167]. Neither adopted. Remaining value is in the endpoint (hub/A-loop stratum), not more CV search |
| — | ~~Kinase-gated priming blend~~ | Done post-TEST (§11.1 (b)): a rank shift reproduced by a label-free priming-site flag; +0.001 on all positives. A real test needs flag-matched negatives |
| 2 | M2 for understudied kinases: a blend toward Cantley below 20 labels (not yet tried) | CV 5–19 bin 0.142 vs Cantley 0.162; TEST has only 13 such kinases (`audit_test.md` §2). The intercept route is done (§11.1 (a)): it trades M2 between kinases |
| — | ~~Per-kinase intercepts for the site-centric score~~ | Done post-TEST (§11.1 (a)): the macro-M2 gain needs a post-hoc balanced variant and costs ≥100-label kinases −0.124 |
| 2 | Y-site assignment calibration: Y probabilities are capped at 0.445 and mildly overconfident | `audit_test.md` §6 |
| 2 | True leave-drug-out retrain (F6 with `exclude_drugs`) | separates the drug effect from the Hijazi-dataset effect in the downstream gap |
| 2 | Prospective or external evaluation for the analyses the batch could not provide (hub-excluded Y, temporal/external calibration, paired T3, stratum CIs), for example on a future PSP release | §7.9; TEST cannot be reopened |
| 3 | Site-uniform v3 negative canary | turns the pooled-only gate into a clean null |
| 3 | R5 revisited: LoRA or last-k-layer fine-tuning, a lower trunk learning rate, protein-grouped early stopping | the trunk overfit by step 100–600 |
| 3 | decryptM potency coherence (EC50 vs Kd); Bayer–Kuster 2025 dose–response data | closer to direct-substrate evidence (DATAINEED 4, 8) |
| 3 | Temporal retrain on 2021-only labels; L3 variant; KolossuS / AF-Multimer pilot | deferred in DESIGN §0 |
| 3 | Seeds 1–4 for B3; 100 bags for B4 | baseline strength (only relevant for a future benchmark) |

---


### 11.1 Post-TEST exploratory results (A34)

These ran after the TEST batch was spent, on **DEV CV only**, with every tunable chosen inside each outer fold (nested). The released model is unchanged. Bootstrap CIs resample substrate clusters (B = 1,000). An independent audit re-ran all three scripts and reproduced every number exactly (`docs/audit/audit_posttest.md`). One caveat is shared by all three: they are second-stage fits on OOF predictions whose base models saw the held-out fold's *training* partner folds (ordinary stacking). Write-ups: `docs/posttest_cv.md` for (a) and (b), `docs/posttest_tyr.md` for (c).

**None of the three is a clean improvement, and none is adopted.**

| Target | Intervention (nested) | Target stratum Δ [95 % CI] | Cost elsewhere | Audit verdict |
|---|---|---|---|---|
| (a) M2, understudied kinases | per-kinase intercepts on `score_site`, **kinase-balanced variant added after seeing the plain variant's first grid point** | M2 5–19-label bin 0.142 → 0.227, +0.085 [0.048, 0.116]; pooled macro M2 0.263 → 0.299 | M2 ≥100-label bin 0.621 → 0.497 (−0.124; 28/31 kinases worse); site-averaged (micro) Hit@5 0.535 → 0.465; unlabelled kinases' share of top-5 slots 14 % → 2.3 %. **Pre-declared plain variant alone: 5–19 bin −0.013 [−0.021, −0.004], pooled −0.006** | confirmed with caveats |
| (b) primed sites | priming-gated blend with priming-aware Cantley (γ 0.2–0.3 per fold) | primed positives 0.893 → 0.911, +0.018 [0.013, 0.024] | unprimed −0.006 [−0.010, −0.001]; all priming-kinase positives +0.001 [−0.003, +0.005]; M3 −0.0017 (CI < 0); overall M1 ±0.000 | confirmed with caveats |
| (c) Tyr specificity | Y-specific fusion weight | Y M1 (38 K_cv kinases) 0.6387 → 0.6347, −0.0039 [−0.0097, +0.0010] | none (ST bit-identical) | confirmed (fails) |
| (c) Tyr specificity | motif repair: own / nearest / cross-family Cantley Tyr row | Y M1 → 0.6429, +0.0042 [−0.0103, +0.0167] | none (ST bit-identical) | confirmed (fails) |

**(a) moves M2 around rather than improving it.** The balanced intercepts take top-5 slots from well-studied kinases and from kinases with no labels at all, which M2 never scores, and give them to the 5–19-label kinases. That raises the kinase-macro M2 and lowers the per-site one. The gain exists only for the variant added after the first result was seen. The pre-declared variant makes the target bin worse. A per-kinase intercept is cheap to apply, so a user who wants rankings balanced across kinases can add one. It would need a pre-registered test on new labels before it could be claimed as an improvement.

**(b) is a rank shift, not new information.** The primed stratum is defined by the same property the priming-aware Cantley arm rewards: a phosphosite at the priming position. A label-free 0/1 indicator for "the priming position is an observed phosphosite", with no Cantley term, reproduces the gain (0.909–0.917). Against negatives that carry the same flag, the gain is +0.002. It would also fail the earlier A32 bar of ≥ +0.02.

**(c) both fail the A34 bar.** γ was selected at the top of its grid (0.5) in all five folds, so "fails" holds for γ ≤ 0.5. The Y-specific weight re-chose w = 0.6 in three folds (Δ = 0 by construction) and was negative in the other two. Motif repair gains +0.012 to +0.027 on the inner folds but only +0.0056 on the held-out folds (3 of 5 positive). With 38 kinases and ~1,150 Y positives, this stratum cannot support selecting even one parameter. The per-family pattern is mixed:
- **Gains:** 18 of 38 kinases, with FGFR +0.060, Abl +0.038 and Jak +0.024.
- **Losses:** 19 kinases. They include two of the families `features_tyr.md` §1.5 predicted to benefit, TAM/Axl −0.053 and Tec −0.031, as well as Ret −0.052, PDGFR −0.038 and Trk −0.035.
- **Untested:** Eph and VEGFR, which have no K_cv kinase.

B0-pop on the Y stratum is 0.592, so the fusion's lead over popularity on Y is marginal on CV as well as on TEST (§7). Tyr remains the main open problem.
