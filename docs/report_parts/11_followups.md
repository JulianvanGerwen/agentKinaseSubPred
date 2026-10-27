## 11. Follow-ups and remaining ideas

The TEST batch, FREEZE and ledger are done (`results/test/LEDGER.md`, 2026-10-27 00:57), and so is the Phosformer fidelity check (`results/baselines/phosformer/fidelity_check.json`); both are removed from this list.

**Any further evaluation must use new data, a new split, or a logged `--override-completed` (A29).** The sealed TEST set has been spent.

| Priority | Idea | Rationale / evidence |
|---|---|---|
| — | ~~Deployment model (DEV+TEST) plus U_obs coverage~~ | Done (G5, §9, `results/release/`) |
| 1 | Tyr specificity: nearest-atlas imputation for anti-correlated matrices (FGFR, Eph, VEGFR, TAM, BTK); a pre-registered A-loop/hub sensitivity analysis; PTB-adaptor recruitment | TEST Y 0.595 vs popularity floor 0.585; T3 Y 0.531 (`audit_test.md` §4); `features_tyr.md` §1.5–1.6. **Post-TEST exploration done (A34, `docs/posttest_tyr.md`): both cheap score-level fixes fail the A34 bar on CV** — Y-specific fusion weight −0.0039 [−0.0097, +0.0010]; nested nearest-atlas motif repair +0.0042 [−0.0103, +0.0167]. Neither adopted. Remaining value is in the endpoint (hub/A-loop stratum), not more CV search |
| 2 | Kinase-gated priming: use primed-hard Cantley only for priming kinases at their own offset, or a per-kinase blend; evaluate against the priming-aware Cantley arm | CV: B1 primed 0.945 vs 0.893; TEST unprimed is a near-tie; A32(a) additive features failed |
| 2 | M2 for understudied kinases: per-kinase prior removal fitted on CV, or a blend toward Cantley below 20 labels | CV 5–19 bin 0.142 vs Cantley 0.162; TEST has only 13 such kinases (`audit_test.md` §2) |
| 2 | Y-site assignment calibration: Y probabilities are capped at 0.445 and mildly overconfident | `audit_test.md` §6 |
| 2 | True leave-drug-out retrain (F6 with `exclude_drugs`) | separates the drug effect from the Hijazi-dataset effect in the downstream gap |
| 2 | Prospective or external evaluation for the analyses the batch could not provide (hub-excluded Y, temporal/external calibration, paired T3, stratum CIs), for example on a future PSP release | §7.9; TEST cannot be reopened |
| 3 | Site-uniform v3 negative canary | turns the pooled-only gate into a clean null |
| 3 | R5 revisited: LoRA or last-k-layer fine-tuning, a lower trunk learning rate, protein-grouped early stopping | the trunk overfit by step 100–600 |
| 3 | decryptM potency coherence (EC50 vs Kd); Bayer–Kuster 2025 dose–response data | closer to direct-substrate evidence (DATAINEED 4, 8) |
| 3 | Temporal retrain on 2021-only labels; L3 variant; KolossuS / AF-Multimer pilot | deferred in DESIGN §0 |
| 3 | Seeds 1–4 for B3; 100 bags for B4 | baseline strength (only relevant for a future benchmark) |

### 11.1 Post-TEST exploratory results (A34)

Run after the TEST batch was spent, on **DEV CV only**, with nested selection; the released model is
unchanged. Full write-up: `docs/posttest_tyr.md`; artefacts: `results/models/posttest_tyr/`.

**Target (c), Tyr specificity.** Two cheap score-level interventions on the frozen R2+R3 fusion's
pooled CV OOF predictions, both selected inside each outer fold:

| | overall M1 | ST M1 | Y M1 (38 K_cv kinases) | Y Δ [95 % CI] | p_boot |
|---|---|---|---|---|---|
| frozen fusion (w = 0.6) | 0.7864 | 0.8208 | 0.6387 [0.5925, 0.6674] | — | — |
| Y-specific fusion weight | 0.7859 | 0.8208 | 0.6347 | −0.0039 [−0.0097, +0.0010] | 0.128 |
| Tyr motif repair (own / nearest / cross-family atlas, nested) | 0.7873 | 0.8208 | 0.6429 | +0.0042 [−0.0103, +0.0167] | 0.549 |

Neither clears the A34 bar (CI above 0), so neither is adopted. ST M1 is unchanged to machine
precision in both. B0-pop on this stratum is 0.592, so the fusion's Y advantage over popularity is
marginal on CV as well as on TEST.

The informative part is the failure mode: the motif repair gains +0.012 to +0.027 on the *inner*
folds in all five outer folds, and only +0.0056 on the held-out folds (3 of 5 positive). With 38
kinases and ~1,150 Y positives the stratum cannot support selecting even one or two parameters —
reported without nesting the same experiment would have claimed a ~+0.019 Y gain. The direction and
the per-family pattern do match the prediction in `features_tyr.md` §1.5 (FGFR +0.060, Abl +0.038,
Jak +0.024), so the idea is not dead; it needs a pre-registered rule tested on new data, not another
CV search.
