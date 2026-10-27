## 11. Follow-ups and remaining ideas

The TEST batch, FREEZE and ledger are done (`results/test/LEDGER.md`, 2026-10-27 00:57), and so is the Phosformer fidelity check (`results/baselines/phosformer/fidelity_check.json`); both are removed from this list.

**Any further evaluation must use new data, a new split, or a logged `--override-completed` (A29).** The sealed TEST set has been spent.

| Priority | Idea | Rationale / evidence |
|---|---|---|
| 1 | Deployment model (DEV+TEST) plus U_obs coverage (§9) | Priority 2 of the brief; in progress |
| 1 | Tyr specificity: nearest-atlas imputation for anti-correlated matrices (FGFR, Eph, VEGFR, TAM, BTK); a pre-registered A-loop/hub sensitivity analysis; PTB-adaptor recruitment | TEST Y 0.595 vs popularity floor 0.585; T3 Y 0.531 (`audit_test.md` §4); `features_tyr.md` §1.5–1.6 |
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
