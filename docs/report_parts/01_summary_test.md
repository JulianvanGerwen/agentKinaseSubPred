<!-- Replaces in REPORT.md §1: the "TEST (placeholders)" table and the "main caveats" list. Sources: results/test/LEDGER.md, results/test/721dc1977366/, docs/audit/audit_test.md -->

**TEST results: single batch, FREEZE 721dc197, completed 2026-10-27 00:57, independently reproduced.**

Standalone values use missing = 'worst', with 95% cluster-bootstrap CIs (B = 1,000). Each Δ is the final model minus the baseline in the matched-universe primary test (B = 10,000, Holm over 9 tests).

| Endpoint | Final | Cantley B1 | B3 Phosformer (our retrain) | B4 SELPHI2.0 (our retrain) |
|---|---|---|---|---|
| M1 | **0.794** [0.774, 0.811] | 0.714 (arm i); Δ **+0.075** [0.046, 0.099] | 0.708; Δ **+0.086** [0.059, 0.105] | 0.705; Δ **+0.089** [0.063, 0.117] |
| M2 Hit@5 | **0.363** [0.324, 0.397] | 0.189 (arm ii); Δ **+0.179** [0.147, 0.213] | 0.163; Δ **+0.199** [0.151, 0.229] | 0.233; Δ **+0.130** [0.091, 0.156] |
| M3 top-1% recall | **0.363** [0.309, 0.401] | 0.161 (arm i); Δ **+0.205** [0.165, 0.235] | 0.104; Δ **+0.259** [0.213, 0.291] | 0.304; Δ +0.059 [0.010, 0.120] (modest) |

Sources: `results/test/LEDGER.md`; `results/test/721dc1977366/{final,B1_i,B1_ii,B3,B4}_ci.tsv`; `superiority.tsv`.

**All 9 primary tests meet the "beats" rule.**
- 8 of them have p < 0.0002 (bootstrap floor), with Holm p ≤ 0.0018.
- The B4 M3 test has p 0.018 (Holm 0.018).

**Unseen family and unseen protein (T3).** M1 is 0.752 (`results/test/LEDGER.md`).

**Runtime.** Training takes 1,137 s and scoring the 25.9 M TEST pairs 317 s (`results/test/LEDGER.md`).

**Read these numbers with the following in mind.** The independent review's caveats (`docs/audit/audit_test.md` §0, §7) govern how they may be stated.

- **M1 shows no CV→TEST optimism.** On the same 117 kinases, CV M1 is 0.792 and TEST is 0.793 (paired Δ +0.001 ± 0.012).
- **TEST M2 (0.363) is not an improvement over CV.** The rise is K_eval composition: on the same 117 kinases, CV M2 is 0.348.
- **The TEST margin over Cantley (+0.075) is partly a weak-baseline draw.** On CV the margin was +0.063 [0.048, 0.074]. On the same kinases Cantley fell by 0.025 from CV to TEST, mostly on Y kinases.

The main caveats:
1. **No demonstrated Tyr-kinase specificity on TEST.**
   - Y kinases score M1 0.595, barely above the site-popularity floor B0-pop (0.585), and the final model beats that floor on only 10 of 23 Y kinases.
   - Cantley is at chance (0.496).
   - On unseen families (T3), Y kinases score 0.531.
   - Y-kinase M1 dropped by 0.072 ± 0.033 from CV.
   - The review rates this NOT TRUSTED as a skill claim (`docs/audit/audit_test.md` §0 #10, §4).
2. **No reliable high-confidence site assignment.**
   - On labelled TEST S/T sites, the ≥0.5 tier is honest: 0.615 claimed vs 0.603 observed (n = 345).
   - There is no ≥0.8 tier, and no Y site can exceed 0.445.
   - Coverage is 5,773 TEST sites at 0.5 and 0 at 0.8 (`final_reliability.tsv`; audit §6).
3. **M2 is still weak for understudied kinases.** TEST includes only 13 kinases with 5–19 labels, and they are the best-labelled of that bin. The CV finding (Hit@5 0.142 against 0.162 for Cantley) is therefore not refuted (audit §2).
4. **B3 and B4 are our faithful retrains on this benchmark.**
   - The B3 fidelity check passed.
   - 66.6% of B4's TEST pairs use SELPHI's NA-default site features.
   - The B4 M3 win is modest: CI lower bound 0.010, with 52 kinases favouring the final model and 27 favouring B4.
5. **Priming.** The primed-site advantage (0.897 vs 0.704) holds only against non-priming-aware Cantley. On unprimed positives the two are near-tied (0.740 vs 0.725). On CV the priming-aware Cantley arm led (0.945); that arm was not in FREEZE.
6. **Strata carry no CIs in the batch, so all stratum statements are descriptive.** The gains on hub sites and external tier B are mostly popularity (B0-pop 0.788 and 0.724).
7. **Downstream, the library is at parity with PSP, not better, and trails PSP on unseen drugs** by 0.04–0.06 (§8).
8. **About 90 configurations were tried on CV.** TEST confirms the audited optimism bound of ≲0.005.
