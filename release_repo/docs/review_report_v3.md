> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Independent review of REPORT.md (v3, 1,113 lines)

**Reviewer:** independent results reviewer. I wrote none of the audited work.
**When:** 2026-10-27, 09:29:28–09:45:59 CET (`date`). Wall time **16.5 min**, single session, ≤8 threads, peak RSS well under 60 GB (largest job was a TSV/JSON scan; no model was run).
**Scope:** REPORT.md against the committed artefacts it cites, plus PROMPT.md for completeness.

**Constraints kept.** No `git` was run. Nothing under `sealed/` and no TEST label file was opened; only metrics/superiority/strata/per-kinase/reliability outputs under `results/test/721dc1977366/` were read. No existing file was modified; this file is new.

**What I checked.** Every quantitative claim in §1, §4.1–4.5, §5, §6.1–6.8, §7.1–7.9, §8, §9, §11.1 and §12 that cites a file, against: `results/test/LEDGER.md`, `results/test/721dc1977366/{*_ci.tsv,*_ci_kinases.tsv,*_strata.tsv,*_ext_tiers.tsv,*_per_kinase.tsv,*_t3*.json,final.json,final_reliability.tsv,summary.{json,md},superiority.tsv,secondary.tsv}`, `results/cv/LEADERBOARD.md`, `results/release/{RELEASE_MANIFEST.json,coverage.tsv,kinase_top200.tsv.gz,site_top5.tsv.gz}`, `results/runs/{r2_deploy,r3_deploy,r2r3_release}*/manifest.json`, `results/freeze_entries/*`, `results/baselines/simple/summary.md`, `results/baselines/phosformer/fidelity_check.json`, `results/models/r2/compare__*.tsv`, `results/splits/split_report.md`, `FREEZE.json`, and `docs/{audit/audit_test.md,audit/audit_posttest.md,audit/audit_cv_r2.md,strata_cv.md,posttest_cv.md,posttest_tyr.md,model_r2.md,model_r3.md,baseline_phosformer.md,baseline_selphi2.md,regimes_labels.md,downstream.md,labels_universe.md,features_*.md}`.

**Headline verdict.** Number fidelity is **very high**. All 9 primary TEST tests, all 17 secondary tests, every §7.1 standalone metric and CI, every §7.5/§7.6 stratum value, the §7.7 reliability table (including the merged tiers and the weighted ST figures), the §9 coverage and release counts, the §4.3 fusion sweep, the §4.4 regime table, the §4.5 calibration table, the §5 baseline tables, the §6.2 ablation, the §6.4–6.5 strata, the §6.8 runtimes and the whole of §11.1 reproduce their sources to the quoted precision. The audits' required wording has largely landed. The defects below are stale pre-TEST/pre-deployment sentences, five broken leaderboard row references, one value/citation clash, and one systematic label error in the strata pair counts. No number that bears a scientific claim is wrong.

---

## (A) Must-fix factual errors

**A1. `REPORT.md:1010` — "The deployment model is still running (§9)."**
Stale and self-contradicting. `REPORT.md:3` says "Deployment model and coverage (§9) done", `REPORT.md:1024` strikes the item through as "Done (G5, §9, `results/release/`)", and §9 reports completed runtimes. Evidence it finished: `results/runs/r2r3_release__20261027-012616__s0/manifest.json` (`status` set, `runtime_s` 167.92), `results/release/RELEASE_MANIFEST.json`, `docs/model_r2.md:583` ("### 7.1 G5 deployment: DONE (01:29)"). Delete the bullet.

**A2. `REPORT.md:180–183` — the "Current state" block in §2.5 is pre-TEST text.**
It reads "The draft FREEZE validates with all 9 primaries" and "FREEZE.json is not yet committed." Both are superseded: the batch ran against a committed FREEZE (`results/test/LEDGER.md` STARTED block records `freeze_sha 721dc1977366…` and `git_head c3e2a1448e…`, and §12 at `REPORT.md:1105` lists `scripts/build_freeze.py ; git commit FREEZE.json ; src/ksp/eval_test.py --real` in that order), and the batch COMPLETED at 00:57:29 with all 9 primaries tested, not a draft. Replace with the post-batch state.

**A3. `REPORT.md:455, 458, 459, 460, 463` — five LEADERBOARD row references do not resolve.**
`REPORT.md:466` itself states the board was regenerated at 01:23 with 135 runs; the row numbers were not updated.

| REPORT line | cites | actual row in `results/cv/LEADERBOARD.md` |
|---|---|---|
| 455 (R2 full, seed 0) | "LEADERBOARD #5" | **#15** `r2_full_cv5__20261026-162445__s0` (line 21) |
| 458 (R3 final distilled) | "LEADERBOARD #11" | **#42** `r3_final_distill_cv_pooled__20261026-173414__s0` (line 48) |
| 459 (B7 simple-strong) | "LEADERBOARD #19" | **#53** `simple_b7_lr_reeval__20261026-160451__s0` (line 59) |
| 460 (B1 Cantley (i)/(ii)) | "LEADERBOARD #26–27" | **#60** (arm ii, line 66) and **#63** (arm i, line 69) |
| 463 (B2 learned PSSM) | "LEADERBOARD #38" | **#74/#75** `simple_b2_pssm*` (lines 80–81) |

All the *values* on those rows are correct; only the pointers are wrong.

**A4. `REPORT.md:407` and `REPORT.md:462` — B3a M4-CV = 0.657, but the cited source says 0.655.**
`REPORT.md:462` sources the B3a row to "LEADERBOARD.md (A24 re-eval)". `results/cv/LEADERBOARD.md:79` gives M4 = **0.655** for `b3a_phosformer_cv_pooled__20261026-222734__s0`. The value **0.657** comes from `docs/baseline_phosformer.md:124` (pre-re-evaluation). Either change the number to 0.655 or change the citation — and note that the same row already handles the M1 clash explicitly ("0.677 … the doc gives 0.679"), so the M4 clash should be handled the same way.

**A5. `REPORT.md:400` — "`model_r2.md` and `RESEARCH_TREE.md` 5.14 quote B1 M1 = 0.719."**
Only `docs/model_r2.md:11` quotes 0.719. `RESEARCH_TREE.md:133` (item 5.14) quotes **0.717** and already flags the other value: "vs B1 0.717/0.171/0.133 (0.719 was a stale pre-15:31 value)". Drop RESEARCH_TREE from the sentence.

**A6. `REPORT.md:69` — wrong file cited for the TEST coverage counts.**
"Coverage is 5,773 TEST sites at 0.5 and 0 at 0.8 (`final_reliability.tsv`; audit §6)." Those counts are **not** in `results/test/721dc1977366/final_reliability.tsv` (which holds only the per-bin reliability rows). They are in `results/test/721dc1977366/summary.json` (`coverage@0.5: 5773`, `coverage@0.5_unlabelled: 5428`, `coverage@0.8: 0`) and are confirmed in `docs/audit/audit_test.md:73`. `REPORT.md:967` repeats the same wrong citation. The numbers themselves are correct.

**A7. `REPORT.md:722–724, 767–768` — the "kinases / pairs" column pairs a K_eval kinase count with a K_any pair count.**
The M1 values quoted in those rows are the K_eval (K_primary) ones, but the pair counts are the K_any ones. From `results/test/721dc1977366/final.json` → `pair_strata` (K_primary vs K_any `n_pos`):

| REPORT line | REPORT says | K_primary (matches the quoted M1) | K_any |
|---|---|---|---|
| 722 in vivo only | 88 / **535** | 88 / **468** | 137 / 535 |
| 723 in vitro only | 96 / **807** | 96 / **721** | 152 / 807 |
| 724 both | 107 / **862** | 107 / **771** | 172 / 862 |
| 767 new since 2021 | 61 / **262** | 61 / **215** | 96 / 262 |
| 768 pre-2021 | 118 / **1,942** | 118 / **1,745** | 234 / 1,942 |

This is inherited from the batch's own `*_strata.tsv`, which collapses the two views into one row — the external-curated row is the only one the report splits correctly (`REPORT.md:771`: "49 / 130 (K_any 75 / 183)"), and the json confirms that split. Note `docs/audit/audit_test.md:282` endorses "262 temporal pairs (K_any)" as the batch count, so the 262 figure may stay, but the column must be labelled "K_eval kinases / K_any pairs" or the K_primary counts used. This is a labelling fix, not a wrong statistic.

---

## (B) Overstatements

**B1. `REPORT.md:418` — "B3a's shortfall on our benchmark therefore reflects the split, labels and negatives, not the implementation."**
The fidelity check (`results/baselines/phosformer/fidelity_check.json`) trained on **60 % of Phosformer's own train split for one epoch** and scored 0.888 on their test split. That rules out a broken or badly-initialised pipeline; it does not rule out under-training on *our* folds, which the report itself concedes two lines earlier and at `REPORT.md:980` ("1 seed and a PI-fixed 1:16 negative ratio", ≤3 epochs per `docs/baseline_phosformer.md:68`). `docs/baseline_phosformer.md:169` uses the weaker phrasing "not … an undertrained or broken reimplementation". Soften to match, e.g. "is not explained by a broken implementation".

**B2. `REPORT.md:328` — "**One global model is best.**"**
Supported only over the 10 regimes in `docs/regimes_labels.md`, all at seed 0 and one R2 configuration. Add "of the regimes tested" (the sub-bullets already give the right evidence).

**B3. `REPORT.md:919` (§9 "Reading") — TEST reliability is attributed to the deployment model.**
"on labelled sites the final model's ≥0.5 tier is calibrated on TEST (S/T: 0.615 claimed vs 0.603 observed, §7.7)". That reliability was measured on the **DEV-trained** TEST model. The deployment model is retrained on DEV+TEST and reuses the DEV-OOF calibrators with no OOF of its own — stated at `REPORT.md:908` and in `results/release/RELEASE_MANIFEST.json` `limitations[0]`, but not repeated where the claim is made. Add the transfer caveat in the bullet.

**B4. `REPORT.md:597` — "The final model is both the most accurate and among the cheapest to score."**
Accuracy: true on all 9 primaries, but the B4 M3 margin is "modest" (Δ +0.059, CI lower bound 0.010) and Y-kinase M1 is at the popularity floor, so "most accurate" should carry the §1 caveats by reference. Cost: the final model infers in 316.9 s against B1 4.4 s, B2 39.7 s and B0-pop 0.4 s (`results/test/LEDGER.md`), so "among the cheapest" holds only against the *learned* baselines. Suggest "the cheapest of the learned models to score".

**B5. `REPORT.md:466` — "LEADERBOARD.md … now includes every run listed here."**
The §6.1 B4 row (matched 0.708 / 0.176 / 0.316 / 0.714) is not on the board at those values; `results/cv/LEADERBOARD.md:85–86` carries only the full-universe `b4_selphi2_{noleak,faithful}` runs at 0.589 / 0.587. Reword to "every CV run listed here, at its full-universe values; the B4 matched figures come from `baseline_selphi2.md`".

**B6. `REPORT.md:186` — "Exploration on CV closed at 01:00 (A28)."**
Read literally this conflicts with §11.1, whose A34 runs are dated 01:34–02:34 (`docs/posttest_cv.md:9`) and 02:34–03:16 (`docs/audit/audit_posttest.md:3`). Add "for model selection"; §11.1 already says the released model is unchanged.

**Checked and found *not* overstated.** The audits' binding wording has landed: "p < 0.0002 (bootstrap floor)" (`REPORT.md:46`, audit §0 #6); M2 0.363 never presented as a CV improvement and 0.348 quoted alongside at `:56`, `:838`, `:971` (audit §0 #7); "no demonstrated Tyr-kinase specificity … Cantley is at chance" at `:60–65`, `:954–961` (audit §0 #10); B4 M3 called "modest" at `:41`, `:74`, `:667`, `:982` (audit §0 #5); "our faithful retrained reimplementations under this benchmark's split, labels and negatives" at `:664`, `:983` (audit §3); priming limited to the non-priming-aware arm at `:75`, `:759`, `:976–978` (audit §0 #11); hub and tier-B gains attributed to popularity at `:76`, `:761`, `:986`; temporal/external ECEs explicitly not called miscalibration at `:799–802` (audit §6); the six not-delivered items listed at `:849–858` and `:987–994` (audit §7 #10); count drift acknowledged at `:776` (audit §7 #9). No adopted-model claim leaks from §11.1: "The released model is unchanged" (`:83`, `:1043`), "none is adopted" (`:1045`), "Neither adopted" (`:1025`), and the one user-facing suggestion at `:1054` is immediately fenced with "would need a pre-registered test on new labels".

---

## (C) Gaps against PROMPT.md

| Brief ask | Status | Where / what is missing |
|---|---|---|
| Beat Cantley, Phosformer, SELPHI2.0 | **answered** | §1 table, §5.2, §5.3, §7.2 (9/9 Holm-significant) |
| "Retrain these benchmarks so that it is a fair comparison" | **partly** | B3 and B4 are retrained on our folds (§5.2, §5.3). The Cantley atlas is never retrained and §5.1 never says why (it is label-free; B2, the learned PSSM on our folds, is the retrained motif analogue). One sentence in §5.1 would close this. |
| In vivo vs in vitro | **answered** | §4.4 (L2_vitro −0.003 M1 / −0.047 M2; L2_vivo −0.029 / −0.090), §6.4, §7.5 evidence strata |
| Species pooling (human/mouse/rat) | **answered** | §4.4 L1 vs L2 (M2 +0.029 at no M1 cost) |
| Other annotation sources, critically assessed | **answered** | §3.2 exclusion table (iKiP-DB, Sugiyama, STRING, Ochoa/PSP-derived), §4.4 L2+ext ablation-only |
| Leakage (cluster by protein; split on sites not pairs) | **answered** | §2.1–2.2, §3.2, §6.7 canaries, §7.8 CV→TEST alarm |
| Graded difficulty / leave-kinase-out | **answered** | §2.2 T1/T2/T3, §6.6, §7.4 |
| Harmonious performance across kinases, incl. understudied | **answered** | kinase-macro endpoints (§2.4), M5 bins (§4.4), §6.4, §7.5, §10.3, §11.1(a) |
| Training regimes (global / per-kinase / in between) | **answered** | §4.4 (10 regimes) |
| Unorthodox / priming kinases | **answered** | §2.4 strata, §6.5, §7.5, §11.1(b) |
| Runtime alongside accuracy | **answered** | §6.8, §7.1 (train/infer columns), §9 |
| Licence and data provenance | **answered** | §3.3, §9 (CC BY-NC-SA, non-commercial), DATA.md (56 SHA-256 rows, 57 licence mentions) |
| Coverage of as many phosphosites as possible | **partly** | §9 covers all 350,736 U_obs sites, but U_obs is ~19 % of the 1,870,491 proteome S/T/Y residues (`docs/features_structure_esm.md:38`). U_all is listed as unbuilt at `REPORT.md:1009`, but neither §1 nor §9 states the denominator. Add one clause. |
| Design reviewed in ≤3 rounds by an external agentic PI and specialist scientists | **missing from REPORT** | `docs/reviews/` holds round 1 (external PI, prior work, statistics, data leakage, algorithms/compute, AI failure modes, kinase biology) and round 2 (4 reviewers), and `RESEARCH_TREE.md:42` records "~80 of 84 round-1 changes accepted". REPORT never mentions the review programme or what it changed. |
| Verification scripts with toy data, positive/negative controls, sanity checks | **partly** | Canaries (§2.2, §6.7) and the B0 floors (§2.4, §5.1) are well covered. The 14-file `tests/` suite and the `--toy` leak test appear only as two command fragments (`REPORT.md:1091`, `:1095`); nothing reports what they verify or that they pass. |
| "Fix the random seeds list" | **partly / disclosed deviation** | `REPORT.md:1078`: "There is no single central seed list in code; seeds are recorded per run manifest." Honest, but it is a deviation from the brief and is not flagged as one. |
| Commit results to git | **not assessable here** | git is PI-only; see C-note in §(A)/§(D) on `be635d1`. |
| LESSONS.md | **not referenced** | `LESSONS.md` and `docs/retro/` exist but are never cited in REPORT. Low priority. |

---

## (D) Minor and wording

**D1. §5.1 and §6.1 give different last digits for the same CV quantities** (two sources: `results/baselines/simple/summary.md` vs `results/cv/LEADERBOARD.md`). The report's own Conventions (`REPORT.md:5`) promise one convention.

| quantity | §5.1 | §6.1 |
|---|---|---|
| B0-miss CV M1 | 0.538 (`:396`) | 0.539 (`:464`) |
| B2 CV M3 | 0.108 (`:394`) | 0.109 (`:463`) |
| B2 CV M4 | 0.639 (`:394`) | 0.638 (`:463`) |
| B7 CV M3 | 0.177 (`:395`) | 0.179 (`:459`) |
| B7 CV M4 | 0.709 (`:395`) | 0.708 (`:459`) |

Also `REPORT.md:456` quotes R2-clean M3/M4 as 0.300 / 0.768 (= `docs/model_r2.md:430`) where `LEADERBOARD.md:24` rounds to 0.299 / 0.767. All are the same underlying values; pick one source per table and say which.

**D2. `REPORT.md:398`** — "Arms (i) and (ii) tie on M1 (difference 9e-7)". `results/baselines/simple/summary.md` (best-arm paired table) gives the CV M1 difference as **1e-06**. The +0.0135 [0.007, 0.020] M2 figure is exact (summary gives −0.013534 [−0.020117, −0.00698]).

**D3. `REPORT.md:826`** — a sentence is left in quotation marks: `"On the same 117 kinases, CV M2 = 0.348 vs TEST 0.363 (paired Δ +0.014 ± 0.021). No leak signal."` It is a verbatim lift from `docs/audit/audit_test.md:39`. Either attribute it or unquote it.

**D4. §11.1 omits one audit-required item.** `docs/audit/audit_posttest.md:97,105` requires reporting that Spearman(τ·b_k, log label count) is **+0.39 overall but −0.79 within K_cv**, "by omission" being the flagged problem. It landed in `docs/posttest_cv.md:4` but not in `REPORT.md:1049–1054`.

**D5. `REPORT.md:910`** — the §9 column header "Coverage@0.5 (no-PSP sites)" is ambiguous: the bold number (29,886) is over all 350,736 sites and the bracketed number (26,762) is the no-PSP count. Rename to "Coverage@0.5 (no-PSP in brackets)".

**D6. `REPORT.md:295–299`** — the fusion-weight table drops the w = 0.4 column present in `docs/model_r2.md:316` (M1 0.784 / M2 0.261 / **M3 0.363**). As printed, M3 appears to peak at the chosen w = 0.6 (0.361); in the full grid it peaks at w = 0.4. Immaterial to the choice (made on M1), but the abridgement flatters it.

**D7. `REPORT.md:908`** — "15,202 training pairs, L2". `docs/model_r2.md:584` adds "including canaries", and notes only the 503 real kinases are released. Worth one clause.

**D8. `REPORT.md:116`** — "chosen by simulated annealing (seed 0) to give per-kinase test fractions of 15–25%" states the target as if achieved. `results/splits/split_report.md`: 187 of 228 kinases with ≥5 positives fall in 15–25 % (223 within the integer band), and 10 kinases have 0 TEST positives.

**D9. §7.6 vs the audit on the temporal margin.** `REPORT.md:775` gives "+0.056" (standalone, 0.766 vs 0.710, `final_strata.tsv`); `docs/audit/audit_test.md:212` gives the **matched** figures 0.771 vs 0.713 (+0.058). Both are correct; say which view is quoted, as §7.5's footnote does for the other strata.

**D10. `REPORT.md:824`** — "most Y kinases have only 3–4 TEST positives". From `final_per_kinase.tsv`, SYK has 11 and LYN 14 of the three named collapses. The sentence is inherited verbatim from `audit_test.md:180`; "many" would be safer.

**D11. `REPORT.md:1039–1041`** — a stray `---` plus two blank lines precede "### 11.1", which is an H3 sitting after the §11 table. Either promote it to `## 11.1` or move it under a §11 sub-heading.

**D12. §12 reproducibility — paths and seeds check out; one commit cannot be verified.**
All 32 scripts/paths named in §12 and §6.8 exist on disk, with one deliberate exception: `scripts/run_all.sh` is absent, which §10.15 (`REPORT.md:1011`) already discloses (promised in `DESIGN.md:345`, §11 item 8). `results/freeze_entries/{final,r2_alone,r3,b3,b4}.json`, `results/baselines/simple/freeze_entries/ALL.json`, `results/baselines/selphi2/freeze_entry_B4.json`, `results/splits/SEAL.json`, `logs/dl/sha256_wave1.txt`, `results/index_manifest.json`, `results/cv/REGISTRY.tsv` and both `results/models/r2/compare__*.tsv` are all present and carry the quoted values. The harness commit `974d520` is real: `results/runs/*/manifest.json` records the full hash `974d52074de936026f0125333a04666fdad22f73`. **`be635d1` (`REPORT.md:1069`, "Repository HEAD at drafting time") appears nowhere on disk outside REPORT.md** and could not be verified without git — the PI should confirm it. Seeds as listed (splits 0; R2/R3 0–4; ablations/B3/B4/regimes 0, L2ext also 1–2; bootstrap 0; downstream 42 / tumour 1) are consistent with the manifests and docs I read.

---

## Verification notes (what reproduced exactly)

For the record, the following were recomputed or matched line-for-line and need no action:

- **§1 and §7.1–7.3:** every point estimate and CI against `*_ci.tsv`, `final_ci_kinases.tsv`, `superiority.tsv`, `secondary.tsv` and `LEDGER.md`; all 9 primary Δ/CI/p/Holm and all 17 secondary Δ/CI; matched-universe figures (386 kinases, 74.13 % of pairs, 116 K_eval, matched M2 0.372); runtimes 1,136.6 / 316.9 s and the 25,922,452-pair pre-flight.
- **§7.4:** T3 M1/M3/M4 for all five configs against `*_t3.json`; Y-kinase T3 0.531 over 23 kinases recomputed from `final_t3_per_kinase.tsv`; "beats B0-pop on only 10 of 23 Y kinases" recomputed from `final_per_kinase.tsv` + `B0_pop_per_kinase.tsv` (exactly 10).
- **§7.5–7.6:** all 19 stratum rows × 5 configs, the 9 Manning-group rows, the M2-by-bin table, both external tiers and the external-all M1/M4 split.
- **§7.7:** the merged [0,0.2) and [0.2,0.5) tiers recomputed from the batch's 0.1-wide bins (0.12/0.14 and 0.34/0.32); the ST ≥0.5 tier 0.6155 claimed / 0.6029 observed over n = 345; ECE 0.0257; temporal/external ≥0.5 accuracies 19/53 = 0.36 and 6/28 = 0.21.
- **§9:** `coverage.tsv` row-for-row; 131,144,410 pairs; sha `cceb32c7…`; release files 100,600 and 1,753,680 data rows (verified by `zcat | wc -l`); all three spot checks reproduced from `site_top5.tsv.gz` (TP53 S15 → ATM 0.69863; AKT1 S473 → PDPK1 / AKT1 / MTOR; EIF2S1 S52 → PKC family); every deployment runtime (309.7 / 985.4 / 446.6 / 11.8 / 169) and the 53 / 235–543 early-stopping rounds.
- **§11.1:** every number in the table and the (a)/(b)/(c) prose against `docs/audit/audit_posttest.md`, including the pre-declared-variant penalty, the micro Hit@5 drop, the 14 % → 2.3 % dark-kinase share, the flag-matched +0.002, the γ-boundary and w_Y caveats, and the 18-up/19-down per-family split with TAM/Axl and Tec named as losers.
- **§2–§6 and §8:** split counts, K_eval/K_any/K_cv, label funnel, evidence mix, all feature-table coverages I sampled, the full ablation curve, the fusion sweep, the regime table, the calibration-in-the-large table, both coverage-before/after-A31 sets, the canary per-fold values, all baseline tables, the Phosformer fidelity numbers, the SELPHI imputation numbers, and every downstream Δ and CI.
