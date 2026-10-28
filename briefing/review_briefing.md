# Independent review of the research briefing (briefing.md, briefing.html, 8 figures)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Reviewer: independent agent; I wrote none of the reviewed files. No git was run, nothing under `sealed/` and no TEST label or site-level file was opened, and no data from `data/` or `external/` is reproduced here. Only this file was written. Line numbers refer to `briefing/briefing.md` as of 2026-10-28 10:29.

**Overall.** Number fidelity is high. I recomputed or traced 60+ numbers (list in §6), and all of them match their sources except the items in §1. The audits' required caveats are almost all present in the body text, but two are missing from "In brief". The main problems are one wrong time in a figure caption, one wrong statement about PU precision in a caption, a few claims stronger than the audits allow, missing definitions for several terms, and a main text about 25% longer than the brief asks for.

---

## 1. Errors (must fix)

**E1. Fig 8 caption, line 270: wrong gap to the gate.** "the test was spent at T+11.6 h, 5.3 h before its own gate". The gate is T+16 h (05:20), so 16 − 11.6 = 4.4 h. The 5.3 h figure is the FREEZE commit (00:04 → 05:20; `docs/retro/final_retro.md` §1.2 item 3).
→ Replace with: "The model stopped improving at T+6 h. FREEZE was committed at T+10.7 h, 5.3 h before its planned gate, and the single test run ended at T+11.6 h."

**E2. Line 13 (In brief): same conflation.** "opened the single sealed test set at the halfway point of a 24-hour budget, 5.3 h before its own schedule". The opening was at 00:10, 5.2 h before the gate; the 5.3 h figure belongs to FREEZE.
→ "The team also froze the model 5.3 h before its own planned gate and used up the single sealed test set by the halfway point of the 24-hour budget. Several rules that existed only as text were not followed."

**E3. Fig 7 caption, line 229: wrong PU statement.** "Precision on labelled sites is a lower bound for, and may differ from, precision on unlabelled sites." Nothing supports a lower-bound relation between labelled-site precision and unlabelled-site precision. The audit says only that it "is not a precision on unlabelled sites" (`docs/audit/audit_test.md` §6).
→ "These precisions apply to sites that already have a curated kinase. They say nothing about precision on unlabelled sites."

**E4. Line 87: overstated access control.** "only the final test code could read them". The independent TEST auditor also read `sealed/test_labels.parquet`, read-only (`audit_test.md` "What I did"). In addition, before it was fixed the single-run guard could be bypassed (line 279).
→ "Their labels were stored with a checksum. A software test checked that only the test code read them; after the run, the independent auditor read them read-only to reproduce the results."

**E5. Line 35 (glossary) and line 87: the ledger does not by itself block a second run.** The harness refuses a second run, but an explicit override exists (A29 `--override-completed`; `final_retro.md` §3.1).
→ Glossary: "**FREEZE and ledger:** a file that fixes every prediction file by checksum before the test, and a log of the start and end of the test run. The test code refuses a second run unless an explicit, logged override is used." Line 87: "A ledger logged the start and end of the test run, and the test code refuses a second run without an explicit override."

**E6. Line 280: the post-TEST audit did not "overrule" anything.** Both verdicts were "CONFIRMED-WITH-CAVEATS": the numbers reproduced, but the interpretation was deflated (`docs/audit/audit_posttest.md`, Verdicts).
→ "**Two 'promising' post-test results downgraded.** The post-test audit found that one M2 gain existed only for a variant added after a first result was seen, and that a label-free 'priming site present' flag reproduced the other gain."

**E7. Line 281: incomplete account of the report review.** The review also flagged six overstatements (B1–B6) and several gaps (`docs/review_report_v3.md` §B–C).
→ "**A report review** found 7 factual errors (stale sentences, broken row references, one mislabelled count column) and 6 overstatements. No number behind a scientific claim was wrong."

**E8. Line 10: misdescribes the tests.** "our predictor ranked the correct kinase better than the three comparison methods". Only M2 ranks kinases. M1 and M3 rank sites for a kinase.
→ "**On that test, our predictor beat the three comparison methods on all three main metrics.**"

**E9. Line 285: "with no recorded reason".** The M3 decision note (`RESEARCH_TREE.md` 00:05) gives only "Nothing else beat the noise floor". The retrospective reads that as the de facto reason ("CV had stopped improving").
→ "FREEZE was committed 5.3 h before its planned gate. The only reason recorded is that nothing else had beaten the noise floor on CV."

**E10. Line 340: wrong scope.** "beat the Cantley atlas and our retrains of two learned predictors on Ser/Thr kinases". The nine pre-registered tests used all kinases. The Ser/Thr result is a stratum with no confidence interval.
→ "…a fused model beat the Cantley atlas and our retrains of two learned predictors on all nine pre-registered tests. The advantage rests on Ser/Thr kinases (TEST M1 0.841); Tyr kinases showed no specificity."

**E11. Line 217: incomplete priming list.** The pre-fixed list also includes PLK1 (line 254 of the paper draft).
→ "GSK3, CK1, CK2α and PLK1 prefer sites near an already phosphorylated residue…"

**E12. Line 233: mixed library sizes.** The +0.176 comparison is top-100 against Cantley top-100, while the preceding numbers are top-50 (`REPORT.md` §8 table).
→ "…It was much better than a library built from the Cantley atlas (top 100 per kinase for both: +0.176 [0.091, 0.263] on perturbations)."

**E13. Line 149: missing harness qualifier.** "+0.019 [0.010, 0.033]" was measured on the pre-A24 harness (line 268 of the paper draft).
→ "Distillation raised its M1 by +0.019 [0.010, 0.033] (measured before the A24 metric fix)."

**E14. Line 126: missing context for 0.796.** The released Phosformer weights could score only 2% of required pairs (LEADERBOARD row 1, coverage 0.020).
→ "The released Phosformer weights scored 0.796 on the 2% of pairs they could score, but…"

**E15. Fig 2 and caption (line 83).** The right panel seals 2 of 5 "families" (40%), but the caption says about 20%. The panels also label homology clusters as "family", which clashes with "kinase family" in T2/T3 (lines 92–95).
→ Redraw with 1 of 5 clusters sealed (or add "proportions not to scale"), and relabel rows "protein cluster A…E". Caption: "Each grey bar is a protein with five sites; two homologous proteins form a cluster."

**E16. Line 182: self-contradictory p-value wording.** "Eight had a bootstrap p-value below 0.0002, which is the smallest value that 10,000 resamples can show."
→ "Eight reached the smallest p-value that 10,000 resamples can show (0.0002: no resample favoured the baseline), so their true p-values are smaller. After Holm correction they are 0.0018."

**E17 (source discrepancy, for the PI).** Line 258 "about 80 of their 84 requested changes were accepted" matches `RESEARCH_TREE.md:47`, but `final_retro.md` §1.1/§3.4 says "84 … accepted". Keep the research-tree value and note the discrepancy.

---

## 2. Overstatements and audit constraints

**O1. Line 13: "Independent agent audits caught real errors before they reached the results."** The A24 defect had already affected 47 scored runs, and one TEST claim (priming) had to be qualified after TEST.
→ "Independent agent audits caught real errors, most of them before the test set was opened."

**O2. Line 340: causal claim about controls.** "both retrained learned methods fell below the label-free atlas once the controls were in place". The manuscript (line 193 of the draft) says the controls were not ablated and that limited retraining budgets may contribute.
→ "…both retrained learned methods fell below the label-free atlas under our controls. We did not remove the controls one at a time, and our limited retraining budgets may also contribute."

**O3. Line 67: "every precision and probability that we report is a lower bound on the true value".** This holds only if every curated pair is a true direct pair, and line 61 says some are indirect.
→ "Cost: if curated pairs are correct, the precisions we report underestimate true precision, because some unlabelled pairs are real. Some curated pairs are indirect, though, so even this is not guaranteed. We cannot measure the false-positive rate on unlabelled sites."

**O4. Line 225: "well calibrated".** The audit says honest for S/T, and Tyr has ECE 0.069 and mild overconfidence.
→ "On 1,458 labelled TEST sites, the assignment probability was well calibrated for S/T sites (overall expected calibration error 0.026; see glossary)…"

**O5. Line 209: "explain part of the atlas's Tyr gap".** This is a CV diagnosis that was linked to the gap, not tested as a cause.
→ "On CV, promiscuous tyrosine 'hub' sites and weak atlas matrices were linked to the atlas's Tyr gap…"

**O6. Line 126: "We traced SELPHI2.0's published 0.97 mainly to random-pair negatives".** One contrast supports this; "traced" is too strong.
→ "SELPHI2.0's published 0.97 probably reflects random-pair negatives: on our cluster split, the same 'positives versus random pairs' measure still gave 0.947."

**O7. Line 274: "Independent checks were the strongest part of the process."** This is the retrospective's judgement.
→ "The retrospective judged the independent checks to be the best part of the process. They caught:"

**O8. In brief, bullet 2 (line 10): two required caveats are missing.** The bullet quotes the margin over Cantley but omits the weak-draw caveat and the M2 composition caveat (`audit_test.md` §0 #7, §7 items 1 and 4).
→ Append: "Part of the margin over the atlas reflects a weak test draw for it (the cross-validation margin was +0.063), and the test rise in the site-level score M2 reflects which kinases were tested, not an improvement."

**Caveat checklist (body text).**

| Caveat | Where | Status |
|---|---|---|
| Tyr not trusted | lines 11, 209, 342 | present |
| No ≥0.8 tier | lines 12, 225, 342 | present |
| Cantley margin partly a weak draw | lines 192, 342 | present; missing from In brief (O8) |
| M2 rise is composition | lines 191, 342 | present; missing from In brief (O8) |
| Baselines are our retrains | lines 10, 193, 342 | present |
| Single TEST spent early | lines 13, 285, 342 | present |
| Post-TEST exploratory, model unchanged | lines 285, 342 | present |
| assignment_prob rank-1 caveat | line 239 | present and accurate (25,532/350,736 = 7.3%; 19 Y sites at 0.55) |
| PSP label bias | lines 49, 61, 342 | present |
| One split, one run | line 342 | present |
| Agent-process failures | §8 | present |
| Strata without CIs | lines 213, 221, 342 | present |

---

## 3. Accessibility

**Jargon used before it is defined, or never defined:**

| Line | Term | Fix |
|---|---|---|
| 10 | "pre-registered" | glossary: "**Pre-registered:** fixed in writing before the data that test it were seen." |
| 12 | "assignment probability" (defined only at line 159) | add "(a calibrated probability that the top-ranked kinase is right; Section 5)" |
| 47 | "leak/leakage" (never defined) | glossary: "**Leakage:** information about the test answers reaching the model during training, for example through the same protein appearing on both sides of a split." |
| 99 | "B0-pop", "B0-prior" | "B0-prior (a control that ranks kinases only by their number of known substrates)"; drop "B0-pop" from the text, since "popularity floor" is enough |
| 155 | "standardised", "nested check", "stacking" | "…standardised scores (each rescaled to mean 0, spread 1)… A check that re-chose the weight inside each fold picked the same weight every time. The alternative, feeding R3's output into R2 as a feature ('stacking'), would have needed 36 neural trainings instead of 6." |
| 159 | "per-kinase scaling (Platt)", "softmax", "isotonic map" | "…first put each kinase's scores on a common scale, then shared the probability across all candidate kinases so that it sums to one, then corrected the top value with a monotone map fitted on labelled sites, separately for S/T and Tyr…" |
| 182 | "matched kinases" | "on the 116 kinases that both methods score" |
| 205 | "Cantley's standard-tier value (0.714)" | "Cantley's overall TEST value (0.714; the atlas uses no labels, so holding out families does not change it)" |
| 225 | "expected calibration error" | glossary: "**Expected calibration error (ECE):** the average gap between claimed and observed precision across probability bins, weighted by the number of sites; 0 is perfect." |
| Fig 4 label; line 95 | "Manning group" | "kinase group (the standard kinase family tree of Manning et al. 2002)" |
| 108 | K for CV not given | add "On CV, the 203 kinases with at least 5 development substrates on at least 2 clusters." |

**Glossary corrections:**
- Line 32, CI: "the range that would contain the true value in 95% of repeats" is imprecise. → "**95% interval (CI):** a range built so that, over many repeats of the study, 95% of such ranges would contain the true value. Here it comes from a **bootstrap**: protein clusters are resampled many times, with the set of kinases held fixed."
- Line 31, calibration: → "adjusting scores so that '0.6' means right about 60% of the time, measured on sites that have a curated kinase."
- Line 22, hard negative: "a site that is phosphorylated by some *other*, unrelated kinase" asserts a fact the labels do not establish. → "a site with a curated label to another kinase outside the kinase's paralog group."
- Line 33, Holm: → "a rule that adjusts p-values when several tests are run together, so that the chance of any false 'win' among them stays below 5%."
- Line 35: see E5.
- Correct as written: AUROC (line 24), macro average (25), PU (20), popularity floor (26; it ranks by literature and MS observation counts, `DESIGN.md:139`), canary (27), distillation (29), late fusion (30).

**Metric explanations.** M1 (line 105; the toy example is correct: 28 of 32 pairs = 0.875), M2 with paralog filtering (106, matches DESIGN A18), M3 (107; the ~60,000 candidates match the per-kinase E-all counts) and M4 (108) are correct. Two fixes are needed: the p-value sentence (E16) and the M2 candidate count in Fig 3 (§5).

---

## 4. Gaps against the brief

**G1. Length.** The main text before §9 is about 6,800 words; the brief asks for 4,500–5,500. Suggested cuts:
- §5 "What information was used" → one sentence per kind of evidence.
- §6 "Popularity and other strata" → merge into the per-kinase section.
- §7 → 6 items.
- §8 "What we would do differently" → 4 items.

**G2. Fusion weights for M2.** The site-level score uses equal weights (score_site = 0.5·R2 + 0.5·R3; `audit_test.md` §1), and the A31a sentence (line 286) cannot be followed without this. Add to line 155: "For ranking kinases at a site (M2), the two scores are weighted equally."

**G3. Human checking.** The files say that no human check of the code or analyses is recorded (lines 11, 207, 314 of the paper draft). Add to "What is not known" (line 262): "Whether the human checked any code or analysis is not recorded."

**G4. Auditor access to sealed labels.** Not mentioned; see E4.

**G5. Credit table (lines 307–336).** Most attributions match LITERATURE.md and the manuscript. Corrections:
- Phosformer-ST (313): "Idea behind R3 distillation" inflates the credit. The external-PI reviewer proposed the idea, and it *follows* Phosformer-ST (line 268 of the paper draft). → "R3's distillation follows its training on atlas percentiles; our split design drew on its kinase and peptide grouping."
- DeepKinZero/InfoNCE (330): the manuscript says R3's head is "similar to" DeepKinZero; InfoNCE was used directly. → "R3's scoring head resembles its zero-shot design; InfoNCE term in R3's loss."
- Sugiyama 2019 (328): → "Planned as a feature but never built; not used as labels, because it hurt activity inference in Müller-Dott et al. 2025."
- CoPheeMap (325): add "its network outputs were excluded as PSP-derived; only its abundance-adjustment approach was followed."
- SIGNOR (335): DATAINEED #3 says the SIGNOR host was unreachable and OmniPath's SIGNOR-derived rows were used. → "SIGNOR (via OmniPath)". The PI should confirm against line 221 of the paper draft.

**G6. Agent-science section vs sources.** The timeline matches. Checked against LEADERBOARD.tsv, file mtimes and the retro:
- first CV run 14:57; last 23:17;
- 62 of 135 runs started at or after 20:00;
- 11.8 h of serial CV time;
- FREEZE mtime 00:04:07; ledger 00:10:01–00:57:29;
- last pre-gap writes 03:15–03:17; retro 09:27–09:34; report review 09:29–09:46.

The counts also match: 34 amendments, 12 in the last 8.6 h, 6 in 40 min, A24 at T+3.7 h after 47 runs, 11 `[~]` markers with 8 stale, and 2 human-resolved versus 10 open data requests (DATAINEED.md). The OOM wording "at least two" is right: the retro names 2, and `docs/regimes_labels.md` §6 records further kills between 20:14 and 21:05. Line 266: "62 started after 20:00" → "at or after 20:00".

**G7. HTML alt text.** Every `<img>` has `alt="Figure N"`. Use the bold caption title as alt text, for example `alt="The study in six steps: labels, fair test, evidence, two models, calibration, outputs"`.

---

## 5. Figures and HTML (minor unless noted)

**HTML.**
- Eight `<img>` tags, each a base64 PNG. Their decoded hashes match `fig1, fig2, fig3, fig6, fig4, fig5, fig8, fig7` in the order of Figures 1–8, which is correct for the renumbering.
- An HTML parser found no unclosed or mismatched elements.
- No `<script>`; the only external URL is the kinase-library GitHub link in the references.
- No personal data (no absolute paths, user, host or session strings).
- Light and dark tokens are present, images sit on a white plate, and the layout uses 16 px gutters with a 46 rem measure.

**Figures.**
- **Fig 1:** readable and accurate. Box 5 "average the two standardised scores" could add "(weighted)".
- **Fig 2:** see E15.
- **Fig 3:**
  - M2 header "top 5 of ~500 kinases": candidates are kinases of the matching residue class, about 346 per positive on TEST (`audit_test.md` §2). → "top 5 of the few hundred candidate kinases".
  - The M3 top-1% band takes about 16% of the bar. → Add "(not to scale)".
- **Fig 4 (fig4_results):** every plotted value matches the `*_ci.tsv` files and the LEADERBOARD official rows:
  - Cantley TEST 0.714 / 0.189 / 0.161; CV 0.717 / 0.171;
  - Phosformer CV 0.677 / 0.132 / 0.104;
  - floor CV 0.552 / 0.020 / 0.142.
  - The caption is correct: Cantley's M3 is blank on CV in the official arm.
- **Fig 5:**
  - Panel B values match the strata files (S/T 0.841, 0.772, 0.736, 0.751, 0.576; Tyr 0.595, 0.496, 0.575, 0.534, 0.585). The panel A counts match my recomputation (3 hollow S/T, 13 hollow Tyr, 2 dual-specificity filled).
  - Colour problem: orange means "Tyr kinase" in A and "Cantley" in B, and blue means "Ser/Thr" in A and "our model" in B, so a reader may misread B. → Use a different hue pair in A (for example dark grey circles and purple triangles), or state the switch in the caption.
  - The x = 1 tick stands for kinases with 0 labels. → Label it "0".
- **Fig 6 (regimes):** all values match `docs/regimes_labels.md` §1.
  - Caption gaps: split-model M2 uses a rank-preserving map onto the global score scale (qmap), and partial pooling equals the global model on M1 but not on M2 (0.199 against 0.192).
  - → "Partial pooling (blue) chose a setting that matches the global model on M1. M2 for split models uses their scores mapped onto the global model's scale."
  - Write "in-vivo-only labels … a set three times smaller than the in-vitro set" (relative to L2 it is four times smaller).
- **Fig 7 (calibration):** points match `final_reliability.tsv`. See E3.
- **Fig 8 (timeline):** event times match sources: audit bars at about 17:00, 21:20–22:04, 01:00–02:00 and 02:34–03:16; problem markers at 16:55, 17:50, 19:55, 21:11, 22:12 and about 23:00.
  - Fix the caption (E1).
  - The short unlabelled green bar at about T+7.3 h belongs to the "test-code review" label below it; place the label next to the bar. The review ran from 20:25 to 23:30, so the bar is too short.

**Other minor items:**
- Line 246: "M1 0.725 against 0.751" → "M1 0.725 against a required 0.751".
- Line 41: give the 95% source in NUMBER_MAP (it is line 21 of the paper draft).
- Line 192: "On kinases scored in both sets" → "On the 115 kinases scored in both sets".

---

## 6. Numbers recomputed or traced (selection)

**TEST primary comparisons** (`superiority.tsv`; all `beats=True`; Holm 9×0.0002 = 0.0018):

| Comparison | Δ [95% CI] | p |
|---|---|---|
| vs B1_i, M1 | +0.0745 [0.0465, 0.0991] | 0.0002 |
| vs B1_i, M3 | +0.2051 [0.1650, 0.2346] | 0.0002 |
| vs B1_ii, M2 | +0.1787 [0.1469, 0.2130] | 0.0002 |
| vs B3, M1 | +0.0861 [0.0594, 0.1051] | 0.0002 |
| vs B3, M2 | +0.1993 [0.1507, 0.2290] | 0.0002 |
| vs B3, M3 | +0.2592 [0.2131, 0.2908] | 0.0002 |
| vs B4, M1 | +0.0894 [0.0630, 0.1172] | 0.0002 |
| vs B4, M2 | +0.1300 [0.0908, 0.1565] | 0.0002 |
| vs B4, M3 | +0.0588 [0.0096, 0.1198] | 0.0178 |

**TEST standalone and strata:**
- Final 0.794 [0.774, 0.811] / 0.363 [0.324, 0.397] / 0.363 [0.309, 0.401]; T3 0.752; 25,922,452 pairs; 12 configs.
- 119 kinases in M1 and 256 in M4 (recomputed from `in_K_primary` / `in_K_any`).
- Recomputed from per-kinase tables:
  - S/T 0.8413 (94 kinases), Tyr 0.5947 (23), Tyr floor 0.5847;
  - kinases beating their own floor: 91/94 S/T, 10/23 Tyr;
  - B4 M3: 52 better / 27 worse / 40 tied;
  - 16 of 23 Tyr kinases between 0.4 and 0.7.
- 97/116 is an audit diagnostic; the standalone count is 100/119.
- Tyr: Cantley 0.496, B3 0.534, B4 0.575.
- Hub 0.863 / floor 0.788. New since 2021: 0.766 − 0.710 = 0.056.
- Priming: primed 0.897 / 0.704; unprimed 0.740 / 0.725.
- By training bin: 5–19 0.759; ≥100 0.795; M4 bin 0 0.807 (16 kinases).

**Calibration** (`final_reliability.tsv`):
- ≥0.5 tier: 236 + 37 + 72 = 345 sites; claimed (weighted) 0.6155; observed 208/345 = 0.603.
- ECE 0.0257; top-1 accuracy 0.3532; n = 1,458.

**Coverage:** 29,886/350,736 = 8.52%; 26,762; Cantley-only 30,258; random 6; 25,532/350,736 = 7.28%.

**CV:**
- Final 0.786 / 0.263 / 0.361; R2 0.776; R3 0.754.
- Δ vs Cantley +0.063 / +0.105 / +0.227.
- Same-kinase comparisons: 0.792 → 0.793 (M1), 0.348 (M2), Tyr 0.667 → 0.595.
- Regimes: global 0.768; per-kinase 0.756 / 0.755; per-group 0.761; 78 data-rich kinases 0.805 vs 0.779; human-only M2 −0.029; L2ext 3-seed mean (+0.0058 − 0.0040 + 0.0022)/3 = +0.0013; in-vitro/in-vivo pair counts 8,987 / 2,808 = 3.2×.
- Ablation: M3 0.125 → 0.236; descriptors 0.759 → 0.773 = +0.014.

**Timeline** (LEADERBOARD.tsv `started` column, file mtimes, ledger):
- 135 runs; 62 started at or after 20:00; 11.82 h total runtime.
- FREEZE 00:04 = T+10.73 h, 5.27 h before 05:20.
- TEST 00:10:01–00:57:29 (T+11.6 h).

**Toy examples:** Fig 3 AUROC 28/32 = 0.875; M3 600/60,000 = 1%, 3/5 = 0.6.

**References:** all 52 entries in the briefing's list occur verbatim in the accompanying paper or its reference list; Platt keeps its "[verify]".
