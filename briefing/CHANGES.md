# Changes after the independent review (`review_briefing.md`), 2026-10-28

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Applied to `briefing.md` and `figures/make_figures.py`. All figures were re-generated, and `briefing.html` was rebuilt; the changed PNGs (Fig 1, 2, 3, 6/regimes, 5/per-kinase, 7/timeline) were re-read. `NUMBER_MAP.tsv` was updated: rows for removed claims were dropped, and rows for new numbers were added. No git was run.

## 1. Errors (E1–E17): all applied after checking the source

- **E1, E2 (gate timing).** FREEZE was committed at 00:04 (T+10.7 h), 5.3 h before the 05:20 gate (T+16 h). The TEST run ended at 00:57 (T+11.6 h), 4.4 h before the gate. Both facts are now stated the same way in "In brief", §8 Timeline, the Fig 8 caption and "What went wrong". Sources: `final_retro.md` §1.2, `results/test/LEDGER.md`.
- **E3.** The Fig 7 caption no longer claims a lower-bound relation to unlabelled sites (`audit_test.md` §6).
- **E4.** The access statement now says that a software test checked the access and that the auditor read the labels read-only. This was verified in `audit_test.md` ("What I did").
- **E5.** The glossary and §2 now mention the logged override, which the manuscript Methods describe.
- **E6.** "Overruled" was replaced. Both verdicts were CONFIRMED-WITH-CAVEATS (`audit_posttest.md`, Verdicts).
- **E7.** The text now says "7 factual errors and 6 overstatements" (`review_report_v3.md` §A, §B). The reviewer's list of error types was not copied, because I did not check each one.
- **E8.** The "In brief" heading now says the predictor beat the three comparison methods on all three main metrics.
- **E9.** The text now says that the only reason recorded for the freeze is that nothing else had beaten the noise floor (`RESEARCH_TREE.md` M3).
- **E10.** The "What the study shows" scope is now all nine tests, with the advantage resting on S/T kinases.
- **E11.** PLK1 was added to the priming list (manuscript Methods).
- **E12.** The Cantley comparison now states top 100 per kinase for both libraries (`REPORT.md` §8 table).
- **E13.** The distillation gain now carries the pre-A24 qualifier.
- **E14.** The 2% coverage of the released Phosformer weights was added.
- **E15.** Fig 2 now seals 1 of 5 clusters, and its rows are labelled "cluster", with matching titles and notes. The caption was updated to match.
- **E16.** The p-value wording was rewritten, and the Holm value 0.0018 was added.
- **E17.** "About 80 of 84" (`RESEARCH_TREE.md`) was kept. The retrospective's "84 accepted" is recorded as a discrepancy in `NUMBER_MAP.tsv`; **the PI should decide**.

## 2. Overstatements (O1–O8): all applied

- **O8.** The weak-draw caveat for Cantley and the M2-composition caveat are now in "In brief".
- **O3.** The PU "lower bound" claim was replaced with the conditional wording. In §10, "all precisions are lower bounds" became "precision on unlabelled sites is unknown".
- **O1, O2, O4, O5, O6, O7.** Applied with the reviewer's wording, lightly shortened.

## 3. Definitions

- **New glossary entries:** pre-registered, leakage, kinase group, expected calibration error.
- **Corrected glossary entries:** hard negative, calibration, 95% CI with bootstrap, Holm, FREEZE/ledger.
- **Defined inline:**
  - assignment probability, at first use in "In brief";
  - B0-prior; "B0-pop" was dropped from the text;
  - standardised, nested check and stacking;
  - the steps of the calibration;
  - "matched kinases" (now "the 116 kinases that both methods score");
  - Cantley's value for unseen families;
  - the CV kinase set (203 kinases).
- **Renamed:** "Manning group" became "kinase group" in the text and in Fig 4 (regimes).

## 4. Gaps

- **G2.** The 0.5/0.5 weighting of the site score (M2) was added to Late fusion and to the A31a item.
- **G3.** §8 "What is not known" now opens with: "No human check of the code or analyses is recorded in the project files." This is supported by the paper draft (Abstract, Limitation 9, The human's role).
- **G5 (credit table).** Each change was checked against the manuscript, LITERATURE.md and DATAINEED.md:
  - Phosformer-ST is now credited as the work that R3 follows, not as the source of the idea;
  - DeepKinZero is credited as "resembles";
  - Sugiyama is now described as planned but never built;
  - CoPheeMap's network outputs are marked as excluded;
  - SIGNOR is now "via OmniPath" (DATAINEED #3).
- **G6.** "62 started at or after 20:00".
- **G7.** Every image now has descriptive alt text, set in the Markdown and passed through to the HTML.

## 5. Figures

- **Fig 1:** box 5 now says "weighted average".
- **Fig 3:**
  - the M2 title now says "the few hundred candidate kinases" (345.7 per positive; `audit_test.md` §2);
  - the M3 bar is marked "not to scale".
- **Fig 5 (per kinase):**
  - panel A uses dark grey for Ser/Thr and violet for Tyr, so its colours no longer clash with the method colours in panel B;
  - the caption says that panel B colours denote methods;
  - the first x tick is labelled "≤1", because kinases with 0 or 1 label share that position;
  - panel B shows three decimals, which avoids a 0.59 against 0.58 rounding mismatch with the text.
- **Fig 4 (regimes):**
  - the labels now read "kinase group" and "Partial pooling (selected setting)";
  - the caption now explains that the split models' M2 uses scores mapped onto the global model's scale;
  - the caption now says that partial pooling matches the global model on M1 only;
  - the caption now says "three times smaller than the in-vitro set".
- **Fig 8 (timeline):**
  - the review label now gives 20:25, the finding and "fixes applied 23:30";
  - a dotted line runs from the review to a tick at 23:30;
  - the caption was fixed (E1).

## 6. Length

The main text before §9 is now 5,300 words of prose, excluding headings, the glossary, figure captions and the "In brief" box; §10 adds about 340. Before the review it was about 5,800 words, counted the same way.

Repetition was cut:
- §1 and the §5 evidence list were condensed;
- popularity and strata were merged into the per-kinase section;
- §7 went from 8 items to 6; the similarity-feature and intercept items were dropped, and both remain in the supplement;
- "What we would do differently" went from 6 items to 4;
- the timeline paragraph and §10 were shortened.

No caveat was removed.

## 7. Not applied, or applied differently

- **Fig 8, length of the test-code review bar.** The reviewer says that the review ran from 20:25 to 23:30. `RESEARCH_TREE.md` 5.12 records the review launch at 20:25, its findings at 20:45, and the fixes applied by 23:30. I kept the bar at 20:25–20:45 and added a dotted extension to a 23:30 tick, instead of drawing one review bar to 23:30. The reason is that CV audit 2 (21:24–22:15) sits in the same lane during that time.
- **O-numbers quoted as reviewer wording.** Some replacements were shortened to keep the length target. The meaning is unchanged.
- **E7 detail.** The reviewer's list of error types ("stale sentences, broken row references, one mislabelled count column") was not copied, because I did not check each one.
