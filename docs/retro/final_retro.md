# Final retro — outside view (fresh agent)

**Author:** a fresh agent with no prior involvement. I read `PROMPT.md`, `RESEARCH_TREE.md`, `LESSONS.md`,
`DESIGN.md` §0/§15–17, `docs/retro/M{1,2,3}_retro.md`, `docs/DELEGATION.md`, `REPORT.md` §1/§6.8/§7/§11,
`docs/audit/audit_test.md` §0, `docs/audit/audit_posttest.md`, and `results/cv/LEADERBOARD.md`
(137 rows, 135 timed runs). I ran no git, read nothing under `sealed/`, and modified no existing file.

**Clock.** T0 = 2026-10-26 13:20 CET, deadline 2026-10-27 13:20. The last artefact written in the repo is
`REPORT.md` / `RESEARCH_TREE.md` at **03:17** (T+13.95 h). `date` at the time of writing reads
**2026-10-27 09:32 CET** (T+20.2 h), so **about 6h 15min of wall clock have passed with no file written**
and 3 h 48 min remain. That gap is itself the largest single finding of this retro and is treated as such in §1.

This is a retro on **how the work was run**. The science is not re-litigated; where I cite a number it is to
date a decision or size a cost.

---

## 1. Where the time actually went

### 1.1 Measured phase budget

| Window | Elapsed | What happened | What it changed in the delivered model |
|---|---|---|---|
| 13:20–14:33 | 1.2 h | 4 literature agents, data profiling, 2 acquisition waves (12.6 GB), DESIGN v0→v2, **11 reviews** (7 in r1 at 13:49–14:14, 4 in r2 at 14:26–14:33), 84 + 22 accepted changes | Everything downstream. Cheap and high-yield |
| 14:33–16:27 | 1.9 h | 8-agent implementation fan-out; splits frozen (`results/splits` 14:51); first CV runs 14:46; **G2 met 16:27, 4.4 h early** | The pipeline |
| 16:27–19:19 | 2.9 h | R2 ablation curve, R3 Cantley distillation, late fusion | **The final model.** `r2_fusion_r3__20261026-191901__s0`, M1 0.786 |
| 19:19–00:04 | 4.8 h | R5, knn_sim, F9, regimes/labels, A32a/b, B3 retrain, B4, strata, downstream, G4 hardening, audits r1/r2, deliverable pinning, FREEZE | **Nothing.** No CV run after 19:19 exceeded M1 0.786 |
| 00:04–00:57 | 0.9 h | The single sealed TEST batch (2,838 s, 12 configs, B = 10,000) | The headline numbers |
| 00:57–03:17 | 2.3 h | `audit_test`, G5 deployment + 131 M U_obs scores, 3 post-TEST experiments, `audit_posttest`, REPORT v2→v3 | **Nothing** (A34 forbids it by construction) |
| 03:17–09:32+ | **6.3 h and counting** | No artefact written | Nothing |

### 1.2 The three numbers that matter

1. **The delivered model was finished at T+6.0 h.** The frozen final config is the 19:19 fusion run.
   The 18:37 run was already at 0.785. Of 135 timed CV runs, **62 (46 %) started at or after 20:00**, and the
   best of those is 0.774 (`regime_L2ext_refit_s0`, 21:47). The remaining 18 hours of the budget produced
   audits, baselines, documentation and negative results — all legitimate, none of it moved the model.
2. **Compute was never the constraint.** The 135 CV runs total **11.8 h of serial wall time** on a 64-core /
   4×4090 node. The whole modelling leaderboard would fit in one evening on one machine. The binding
   constraints were agent coordination, review latency and the PI's own decision throughput.
3. **The single sealed TEST was spent at T+11.6 h of 24** (COMPLETED 00:57:29), against a planned gate of
   **G4 = T+16 h (05:20)** and the PI's own A28 rule that *"exploration on CV closes at 01:00"*. FREEZE was
   committed at 00:04 — **5.3 h before its gate and 56 min before the PI's own exploration deadline.**
   PROMPT.md says "the goal is the best possible outcome, not the fastest completion" and "keep going until
   you genuinely run out of promising ideas or close to the time limit". The project did the opposite: it
   closed the only irreversible door at the halfway mark and then spent 12.4 h unable to change anything.

### 1.3 Decisions that cost the most for the least return

Ranked by (hours consumed) × (probability the hours could have bought something):

1. **Freezing and spending TEST 5.3 h early (≈ 12 h of budget made unusable).** The direct, visible cost is
   A34: three experiments on the project's own named weaknesses (understudied-kinase M2, primed sites, Tyr)
   ran at **02:18–02:23**, six hours after they could have been adopted. All three failed, so no accuracy was
   lost — but that is luck, not design. The option value was thrown away deliberately. No second sealed shard
   was reserved at split time, although 11,609 clusters would have supported DEV / TEST-A / TEST-B easily.
2. **The regime + label-variant sweep: 55 of 135 runs, 4.53 h of 11.8 h serial CPU (38 %), 121 log files,
   zero adoptions.** This was a *PROMPT-mandated* question ("think intelligently about the training regime …
   test different things"), yet it was delegated at **19:55 (T+6.6 h)** and closed at **00:03**, producing
   amendment A33 **one minute before FREEZE**. It needed no new features — only the existing R2 pipeline — so
   it could have run in the 14:33 fan-out. Late scheduling converted a cheap, required experiment into the
   critical path of the freeze.
3. **The A24 pooled-OOF defect, found at 17:01 by an audit rather than by a control run.** The `B0-prior`
   control read M1 = 0.40 instead of 0.50. **47 of 135 leaderboard runs had already been scored** with the
   broken metric; every doc written before 17:01 became stale (M2 retro's own words), a re-evaluation script
   had to be written (`scripts/reeval_cv_runs.py`, 18:54) and re-run at 01:23, and the REPORT now has to carry
   "pre-A24 / A24" dual values throughout. The control that caught it existed from G1 and could have been run
   against the harness before any model.
4. **Seed and knn sweeps past the noise floor.** `r2_final*` and `r2_final_knn*` account for **74 of 404 run
   directories**; knn_sim was not adopted (+0.003 M1, −0.018 M2) and F9 was not adopted (+0.002). The M2 retro
   flagged exactly this at 18:00 ("chasing marginal DEV gains"); the sweeps continued for another 4 hours.
5. **R5 (ESM cross-encoder).** Correctly gated and correctly dropped at 19:53 — but DESIGN §10's own compute
   estimate (3.6–14 GPU-h/fold) already said the real thing was infeasible, and the ColBERT fallback was
   invented after the agent had started. The gate saved this from being a disaster; the scoping should have
   caught it earlier.
6. **B3 Phosformer's 4.3 GPU-h MLM pretraining delivered at 23:30, 34 min before FREEZE.** PROMPT requires a
   fair retrain, so the work was mandatory and the fidelity check (0.888 vs 0.885) was worth having. But the
   highest-latency mandatory deliverable slipped twice (ETA 19:15 → 22:30 → 23:30) and came within half an
   hour of triggering A29's "NOT TESTED, p = 1" contingency on a primary test.

---

## 2. Which of the 14 LESSONS changed behaviour, and which were written then ignored

### Changed behaviour (evidence of later effect)

| # | Rule | Evidence it bit |
|---|---|---|
| 1 | Reviewers get read-only data + a measured number | r1 produced the three findings the whole design rests on (68 % co-labelled → cluster cross-fitting; >50 % paralog hard negatives; STRING 34 % vs 4 %). r2 then found the coverage metric was *vacuous* (A17) and iKiP at 61 % overlap (A11) |
| 3 | Label-derived / "unbiased" sources are leaky until cross-fitted | Produced A1 (single `crossfit.py` driver with `LeakageError`), A11 forbidden-source list, `tests/test_feature_overlap.py`, and the F6-clean alarm run (0.773 = full → no circularity) |
| 4 | Null controls must match positives on every label-free axis | Canary v1 (0.76 from site density) was replaced by v2 **before** it was used to gate anything. Directly reused post-TEST: `audit_posttest` killed finding (b) with a label-free priming flag |
| 6 | Pooled-OOF needs fold-invariant scores | A24 applied to all methods; whole leaderboard regenerated as `metrics_a24.json` |
| 8 | Strongest comparator; label contaminated references | B3-pub is tagged `CONTAMINATED (trained on PSP; contrast only)` in the leaderboard; A23 explicitly picks the *stronger* Cantley arm per endpoint. **But see §3.3 — it failed on the priming arm** |
| 9 | An independent auditor recomputes headline numbers | Four audits (r1, r2, test, post-test). They caught the A24 artefact, stale fold-0/3 stores, the overconfident ≥0.8 tier, a bypassable single-batch guard, and the post-hoc "balanced" variant in A34(a). This is the single best thing the project did |
| M1 #10 (dropped later) | ≤4 reviewers with distinct remits | r1 used 7, r2 used 4. The rule worked and was then deleted from LESSONS |

### Written and then ignored

| # | Rule | Where it was broken afterwards |
|---|---|---|
| 12 | "status per event; delete status columns that aren't maintained" | Flagged in **all three** retros (M1 §1.3, M2 §1.5, M3 Stop #5) and still broken now: `docs/DELEGATION.md` (last touched 00:05) says `core-infra (G4) … deploy loader running (00:05)` 9 hours later; `RESEARCH_TREE.md` carries **11 `[~]` markers of which 8 are demonstrably done** (0 Setup, 1 Literature, 2 Profiling, 3 DESIGN, 4.8/5.4 B3, 5.2 R2, 5.13 regimes). Three retros told the PI to either update or delete the column. Neither happened |
| 7 | "Write the config log, selection rule and adoption threshold *before* running" | Honoured for A32 (rules at 22:16, experiments at 22:46) and the R5 gate. Then broken at 02:23 in exactly the way it warns about: `audit_posttest` found A34(a)'s **kinase-balanced variant was added after seeing the plain variant's first grid point**, and the pre-declared variant is *worse* on its own target stratum (−0.013 [−0.021, −0.004]). Eight hours after the rule, with the rule in the file |
| 5 / 14 | "Freeze the harness"; "After FREEZE, no protocol amendments" | A28 froze `eval/labels/splits/crossfit/canary` at 17:35, but `eval_test.py` — the code that touches the sealed set — kept changing until 23:30. LESSON 14 was written at **00:05, one minute after FREEZE.json was committed**, and A34 (a 35th amendment) landed at **01:31**. In substance the model was not reselected, which is the thing that mattered; in form the rule was written after the behaviour it was meant to govern, every time |
| 13 | "Audit PROMPT.md asks at every milestone" + "≤60 GB per agent" | The rule exists *because* the regimes ask was discovered at 20:00. The memory half was added after the 17:50 OOM and there was a **second OOM at ~23:00**. A rule that is created by a failure and then immediately suffers the same failure has not changed behaviour |
| 11 | "Fan out implementers only against frozen interfaces" | r2 review (finishing 14:33) ran in parallel with 6 implementers started ~14:00; A1 (one cross-fit driver) and A2 (`store.py` with index hashes) are interface rewrites issued mid-implementation. M1 predicted the rework at 14:15 and the project did it anyway |
| 2 | "Key by UniProt accession" | Not a process rule at all — a project fact. It never governed a later decision; it belongs in DESIGN, not LESSONS |

**Pattern.** The rules that changed behaviour are the ones that were *compiled into an artefact*: a
`LeakageError` in `crossfit.py`, a canary gate in the eval harness, a FREEZE.json SHA pin, a standing auditor.
The rules that were ignored are the ones that stayed prose in a markdown file and depended on the PI
remembering them under time pressure. Every single one of the four "written then ignored" rules is of the
second kind. This is the most transferable finding in this document.

---

## 3. What the PI got wrong (outside view)

### 3.1 Optimising for finishing instead of for the result

Gates were beaten by 4.4 h (G2), 5.3 h (G4), 7.8 h (G5). Each early completion was recorded as a success in
the tree. None of them was converted into more science; they were converted into *idle budget*. The decisive
error is that FREEZE at 00:04 was not a deadline-driven decision — the deadline was 05:20 — it was a
momentum-driven one. An outside observer reading PROMPT.md ("best possible outcome, not fastest completion")
and then the ledger would ask why the irreversible step happened 5.3 h early with four known, named, unsolved
weaknesses on the board (Tyr, M2 for 5–19-label kinases, primed sites, unseen-drug downstream gap). The
honest answer in the record is that CV had stopped improving — but that is an argument for *changing what you
try*, not for closing the test set.

Compounding this: **no second sealed shard was reserved.** With 11,609 substrate clusters, holding out a
10 % TEST-B at G1 would have cost almost nothing in training data and would have made the entire A34
exploration adoptable. The design anticipated the single-shot problem (A29's `--override-completed`,
STARTED/COMPLETED ledger) but solved it with process discipline rather than with a spare shard.

### 3.2 The amendment rate before the freeze

34 numbered amendments. **A23–A34 (12, i.e. 35 % of all amendments) landed in the last 8.6 h before/around
the freeze**, at ~1.4/h, and **six of them (A23–A28) in the 40 minutes between 16:55 and 17:35**. Every one of
them was issued *after seeing a result*: A23 after the baselines report, A24 after the audit, A26 after the Tyr
doc, A27 after audit r1, A29 after the code review, A30 after the B1 TEST prep, A31 after audit r2, A32 after
the strata doc, A33 after the regimes doc, A34 after TEST.

Most are genuinely conservative (A30 can only help the baseline; A24 is a correctness fix applied to all
methods; A23 strengthens the comparator). That is not the problem. The problem is that **the protocol was
still being written while the results were arriving**, which means no reader — including the PI — can
distinguish "we fixed a bug" from "we adjusted after seeing a number" without trusting the PI's own
annotation. M2's retro said so at 18:00; A29–A34 landed anyway.

The specific case that an outside observer would flag hardest is **A31(a)**: M2, the endpoint with the largest
claimed margin (TEST +0.179), had its *scoring column* changed at 22:12 from the designed calibrated
assignment probability (0.205) to the raw fused score (0.263) — a +0.058 swing, decided after seeing both.
The nested check (0.259) and the CV-only selection are the right mitigations and are documented. It is still
an endpoint definition changed late, in the favourable direction, on the endpoint that carries the headline.
M2's definition moved four times in total (A18, A20, A27, A31a).

### 3.3 Review that happened too late

- **The TEST harness code review was launched at 20:25 and its fixes landed at 23:30** — 34 min before
  FREEZE, 94 min before an irreversible read of the sealed set. It found that the single-batch guard was
  *bypassable* and that a crash would leave no trace. These are exactly the bugs that destroy a sealed test
  set. Code on the irreversible path should be reviewed at G2, not at G4 − 5 h.
- **The priming comparator.** `docs/strata_cv.md` (21:11) and A32 (22:16) both record that the
  **priming-aware Cantley arm beats the final model on primed positives, 0.945 vs 0.893**. FREEZE was
  committed 1h48m later **without that arm**. The consequence is in `audit_test.md` verdict #11: the TEST
  primed-site claim (0.897 vs 0.704) holds only against the non-priming-aware arm, and "TEST cannot say we
  beat priming-aware Cantley". This is a direct, documented, three-hours-of-warning violation of the
  project's own LESSON 8.
- **A30** (10 dual-specificity kinases with only one atlas matrix) is a baseline-fairness bug found at 21:50,
  i.e. during TEST prep rather than when B1 was built at ~15:00.
- **A24** (§1.3 item 3): the metric defect was found by an audit at 17:01, after 47 runs.

### 3.4 Delegation that misfired

- **The regimes agent** (§1.3 item 2): a PROMPT-mandated question, 24 cores, delegated at T+6.6 h, delivering
  38 % of all CV compute and 0 adoptions, with its decision note written 1 minute before the freeze.
- **Review volume was front-loaded and under-filtered.** 11 reviewers produced ~190 KB in 73 minutes;
  84 changes were accepted. The PI then carried §15 (change map) + §16 (A1–A22) + §17 (A23–A34) as a standing
  maintenance burden for the next 13 hours. M1's own retro said the reviewers overlapped (prior_work 407
  lines vs algorithms 429) and proposed a cap; the cap was applied to r2 and then deleted from LESSONS.
- **Two OOM kills** (17:50, ~23:00) from budgeting against the node's 375 GB rather than the 128 GiB cgroup.
  The second happened after the lesson was written.
- **The status table was never owned.** `docs/DELEGATION.md` is PI-owned, has a Status column, and was stale
  at M1, M2, M3 and is stale now. Three retros is enough evidence that this particular PI will not maintain
  it; it should have been deleted at M1 as M1 itself recommended.

### 3.5 What the PI got right (so the above is read fairly)

The standing-auditor pattern, the sealed set + FREEZE SHA pins + STARTED/COMPLETED ledger, pre-registered
adoption gates that were actually honoured when they failed (R5, A32a, A32b, A34 a/b/c — **six negative
results accepted without relitigating them**), the fidelity check that makes B3's 0.679 credible, paired CIs
against the stronger comparator, and a REPORT §1 that leads with its own caveats (Y not trusted, M2 rise is
composition not improvement, margin partly a weak-baseline draw, Coverage@0.8 = 0). Most 24-hour projects do
none of this. The criticisms above are about scheduling and about rules that stayed on paper, not about
integrity — the integrity machinery is the strongest part of the work.

---

## 4. What I would do differently, in priority order

1. **Reserve two sealed shards at G1 and never spend the first before its gate.** DEV / TEST-A / TEST-B from
   the same 11,609-cluster split. Cost: ~10 % of training data. Benefit: the A34 experiments, and anything
   else found in the last third, remain adoptable instead of being exploratory footnotes. Put the gate time
   in the ledger pre-flight so the harness refuses to run early.
2. **Turn PROMPT.md into an owned checklist in the first 30 minutes** and re-audit it at every milestone
   (this *is* LESSON 13; it was written too late to help itself). Line-by-line: training regime comparison,
   per-kinase vs global, label-source variants, unorthodox/priming kinases, coverage, runtime reporting.
   Each gets an owner in the 14:33 fan-out. The regimes sweep then runs at T+2 h on idle cores instead of
   T+6.6 h on the critical path.
3. **Run the full control suite through the real eval harness before the first model run.** B0-prior, B0-pop,
   B0-miss, shuffled scores, label-free canaries, positive-control planted motif. A24 would have surfaced at
   T+1.5 h instead of T+3.7 h, 47 runs would not have needed re-scoring, and no document would carry dual
   pre-/post-A24 numbers.
4. **Fix the comparator set from the strata table, not from the baseline document, and freeze it with the
   model.** Rule: *any* arm that beats us on *any* pre-registered stratum enters FREEZE. The priming-aware
   Cantley arm was known to beat us at 21:11 and was still omitted at 00:04.
5. **Review the irreversible path one gate early.** `eval_test.py`, the seal, the ledger and the pin
   verification get their independent code review at G2 (T+7.5 h), not at G4 − 5 h. Everything else can be
   reviewed late; this cannot.
6. **Batch amendments and make them cost something.** After the harness freeze, amendments land only in
   scheduled windows, each naming (a) the numbers it invalidates, (b) who re-runs them, (c) whether it was
   prompted by a result. Six amendments in 40 minutes (16:55–17:35) is a symptom that the design was being
   finished during execution.
7. **Stop exploring when the search stops paying, and spend the freed time on the named weaknesses.** CV M1
   plateaued at 19:19; the project kept sweeping seeds and label variants until 00:03. The Tyr / understudied-M2
   / priming experiments that eventually ran at 02:18 should have run at 20:00, under A32-style pre-registered
   rules, when adopting them was still legal.
8. **Delete the DELEGATION Status column and keep one append-only event log** written by the PI at each
   delegation and each completion. Three retros proved the column will not be maintained.
9. **Budget memory against the cgroup, not the node, in the agent brief template** — not as a lesson learned
   after the first OOM.

---

## 5. Proposed LESSONS.md for the next project (6 rules)

The current 14 are mostly true and partly inert. The test I applied: *did this rule change a later decision,
and can it be compiled into an artefact rather than remembered?* Rules that fail both are dropped or folded in.

> **LESSONS — process rules (next project)**
>
> 1. **Reserve a second sealed test shard, and never open the first before its scheduled gate.** An
>    irreversible evaluation ends your ability to improve: TEST was spent at T+11.6 h of 24, 5.3 h before its
>    own gate, which made 12.4 h of budget unusable and forced three experiments on known weaknesses
>    (`docs/posttest_*.md`, 02:18) to be exploratory-only. After the freeze: no reselection, no new
>    amendments.
> 2. **Run the null/floor control suite through the real harness before the first model run.** `B0-prior`
>    reading 0.40 instead of 0.50 exposed the pooled-OOF defect (A24) only at T+3.7 h, after 47 runs had been
>    scored and had to be re-evaluated. Controls must match positives on every label-free axis — canary v1
>    scored AUROC 0.76 on site density alone, and a label-free priming flag later reproduced an entire
>    "finding".
> 3. **Treat every label-derived feature and every "unbiased" resource as leaky until it is cross-fitted by
>    cluster and overlap-tested.** 68 % of KS pairs have a same-protein co-labelled site; STRING hits 34 % of
>    PSP pairs vs ~4 % of controls; the iKiP combined GMT contains 61 % of PSP pairs. Enforce it in code (one
>    cross-fit driver that raises, plus an automatic overlap test), not in prose.
> 4. **Before running, write down the selection rule, the adoption threshold, the number of selectable
>    configs, and the comparator set — and put in the comparator set any arm that beats you on any
>    pre-registered stratum.** Honoured, it killed six candidates cleanly (R5, knn, F9, A32a/b, A34c).
>    Broken, it produced the one TEST claim the audit had to retract (priming-aware Cantley, known at 21:11,
>    omitted from FREEZE at 00:04) and the one post-TEST "gain" that exists only in a variant invented after
>    seeing the first grid point.
> 5. **Keep a standing independent auditor who recomputes every headline number from the saved predictions,
>    and review the code on the irreversible path one gate earlier than the rest.** Four audits caught the
>    A24 artefact, stale fold stores, an overconfident ≥0.8 tier, and a bypassable single-batch guard — the
>    last of these only 3.5 h before the single sealed read. Reviewers and auditors get read-only data access
>    and must attach a measured number to every change they request (6 of 7 r1 reviewers did).
> 6. **Extract the brief's explicit asks into an owned checklist at T+0 and re-audit it at every milestone;
>    one owner per path, PI-only git, and delete any status field you will not update.** The mandated
>    training-regime comparison was discovered at T+6.6 h, consumed 38 % of all CV compute, and closed one
>    minute before FREEZE with zero adoptions. The DELEGATION status column was stale at M1, M2, M3 and is
>    stale now — it should have been deleted at M1.

**Dropped from the current 14, with reasons:**

| Current | Disposition |
|---|---|
| 2 (UniProt accession) | **Drop.** A project fact, not a process rule; it governed no later decision. Belongs in `DESIGN.md` §2 |
| 4 (null controls match positives) | **Merge into new 2** |
| 6 (pooled-OOF fold invariance) | **Merge into new 2.** It is one instance of "run the controls through the harness first"; stated generally it keeps its force without being kinase-specific |
| 5 (freeze harness, version numbers) | **Merge into new 1** (freeze discipline) and new 5 (every number cites its harness version) |
| 8 (strongest comparator; label contaminated refs) | **Merge into new 4.** Keeping it separate is what let the priming arm fall through: the rule was in LESSONS but not in the FREEZE checklist |
| 10 (verify artefacts by mtime/hash) | **Merge into new 5.** It is what an auditor does |
| 11 (fan out against frozen interfaces) | **Merge into new 6** (owners + checklist). Kept as a sentence, not a rule: its evidence (M1) is a prediction of rework, not a measurement of it |
| 14 (no post-FREEZE amendments) | **Merge into new 1.** Separating "freeze" from "when to freeze" is what produced a rule written one minute after the freeze it was meant to govern |
| 1, 3, 7, 9, 12, 13 | Survive, in new 5, 3, 4, 5, 6, 6 respectively |

---

## 6. Runtime of this retro

Start 09:27:39 CET, finish 09:34:27 CET, 2026-10-27. Wall time **6 min 48 s**, 15 tool calls, no model runs,
no git, no sealed access, one new file written (`docs/retro/final_retro.md`).
