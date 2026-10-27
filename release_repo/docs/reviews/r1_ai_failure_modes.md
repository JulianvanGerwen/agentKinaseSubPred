> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Review r1: AI-agent failure modes (DESIGN.md v0, §11 focus)

Reviewer role: how autonomous multi-agent research projects fail. Scope: whether §11 is sufficient and enforceable for a ~21 h run with parallel agents sharing one NFS working tree, 64 cores and 4 GPUs.

## Verdict
§11 lists the right *intentions*: toy tests, controls, independent review, a locked test set and seeds. Almost none of them is *enforced by code*, though. Each one depends on an agent remembering a rule. Specific problems:

- "Test labels are only read by `eval_test.py --final`" is a convention, not a barrier. The label table that every feature builder needs also holds the test labels.
- "Every invocation is logged to RESEARCH_TREE.md" is a manual step.
- "Locked test: one final evaluation *per final model*" lets someone pick a model on TEST by declaring several "final models".
- The label-permutation negative control, as specified, **cannot detect** the most likely leak: label-derived features cached from all labels. The permuted held-out labels are uncorrelated with real-label features too, so the control reads ≈0.5 even when features leak.
- No module owns the IDs, and the project already has inconsistent ones (evidence below).
- Nothing stops an agent from computing metrics its own way, dropping kinases from a macro-average, averaging over 3 of 5 folds, or reporting a number that no file supports.

Quick evidence gathered for this review (read-only checks):
- **Three different ID schemes are already in use.**
  - `docs/data_profile.md` defines a site as "(human UniProt accession **incl. isoform suffix**, position)", while DESIGN §2 says canonical. 892 of the 14,711 human-native PSP rows in `cache/psp_human_mapped_rows.parquet` have an isoform-suffixed accession.
  - The compendium uses `ACC_POS_M1` IDs.
  - In kinase-library `kinome_information.tsv`, MATRIX_NAME ≠ GENE_NAME for **174 of 404** kinases (ERK2 vs MAPK1, and so on).
  - Any agent that joins on its own idea of an ID will silently lose rows.
- **The size of the macro-average denominator depends on the split.** 150 kinases have ≥15 unique human-native sites, so about 3 test positives at a 20% split. 324 kinases have ≥3 sites. Which kinases reach "≥3 test positives" depends on cluster lumpiness, so the set has to be frozen. Cantley covers 148 of these 150 by gene (missing: BCR/ABL, LRRK1). The coverage of Phosformer and NetworKIN is unknown, so a "drop kinases the baseline can't score" rule would change the denominator between rows of the same table.
- **Leak-detection power of a random-label "canary kinase".** I simulated this (`<agent-scratchpad>/canary_sim.py`) using real S/T 15-mers from the compendium and a per-kinase PSSM like F2a. I compared a PSSM built on all labels (the leaky case) with 5-fold out-of-fold scoring (the clean case):

  | Canary positives | Leaky AUROC | Clean AUROC |
  |---|---|---|
  | n=10 | 0.994 | 0.500 ± 0.112 |
  | n=30 | 0.977 | 0.501 ± 0.073 |
  | n=100 | 0.870 | 0.494 ± 0.037 |

  The separation is huge, so canaries are a cheap, decisive end-to-end leak alarm.
  - Side finding: when negatives were scored by a *different* model than the held-out positives (a full-data PSSM instead of the fold PSSM), clean AUROC fell to 0.45–0.48. Mixing models within one per-kinase AUROC biases it, and the harness must forbid that (AIFM-5).
- **Silent fallbacks already exist.**
  - `scripts/profile_compendium.py:188`: `except FileNotFoundError: pass` silently drops three summary fields.
  - `io_psp.py` uses `errors="coerce"` on positions without counting the NAs.
  - DATA.md records an unverified assumption: DepMap gene_effect "assumed equal" to CRISPRGeneEffect.
  - These are fine for profiling, but the same patterns in feature code would be invisible.
- **Repo and infrastructure hazards.**
  - `.gitignore` excludes `*.parquet`, so any split or label artefact written as parquet is silently uncommitted.
  - `<project root>` is NFS (`<nfs-server>:/<share>/...`), where `flock` is unreliable and concurrent `git` from several agents will collide on `index.lock`.
  - `agents/` exists but is empty, so there is no delegation registry yet.

Estimated cost of everything proposed: **about 3 agent-hours**, mostly inside work already planned in §12 step 1–2 (labels, splits, harness), and **under 1 h on the critical path** if one agent builds it. That is about 5% of the budget, and it buys protection against the failures that would invalidate the headline result.

---

## Proposed changes (prioritised)

### AIFM-1 (critical): canonical modules that every agent must import
**Section:** new §11.1 (and references in §2–§4, §6).

**Change:** add this text:
> All agents use these modules and must not re-implement them:
> - `ksp.ids`: `site_id(acc, pos)` returns `"P12345_S123"`. Only the canonical accession is used. The residue letter is part of the ID and is asserted against the UniProt 2026_03 sequence.
> - `ksp.ids.kinase_id(name, source)` maps a name to the HGNC symbol through alias columns in `kinase_list.tsv`. Its sources are PSP GENE/KINASE, Cantley MATRIX_NAME, Hijazi row names, Phosformer and NetworKIN names, and OmniPath. An unknown name **raises**. Bulk mapping returns `(mapped, unmapped_report)` and never drops rows silently.
> - `ksp.ids.map_window(15mer)` is the single 15-mer to canonical mapper, used for PSP, the compendium, CPTAC (RefSeq-based) and Ochoa.
> - `ksp.splits`: loads `results/splits/{substrate_clusters,outer_split,cv_folds,lko_folds}.tsv` and exposes `fold_of(site_id)` and `is_test(site_id)`.
> - `ksp.labels`: the only reader of label caches (see AIFM-2).
> - `ksp.eval`: the only metric code (see AIFM-5).

Add a pytest guard (`tests/test_guards.py`). It greps `src/` and `scripts/` and fails if any file other than `ksp/labels.py` and `ksp/eval_test.py` references the label-cache paths, or reads PSP `Kinase_Substrate_Dataset` outside `io_psp`/`labels`.

**Why:** the project already has three incompatible ID schemes and two definitions of "site". Inconsistent IDs across agents is the most common silent row-loss and mismatch failure in parallel projects. The alias table also makes baseline kinase coverage explicit, which AIFM-6 needs.

### AIFM-2 (critical): physically seal TEST and replace "one run per final model" with a freeze-then-batch protocol
**Section:** §4.2 and §11 "Locked test".

**Change:**
1. **Separate files.** The split script writes `cache/labels/dev_labels.parquet` and `cache/sealed/test_labels.parquet`, plus the external test set, as separate files.
   - `ksp.labels.dev()` asserts that no returned site has `is_test=True`.
   - Only `eval_test.py` knows the sealed path; the AIFM-1 guard enforces this.
   - The feature pipeline receives a `LabelView(train_site_ids)` object, never a path.
2. **FREEZE.json.** Before any test evaluation, commit `results/test/FREEZE.json`. It lists *every* configuration to be scored on TEST: the primary model, B1–B6, the label-variant and regime ablations, T3 LKO, the external test and the benchmarKIN variant. For each it records the run_id, git commit and the sha256 of the TEST prediction file. The primary model is named there.
3. **Single batch.** `eval_test.py --final` refuses to run without a committed FREEZE.json, or if any prediction sha256 differs. It scores all listed configurations in **one batch** and appends to an append-only, committed `results/test/LEDGER.tsv` (timestamp, git hash, FREEZE sha, caller). A re-run needs `--rerun-reason` and appends a new row; it never overwrites.
4. **No development after unsealing.** After the ledger has its first row, no model or feature changes are allowed except the DEV+TEST deployment retrain.
   - A CV>test discrepancy (see AIFM-12 for the threshold) is investigated and reported.
   - Any fix that comes out of it is reported as **post-hoc** next to the unchanged headline, never in its place.

**Why:** "one run per final model" allows model selection on TEST. Declaring several "final" candidates and reporting the best is the agent version of p-hacking. Today's barrier is a sentence in a brief, and agents writing a quick sanity check routinely load whatever label file is at hand.

### AIFM-3 (critical): canary kinases as the end-to-end leak alarm, plus a planted-signal positive canary
**Section:** §11 "Negative controls" and "Positive controls".

**Change:** replace "label-permutation within kinase → ≈0.5" with:
> **Canary run**, required for every model family that is promoted (R2, R3, the stack, R4, and the B3/B4 retrains).
> - Add 12 pseudo-kinases `CANARY_01..12` to the DEV labels. Each one's positives are random U_obs DEV sites of matching residue class: 4 canaries each at n=10, 30 and 100; 10 S/T and 2 Y.
> - Each canary "shadows" a real kinase's descriptors, so that every kinase-level feature path (F1, F3, F4, F5, F6b/c, F7) is exercised.
> - Run the **full** pipeline: features, negatives, training, stacking, calibration.
> - Pass condition: mean canary CV AUROC ≤ 0.57 (about 3 sd above 0.5 for this mix, by simulation), and no canary with n=100 above 0.62. A failure blocks promotion.
> - Canaries are excluded from all real metrics. They are run as a separate control run, so they never enter real kinases' hard-negative pools.
>
> **Planted positive canary.** A synthetic kinase gets 60 DEV positives that carry a planted motif, a planted PPI edge and a planted co-regulation profile. It must reach CV AUROC ≥ 0.9. This catches broken feature joins (off-by-one positions, index misalignment after merges) that would otherwise just look like a weak model.
>
> Keep label permutation, but specify that labels are permuted **before** feature computation.

**Why:** simulation shows a leaky label-derived feature gives a canary AUROC of 0.87–0.99, against 0.50 for a clean one. A stale feature cache, early stopping on the evaluation fold, a non-nested stack or calibration on held-out data would all surface here. The existing toy planted-signal unit test checks one function; the canary checks the assembled pipeline, which is where convenience code leaks. Cost: one extra CV run per promoted model family.

### AIFM-4 (critical): rules for leakage introduced by convenience code, with an assert for each
**Section:** §7 header and §8.

**Change:** add this text:
> **Training-row features must be computed leave-*cluster*-out.** F2a PSSM, F3 "proximity to the kinase's training substrates", F6a co-regulation, F6d activity correlation and F7 n_labels must, for a training row, exclude every training substrate in the same substrate cluster, not only the site itself.
> - Assert: a planted toy in which a training substrate's own protein is in the kinase's substrate set must give proximity equal to the value without it.
>
> **No early stopping or HPO on the evaluation fold.** LightGBM and neural early stopping use an inner split of the training folds.
> - Assert: the `ksp.eval` OOF writer checks that each prediction's `model_id` was trained with `fold ≠ eval_fold`.
>
> **No preprocessing fitted on held-out data.** Scalers, imputers, PCA (F7) and calibrators are fitted on training rows only.
>
> **Stacking (R3 → R2):** the R3 score for evaluation-fold rows comes only from R3 trained without that fold. This is asserted through `model_id` lineage in the manifest.
>
> **Scoring consistency:** for any per-kinase AUROC, positives and negatives in a fold are scored by the same model.

**Why:**
- In F3, a training substrate's own protein is in the reference set, so its proximity is trivially maximal on training rows but never on evaluation rows. The model then over-trusts the feature and degrades silently.
- Early stopping on the evaluation fold is the most common source of optimistic CV in agent-written LightGBM code.
- Mixing models within one AUROC biased it by −0.02 to −0.05 in my simulation.

### AIFM-5 (major): one evaluation harness, a run directory and manifest, and PI verification by recomputation
**Section:** §6 last bullet, §11 "Reproducibility".

**Change:**
> Every experiment runs as `python -m ksp.run --config configs/<name>.yaml` and writes `runs/<YYYYMMDD-HHMM>_<workstream>_<name>_<git7>/` containing:
> - `manifest.json` with:
>   - git commit, a dirty flag (refuse to run if the tree is dirty, or save `git diff` as `patch.diff`), the config, seeds;
>   - sha256 of locked artefacts and inputs;
>   - start/end time, wall-clock, peak RAM, GPU id, n_jobs;
>   - **fallback counters** and **feature coverage** (AIFM-7).
> - OOF and test predictions, saved to `cache/runs/<run_id>/preds.parquet` with their sha256 in the manifest. Columns: `kinase_id, site_id, fold, seed, model_id, score`.
> - `metrics.json` and `per_kinase.tsv`, produced **only** by `ksp.eval.score(preds_path)`.
>
> `ksp.eval` asserts:
> - completeness: every (eval kinase × expected site) pair is present for all 5 folds and all seeds;
> - no NaN or inf scores;
> - unique site IDs, residue-class consistency, and canaries excluded;
> - it stamps `harness_sha256` and `locks_sha256`.
>
> `scripts/make_tables.py` refuses to put rows with different harness or lock hashes, or `complete=false`, in the same table.
>
> **PI rule:** a number enters RESEARCH_TREE only after `scripts/verify_run.py <run_id>` passes. That script recomputes metrics from the saved predictions with the harness (seconds), checks the sha256 values, and checks that the git commit exists. The PI does not trust metrics stated in an agent's summary.

**Why:** this catches hallucinated or hand-edited numbers, partial runs presented as complete ("5-fold" when 2 folds OOM'd), runs without a commit, and metric code that differs between agents. Recomputing from predictions is cheap, so verification does not depend on trusting the producer.

### AIFM-6 (major): freeze the evaluated kinase sets and the evaluation labels; give missing kinases a fixed handling rule
**Section:** §6 and §3 "Label variants".

**Change:**
> The split script writes and locks `results/splits/eval_kinases_{cv,test,lko}.tsv`. These are the kinases with ≥3 human-native positives in the respective held-out set; expect about 150 for TEST.
> - Every macro metric uses exactly this set.
> - A method that cannot score kinase *k* gets a constant score for *k*, so its AUROC is 0.5 and its AUPRC equals prevalence, and it is counted. Its coverage, as n_scored over n_eval, is reported next to it.
> - A secondary table restricted to the intersection of covered kinases is reported for each baseline.
>
> **CV evaluation labels are always human-native (L1-type) pairs, for every label variant**, exactly as for TEST. L2, L3 and L2+ext change only the training labels.
> - Assert: `eval_labels_sha256` is identical across all label-variant runs.

**Why:**
- If baselines drop unscored kinases, the denominator changes per row, and the method with the narrowest coverage looks best.
- §3 fixes test labels, but it says nothing about how CV ablations are scored. If L3 were also evaluated on L3 labels, it would "win" just by having more and easier evaluation positives. This is metric drift that happens without anyone deciding it.

### AIFM-7 (major): no silent fallbacks; every fallback is counted, and feature coverage must fall in expected ranges
**Section:** §7 header and §11.

**Change:**
> - Missing inputs **raise**. No `except …: pass`, and no `fillna(0)` or mean-imputation in feature code; NA plus a missingness flag is the only allowed representation.
> - Every deliberate fallback calls `ksp.log.fallback("<name>", n)`, and the counts go into the manifest. Examples: nearest-atlas Cantley for non-atlas kinases, ortholog mapping by symbol, sequence-based site mapping, NA from pairwise-complete minimum-n.
> - Each feature module declares an expected coverage range in `docs/feature_provenance.md`, for example AFDB pLDDT ≥95% of U_obs sites, F6a non-NA fraction from a pilot ±10%. A run outside the range fails.
> - Known traps that need explicit handling and a unit test:
>   - ESM-2's 1,022-token limit: embed a ±500-aa window around the site and never truncate from the N-terminus.
>   - AFDB splits human proteins >2,700 aa into fragments: map by fragment offset.
>   - CPTAC RefSeq positions: use `map_window`.
>   - `errors="coerce"`: count the NAs produced.

**Why:** an agent that hits a missing file or an edge case under time pressure writes a fallback and reports success. The model still trains, so nothing downstream fails. Coverage ranges turn "silently 30% NA" into a hard error.

### AIFM-8 (major): baseline fidelity gate and an equal tuning budget
**Section:** §8 R0.

**Change:**
> Before a re-trained baseline (B3 Phosformer, B4 SELPHI2.0) is run on our splits, it must reproduce its published number on its published protocol within a stated tolerance: Phosformer random-pair split AUROC ≈0.88 ± 0.03, and SELPHI2.0's pair-CV metric ± 0.03.
> - If it fails, the baseline is labelled "approximate re-implementation" in every table, and the gap is reported.
> - Baselines get the same HPO budget as our models (the same number of configurations and the same CV), recorded in the manifest.
> - The baseline agent must not be the agent developing R2/R3.

**Why:** straw-man baselines are a characteristic agent failure. A re-implementation that silently underperforms (wrong tokenisation, no pretraining, default hyper-parameters) makes "we beat Phosformer" meaningless. The fidelity check costs about 1 GPU-h.

### AIFM-9 (major): delegation, review and shared-document protocol
**Section:** §11 "Independent review", plus a new "Agent operations" paragraph.

**Change:**
> **Delegation.** The PI keeps `agents/REGISTRY.tsv` with these columns: agent_id, model, task, inputs (with sha256), output paths, acceptance checks (shell commands that must pass), wall-clock budget, GPUs/cores, status and reviewer.
>
> **Task briefs** use a template in `docs/AGENT_BRIEF.md` with these parts:
> - goal;
> - inputs;
> - outputs;
> - **executable acceptance criteria**: pytest targets, `verify_run.py`, canary pass;
> - budget;
> - forbidden actions: reading sealed paths; editing `ksp/eval.py`, `ksp/splits.py`, locked artefacts or existing asserts in `tests/`; editing RESEARCH_TREE, LESSONS or DESIGN.
>
> **Agent output.** Agents write `agents/<agent_id>/SUMMARY.md`. Its first line is the harness's `headline()` string (primary and co-primary metrics, n_kinases, complete flag, run_id, harness hash), and it has an explicit "not done / partial" section.
>
> **Review.**
> - The reviewer receives the run_id and the code diff, not the producer's summary, so the review is not anchored on the producer's claims.
> - The reviewer re-runs one fold × one seed from a clean checkout and must match within 0.005, or within the seed sd for GPU models.
> - The reviewer inspects `git diff` of `tests/` for weakened or removed asserts, checks the fallback counters and coverage, and writes `runs/<id>/REVIEW.md` using `docs/REVIEW_CHECKLIST.md`.
> - [x] and [!] in the tree require a REVIEW.md that passes.
>
> **Shared documents.** Only the PI edits RESEARCH_TREE, LESSONS and DESIGN. DATA.md rows are added only via `scripts/make_data_rows.py`.

**Why:** the current text ("reviewed by an agent that did not produce it") does not say what is reviewed, so reviews become a re-reading of the summary. Agents weakening a failing test to get green is a known failure, and so are lost updates to shared markdown. The registry plus the tree makes it possible to resume after a PI context reset.

### AIFM-10 (major): git, concurrency and resource discipline on a shared NFS tree
**Section:** §11 "Reproducibility" and §10.

**Change:**
> **Git.**
> - Agents never run `git add -A`, `git stash`, `git checkout -- .` or `git reset`.
> - They commit only their own paths via `scripts/commit_paths.sh <paths>`, which takes a lock with `mkdir`, because `flock` is unreliable on NFS, and retries.
> - Commit within 10 min of a completed run.
> - Locked artefacts are TSV or JSON, because `*.parquet` is gitignored.
> - A guard (`scripts/check_no_restricted.py`, run by the commit helper) refuses staged files that contain more than 50 PSP 15-mer windows, so no restricted data is committed.
>
> **Caches.** Each cache artefact has a single owner (listed in REGISTRY). Writes are atomic (write to `.tmp`, then `os.replace`), and file names carry a content or config hash.
>
> **Resources.**
> - Every agent gets explicit `CUDA_VISIBLE_DEVICES` and `n_jobs`/`OMP_NUM_THREADS`, at most 16 by default; LightGBM's default would use all 64 cores in every agent at once.
> - The RAM cap is noted in the brief.

**Why:** with 3–5 agents in one working tree, one agent's `git add -A` commits another's half-written code, and concurrent writes to shared caches corrupt them. Oversubscribed threads turn 5-minute jobs into hour-long ones, which then hit time gates.

### AIFM-11 (major): long-run hygiene and hard time gates
**Section:** §10 and §12.

**Change:**
> **Before a run.** Any job expected to take more than 15 min first runs a timing pilot (1 fold, about 5% of rows) and records the extrapolated time in REGISTRY. It then runs in the background under `timeout` set to 2× the estimate.
>
> **During a run.**
> - It writes per-fold checkpoints and can resume.
> - It writes `runs/<id>/status.json` as a heartbeat (fold, step, ETA) every 5 min.
> - `scripts/runs_status.py` lists stale runs (no heartbeat for more than 15 min).
>
> **Time gates** (T0 = 2026-10-26 13:20):
>
> | Time | Gate |
> |---|---|
> | T+15 h | no new feature families |
> | T+17 h | candidate configurations chosen, canary runs passed |
> | T+18 h | FREEZE.json committed and test batch run |
> | T+19.5 h | coverage/deployment inference started (the ~1 h, 130M-pair job) |
> | T+21 h | REPORT draft from generated tables |
>
> Work not finished at a gate is recorded as `[ ] not run` with its reason. It is never presented as a result.

**Why:** runs that die overnight without anyone noticing, and finishing too late to freeze, are the most common way agent projects end with a partial headline. The 130M-pair inference and the test batch must not compete with late experiments.

### AIFM-12 (minor): quantified alarm thresholds and number provenance in the report
**Section:** §11 "Leakage alarms" and §12 step 10.

**Change:**
> **Alarm thresholds.**
> - A CV-to-test drop in macro AUROC (E-hard) greater than 0.05, or outside the paired bootstrap 95% CI, triggers the AIFM-2 investigation.
> - A single feature family holding more than 50% of LightGBM gain, or any one kinase-level feature with |ρ| > 0.9 with the label within a kinase, triggers a leakage audit before promotion.
>
> **Coverage at precision.** The claims of ≥0.5 and ≥0.8 precision are checked on TEST: estimated precision from the CV calibration against observed test precision per bin, reported in the same table.
>
> **Number provenance.**
> - REPORT.md tables are generated by `scripts/make_report_tables.py` from committed `metrics.json`.
> - Every [!] line in RESEARCH_TREE carries `(run: <run_id>)` or `(src: <file>)`.
> - `scripts/check_numbers.py` cross-checks numbers in REPORT against the cited JSON or TSV.

**Why:** "investigated" and "implausible dominance" are not decidable rules for an agent. And model-based coverage claims are exactly the kind of unverified number that propagates into a report.

---

## Build order and cost (all fits in §12 step 1–2)
| Item | Who | Cost |
|---|---|---|
| `ksp.ids` (kinase alias table plus site and window mapper), `ksp.splits`, `ksp.labels` with sealed TEST | splits/labels agent | +45 min beyond what is already planned |
| `ksp.eval` harness, manifest, `verify_run.py`, `make_tables.py` | harness agent (in parallel) | about 1.5 h |
| guard tests, `commit_paths.sh`, restricted-data check, REGISTRY and brief template | PI or a small agent | about 30 min |
| canary generator and planted positive canary | harness agent | about 30 min, plus one CV run per promoted model family |
| FREEZE.json and LEDGER mechanics in `eval_test.py` | harness agent | about 20 min |

## Risks I could not resolve
- **Same Unix user for all agents.** The sealed TEST file cannot be protected by permissions, so sealing relies on the guard test and the brief. An agent that runs ad-hoc Python outside `src/` and `scripts/` (for example `python -c`) bypasses the grep guard. The only real mitigations are reviewer diff inspection and the canary run.
- **GPU nondeterminism.** R3 and B3 will not reproduce bit-exactly. The review tolerance must use the seed sd, and that needs at least 3 seeds before review.
- **Reviewer independence is partial.** Reviewers are the same model family with the same priors, so they may share blind spots, for example about which PPI sources are "unbiased". Canaries and recomputation are the model-independent checks.
- **Too many FREEZE entries.** Listing all ablations makes the test batch large. Agents may also be tempted to exclude ablations that look bad on CV. FREEZE.json should be drafted at T+17 h from the pre-declared §12 list, not from CV results.
- **NFS behaviour of `os.replace` and git locks across hosts.** I have not tested this. If agents run on different hosts against the same mount, have a single committer (the PI) instead of `commit_paths.sh`.
