# Core infrastructure (splits, labels API, eval harness, canaries, manifests, TEST seal)

Owner: core-infra. Frozen at G1 (split frozen 2026-10-26 14:34 CET, seed 0; label tables relabelled 14:50 after the
v2 label amendments, with the split unchanged). Every number below comes from `results/splits/split_numbers.json`
and `split_report.md`.

```
/workspace/.venv/bin/python -I scripts/build_kinase_list.py      # residue-class overrides (config/kinase_curation.tsv)
/workspace/.venv/bin/python -I scripts/build_labels.py           # -> cache/labels/ks_pairs.parquet (pre-split, full)
/workspace/.venv/bin/python -I scripts/build_splits.py           # FREEZE (refuses if SEAL.json exists; --force = PI only)
/workspace/.venv/bin/python -I scripts/build_splits.py --relabel # after a label change: same split, new tables, re-seal
/workspace/.venv/bin/python -m pytest -q tests/                  # 63 tests, ~45 s
```

## 1. Residue classes (DESIGN §2) and kinase list

`config/kinase_curation.tsv` has 36 `residue_class` rows:
- **ST:** LMTK1–3, STYK1; DYRK1A/1B/2/3/4, CLK1–4, HIPK1–4, GSK3A/B. For DYRK, CLK, HIPK and GSK3 the Y labels are cis-autophosphorylation. GSK3 and HIPK are written explicitly even though ST is already their default.
- **STY:** MAP2K1–7, TTK, TESK1/2, PKMYT1, PKDCC, NEK10.
- **Y:** WEE1/2.
- **ERN1 → ST.** IRE1α is a Ser/Thr kinase. Its only label is O75460_Y628 (autophosphorylation, in vitro), which is now dropped as out of class.
- **MAP3K21 → ST.** MLK4 is a Ser/Thr MAP3K; it keeps its EXT labels MAP2K1 S218/S222. Its PSP label Q9NRM6_Y447 is dropped.

Class totals for the 503 included kinases: ST 402, Y 88, STY 13. Inclusion and order are unchanged.

`cache/index/kinases.parquet` was rebuilt with `scripts/build_index.py`:
- Only the `residue_class` column changed (15 kinases). New SHA-256: `e278ec9b…`.
- sites.parquet and proteins.parquet are byte-identical.
- Any kinase-axis matrix saved through `ksp.store` before 14:16 must be re-saved.

Regulatory-subunit rows (CSNK2B, PRKAB1, PRKAG2, PIK3R1, PHKA1, PHKB) are excluded and have 0 L2 labels (verified). Under A7 they become mask-only rows for their catalytic kinases.

Strict name resolver (A4): `ksp.kinases.resolve_kinase_names(names, source, accs=None, on_unknown='collect'|'raise')`.
- It returns gene, accession, `how`, family members, and the unresolved names.
- Per-source homonym rules are in `config/kinase_homonyms.tsv`. "PDK1" resolves to PDPK1 for hprd/signor/cantley/kinase_library/kinhib, is dropped for uniprot/other, and is PDK1 for PSP rows (which resolve by accession). It also covers RSK*/MSK*, PKN/CK1/CK2/PKA/PKC/AMPK family names, MST1/2, PKD1/2 and CDK11.
- UniProt EXT labels are now resolved through it: 2,918 → 2,914 pairs.

## 2. Labels API (`ksp.labels`)

| function | what |
|---|---|
| `load_dev(variant='L2', with_canaries=False, ums_only=None)` | **The only modelling label reader.** Reads DEV proteins only. Variants: L1/L2/L3/L2_vitro/L2_noauto/EXT/EXT_not_in_psp, `eval` (evaluation positives) and `all` (every row, incl. mask-only, for masks). U_ms-only by default, except for `all` (A6). |
| `evaluable_positive(df)` | PSP ∧ L1 ∧ ¬autophos. This is the same set for TEST and CV whatever the training variant (AIFM2-7). |
| `known_positive_sets(pairs)` | KP(k) over all sources, species and flags, plus identical-15-mer U_obs sites. Sensitivity-only rows are excluded unless `include_sensitivity=True`. |
| `paralog_sets()` | Par(k): ≥60% kinase-domain identity, OR same Manning subfamily, OR (no subfamily) same family with ≥50% identity (A8). Symmetric; mean size 2.1. |
| `kinase_domain_identity()` | KinCore MSA identity (identical / both-residue columns, ≥50% overlap of the shorter domain; max over multi-domain entries). The 28 atypicals missing from KinCore use a within-family Biopython global alignment of UniProt domains. |
| `NegativeMasks(pairs)` | `.excluded(k)` = KP(k) ∪ Hamming≤2 15-mers of KP(k) sites ∪ sites labelled to Par(k). `.hard_candidates(k)` = sites with a human-native non-autophos label to a kinase ∉ {k}∪Par(k) (PI r2). |
| `hamming_neighbors(q, ref, 2)` | Exact pigeonhole search; tested against brute force. |
| `new_since_2021`, `external_test_filter`, `ext_test_tiers`, `mask_only_pairs` | A10, DL-5+A9 and A7 respectively; details below. |

Notes on those last four:
- **Temporal (A10):** a pair is new if it is absent from the 2021 kinase-library PSP file in every species, matched by acc+window, acc+residue/position, or kinase symbol+SITE_GRP_ID.
- **External test (A9):** OmniPath rows that list PhosphoSite are dropped, pairs are capped at 20 per PMID, and tier A = SIGNOR/HPRD/UniProt-ECO. Flags are `ext_test` and `ext_tier`.
- **Mask-only rows (A7):** they come from `config/complex_subunits.tsv` (cyclins, AMPKβ/γ, PKA-R, TFIIH/CAK, mTORC, Ku, NEMO) and UniProt family notes, 3,159 rows in total. The 20,652 OmniPath text-mining/KEA rows are `mask_sensitivity_only`.

Columns added at freeze: `in_index`, `in_U_ms`, `new_since_2021`, `ext_test`, `ext_tier`, `mask_only`, `mask_sensitivity_only`, `is_canary`.

## 3. Splits (`ksp.splits`, `scripts/build_splits.py`)

**Substrate clusters.**
- Built from an mmseqs2 easy-search all-vs-all (`-s 7.5 --max-seqs 2000 -e 1e-3 --min-seq-id 0.4 -c 0.5 --cov-mode 0`, 12 threads, 109 s).
- 167,275 edges plus 10,477 identical-15-mer edges, connected into components.
- Result: **11,609 clusters** (8,874 singletons). The largest has 455 proteins (the kinases) and holds 2.35% of labels. That is below the 3% threshold, so nothing was re-thresholded.

**TEST.** Chosen by simulated annealing with soft penalties (seed 0; 400k moves, 11 s):
- 2,353 clusters (20.3%), **3,957 proteins**, 69,430 U_ms sites, 2,952 L2 pairs and **2,205 evaluation positives**.
- Label, U_ms and S/T/Y fractions are all between 19.96% and 20.04%.
- Per-kinase band, for the 228 kinases with ≥5 positives:
  - 187 fall in 15–25%;
  - 223 fall in the integer band [⌊.15n⌋, ⌈.25n⌉];
  - 92 of 94 kinases with ≥20 positives fall in 15–25%.
- The misses are structural, because all their positives sit in 1–3 clusters (e.g. OXSR1 0/15, BMPR1B 0/6).

**CV.** 5-fold annealing over DEV clusters (800k moves). Each fold holds 19.98–20.03% of labels.

**LKO.** 120 kinase clusters (≥50% domain identity ∪ Manning family) in 5 folds of 97–103 kinases, with no family split.

**Kinase sets.**

| | n | AGC | Atyp | CAMK | CK1 | CMGC | Other | STE | TK | TKL | ST/Y/STY |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K_eval | **119** | 24 | 4 | 15 | 4 | 23 | 15 | 8 | 23 | 3 | 94/23/2 |
| K_any | **256** | 41 | 11 | 37 | 7 | 41 | 30 | 26 | 48 | 15 | 205/45/6 |
| K_cv | **203** | 33 | 5 | 28 | 6 | 34 | 25 | 23 | 38 | 11 | 161/38/4 |

K_eval is 119, below the ≈135–160 expected in DESIGN. The cause is the U_ms restriction (A6), which drops 2,431 of 13,392 evaluation pairs.

**Dropped (A6).** 956 L2 pairs have no index row, and 2,934 L2 pairs plus 2,431 evaluation pairs are not in U_ms.

**Readers.** `load_clusters`, `load_test_proteins`, `load_cv_folds`, `load_kinase_folds`, `load_eval_kinase_sets`, `site_cluster`.

**Helpers.**
- `inner_folds(train_sites, k=5, seed)`: cluster-grouped and balanced.
- `train_valid_split_inner`: early stopping inside training folds only (AIFM2-7).
- `cv_eval_rows(fold)`.

**Frozen files.** In results/splits/: `substrate_clusters.tsv`, `test_proteins.tsv`, `cv_folds.tsv`, `kinase_folds.tsv`, `eval_kinase_sets.tsv`, `canary_ppi_partners.tsv`, `split_report.md`, `split_numbers.json` and `SEAL.json` (SHA-256 of every sealed file, plus history).

**Seal (A3).**
- Moved into `sealed/`: `ks_pairs_full.parquet`, `ks_unmapped.parquet` and `kinase_label_counts.tsv`, plus the new `test_labels.parquet` (TEST rows and TEST canaries).
- DEV tables are `cache/labels/dev_pairs.parquet` and `canary_dev.parquet`, plus public `results/labels/kinase_label_counts_dev.tsv`.
- `tests/test_labels_universe.py` now skips its label checks because `ks_pairs.parquet` is sealed.

## 4. Cross-fitting (`ksp.crossfit`; DESIGN §7, §16 A1)

- **Contract:** subclass `LabelDerivedFeature`, decorate it with `@register`, and implement `fit(train_pairs, ctx)` and `transform(query, ctx)`. `transform` returns len(query) rows in the same order.
- **Driver:** `crossfit(factory, train_labels, train_rows, eval_rows, k_inner=5, seed, eval_mode='inner')`.
  - Training rows get an inner 5-fold over substrate clusters.
  - By default, eval rows go to a randomly assigned inner-fold fit (A1). `eval_mode='all'` is reserved for a PI decision.
  - `LeakageError` is raised if any fitted label lies on a query protein or cluster.
- **Outer folds:** `outer_split(fold)`.
- **Toy test:** in the planted-signal test, the "is labelled" feature gives AUROC 1.0 without cross-fitting and exactly 0.5 under the driver.

## 5. Evaluation (`ksp.eval`; the only metric code)

**API.** `prepare('cv', fold=None|f)`, then `evaluate(preds, data, missing='worst'|'drop'|'error')`. The result is an `EvalResult` with `.endpoints`, `.per_kinase`, `.site_level`, `.pair_strata`, `.missing`, `.canary` and `.to_json()`; `strata(res)` gives the stratum table.

**Input.**
- Predictions are (kinase_acc, site_id, score), or (kinase_idx, site_idx, score) for very large outputs. Canaries are appended after kinase index 502.
- An optional `score_site` column holds the cross-kinase score used for M2 (A20).
- Unknown IDs, duplicate pairs and NaN scores raise. Missing pairs are counted per kinase.

**Endpoints.**
- M1 and M4: E-hard AUROC.
- M2: filtered site-centric Hit@5, which drops other labelled kinases and Par(k) (A18). Expected values under random tie-breaking. Hit@1/10, MRR and family Hit@5 are secondary.
- M3: top-1% recall in E-all, with fractional tie credit.
- Secondary: pAUROC at FPR≤0.01 (McClish-standardised, matches sklearn), AP and lift, and E-all AUROC.
- CV uses K_cv (so M5 = M1-CV).

**Strata.**
- Per kinase: Manning group, residue class, training bin (0/1–4/5–19/20–99/≥100 DEV L2), priming, PIKK.
- Per pair: evidence class, new-since-2021 versus pre-2021, external curated.

**Statistics.**
- `bootstrap(res, B=1000|10000, resample='clusters'|'kinases')` gives cluster-bootstrap CIs. M2 and M3 keep each kinase's top-1% cut-off fixed, which is an approximation.
- `paired_bootstrap(A, B)` refuses unmatched universes and returns bootstrap p-values; `holm(p)` applies Holm correction.
- `matched_universe`, `coverage` and `coverage_at` (A17).

**Pooled CV (A24, DESIGN §17).**
- With `split='cv'` and `fold=None`, scores are by default replaced by their mid-rank percentile ((rank − 0.5)/n) within (CV fold, kinase), over universe sites. This is `evaluate(..., fold_normalise=None)`, which turns ON for pooled CV and is a no-op for a single fold or TEST.
- It is applied before M1/M3/M4/M5. M2 uses the supplied `score_site`, or the raw score when none is given.
- `endpoints` gains `fold_normalised`, `M1_per_fold` (macro M1 per fold) and `M1_fold_mean`.
- Leave-kinase-out pooled predictions come from one model per kinase, so pass `fold_normalise=False` for them.
- Check: B0-prior with fold-training counts had pooled M1 0.400 raw; with the option it is 0.500, and per-fold M1 is 0.5.
- Pooled evaluation of 70M rows takes 79 s.

**Speed.** CV fold 0 has 21.3M (k,s) rows, i.e. 516 kinases × 55.6k sites; one evaluate takes 27 s and B=200 bootstrap takes 1.3 s. TEST (≈70k sites) should take about 35 s.

## 6. Canaries (`ksp.canary`, **v2**, PI decision after the R3 audit)

**Negative canaries (12; 6 ST, 4 Y, 2 STY).**
- Labels are **real DEV-labelled sites**: evaluation-positive sites of real kinases (human-native PSP, non-autophos, in U_ms), pooled across kinases.
- Sampling is protein-clustered and site-uniform. Proteins are drawn without replacement with weight n_p/E[min(K,n_p)], and then min(K,n_p) of their labelled sites are taken, with K ~ U{1..4}.
- The v1 sizes are kept (10–70).
- About 70% of canary sites have a same-protein sibling, against ~66% for real labels.
- KP mask = the canary's own sites only.
- The pool is DEV-only, so negative canaries have **no TEST labels** (TEST share 0) and are CV-only. The split did not change.

**Positive canary.** Unchanged from v1: 138 R-x-x-[ST]-[FILMV] sites on 120 planted PPI partners (`canary_ppi_partners()`), TEST share 0.196.

**Acceptance** (`eval.canary_check`).
- The gate is mean negative-canary AUROC on **E-hard** ≤ 0.57. E-all is reported but not gated, because generic site-likeness legitimately lifts it.
- The positive canary must reach ≥ 0.90 on both E-hard and E-all.

**Bias probes** (pooled DEV, mean negative E-hard AUROC; all are regression-tested on fold 0).

| score | v1 | v2 |
|---|---|---|
| random | — | 0.495 |
| −(#U_ms sites on protein) | 0.76 (R3 audit) | 0.489 |
| #DEV-labelled sites on protein | — | 0.453 (E-all 0.894, not gated) |
| "is a labelled site" | — | ties, = random |

An intermediate v2 draft that drew a labelled site uniformly and then took its protein over-weighted label-rich proteins: the #labelled score reached 0.61. It was fixed before release.

**Provenance.**
- v1 (random U_ms sites, used at the 14:34 freeze) is kept in `cache/labels/canary_dev_v1.parquet`. The `generate_v1` function stays for provenance.
- `sealed/test_labels.parquet` lost the v1 negative-canary TEST rows, so the seal was rehashed and SEAL.json `history` records the previous hashes.

**Use.** `load_dev(..., with_canaries=True)` adds the canaries, and `kinases_with_canaries()` gives the extended kinase list (canaries are kinase rows 503–515). Canaries never enter endpoints.

## 7. Manifests and TEST

- **`ksp.manifest.RunManifest(name, config, seed, inputs)`** is a context manager. It records the git hash (read-only), input SHA-256, runtime, n_rows, fallback counters and status. Predictions go to `cache/runs/<id>/predictions.parquet`, metrics to `results/runs/<id>/metrics.json`, and existing runs are never overwritten.
- **`ksp.eval_test` (G4; see `docs/G4_RUNBOOK.md`):**
  - **CLI:** `--real`, `--dry-run` or `--check-predictions`; with none of them it prints usage and exits.
  - **Guards:**
    - FREEZE.json must be git-tracked and unmodified.
    - Schema: exactly one `final`; the 9 pre-declared (baseline_id, endpoint) primaries, or the `unavailable_primary` contingency (p = 1, m stays 9); `missing` must be `worst`; runtime fields are required.
    - Every input file's SHA-256 is pinned in FREEZE.json.
    - Real mode needs B = 10000 and harness_commit = 974d520, with the frozen files unchanged (`git diff --quiet`).
    - Every prediction file passes pre-flight: exhaustive per kinase, finite scores, no canaries, score_site any finite value, assignment_prob in [0, 1], and T3 covers all of K_any.
    - The seal check passes.
  - **Ledger:** machine-readable STARTED (written before sealed/ is read; FREEZE sha, input hashes, git HEAD, SHA-256 of every `src/ksp/*.py`) and COMPLETED blocks.
    - Any COMPLETED TEST block refuses every further batch, unless `--override-completed "<reason>"` (logged).
    - An unfinished STARTED block can only be rerun identically, logged as attempt n.
    - Outputs go to `<sha12>.partial/` and are renamed atomically.
  - **Batch:**
    - Per config: cluster and kinase-resampling CIs (B_secondary, default 1000), per-kinase and site-level tables, strata (from ksp.strata), external-curated ext_tier A/B, A17 reliability with ECE plus Coverage@τ (labelled, external and temporal sites; configs with assignment_prob), and T3 with missing-kinase counts.
    - The 9 primary tests at B = 10000 on matched universes, reporting matched-kinase counts, the fraction of pairs kept and the M2 score column. Holm over the primary 9, and separately over secondary tests.
  - **Rehearsal:** `--dry-run` uses DEV fold 0 as a fake TEST, with every path under `cache/eval_test_dryrun/`.

## 8. Tests (`tests/`, pytest added to requirements.txt)

### 8.0 Deployment labels (G5; `src/ksp/deploy.py`, new module, so labels.py stays frozen under A28)

**Call.** `ksp.deploy.load_deploy_labels(variant="L2", with_canaries=False, ums_only=None)` returns DEV + TEST labels for the deployment model (DESIGN §8.6: the frozen config retrained on DEV+TEST for the released predictions).

**Semantics:** the same as `load_dev`.
- Variants: L1/L2/L3/L2_vitro/L2_noauto/EXT/EXT_not_in_psp, `eval`, and `all` (every row, including mask-only rows, for KP and `NegativeMasks`).
- U_ms-only by default, except for `all`.
- The columns are the same, including `is_canary`.
- DEV and TEST must be protein-disjoint; this is asserted.

**Refusals.** It raises `DeployError` unless:
- FREEZE.json is git-tracked and unmodified;
- `results/test/LEDGER.md` has a **COMPLETED TEST** block for exactly that FREEZE sha (parsed with the `eval_test` ledger markers; STARTED-only, DRYRUN, or another sha do not count);
- the sealed files match SEAL.json.

**Canaries are off by default.** Released predictions must not contain canary pseudo-kinases, and after DEV+TEST retraining there is no held-out data on which a canary gate could be evaluated.
- `with_canaries=True` adds the DEV canaries (v2) and the sealed TEST canary rows (positive canary only), for pipelines whose code path expects canary kinases.
- Their CANARY_* outputs must be dropped before release.

**Access and tests.** `deploy.py` is on the sealed-access allowlist (`tests/test_no_test_access.py`). The tests (`tests/test_deploy.py`, using a toy sealed set on public TEST proteins) check four things:
1. It refuses with an empty ledger, with STARTED only, with COMPLETED for another sha, with DRYRUN, and with an untracked FREEZE.
2. Once COMPLETED exists, it returns a strict superset of DEV for L2, all and eval, with identical columns, and is canary-free by default.
3. It refuses a tampered seal.
4. It refuses a DEV/TEST protein overlap.


- **test_metrics:** hand-computed toys and sklearn cross-checks.
- **test_splits:** no protein, cluster or identical-15-mer overlap between TEST and DEV or across CV folds; mmseqs edges stay inside clusters; constraints hold.
- **test_masks:** Hamming search vs brute force; KP, paralog and E-hard semantics.
- **test_crossfit:** the planted leak.
- **test_canary.**
- **test_eval_harness:** ID errors, missing-prediction policies, oracle and random behaviour, A18/A20, canaries, bootstrap, TEST refusal, manifest.
- **test_feature_overlap:** A11 helper `check_feature_overlap(pairs)`, which fails above 3× chance overlap with DEV PSP pairs.
- **test_no_test_access:** grep over the PI's pattern list, plus "feature modules never load labels".

**Grep-guard violators, for their owners.**
- `scripts/features/phospho_study_lists.py` reads OmniPath enzsub and iPTMnet directly for PMIDs only (the A12 exclusion list). The PI has added it to ALLOW.
- Docstring-only, listed in `KNOWN_VIOLATORS`: `src/ksp/sites.py`, `src/ksp/ids.py`, `src/ksp/features/ppi.py`, `scripts/make_data_rows.py`, `scripts/profile_psp.py`.
- `scripts/validate_structure_esm_features.py` does not match any pattern.

## 9. Biological controls

`config/bio_controls.tsv` has 44 canonical pairs, all in PSP L2 and none autophosphorylation, covering 9 groups with 10 Y sites.
- Numbering follows the human canonical sequence, e.g. CREB1 S133 → S119, H3 S10 → S11, EIF2S1 S51 → S52, PRKAA1 T172 → T183, RPS6KB1 T389 → T412, LAT Y191 → Y220.
- `in_test_split` is set for 8 of 44.
- The panel is not used for tuning.

## 10. Assumptions and open points

- Coverage for mmseqs edges means coverage of both query and target ≥ 50%. Identity is mmseqs fident over the alignment.
- In CV, masks use DEV labels only. In TEST, masks use DEV ∪ TEST labels, because Hamming neighbours can cross splits.
- M3/M2 bootstrap keeps each kinase's top-1% threshold fixed.
- Pooled out-of-fold CV evaluation mixes fold models within a kinase. Per-fold evaluation is available through `prepare('cv', fold=f)`.

## 11. API stability (R2, R3, B3, B4 and baseline agents)

I expect no breaking changes. Every public signature above is frozen, including load_dev, NegativeMasks, paralog_sets, splits readers, inner_folds, crossfit, prepare/evaluate/bootstrap/paired_bootstrap, RunManifest and canary helpers. Possible changes are additive only, with defaults that keep current behaviour:
- `evaluate(..., fold_normalise=None)` was added (A24). It changes pooled-CV M1/M3/M4/M5 relative to runs evaluated before 17:00, so re-evaluate those from their saved predictions. Further optional arguments may follow, always with defaults.
- `EvalResult.endpoints` may gain extra keys, e.g. bootstrap summaries. Existing keys keep their names.
- If the PI adopts it, a shared `ksp.labels.sample_negatives(...)` generator would be added (A15/A22). Nothing existing would be renamed.
- Label tables change only through `build_splits.py --relabel`. That mode keeps the split, evaluation positives and canaries fixed and fails otherwise, so DEV evaluation sets cannot shift silently.

## 12. CV leaderboard (`scripts/reeval_cv_runs.py`)

The script re-evaluates every pooled-CV run that has `cache/runs/<id>/predictions.parquet`, using the current `ksp.eval`:
- A24 fold rank-normalisation is on, except for LKO pooled runs (`fold_normalise=False`).
- Missing pairs use `missing='worst'`.
- It writes only `results/runs/<id>/metrics_a24.json`; existing metrics files are kept.
- It then builds `results/cv/LEADERBOARD.{md,tsv}`.

The script is re-runnable: already re-evaluated runs are skipped unless `--force`. Other options are `--only REGEX`, `--list` and `--board-only`.

How the board treats special cases:
- **Partial coverage** (scored fraction of the required pairs < 0.95, e.g. B3-pub scored only positives + E-hard, B4 only the U_B4 sites): M3 is blanked, because unscored E-all pairs ranked last inflate E-all metrics. Compare such runs with `ksp.eval.matched_universe`.
- **Canary gate:** shown only for runs that scored canaries. A run trained before the v2 canary labels is marked "(v1-trained)".
- **Official baselines (A23, latest run of each arm):** B1 arm (i) for M1/M3, B1 arm (ii) for M2, B7, and the B0-prior/pop/miss floors.
- **R3 pooled runs** have no saved pooled predictions, so they cannot be re-evaluated until the R3 agent saves `predictions.parquet` for the pooled runs.
