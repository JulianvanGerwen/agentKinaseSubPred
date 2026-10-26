# Review r2: AI-agent failure modes (DESIGN v1 §0, §4.3, §11, plus the code written so far)

Reviewer role: how autonomous multi-agent research projects fail. I compared v1 against r1 (`r1_ai_failure_modes.md`), and I read every file in `src/ksp/` and `scripts/` as it stood at about 14:20–14:25. That includes the in-progress `splits.py`, `features/ppi*.py` and `baselines/phosformer/`. All checks were read-only.

## Verdict
v1 adopted most of r1 on paper: canonical modules, the seal, FREEZE/LEDGER, canaries, manifests, PI-only git and the delegation registry. The code that now exists shows the expected gap between the design and what agents are building:

- The canonical index uses a **different row order** from the feature caches that are already built.
- The **all-labels table stays freely readable**, and one script already uses all PSP KS labels in feature validation.
- **18% of human-native positives sit outside the evaluation negative universe.**
- A **third ID mapper** has appeared, and kinase-name homonyms leak into the external labels.
- No artefact records the versions of the inputs it was built from, while upstream tables were being rebuilt during this review.
- Several silent fallbacks, and tests that pass by skipping.

None of this invalidates work done so far. I checked that ESM covers 350,736 of 350,736 index sites and that the ESM kinase domains match `cache/index/kinases.parquet`. But each item becomes invisible once G2 feature stores and models exist. All eight fixes fit in the G1 core-infra work, at about 2–3 agent-hours in total.

---

## Proposed changes (prioritised)

### AIFM2-1 (critical): a feature-store contract, because row order already differs between index and caches
**Sections:** §7 header, §11.1.

**Evidence:**
- `scripts/build_index.py` sorts `cache/index/sites.parquet` by the `site_id` *string*, giving `A0A024RBG1_S149, A0A024RBG1_Y95, …`.
- `features/esm_embed.py:merge` (line 161) sorts by `(acc, pos)`. `structure.py` writes in sequence order.
- I compared the two orders: **328,509 of 350,736 rows (94%) sit at different positions.** A feature agent that `np.load`s an ESM or structure array and indexes it by `ksp.index.site_pos()` gets permuted features. Nothing errors, and the model is just weaker. The positive canary may catch this, but only at G2.
- `ksp/index.py` loaders do not check `results/index_manifest.json`. Its docstring says canaries are "appended later" to `kinases.parquet`, which would turn the 503-row kinase axis into 515 rows after some matrices have been built.

**Change:** add `ksp/store.py`, owned by core-infra, as the only way to write and read feature matrices.
- `write(name, matrix, row_ids, col_ids=None, meta)` reindexes to the canonical order itself. It stores `site_id`/`kinase_acc` arrays next to the matrix, and stamps `index_sha256` (from `index_manifest.json`) plus `inputs_sha256` (see AIFM2-5).
- `load(name, site_ids=None)` asserts that the stored IDs equal the current index and that `index_sha256` matches. Otherwise it raises.
- ESM and structure caches are exposed through thin adapters that gather by `site_id`, never by position.
- Canary pseudo-kinases are never added to `cache/index/kinases.parquet`. `ksp.canary` keeps its own table that maps each canary to the row of the real kinase whose descriptors it borrows.
- Add a test that loads every registered store and checks `ids == index`.

### AIFM2-2 (critical): close the seal around the all-labels table, not only `sealed/`
**Sections:** §4.3, §11.1.

**Evidence:**
- §4.3's guard only looks for reads of `sealed/`. But `cache/labels/ks_pairs.parquet` holds **all** labels, future TEST labels included (31,605 pairs), and anyone can read it.
- The same goes for `results/labels/kinase_label_counts.tsv` (all-label counts per kinase) and for the pre-canonical legacy caches `cache/psp_human_sites.parquet` (which has `ks_any`) and `cache/psp_human_mapped_rows.parquet` (isoform-suffixed accessions).
- This has already happened:
  - `scripts/validate_structure_esm_features.py:35–71` reads `ks_any` from the legacy cache, with its own window check instead of `ksp.ids`. It writes "KS-annotated (PSP)" structure statistics to `results/features/phospho_vs_nonphospho_structure.tsv`. That means feature-validation evidence was built from all labels, including the future TEST set.
  - `scripts/audit_ppi_ks_leakage.py:27–40` parses the raw PSP file with `csv`, bypassing `io_psp`.
- External files that contain KS labels are not covered either: `external/omnipath/omnipath_enzsub.tsv`, `uniprot/ptm_tsv`, KinaseCanvas, the SELPHI2 feature/prediction tables, and Phosformer `hf_dataset_curated/{curated,splits}`.

**Change:**
1. `build_splits.py`, at G1:
   - moves `ks_pairs.parquet` to `sealed/ks_pairs_all.parquet`;
   - writes `cache/labels/dev_pairs.parquet`, with TEST proteins physically removed;
   - regenerates `kinase_label_counts.tsv` with separate DEV and TEST columns.
   - Every n_k-dependent quantity uses **within-fold training counts**: α-weights n_k^−α in §8.2, the calibration tier "≥30 DEV positives" in §8.5, and K_cv.
2. Move the legacy caches to `cache/legacy_pre_canonical/`, with a README that says "do not use".
3. `tests/test_no_test_access.py` greps `src/` and `scripts/` for the patterns `ks_pairs|Kinase_Substrate_Dataset|psp_human_sites|psp_human_mapped_rows|omnipath_enzsub|ptm_tsv|KinaseCanvas|selphi2/.*(feature|prediction)|hf_dataset_curated/(curated|splits)|kinase_label_counts|sealed/`.
   - The allow-list is `io_psp, labels, splits, eval_test, build_labels, build_site_universe, build_kinase_list, profile_*`, and the test fails on any other hit.
   - Re-point the validation script at `ksp.labels.load_dev()` once splits exist, or drop its KS class.

### AIFM2-3 (critical): 18% of human-native positives lie outside the negative universe, and 5% of L2 positives lie outside the index
**Sections:** §2 (U_ms), §5 (E-all, training negatives), §6.5 (B0-miss).

**Evidence** (`cache/labels/ks_pairs.parquet` × `cache/index/sites.parquet`):
- **2,431 of 13,392 L1 non-autophosphorylation pairs (18.2%; 1,929 sites) are not in U_ms.** These are PSP low-throughput-only sites. For comparison, only **0.4%** of unlabelled U_obs sites are outside U_ms (21% of positive sites against 0.4%).
- E-all negatives are U_ms by definition. If positives are not restricted the same way, any model that can spot "no MS evidence anywhere" gets a free boost. That evidence is: every F6 block NA, no Ochoa row, and the R3 phospho-channel bits.
- Training negatives drawn from U_ms teach the same shortcut. That shortcut is useless for deployment, which scores MS-observed sites.
- **956 L2 pairs (5.1%) are not in U_obs at all**, so they have no row in `cache/index/sites.parquet`. Every feature join through `site_pos()` will either throw KeyError or, more likely, be wrapped in `.dropna()` and drop them silently.

**Change:**
- §5/§6: "E-all evaluates positives ∩ U_ms against negatives in U_ms." The number of excluded positives is reported per kinase. E-hard is unchanged, because its negatives share the low-throughput enrichment.
- §5 training: E-all-type training negatives are sampled from U_obs, not only U_ms. The training-set builder logs P(all MS-derived features NA | pos) against the same quantity for negatives, and fails if the gap is above 0.05 after negative sampling.
- B0-miss gets an explicit `in_U_ms` indicator.
- The PI decides **now, before stores are built**, what happens to the 956 out-of-index L2 pairs. Either they are added to the index as `in_U_obs=False` rows (recommended: they are real training positives and F0/F1/F2 exist for them), or they are dropped with the count in `label_funnel.tsv`. Rebuilding the index after G1 would invalidate every store (AIFM2-1).

### AIFM2-4 (major): one ID-mapping authority, strict bulk resolvers, and a homonym filter for external labels
**Sections:** §2, §3 (DL-5 filters), §11.1.

**Evidence:**
- **Three independent name→accession mappers now exist, with different rules:**
  - `ksp.ids.CanonicalStore` gene index: UniProt GN primary, then synonyms that never override a primary;
  - `ksp.kinases.KinaseResolver`: own symbol, then family names, then unambiguous synonyms;
  - the new `ksp/features/ppi_idmap.IdMapper`: HGNC approved first, then UniProt primary, unique previous, unique alias, plus its own Entrez/RefSeq/ENSP tables.
- On failure, `IdMapper.from_*` returns `set()`, and `ppi.load_all_edges(log=None)` builds a throw-away `MapLog`, so mapping losses can go unrecorded. `KinaseResolver.resolve` returns `(None, 'unknown')`, although §2 and §11.1 promise that unknown names raise.
- **Homonyms already leak into labels.** OmniPath rows from HPRD/SIGNOR named "PDK1" are mapped by accession to **PDK1, the pyruvate dehydrogenase kinase**. They cover AKT1 T308 and S473, RPS6KB1 T390, PKN2 T816/T951 and ITGB3 T779, which are all PDPK1 sites. "PDK2" is mapped to AKT1 S473, which is the historical name for mTORC2 activity.
- I checked every EXT pair whose site is PSP-labelled to a different kinase whose synonym set contains the EXT kinase's symbol. **11 of 3,891 `EXT_not_in_psp` pairs** collide this way: PDK1/PDPK1 4, CSNK2A1/A2 3, PRKACA/PRKCA 3, PAK5/PAK6 1. All 11 go into the external test, L2+ext and KP.

**Change:**
- `ksp.kinases.resolve_many(names, source, max_unmapped_frac)` returns `(mapped, report)`. It writes `results/<module>/idmap_<source>.tsv` and raises above the declared fraction.
  - It is required for Cantley matrices, Hijazi/decryptM/kinobeads/Davis drug targets, Phosformer, NetworKIN, SELPHI2 and Sugiyama kinase names.
- `ppi_idmap` must take symbol→accession from `ksp.sites.default_store()`. Alternatively, a test asserts that the two agree for every HGNC symbol and differ in fewer than 0.5% of cases. `MapLog` becomes required, is written by every build, and its rates are copied into the run manifest.
- DL-5 gains a homonym filter that drops EXT pairs matching the collision rule above, with counts. Add AKT1 T308→PDPK1 (and *not* PDK1) to `config/bio_controls.tsv`.

### AIFM2-5 (major): lock artefact versions and record inputs; rebuilds happened mid-review
**Sections:** §0 G1, §11.2.

**Evidence:**
- During this review, the following were rewritten within 35 s of each other while other agents were running: `kinase_list.tsv`/`kinase_synonyms.tsv` (14:20:05), `ks_pairs.parquet` (14:20:20) and `cache/index/*` (14:20:40). `config/kinase_curation.tsv` changed at 14:19:57; DELEGATION assigns it to core-infra.
- `cache/esm2_650m/kinases.parquet` was built at 14:17 from the *previous* `kinase_list.tsv`. It happens to be consistent, which I checked.
- No artefact records the inputs it was built from:
  - `index_manifest.json` hashes only its outputs;
  - `site_structure.meta.json` and ESM `timing_*.json` record no kinase-list, universe or index hash;
  - writes are not atomic, so a reader can see a half-written parquet.

**Change:**
- At G1, the PI commits `results/LOCKS.json` with the sha256 of the curation, kinase_list, synonyms, universe, `index/*`, `dev_pairs`, `sealed/*` and `splits/*`.
- Every cache writer stores a `meta.json` with `inputs_sha256` and writes atomically (`.tmp` followed by `os.replace`).
- `ksp.index`, `ksp.labels.load_dev` and `ksp.store.load` verify against LOCKS.json and raise on any mismatch.
- After G1, a locked artefact can be rebuilt only by the PI, with a version bump and a listed set of dependent stores to rebuild. This goes in DELEGATION.md.

### AIFM2-6 (major): remove the silent fallbacks and silent-pass tests already in the code
**Section:** §11.1, "no silent fallbacks", made enforceable.

**Found:**
- **Tests that pass by skipping.** `tests/test_labels_universe.py:18–22` `_skip()` prints "SKIP" and returns, so a missing artefact gives a green run. The tests are run as scripts, not under pytest.
- **ESM kinase descriptors.** `features/esm_kinases.py`:
  - `kinase_set()` silently falls back to a KW-0418 kinase set when `kinase_list.tsv` is missing;
  - lines 305–308 print WARNING and then drop kinases;
  - line 316 silently discards domains that fall outside the sequence.
- **Structure features.** `features/structure.py`:
  - lines 374–375 fall back silently to a different calibration-file path;
  - lines 400–401 and 406 turn NaN into 0 for `idr`, `high_acc` and `disorder_plddt` (`(NaN <= t).astype(int8)`), which labels missing values as "ordered";
  - line 243 uses a default atom radius.
- **Site loaders.** `sites.py:82–85` silently omits compendium counts when the cache is missing, and `sites.py:119` returns a possibly stale Ochoa cache.
- **Kinase list.** `kinases.py:90` silently skips missing Cantley tables, so the Cantley synonyms disappear.
- **PSP parsing.** `io_psp.py:88–89, 299–307`: `errors="coerce"` and `fillna(0)` on the LT/MS counts that **define U_ms**, and the NAs are not counted.
- **Label build.** `build_labels.py:84`: `kinase_included.fillna(False)`.
- **PPI.** `ppi.Graph.rwr_power` returns without an error when it does not converge within `max_iter`.

**Change:**
- Each of these either raises or calls `ksp.log.fallback(name, n)`, with the counters going into the manifest.
- Tests run under `pytest` with `KSP_STRICT=1`, which turns a missing artefact into a failure. `scripts/run_tests.sh` becomes the acceptance check in every agent brief.
- A lint test fails on `except Exception`, `fillna(0)`, `errors="coerce"` or an `if path.exists()` branch inside `src/ksp/features/` unless the same line carries a `# FALLBACK:<name>` tag.

### AIFM2-7 (major): r1 safeguards missing from v1 that the G2/G3 code needs
**Sections:** §3, §6, §8.1, §11.6.

v1 does not contain these r1 items. Each is cheap, and each is the kind of drift that agents introduce without deciding to.

1. **CV evaluation labels are fixed across label variants.** Every L1/L2/L3/L2_vitro/L2+ext ablation is scored on the same human-native, in-class, non-autophosphorylation DEV positives, and `eval_labels_sha256` must be identical across those runs. §3 fixes this only for TEST.
2. **Training hygiene asserts:**
   - no early stopping or HPO on the evaluation fold (use an inner split);
   - anything supervised (calibration, target encodings, label-aware imputation) is fit on training rows only;
   - label-free transforms (the ESM PCA-64 and F7 PCA) may be fit globally, but this must be stated.
   - `ksp.eval` checks `model_id` lineage: `fold_trained ≠ fold_eval`.
3. **Baseline fidelity gate (§8.1):**
   - B3 must reproduce Phosformer's published protocol AUROC within ±0.03, and B4 its published pair-CV metric within ±0.03. Otherwise the table says "approximate re-implementation" and reports the gap.
   - Baselines get the same HPO budget as R2, recorded in the manifest.
4. **Quantified alarms (§11.6):**
   - a CV→TEST drop in M1 greater than 0.05, or outside the paired-bootstrap 95% CI;
   - one feature family holding more than 50% of LightGBM gain;
   - any kinase-level feature with |ρ| > 0.9 with the label within a kinase.
   - Each of these triggers an audit before promotion. "A large gap" is not decidable by an agent.

### AIFM2-8 (major): the Phosformer MLM corpus is selected by labels and contains TEST proteins
**Section:** §8.1 B3.

**Evidence:** `baselines/phosformer/corpus.py` uses the released corpus unchanged. That corpus is 295,320 kinase domains plus **5,338 full-length substrate proteins**, and those proteins were chosen *because* they are curated KS substrates, mostly from PSP. MLM on them is unsupervised, but the selection carries label information. Pretraining on our future TEST substrates gives B3 a "known substrate protein" familiarity prior. That prior inflates its E-all scores, and it would also inflate R3 if R3 ever reused these weights. A pilot already exists in `cache/phosformer/pilot`, and §12 starts pretraining at G1.

**Change:**
- Pretraining starts **after** splits are frozen.
- The corpus builder drops every substrate sequence that maps (≥40% identity and ≥50% coverage, the §4.1 rule) to a TEST substrate cluster. It records the count and the corpus sha in `meta.json` and lists this as a documented deviation from paper fidelity.
- The released Phosformer weights stay "as-published, contaminated" as v1 already says.

---

## Risks I could not resolve
- **The seal is still a convention.** All agents run as the same Unix user, and ad-hoc `python -c` or notebooks bypass grep guards. Canaries (§11.3) and the reviewer's diff inspection remain the only model-independent check.
- **The machine is shared beyond this project.** `ps` shows other users' long-running jobs (engiopt, GAOT and others) on this host, so the "64 cores / 4 GPUs" budget and the §10 timings may not hold overnight. Time pilots should be re-run right before G2/G3 jobs, and the gates should keep slack.
- **The paralog-coverage rule in `splits.py`.** It requires `min(qcov, tcov) ≥ 0.5` (`--cov-mode 0`), so a short protein that shares one domain with a long paralog is not linked. I did not quantify how many cross-split near-duplicate kinase-domain or SH2-domain substrates this leaves. The identical-15-mer edges cover only exact window copies. This is in the data-leakage reviewer's area; flagged here because the code is being written now.
- **Homonyms outside the name-collision rule.** The AIFM2-4 check catches only collisions where the other kinase has a PSP label at the same site. A homonym error at a site with no PSP label (for example a "PDK1" substrate outside PSP) remains undetectable automatically.
- **Who edited `config/kinase_curation.tsv` at 14:19:57** cannot be determined, since only the PI runs git and agents share a user. Ownership is enforceable only through the registry and review.
