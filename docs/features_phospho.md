# F6 phosphoproteomics features (feat-phospho, 2026-10-26)

DESIGN §7 F6, §6.7, §16 (A6, A11, A12). Code: `src/ksp/features/phospho.py` (blocks, GPU masked correlation, F6a +
crossfit family), `phospho_inhib.py` (F6b), `phospho_cptac.py` (CPTAC residualisation, A-loop, F6c).
Scripts (all `CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/<x>`; GPU 1, ≤10 cores):

```
phospho_build_blocks.py   # 1. blocks (compendium 30 s + CPTAC 10 s)       -> cache/features/phospho/blocks/
phospho_study_lists.py    # 6. benchmarKIN + KS-label-source study lists     -> results/features/phospho/exclude_lists.json
phospho_f6b.py            # 2. F6b (20 s)   [--variant clean_bench --exclude-pmids ... | --exclude-drugs ...]
phospho_f6c.py            # 3. F6c (16 s)   [--variant v --exclude-cptac all|BRCA,...]
phospho_selftest.py       # kernel + exclusion + crossfit-contract tests (synthetic; dummy cluster map)
phospho_sanity.py         # 5. NON-OFFICIAL single-feature AUROC + missingness (40 s incl. 5-fold F6a cross-fit)
```
All kinase × site matrices are written/read through `ksp.store` (`phospho/<name>[_variant]`, float16,
503 × 350,736, canonical index hash checked). NaN = not computable; there are **no count / coverage features** (A6).

## 1. Data blocks
**Compendium** (`data/01–03`, via `cache/compendium_wide.parquet`). 187,200 of 189,962 peptide UIDs map to index sites
(mapping_log); 179,238 sites (125 sites merge >1 UID → nan-mean). Each condition column is divided by its robust SD
(1.4826·MAD; `conditions.parquet:scale`, median 0.6). Blocks = one per PMID (≥ 8 conditions); decryptM (PMID 36926954,
790 contrasts) is split per drug series (Rituximab 306, Dasatinib 42, …; series < 8 → `36926954:other`); the ~150
studies with < 8 conditions are pooled into `pooled_small` (351 conditions). 66 compendium blocks.
Condition metadata (`cache/features/phospho/conditions.parquet`, `results/features/phospho/conditions_meta.tsv`):
PMID, block, perturbation type (keyword heuristic `classify_condition` from scripts/profile_compendium.py), treatment
time in minutes (`Timepoint`; unit-less numbers = minutes, ASSUMPTION; fallback: a single time in the Condition text),
drug names/doses (F6b table below).

**CPTAC** (BCM v1.2 LinkedOmics S3, tumour only, 10 cohorts, 83–122 tumours). Per cohort and site:
`site ~ 1 + host protein (same ENSG, same tumours) + confounders`, batched OLS on GPU (ridge 1e-3), ≥ 10 usable tumours.
Confounders available: **ESTIMATE stromal and immune score (RNA) and WES purity** from the PDC pan-cancer
`Clinical_meta_data_v1` (S3 `*_meta.txt` has none); median-imputed + a missingness indicator when partly missing
(BRCA, LUAD, OV, UCEC). **No TMT plex/batch and no ischemia time are in either metadata** → not regressed (limitation).
The BCM files contain `Inf` entries (551 in BRCA) → treated as missing. Sites without host-protein data are dropped so
every CPTAC value is protein-residualised: 47.8k of 49.6k (BRCA) … kept after ≥ 10 obs: 19.6k (COAD) – 47.3k (LSCC);
`results/features/phospho/cptac_residualisation.tsv`. Kinase protein abundance per cohort is residualised on the
confounders only (`cptac_kinase_protein.npz`; 253–379 kinases per cohort).

## 2. F6b inhibitor response (not label-derived)
**Inhibition matrix A (condition × kinase).** Hijazi (PMID 31959955, 60 min): provided binary `kinhib` matrix
(rows resolved with `ksp.kinases` resolver: 400/420 kinase rows used; unresolved = lipid kinases PI4KB, PIK3C*,
PIP4K2*, PIP5K1*, PIKFYVE, BCR/ABL, PRKY, Q6ZSR9, PDPK2P). Other conditions: drug names found in Condition **and**
Control text; `A_row = occ(condition drugs) − occ(control drugs)` so 'DMSO vs drug' rows are negative and
'drug + stimulus vs stimulus' keeps only the drug. `occ = dose/(dose+Kd)`; combinations `1 − Π(1 − occ)`.
Kd source priority: **Klaeger 2017** Table S3 Kd_app (min over the 4-line, BT-474 and Caki-1 sheets; `n.i.` = measured,
not bound → 0, no fallback; `n.d.` → next source) > Reinecke 2024 Kinobeads (tool compounds) > CurveCurator re-fit
(only drugs absent from Klaeger) > Davis 2011 Kd (non-mutant constructs; blank = > 10 µM → 0) > Hijazi kinhib binary
(listed target → occ 1). Drug identity = union-find over names + synonyms of all sources + the curated
`results/features/phospho/drug_aliases.tsv` (56 aliases: same-compound aliases, and class tokens such as `MEKi`, `p38i`,
`ATRi` with designated targets and assumed potency; rapalogs/AZD8055/dactolisib/tideglusib/Torin1 → MTOR/GSK3B with
literature IC50s because Kinobeads do not capture them — **ASSUMPTION**, values from memory). Doses parsed
(pM/nM/µM/mM, sci notation, decryptM `A+B_1to2_10nM`, PMID 32034161 `_0-1` = 0.1 µM); unknown dose → 1 µM (flag).
Time ≤ 120 min where known (40 conditions dropped); 67 used conditions have unknown time (kept, flagged).
Every condition with its parsed drugs/doses/drop reason: `f6b_conditions.tsv`; the unmapped/dropped drug-like ones:
`f6b_unmapped_or_dropped_drug_conditions.tsv` (unmapped examples: AZ123, 1-NM-PP1 (kinase not named), noscapine,
genistein, PRL3i, p97 inhibitors, RG7356 antibody, Olaparib; SC-514 and omipalisib (no protein-kinase occ ≥ 0.5)).
**Used: 692 conditions** (Hijazi 165, affinity-mapped 527, incl. decryptM dose series; low doses kept as near-zero
anchors), **384 exposed kinases** (Σ|A| ≥ 1).

**Per-site regression** (`ridge_deconvolution`, GPU, 76,380 sites with ≥ 3 observed F6b conditions): scaled log2FC
~ A·β, ridge λ, no intercept; `f6b_ridge_beta` (β < 0 = site decreases when k is inhibited), `f6b_ridge_z` (β / ridge
posterior se, σ² moderated to the median, d0 = 4); NaN where k is not exposed in the site's observed conditions.
λ was chosen on the sanity AUROC (configurations tried: 1, 10, 100, 300, 1000 → macro AUROC of −z 0.539, 0.560,
0.579, 0.590, **0.606**; `results/features/phospho/tuning/`): strong shrinkage wins, i.e. with ~5–50 observed
conditions per site the deconvolution is close to its marginal (Aᵀy) limit. `f6b_contrast` = mean y over observed
conditions with A_k ≥ 0.5 − mean over |A_k| < 0.1 (≥ 2 / ≥ 3 conditions).
**Non-separable kinases** (KB-5): cosine of inhibition profiles ≥ 0.95 → 191 pairs, 109 kinases
(`f6b_nonseparable_pairs.tsv`, `f6b_kinases.tsv`); kinase-level feature `phospho/f6b_kinase_maxcos` (503-vector). Many
are kinases hit only by one promiscuous Hijazi drug (e.g. SGK1/SGK2/HASPIN/ADCK1/CSNK2A3: one inhibitor each).

## 3. F6c CPTAC kinase state (not label-derived)
Per cohort Pearson r (≥ 20 paired tumours) of the residualised site with (i) kinase protein abundance
(`f6c_prot`, 427 kinases) and (ii) the kinase A-loop activity profile (`f6c_aloop`, 99 kinases with an observed A-loop
site in ≥ 1 cohort; 31–71 per cohort). Fisher z averaged over cohorts with weight n−3. A-loop = residues strictly between
the DFG (first `DFG` after 35 % of the UniProt kinase domain, fallbacks `D[FLWY]G`, `D.G`) and the first `APE`
(fallbacks `[ASPTGC]PE`, `.PE`, `.P.`) 12–45 residues later; 452 kinases get a range (257 strict DFG…APE), checked on
MAPK1 (T185/Y187), AKT1 (T308), PRKACA (T198), SRC (Y419), CDK1 (T161), AURKA (T288). Profile = mean of the z-scored
residualised A-loop sites (`f6c_aloop_sites.tsv`). Sites on the kinase's own protein → NaN.

## 4. F6a co-regulation with training substrates (label-derived — caller cross-fits)
`coreg_features(train_pairs, query_site_idx, blocks, cluster_of_site=..., exclude_same_protein=True,
exclude_pmids=..., exclude_cptac=..., m=3, min_shared=8, min_blocks=2)` → {'top', 'med'} (n_kinases × n_query).
Per block: pairwise-complete Pearson (6 masked fp32 matmuls on GPU, ≥ 8 shared columns) of each query with every
training substrate; Fisher z; substrates on the query's protein or substrate cluster and the query itself are excluded;
per kinase: mean of the top-3 z and median z (padded kinase buckets, batched topk/nanmedian). Across blocks: weighted
mean (weight = columns in which the query is observed in that block), NaN unless ≥ 2 blocks. Compendium and CPTAC are
computed separately. **Speed: all 503 kinases × 350,736 sites with ~18k DEV pairs: 12.7 s (compendium) + 4.0 s (CPTAC)
on one RTX 4090.** Crossfit contract: `ksp.features.phospho.CoregF6a` (registered `f6a_coreg`;
fit(train_pairs, ctx) / transform(query, ctx) → columns `f6a_{comp,cptac}_{top,med}`; ctx.clusters used for the
cluster exclusion); variants with `make_f6a(name, exclude_pmids=..., exclude_cptac=..., exclude_drugs=...)`.
Self-test (`phospho_selftest.py`): masked_corr equals pandas pairwise Pearson (max diff 2e-7); with a dummy cluster
map the same-cluster perfect correlate is excluded (r 0.999 → 0.835 = next substrate); query-only kinase → NaN;
1 block → NaN; exclude_pmids works. **F6d** (kinase-activity profile) is not built: an honest version needs a per-query
leave-protein/cluster-out profile, and it overlaps F6a (DESIGN marks it low priority).

## 5. Sanity (NON-OFFICIAL; `sanity_summary.tsv`, per kinase `sanity_auroc.tsv`)
DEV L2 positives (load_dev, U_ms) vs 3,000 random U_ms DEV sites of the kinase's residue class without a label of k;
kinases with ≥ 5 positives; AUROC among non-missing pairs (macro). F6a cross-fitted with a temporary
protein-grouped 5-fold split (no query ever sees its own protein's labels).

| feature | kinases | macro AUROC (non-missing) | coverage pos / neg |
|---|---|---|---|
| f6a_comp_med | 152 | **0.633** | 0.39 / 0.14 |
| f6a_comp_top | 152 | 0.607 | 0.39 / 0.14 |
| f6a_cptac_med | 127 | 0.611 | 0.30 / 0.14 |
| f6a_cptac_top | 127 | 0.516 | 0.30 / 0.14 |
| f6b_ridge_z (−z) | 111 | 0.606 | 0.24 / 0.06 |
| f6b_ridge_beta (−β) | 111 | 0.604 | 0.24 / 0.06 |
| f6b_contrast (−) | 101 | 0.598 | 0.21 / 0.05 |
| f6c_aloop | 35 | 0.555 | 0.08 / 0.03 |
| f6c_prot | 114 | 0.532 | 0.31 / 0.16 |

Top-3 is worse than median for CPTAC (with many substrates the top-3 of noisy z is high for every query).
**Ascertainment (`sanity_missingness.tsv`):** a labelled DEV site has a non-missing value for some kinase 2.1–2.5×
more often than an unlabelled U_ms site (e.g. F6b 47 % vs 21 %, F6a-comp 45 % vs 18 %, F6c 43 % vs 19 %); within kinase,
positives are covered 2–4× more often than random sites. Missingness alone would give per-kinase AUROC ≈ 0.5 +
(cov_pos − cov_neg)/2 ≈ 0.58–0.63 — the B0-miss floor (DESIGN §6.5) must be reported next to any F6 gain.

## 6. Clean variants (DESIGN §6.7, §16 A12) and leave-drug-out
Every builder takes `exclude_pmids` / `exclude_cptac` (+ `exclude_drugs` for F6a/F6b and `load_blocks`).
`results/features/phospho/exclude_lists.json`:
* `clean_bench`: the 22 benchmarKIN perturbation PMIDs present in the compendium (24 Hernandez-Armenta PMIDs, 22 in
  the compendium, 276 conditions, incl. Hijazi 31959955) + all CPTAC (the tumour benchmark uses BRCA, CCRCC, GBM,
  HNSCC, LSCC, LUAD, UCEC). Built: `phospho/f6b_*_clean_bench` (483 conditions, 373 kinases; Hijazi removed); F6c clean
  = all NaN (no CPTAC); F6a: `make_f6a(..., exclude_pmids=..., exclude_cptac=True)`. Matching is by PMID only (neither
  the benchmark metadata nor the compendium has PXDs). `benchmarkin_studies.tsv` lists them.
* `f6_clean` (A12): 18 compendium PMIDs (55 conditions) cited as evidence of KS records in OmniPath enzsub references
  (HPRD/SIGNOR/KEA/ProtMapper/dbPTM) or iPTMnet phosphorylation rows with an enzyme; 10 PMIDs are cited by PSP
  Regulatory_sites. `label_source_studies.tsv`. Built: `phospho/f6b_*_f6_clean`. Limitation: PSP KS rows carry no
  PMIDs, so PSP-curated HTP evidence (e.g. from decryptM/Hijazi) cannot be detected this way; neither decryptM nor Hijazi
  is cited by these sources.
* `clean_bench_perturbation` / `clean_bench_tumour`: the two halves separately.

## 7. Assumptions / limitations
- Kinobeads Kd_app (lysate) ≠ cellular occupancy; allosteric/MTOR/GSK3 drugs use literature potencies (aliases file).
- Class tokens (`MEKi`, `JNKi`, …) assume a potent, selective drug (occ 0.9).
- decryptM dose series share their control → correlated errors in the ridge and in F6a correlations within a series.
- Compendium per-condition scaling uses MAD over all sites; qPTM-sourced contrasts may be significant-only.
- No plex/ischemia regression for CPTAC (not available); purity/immune/stromal only.
- Inhibitor response ≠ direct phosphorylation (downstream kinases respond too) — F6b is evidence, not a label.
- iKiP-DB GMT/activity tables are not used (A11).
