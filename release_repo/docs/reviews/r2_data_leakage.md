> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Review r2: data and leakage (DESIGN v1 §2, §3, §5, §7, plus the implemented labels and universe)

Reviewer role: data engineer / curator. I reviewed DESIGN.md v1 (14:30) together with the current code: `src/ksp/labels.py` (including the new core-infra KP/Par/NegativeMasks section), `src/ksp/crossfit.py`, `src/ksp/features/{ppi_label,motif_cantley,phospho}.py`, `results/labels/*`, `docs/labels_universe.md` and the external files named in the brief.

Round-1 changes DL-1 to DL-12 are mostly adopted, and several are already implemented correctly:
- connected-component clusters;
- KP with identical-window and Hamming ≤2 extension;
- Par(k) from the KinCore identity matrix plus the Manning subfamily;
- `crossfit._check_disjoint`;
- `phospho.py` `exclude_pmids`;
- the Cantley reference.

The items below are new findings or gaps in the implementation.

## 0. Checks run (read-only, ≤4 cores)

Scripts are in `<agent-scratchpad>/r2dl/scripts/` (`c1_kp.py`, `c1b_par.py`, `c2_temporal.py`, `c3_ext.py`, `c4_sug.py`, `c5_ikip.py`). Port them to `scripts/leakage/` if adopted, because the scratchpad is session-scoped.

| # | Check | Result |
|---|---|---|
| R1 | **iKiP-DB GMT** (`external/ikipdb_zenodo_21469258/Table7_iKiPDB_Modified_Combined_PSP.gmt`) | 369 sets and 70,413 members. Of the 50,437 members that resolve to sites, **11,969 are PSP pairs for the same kinase** (23.7%; 252 kinases). That is about 61% of all mapped PSP pairs. The file is iKiP merged with PSP, e.g. the CDK8 set contains SMAD1/3 linker sites. |
| R2 | **Sugiyama 2019 vs PSP** (Table S2, 198,536 rows) | 315 of 385 kinase preps resolve to an unambiguous symbol. Among PSP human-native, non-autophosphorylation pairs whose site Sugiyama detects (n=590), the same-kinase Sugiyama hit rate is: in-vitro-only 18.8%, in-vivo-only 22.4%, both 22.7%. Chance is 3–4%. If PSP had curated Sugiyama, the in-vitro-only rate would be far higher. Only 126 of 129,547 resolved Sugiyama pairs are PSP L1 pairs. **So F8-Sugiyama shows no PSP curation (no leak)**, but it carries real signal (about 6× chance). |
| R3 | Sugiyama name traps | The 70 unresolved preps are:<br>• complexes (`CDK2/CycA2`, `CDK7/CycH/MAT1`, `AMPKa1/b1/g1`);<br>• mutants (`ABL1[T315I]`, `BRAF[V600E]`, `ALK[F1174L]`);<br>• `ALK1` (= ACVRL1, not ALK) and `CGK2` (= PRKG2).<br>Separately, acc+position maps 77.8% of rows directly to U_obs. |
| R4 | KP gaps (`known_positive_sets` reads only `ks_pairs` rows; E-hard counted with the implemented `paralog_sets`) | **A. UniProt family-name ECO:0000269 notes** (585 notes; dropped as `kinase_family_name`): 1,518 member pairs not in KP; **613 are currently E-hard negatives** (64 kinases). Examples: PRKACB 60 and PRKACG 60 (PRKACB has only 11 labels), PRKCI 62, PRKCZ 56, PRKCQ 34.<br>**B. OmniPath complex-subunit-only rows** (cyclins, AMPK β/γ, CK2β, CAK) mapped to the catalytic subunit: 127 not in KP, 44 E-hard (CDK1 18).<br>**C. PSP regulatory-subunit rows:** 29 not in KP, 4 E-hard.<br>**D2. OmniPath rows that are only text mining, KEA, PhosphoNetworks or aggregator:** 8,303 not in KP, 3,063 E-hard (301 kinases; GSK3B 201, RPS6KA3 145).<br>**F. Identical-15-mer copies:** 1,175 not in raw pairs. These are already covered by the Hamming-0 block in `NegativeMasks` (no action). |
| R5 | External test (EXT, not in PSP, in class, not autophosphorylation: 3,787 rows) | Near-copy flags:<br>• ±5 aa from a PSP site of the same kinase on the same protein: 4.9%;<br>• same site PSP-labelled to a same-subfamily paralog: 24.5% (same family: 30.3%);<br>• **OmniPath itself lists PhosphoSite as a source but our PSP match failed: 23.5% (889)**. DESIGN v1 does not filter this case.<br>Residual after the three filters: 2,049 pairs on 235 kinases. Of these:<br>• **FAM20C accounts for 183, all from a single high-throughput MS paper (PMID 26091039)**, and the top-10 PMIDs give 333 pairs;<br>• 29% have no PMID (NCI-PID/Reactome via ProtMapper);<br>• 17% cite a PMID that PSP Regulatory_sites cites for the same protein (PSP read the paper but recorded no KS pair).<br>On a simulated 20% protein TEST (5 seeds): A+B 330–487 pairs with 37–50 kinases having ≥3; **tier A (SIGNOR or UniProt) 221–276 pairs with 24–32 kinases having ≥3**. |
| R6 | Temporal file `Kinase_Substrate_Dataset_count_07_2021.txt` (21,402 rows, release "070821") | It is **usable**. It has the same 3 preamble lines and the same columns, except for `IN_VIVO/IN_VITRO_REF_COUNT` in place of the `_RXN` flags and no DOMAIN. Through the production pipeline it gives 17,647 mapped pairs. Results:<br>• 2025 L1 non-autophosphorylation pairs: 13,392. Of these, 11,891 are in the 2021 file, **1,501 are new (238 kinases)**, and 63 more were in 2021 only as rodent pairs;<br>• new pairs by type: 683 new site on a known substrate protein, 386 known site with a new kinase, 363 new protein;<br>• evidence mix of new pairs: in vivo only 40% (vs 23% for older pairs);<br>• only 13 native 2021 pairs disappear in 2025;<br>• simulated 20% TEST: 275–338 new pairs, with **30–41 kinases having ≥3**. |
| R7 | Cantley "as published" | kinase-library 1.8.0 defaults to `_default_phosprot_name='ochoa'` (82,735 S/T and 5,431 Y Ochoa sites). The PSP 07-2021 file appears **only in docstring examples**, so r1 Q12 was wrong for this version. `motif_cantley.py` uses the Ochoa parquet. **This is label-free; no change needed.** |
| R8 | Ascertainment via the U_ms flag | U_ms covers 347,354 of 350,736 U_obs sites. Of the 3,382 non-U_ms sites, **2,135 (63%) are L2-labelled**. Non-U_ms sites are 18.2% of positive pairs and 24.6% of the E-hard pool. Per kinase (≥20 positives), the non-U_ms fraction ranges 5–30% (10th–90th percentile). |
| R9 | Other sources | **Buljan:** S2 is a plain AP-MS table (no curated columns). S1 has a "number of publications (2018)" column.<br>**OpenCell:** interactions and imaging grades (no curated columns).<br>**Hein 2025:** cluster membership is data-driven (Leiden), but the zip also ships `Table_S6_annotations_literature.xlsx` and `Table_SX_annotation_including_imaging.xlsx`.<br>**COMPARTMENTS:** an `integrated_full` file (knowledge + text mining) is on disk next to the experiments file. |
| R10 | Implementation consistency | `ks_pairs` residue classes equal `kinase_list` (0 mismatches), and the label summaries match. **`docs/labels_universe.md` is stale:** it says ST 391 / Y 92 / STY 20 with DYRK/CLK as STY, and that LMTK1–3 lose all labels. The current classes are 402 / 88 / 13. `results/leakage/label_source_studies.tsv` (DL-6) **does not exist yet**. |

## 1. Proposed changes (prioritised, 8)

### DL2-1 — critical — §7 F8, F6c/F6d, `docs/feature_provenance.md`
**Change.** State in F8 that the source is **Sugiyama 2019 Table S2 only** (MOESM3, with the SIDIC/PTMscore flags). Then:
- **Forbid `Table7_iKiPDB_Modified_Combined_PSP.gmt`** as a feature or library everywhere.
- **Forbid `Table5_Cscore.xlsx` and `Kinase_Scores_Long_Weighted_*`** as kinase-activity inputs to F6c/F6d, because they are kinase activities inferred with a PSP-containing library. If iKiP activity profiles are wanted, recompute them ourselves inside `crossfit` from DEV labels.
- Add `tests/test_feature_sources.py`. For every registered non-label-derived feature source that yields (kinase, site) pairs, it computes the overlap with `ks_pairs` and fails at >2× the Sugiyama rate (≈1% of the source's pairs).

Sugiyama name mapping rules:
- complexes map to the catalytic subunit (`CDK2/CycA2` → CDK2);
- mutant preps are dropped or kept as a separate column, never merged with wild type;
- `ALK1` → ACVRL1 and `CGK2` → PRKG2 (hand-reviewed);
- sites are mapped by peptide/15-mer, not by 2019 positions.

**Why.** R1: the GMT contains about 61% of PSP pairs, so this would be a direct label-in-feature leak for 252 kinases. The data_acquisition note lists it next to Sugiyama as "iKiP-DB", so an agent could easily load it. R2: Sugiyama itself is clean, which resolves r1 risk 4.

### DL2-2 — major — §3 "Known positive set", `src/ksp/labels.py`
**Change.** Extend KP with **mask-only rows**: new provenances in `ks_pairs` whose label variant flags are all False, so they reach `known_positive_sets` through `load_dev('all')` but are never labels or test positives:
- (a) UniProt ECO:0000269 family-name notes, expanded to all catalytic members in `kinase_curation.tsv`;
- (b) OmniPath complex-subunit rows mapped to their catalytic subunit (an explicit `config/complex_catalytic.tsv`: AMPKβ/γ → PRKAA1/2, CK2β → CSNK2A1/2, CCNH/MNAT1 → CDK7, CDK5R1 → CDK5, cyclin → CDK);
- (c) PSP rows of excluded regulatory subunits, plus BCR/ABL → ABL1 and NPM/ALK → ALK.

Text-mining, KEA and PhosphoNetworks-only OmniPath rows go into a **KP-broad** mask that is used only in a pre-declared sensitivity analysis of M1/M2.

**Why.** R4: about 660 (k, s) pairs with curated literature evidence are currently scored as *hard* negatives. They are concentrated in PKA/PKC/CDK paralogs, where per-kinase label counts are small (PRKACB: 11 positives vs 60 spurious hard negatives). The cost is negligible: the mask removes <0.5% of negatives.

### DL2-3 — major — §3 "External curated test"
**Change.** Pre-register these filters and reporting rules:
- (i) Also drop rows where **OmniPath lists any PhosphoSite source** (`omnipath_lists_psp`).
- (ii) The paralog filter uses the implemented `paralog_sets()` (≥60% identity OR subfamily), the same Par(k) as the negative masks.
- (iii) Single-kinase names only (already stated).
- (iv) **Cap each (kinase, PMID) at 10 pairs**, or weight pairs by 1/n per PMID. Report with and without FAM20C.
- (v) Tier A = SIGNOR or UniProt ECO:0000269; tier B = HPRD, phosphoELM, or NCI-PID/Reactome via ProtMapper (no PMID). The primary external metric is tier A, macro over kinases with ≥3 positives.
- (vi) Stratum "PSP-read paper": the PMID is cited by PSP Regulatory_sites for the same protein.

Expected sizes: tier A gives 221–276 pairs and 24–32 evaluable kinases; A+B gives 330–487 pairs and 37–50 kinases.

**Why.** R5:
- 23.5% of the "not in PSP" rows are PSP-derived according to OmniPath itself.
- One high-throughput overexpression/KO MS paper supplies 9% of the residual pairs and would dominate FAM20C.
- Pathway-database rows without PMIDs are often inferred.
- Pairs where PSP read the paper but recorded no KS pair may be indirect ones that PSP rejected.

### DL2-4 — major — §3 "Temporal test", §6.4
**Change.** Make the temporal stratum concrete:
- Parse the 07-2021 file through `io_psp.load_psp` after converting `REF_COUNT>0` to the `_RXN` flags, as in `c2_temporal.py`.
- "New since 2021" = an L1, non-autophosphorylation, in-class 2025 pair that is **absent from the 2021 file in every species**. This excludes the 63 pairs present in 2021 as rodent rows.
- Report it by sub-type: new protein / new site on a known substrate / known site with a new kinase.
- Report F6 gains on this stratum with the **F6-clean** variant as well (see DL2-6), because it is 40% in-vivo-only.
- Use the stratum as the "reduced-contamination" comparison for as-published B3/B6 **only after documenting** that their PSP training snapshot predates 07-2021. Otherwise label it contaminated.

**Why.** R6: the file is usable and gives a 30–41-kinase evaluable stratum. Without the any-species rule, rodent-transferred training knowledge would count as "discovery". The in-vivo skew means that F6's own label-source circularity is strongest exactly here.

### DL2-5 — major — §7 cross-fitting rule, §8.4 stacking, `src/ksp/crossfit.py`
**Verdict.** The rule plus `_check_disjoint` (protein and cluster, all kinases' labels on held-out clusters removed) **is sufficient against split leakage for F3-label and F6a**. Three gaps remain.

**(a) Train/eval feature-distribution mismatch.** Eval rows are fit on 100% of outer-training labels, while training rows are fit on about 80%. Count-type and extreme-value outputs therefore shift by about 1.25× between training and evaluation: `nb_n_used`, `nb_n_adj`, `nb_max`, `nb_top3`, and the F6a number of contributing substrates and top-m mean. The shift is largest for low-n kinases, and it distorts the cross-kinase comparisons that M2 relies on.
- *Change:* eval rows get the fit of a **randomly assigned inner fold**: `eval_mode="matched"` as the default, and `"full"` as an ablation.
- Add a test that the per-feature KS statistic between training-row and eval-row negatives is <0.05 on the toy and the real data.

**(b) Labels may enter only through `fit()`.**
- Add `tests/test_no_label_reads.py`. It greps `src/ksp/features`, `src/ksp/baselines` and `scripts/features` for `ks_pairs`, `dev_pairs`, `load_dev(`, `FULL_PAIRS`, `Kinase_Substrate` and `.gmt`.
- All cross-kinase pooling happens inside `fit()`: F2a empirical-Bayes family/group priors, F6d activity profiles, and any per-kinase standardisation by training positives.

**(c) R3/F2b stacking.** It must run through `crossfit` (5 outer × (5 inner + 1) = 30 R3 fits). If that exceeds the G3 budget, **drop stacking** and blend R2 and R3 at score level, with weights tuned on outer OOF. Never stack in-sample R3 predictions.

### DL2-6 — major — §6.7, §13 (label-source circularity, which cross-fitting cannot fix)
**Change.** Before G2, build `results/leakage/label_source_studies.tsv`. It lists compendium PMIDs/PXDs cited by PSP Regulatory_sites, OmniPath curated refs (`ext_refs`) and UniProt ECO:0000269 PMIDs. This is DL-6, accepted but not yet built.

Run **F6-clean** (F6a/F6b recomputed with `exclude_pmids` = that list, as `phospho.py` already supports) as a pre-declared sensitivity analysis on:
- M1/M2;
- the in-vivo-only test stratum;
- the temporal stratum.

A gain from F6 that disappears under F6-clean is reported as circular.

**Why.** Kinase-targeted studies in the compendium (Torin1, BI-D1870, AMPK-DKO and others) are PSP in-vivo evidence for the very drops that F6b/F6a detect. Cross-fitting removes protein/cluster leakage, not study-level leakage.

### DL2-7 — minor — §7 F0/F1/F6 missingness, §8.4 R3 phospho channel, §9
**Change.** Forbid any **centre-site** ascertainment indicator as a feature: `in_U_ms`, per-source presence (compendium/CPTAC/Ochoa/PSP-MS), or the PSP LT/MS counts. The R3 phospho channel and the priming features mark **flanking positions only**, never position 0. F6 NaN patterns stay as they are and are audited by B0-miss.

**Why.** R8: 63% of non-U_ms sites are positives, and the per-kinase share of such positives varies from 5% to 30%. A centre bit would let models learn kinase × ascertainment rather than specificity, even on E-hard.

### DL2-8 — minor — §7 F3/F4 exclusions, `docs/feature_provenance.md`; `docs/labels_universe.md`
**Change.** Add explicit forbids:
- COMPARTMENTS `integrated` channel (it contains knowledge and text mining);
- Hein 2025 `Table_S6_annotations_literature.xlsx` and `Table_SX_annotation_including_imaging.xlsx` (only data-driven cluster membership and enrichment tables are allowed);
- the HPA `Reliability` field (its "Supported" grade means agreement with the literature);
- the Buljan S1 publication-count column.

Also:
- The Cantley nearest-atlas fallback is chosen by kinase-domain identity only, never by CV performance.
- Refresh `docs/labels_universe.md` §2/§4 to the current residue classes (402/88/13; DYRK/CLK/HIPK/LMTK/STYK1 → ST).

**Why.** R9 and R10: each of these files is on disk next to a legitimate file, and the stale doc contradicts the code that agents read.

## 2. Checked, no change needed
- **F1 Cantley "as published"** uses the Ochoa reference and is label-free (R7). Arm (ii), recalibrated over U_ms, is also label-free.
- **Buljan AP-MS interactome, OpenCell interactions and localisation, Hein 2025 cluster membership, BioPlex, HuRI (non-Lit-BM), Human Cell Map:** no KS curation in content. Coverage bias is handled by B0-miss.
- **Implemented core pieces:**
  - `known_positive_sets` with identical windows, and `NegativeMasks` (KP + Hamming ≤2 + Par(k) labels);
  - `paralog_sets` (symmetrised; atypical kinases fall back to the Manning family);
  - `crossfit._check_disjoint`, which removes all kinases' labels on held-out clusters. This also covers paralog/family-prior leakage at the same site, provided pooling happens inside `fit()`;
  - `ppi_label.py`, which excludes s = k (self) and same-cluster seeds;
  - the canonical site index, which holds no label columns.

## 3. Risks I could not resolve
1. PSP has no per-pair PMIDs, so DL2-6 can only bound the label-source overlap from below.
2. 548 OmniPath rows with residue mismatches (numbering errors) are neither labels nor KP. Masking same-residue sites within ±5 would fix this, but its precision is unknown.
3. KP-broad precision (text mining / KEA) is unknown, so the sensitivity analysis decides whether it matters (largest effect: GSK3B at about 2% of its E-hard negatives).
4. The PSP snapshot dates used to train the Phosformer released weights and KolossuS are not verified (this matters for DL2-4).
5. About 22% of Sugiyama rows do not land on U_obs by acc+position. Some are genuinely unobserved lysate sites, and some are numbering drift. F8 coverage depends on peptide-based mapping.
