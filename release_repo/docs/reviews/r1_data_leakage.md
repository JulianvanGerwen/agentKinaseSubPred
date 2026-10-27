> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Review r1: data and leakage (labels, universe, clustering, leakage paths, external test)

Reviewer role: data engineer / curator (PSP, UniProt, OmniPath, CPTAC, compendium). Reviewed DESIGN.md v0 (2026-10-26).

## 0. Quick checks I ran (read-only, 4 cores)

The scripts are in `<agent-scratchpad>/chk/` (c1–c5) plus inline one-liners. The PI should port the ones it adopts to `scripts/leakage/`, because the scratchpad is session-scoped. Inputs: `cache/psp_*.parquet`, `data/0[1-3]*`, `external/omnipath`, UniProt 2026_03, `<mmseqs install dir>` (v 8cc5ce3…).

| # | Check | Result |
|---|---|---|
| Q1 | Same-kinase sibling sites | 68.3% of HMR KS pairs (13,306/19,477) have ≥1 other site **on the same protein labelled to the same kinase**. |
| Q2 | Co-regulation of same-protein sites (compendium, sites ≥50 conditions) | Same-protein site pairs: median Pearson r 0.23, and 20.3% have r>0.5. Random pairs: median 0.03, and 3.6% have r>0.5. |
| Q3 | Self pairs (kinase acc = substrate acc after human mapping) | 1,124 self pairs. The PSP-organism `autophos` flag misses 53 of them and wrongly flags 3. **435 (39%) lie inside the kinase's own UniProt "Protein kinase" domain**, so the peptide is a verbatim substring of the kinase-domain sequence that R3/B3 encode. 47 kinases with 1–19 labels have ≥50% self labels. 29 kinases with ≥5 labels have <5 non-self labels. |
| Q4 | Identical-window copies | 818 (k, site) pairs in U_obs have a 15-mer identical to a site labelled to k on another protein (14-3-3 isoforms, 4E-BP1/2, PRKAA1/2, H3 variants, MYL9/12B, HSP90AA1/AB1), but are not labelled to k themselves. |
| Q5 | Paralog label sharing (HMR) | MAPK3 shares 398 of its 563 sites with MAPK1. GSK3A shares 75/112 with GSK3B. PRKAA2 shares 51/87 with PRKAA1. AKT3 shares 20/27 with AKT1. CAMK2B shares 30/61 with CAMK2A. PRKACB has only 11 sites against 1,045 for PRKACA, so "PKA" in PSP is effectively PRKACA. |
| Q6 | Ascertainment via feature availability | 26% of HMR KS sites (3,125/12,132) have **no MS evidence** in PSP. Among canonical PSP sites with neither compendium nor PSP-MS evidence, **68% are KS-labelled** (base rate 4.9%), and this flag alone gives a pooled AUROC of 0.606. Compendium coverage count among quantified sites gives AUROC 0.635. Coverage of KS-substrate proteins vs other U_obs proteins: OpenCell target 15.1% vs 4.4%; Human Cell Map prey 39.5% vs 19.0%; BioPlex 293T node 79.9% vs 61.1%. |
| Q7 | Clustering (U_obs = PSP-human ∪ compendium on Swiss-Prot canonical: 303,123 sites, 18,478 proteins) | `easy-cluster` 40%/50% cov-mode 0 gives 11,937 clusters, but **53,894 of the 163,247 protein pairs that meet the 40%/50% criterion are split across clusters** (set-cover clustering is not transitive). Adding identical-window union-find gives a largest group of 417 proteins. Connected components of all 40%/50% hits plus windows: largest group 451 proteins, holding 2.7% of label rows. At 30% identity (coverage on either side) the largest component is 2,612 proteins with 24% of label rows, which is unusable. After window merging, labelled site pairs in different groups with ≥11/15 identical windows: 74 (20 share a kinase); at ≥12/15: 21. |
| Q8 | OmniPath "curated, no PSP source" rows | 6,003 rows. 1,928 have a non-kinase/complex enzyme (OmniPath expands complexes: CCNB1 93, PRKAB1 91, PRKAG1 74, CSNK2B 65; phosphatases PPP2CA 82, PTPN1 60). The 4,075 rows with a PSP kinase break down as: 355 within ±5 aa of a PSP site for the same kinase and protein; 2,169 on a site PSP labels to another kinase, **1,062 of them to a same-family paralog**; 606 at a distant site on a protein PSP already links to the same kinase (only 63% residue-match the canonical sequence, so mostly numbering errors); 945 fully new. In addition, 347 (site, reference) groups are assigned to ≥2 paralogs (PRKAA1+PRKAA2, CHUK/IKBKE/TBK1), which means family-level curation. After the filters in DL-5, about 2,300 candidate pairs remain (before removing regulatory subunits). On a random 20% of design groups (5 seeds) that gives 374–644 test pairs, with 44–59 kinases having ≥3 test positives. The set is skewed: SRC 136, MAPK1 90, MAPK3 72, PRKACA 69. |
| Q9 | Label-source studies inside the compendium | The PSP Regulatory_sites PMIDs include **10 compendium studies** (29 sites, 19 of them KS-labelled). These are kinase-targeted perturbations: Torin1 21659604, BI-D1870/RSK 25002506, AMPK-DKO 30631047, p38i+UV 29523821, VE-821 30001349, FGF7 34086370, IL4 34731634, and others. OmniPath KS evidence cites 17 compendium PMIDs (1,692 rows, 579 PSP H×H pairs; mostly HPRD/KEA citing Olsen 2006/2010, Dephoure 2008, Daub 2008, Mayya 2009). The PSP KS file has no PMIDs and OmniPath strips PSP references, so the true extent cannot be measured. |
| Q10 | Kinase identifiers | Symbol joins have several traps: PSP `PRPF4B` is `PRP4K` in UniProt 2026_03, so it is absent from pkinfam under the PSP symbol. Kinase-library `PDK1` is PDPK1, while HGNC `PDK1` is PDHK1 (both are in K). Other ambiguous symbols in UniProt synonyms: PKD1/PKD2 (PRKD vs polycystins), MST1 (STK4 vs MST1), TAK1 (MAP3K7 vs NR2C2), NIK (MAP3K14 vs MAP4K4), PAK5/PAK7/PAK6, MAK (also an ALPK3 synonym), PASK (STK39 synonym). In kinhib, `ADCK3` and `CABC1` both map to COQ8A. |
| Q11 | Site ID mapping | 12,925 compendium sites are not exact on Swiss-Prot canonical, 11,581 of them because the accession is not a Swiss-Prot canonical. 10,696 are rescued by a unique 15-mer, 376 hit several locations, and 1,853 have no hit. 605 PSP KS sites (1,059 rows) have their 15-mer at >1 canonical location. |
| Q12 | Kinase-library percentile reference | `kinase_library.core.percentile` defaults to a **PSP 07-2021 high-confidence site set** (LT_LIT ≥1 or MS ≥5). |
| Q13 | PhosphoNetworks vs PSP | PhosphoNetworks site calls coincide with 618 PSP H×H pairs, with an evidence mix like the PSP baseline. So PSP did not import them wholesale. |

## 1. Proposed changes (prioritised)

### DL-1 — critical — §7 F6(a), F6(d), F2, F3 "fold-aware"
**Change.**
- Replace "leave-site-out" with **leave-substrate-group-out**. Every label-derived feature for a site *s* is computed from training positives outside *s*'s substrate group (the DL-9 group).
  - For validation and test rows: only from training folds.
  - For training rows: by an inner K-fold over the training folds' groups (nested out-of-fold), so that the train-time and test-time feature distributions match.
  - This applies to F2a PSSMs, F2b/R3 stacked scores, F3 "proximity to training substrates" and F6a/F6d.
- Add a toy test in which two sites on the same protein are planted with identical profiles. The test must fail if the feature for one sees the other.

**Why.**
- Q1: 68% of KS pairs have a same-kinase sibling on the same protein.
- Q2: same-protein sites are strongly co-regulated (protein abundance, shared peptides; median r 0.23 vs 0.03).
- Under leave-site-out, F6a for a validation site would correlate it with its own protein's sibling, which is a positive. That is a near-direct label leak in CV, and it would make CV look much better than TEST, where siblings are not in DEV.
- The §7 header says "inside each fold", but F6a/F6d say "leave-site-out". The ambiguity must be resolved in favour of group-out.

### DL-2 — critical — §5 (training negatives, E-all, E-hard), §3
**Change.** Define a per-kinase **"known-positive mask" = union of all label sources**: L3 (all species), rodent-only transfers, UniProt ECO:0000269, the OmniPath curated rows, family/complex-level labels (DL-8), and identical-window copies of any positive (Q4). Masked (k, s) are *excluded* from negatives in training, E-all and E-hard, and from test positives when they are not human-native.

Additionally, make **E-hard paralog-masked** for the primary endpoint: for *k*, drop E-hard negatives that are positives of a kinase with ≥80% kinase-domain identity to *k*. Report the unmasked version as secondary.

**Why.**
- The design currently treats rodent-only and external-only pairs on TEST proteins as negatives, because test labels are human-native PSP only.
- The 818 identical-window copies would be "hard negatives" with sequence identical to a positive.
- E-hard, the primary endpoint, draws negatives exactly from annotated sites, which are dense with paralog labels (Q5). For example, the 371 MAPK1-only sites (about 74 on TEST proteins) would all be MAPK3 E-hard negatives, and most are probably true MAPK3 substrates. That penalises correct predictions and biases the primary metric against CMGC/AGC families.

### DL-3 — critical — §3 "Autophosphorylation", §6
**Change.**
1. Define self pairs by accession after human mapping (kinase UniProt acc = substrate acc), not by the PSP same-organism symbol flag.
2. **Exclude self pairs from the primary and co-primary endpoints**: both self positives and *k*'s own unlabelled sites as negatives.
3. Report autophosphorylation as a separate stratum. Keep self positives in training (default L2), with L2-noauto as the ablation.
4. For R3/B3 (kinase-sequence encoders), also report T1 with self pairs whose site lies in the kinase's own domain removed from *training*. This shows whether the encoder learned substring matching.

**Why** (Q3).
- 39% of self sites lie inside the kinase domain that R3/B3 encode, so they are trivially solvable by substring matching.
- Pair-level protein features (PPI self, co-expression = 1, co-localisation = 1) are trivially maximal for self pairs.
- Because the macro metric weights kinases equally, the 29–47 low-label kinases whose labels are mostly self would have their scores inflated exactly in the understudied stratum the project cares about. Making this a sensitivity analysis, as currently written, gets the default backwards.

### DL-4 — major — §5 "Evaluation negatives", §6 coverage, §7 F4
**Change.**
1. **E-all is computed over MS-observed sites only** ("E-all-MS"): sites in the compendium, CPTAC, Ochoa 2020, or PSP with MS_LIT+MS_CST ≥1. This applies to both positives and negatives. LT-only positives are scored under E-hard only.
2. Remove substrate-level missingness flags (F4 "plus missingness flags"). Allow only pair-level NA.
3. Extend the popularity diagnostic to a **"missingness-pattern model"**: indicators of which feature families are non-NA, plus compendium/CPTAC coverage counts and PPI-screen node presence. Report it under E-all-MS and E-hard. If it exceeds 0.60 macro AUROC on E-all-MS, the precision/coverage estimates (Priority 2) must use E-hard-calibrated or coverage-stratified calibration.
4. Coverage estimates use E-all-MS calibration.

**Why** (Q6).
- 26% of positives are LT-only, and 68% of MS-invisible sites are positives.
- LightGBM routes NA, so the not-numeric-count rule in §7 F6 does not prevent the model from learning "nothing measured ⇒ studied by LT ⇒ positive".
- The PPI and localisation screens are leak-free in content, but cover labelled substrate proteins 1.3–3.4× more often, through abundance and study bias.

### DL-5 — major — §3 "Independent external test"
**Change.** Pre-register these filters for the external test:
- (a) Residue must match the Swiss-Prot canonical sequence at the stated position. OmniPath has no windows; where windows are available (SIGNOR/UniProt), map by 15-mer.
- (b) Drop rows within ±5 aa of a PSP site for the same kinase and protein, since these are numbering duplicates.
- (c) Drop rows whose site PSP labels to a same-family paralog of the external kinase, or evaluate them only at family level (DL-8).
- (d) Drop regulatory subunits, cyclins, phosphatases and complex-expanded rows. Treat (site, PMID) groups assigned to ≥2 paralogs as family-level.
- (e) Exclude rows that are text-mining-only (Sparser/REACH/RLIMS-P/BEL), KEA-only or PhosphoNetworks-only.
- (f) Use two tiers: **tier A** = SIGNOR ∪ UniProt ECO:0000269 (direct curator calls); **tier B** = HPRD/NCI-PID/Reactome/Phospho.ELM/dbPTM/Li2012-only (pathway-level, often motif-assigned or indirect).
- (g) Use macro per-kinase metrics, with ≥3 positives per kinase.

Expected size: about 40–55 kinases are evaluable, which is enough for a secondary macro metric but not for strata. Report the pooled result too, with a note on the SRC/MAPK skew.

**Why** (Q8). Roughly 35% of the "non-PSP" curated rows with a PSP kinase are PSP near-copies, either numbering or paralog/family re-assignments. They are not independent evidence and would reward models that learned PSP's family conventions.

Also note: OmniPath does not carry SIGNOR's `DIRECT` flag, so the lit-review filter `DIRECT=YES` cannot be applied until the user supplies the PhosphoSIGNOR TSV. Re-request it.

### DL-6 — major — §7 F6(b,d), §6 "Downstream"
**Change.**
1. Build `results/leakage/label_source_studies.tsv`. It lists compendium PMIDs cited as evidence by PSP Regulatory_sites, OmniPath (any source) and UniProt MOD_RES references. Q9 gives 10 + 17 PMIDs.
2. Pre-register an **"F6-clean" sensitivity**: F6 is recomputed with those studies' conditions removed.
3. Report F6's gain separately on test labels with in-vitro support (both + in-vitro-only) and on in-vivo-only labels. If the gain is confined to in-vivo-only labels, F6 is learning PSP's in-vivo curation rule ("site drops on inhibitor X"), which includes indirect targets, not direct targeting.
4. For benchmarKIN, beyond excluding compendium studies by PMID/PXD, the **tumour (CPTAC) benchmark variant must drop all CPTAC-derived features** (F5 CPTAC co-abundance, F6a CPTAC blocks, F6c). The tumour gold standard is defined from CPTAC kinase protein and activating-site levels, which F6c encodes directly. Please confirm against Müller-Dott 2025 Methods.

**Why.** Q9 shows that PSP curated kinase-targeted perturbation studies that are in the compendium. Without PSP KS PMIDs, the overlap can only be bounded from below.

### DL-7 — major — §2 Kinases; all joins
**Change.**
- The **kinase key is the UniProt accession**. Every source is joined through an explicit symbol→accession table, built from current HGNC/UniProt primary symbols, then previous symbols, then synonyms. Only unambiguous synonyms are used, plus a hand-reviewed collision list: PDK1/PDPK1, PKD1/2 vs PRKD1/2, MST1/2, TAK1, NIK, PAK5/6/7, MAK, PASK, CDC2L1/2, ADCK3/CABC1 → COQ8A, PRPF4B → PRP4K.
- Join Cantley matrices by the `UNIPROT_ID` column of `kinome_information.tsv`, never by MATRIX_NAME.
- Verify mouse kinases by pkinfam's mouse accession column, which pairs each human entry with its mouse ortholog, rather than by the assumed same-symbol 1:1 rule (data_profile §1.1 ASSUMPTION).
- Write unresolved symbols to `kinase_list.tsv` with a reason.

**Why** (Q10). One wrong join, for example Cantley "PDK1" attached to PDK1 (PDHK1), silently corrupts a baseline and a feature for a kinase with its own labels.

### DL-8 — major — §2 exclusions, §3, §6 stratification
**Change.**
- Do not remap PSP rows whose kinase is a regulatory subunit onto a catalytic subunit: drop them. These are CSNK2B (15 pairs; 10 duplicated on CSNK2A1/2), PRKAB1 (12) and PRKAG2 (9), the last two all duplicated on PRKAA1/2.
- In external sources, treat complex or family assignments (AMPK, CK2 holoenzyme, IKK, "PKA", "PKC", CDK–cyclin, mTORC1/2) as **family-level labels**. They are used only in a family-level metric and in the DL-2 mask, never as kinase-level positives.
- Add a pre-declared **family-level metric**: credit the maximum over catalytic paralogs. Report it for PKA (PRKACA/B/G, PRKX), CK2α (CSNK2A1/2/3), AMPKα (PRKAA1/2), AKT1–3, GSK3A/B, ERK1/2, CAMK2A–G, PKC conventional, novel and atypical.
- In REPORT.md, state that PRKACA's 1,045 vs PRKACB's 11 labels reflect naming convention, not specificity.

### DL-9 — major — §4.1 substrate grouping
**Change.**
- Keep **40% identity / 50% coverage (cov-mode 0)**, but take **connected components** of an all-vs-all `mmseqs easy-search` (e ≤ 1e-3; edges with fident ≥0.4 and qcov, tcov ≥0.5) instead of `easy-cluster` set-cover representatives.
- Union these with identical-window edges, and also with **≥12/15-identical same-residue windows** (only 21 extra cross-group labelled pairs, so the cost is negligible).
- Do *not* lower to 30%, because a giant component then holds 24% of labels and breaks the 20% stratified test.
- Orthologous rodent sites inherit the human protein's group (already implied; state it).
- Add the test "no edge crosses TEST/DEV" on the edge list itself.

**Why** (Q7). A third of the qualifying pairs fall into different set-cover clusters. Components cost about 3 min on 4 cores, and the largest group is still small (451 proteins, 2.7% of labels).

### DL-10 — minor — §7 F1, §8 B1, §6 site-centric metric, §9 priming
**Change.**
- Recompute Cantley percentiles against a label-free reference: all Swiss-Prot canonical S/T (and Y) residues, or MS-observed U_obs sites. Store raw score, label-free percentile and published percentile. Within-kinase ranking is unaffected, but cross-kinase percentiles, which drive site-centric Hit@k, B1 as used across kinases, and the percentile features in the global model, currently depend on a PSP-curated reference (Q12).
- Compute priming / "flanking known phosphosite" features from MS-observed sites only, not PSP LT-only sites (66% of LT≥1 sites are KS-labelled, data_profile §2).

### DL-11 — minor — §2 Sites
**Change.**
- Rescue non-canonical/TrEMBL compendium (and CPTAC/Ochoa) sites by unique 15-mer to Swiss-Prot canonical (+10,696 compendium sites). Report the 376 multi-hit and 1,853 no-hit sites.
- For multi-location 15-mers, keep the label on the curated accession, put all copies in one group (DL-9) and mask the copies (DL-2).
- Define "human-native" test pairs explicitly as substrate organism = human with kinase organism ∈ {human}. The 797 human-substrate rows with a non-human kinase go into L2, not into test.
- Report the number of isoform-only KS labels dropped.

### DL-12 — minor — §7 provenance, F4, F8
**Change.** Add an explicit exclusion and stratification list to `docs/feature_provenance.md`:
- Exclude the COMPARTMENTS *knowledge* channel (it imports UniProt/GO CC, including IDA from KS papers), or ablate it.
- No GO BP/MF or GO-similarity features.
- No UniProt MOD_RES/"Site" features, and no ELM *instances* (class regexes only).
- PLMs: ESM-2 and ESM-C are fine, since they are unsupervised on UniRef. Exclude ESM3 (conditioned on function/keyword tokens) and any PTM-fine-tuned PLM.
- AlphaFold DB monomer features are fine. For any AF-Multimer/AF3/KinoPlex comparator, flag test pairs with PDB kinase–substrate complexes.
- F8 (Sugiyama, PhosphoNetworks protein arrays): report test metrics stratified by "pair is an F8 hit". Sugiyama's PMID cannot be checked against PSP, and PhosphoNetworks site calls overlap 618 PSP pairs.
- In LKO T2/T3, compute the "hard-negative pool" (annotated sites) without the held-out kinases' labels.

## 2. Things I checked that need no change
- Pooling with SGID plus sequence fallback (99.5% agreement with SGID) and the human-native-only test definition are sound.
- Holding rodent transfers out of the test set is correct.
- Excluding STRING/BioGRID biochemical/curated PPIs is well supported by `results/leakage/` (BioGRID biochemical activity is 23.8× enriched).
- BioPlex, HuRI (without Lit-BM), OpenCell, Human Cell Map and Buljan contain no KS curation, so they are leak-free in content. Their coverage bias is handled by DL-4.
- Cantley PSSMs are not fit to PSP, so they are leak-free for per-kinase ranking. Only the percentile reference is PSP-derived (DL-10).
- The PSP neighbour-lower-case density is a weak popularity proxy (AUROC 0.536), and the number of PSP sites per protein is uninformative (0.513).

## 3. Risks I could not resolve
1. The PSP KS file has no PMIDs, so the overlap between PSP in-vivo KS labels and compendium/CPTAC studies is only a lower bound (Q9). Fetching PSP per-pair references would need site-page scraping, which is probably not allowed under PSP's terms.
2. SIGNOR `DIRECT` and `TAX_ID` are not in OmniPath. Until PhosphoSIGNOR is supplied, tier A of the external test mixes direct and indirect relations.
3. HPRD/NCI-PID often assign kinases by motif or pathway context. Tier B may reward motif models, so treat it as a sanity check, not as evidence of superiority.
4. Whether Sugiyama 2019 or Newman 2013 were curated into PSP in-vitro rows is unknown. If they were, F8 becomes label-in-feature for those pairs.
5. The benchmarKIN tumour gold-standard definition needs confirming before DL-6.4 is finalised.
6. Paralog masking (DL-2) relies on a domain-identity threshold (80%). Near 80% (for example PKC isoforms at about 70%) it is a judgement call. Report sensitivity at 70% and 90%.
