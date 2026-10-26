# Data profile — provided data + UniProt reference (2026-10-26)

All numbers are produced by scripts (re-run order below); tables in `results/data_profile/` (aggregates only —
row-level PSP-derived tables go to `cache/`, gitignored, because PSP is CC BY-NC-SA).

```
python -I scripts/download_uniprot.py                 # external/uniprot/ + MANIFEST.json (UniProt 2026_03)
python scripts/build_compendium_matrix.py             # cache/compendium_wide.parquet (+ long parquet, info json)
python scripts/profile_psp.py                         # results/data_profile/psp_*.tsv, psp_summary.json; cache/psp_*.parquet
python scripts/verify_psp_sequence_mapping.py         # positive control for ortholog sequence mapping
python scripts/profile_compendium.py                  # compendium_*.tsv, kinhib_*.tsv, proteome_*.tsv, compendium_summary.json
python tests/test_io_psp.py                           # toy-data tests for the loader
```
Loaders: `src/ksp/io_psp.py` (`load_psp`, `filter_psp`, `map_kinases_to_human`, `map_substrates_to_human`,
`map_sites_by_sequence`, `check_windows_against_sequences`, `load_psp_sites`, `site_evidence_mask`,
`load_psp_regulatory`, `human_sites_from_site_dataset`), `src/ksp/uniprot.py` (`read_fasta`, `read_annotations`, `parse_pkinfam`).

Definitions used throughout: **site** = (human UniProt accession incl. isoform suffix, position); **KS pair** =
(human kinase gene symbol, site); **HMR** = kinase and substrate organism ∈ {human, mouse, rat}, both mapped to human.

---------------------------------------------------------------------------------------------------------
## 1. PhosphoSitePlus Kinase_Substrate_Dataset (release 051925)

| | value |
|---|---|
| rows | 23,904 (all 15-mers length 15, centre residue lower-case and equal to SUB_MOD_RSD in 100%) |
| human kinase × human substrate rows | 14,724 |
| non-S/T/Y rows (junk, dropped) | 17 (H 8, K 5, R 2, D 1, I 1 — NME1/2 His, histone/“K” sites, TP53 R249, …) → `psp_nonsty_rows.tsv` |
| exact duplicate rows | 0 by accession; 45 rows duplicate (kinase gene, organism, site) via different kinase isoform accessions |

**Organisms** (`psp_organism_crosstab.tsv`): kinase organism human 19,062 / mouse 2,805 / rat 1,419 / cow 405 / other 213;
substrate organism human 15,521 / mouse 4,887 / rat 2,585 / cow 381 / rabbit 200 / other 330. Cross-species rows are
common (e.g. human kinase on mouse/rat substrate — typical for in-vitro assays with rodent proteins).

**Evidence** (S/T/Y rows): in vitro only 10,453 (43.8%), both 6,749, in vivo only 6,685. Every row has ≥1 flag.
Human×human unique pairs: in vitro only 5,854, both 5,170, in vivo only 3,657. `psp_evidence.tsv`.

### 1.1 Kinases → human
* GENE is the kinase gene symbol (organism-specific case, e.g. `Dyrk2`); 10 rows have empty GENE and KIN_ORGANISM
  "no information" → mapped via the PSP kinase name (KINASE → human GENE, unambiguous names only).
* Non-human-kinase rows: 4,840. 4,761 (98.4%) match a human PSP kinase by upper-cased symbol; +37 rows (14 symbols:
  CAMK1D, EPHA3, EPHB6, MUSK, MYO3B, PHKG2, RIPK4, SPEG, STK36, TESK2, TIE1, TSSK2, TSSK3, WEE2) match a UniProt human
  protein kinase; +10 via kinase name → **4,808 / 4,840 (99.3%) mapped**. Unmapped: ALDOB, CAD, CKM, ENPP1, GUCY2E,
  NOTCH2, PHKB, PIKFYVE, viral (UL97, PYV…). ASSUMPTION: same symbol ⇒ 1:1 ortholog (not sequence-checked).

### 1.2 Substrates → human (site groups)
5,549 unique non-human substrate sites (S/T/Y). Mapping cascade (`psp_substrate_site_human_mapping.tsv`):

| method | sites |
|---|---|
| SITE_GRP_ID → human site present in the KS dataset | 1,898 |
| SITE_GRP_ID → human site present only in Phosphorylation_site_dataset | 1,890 |
| fallback: ortholog sequence (same symbol, UniProt Swiss-Prot canonical, ungapped 15-mer ≥11/15 identical, same residue) | 1,259 |
| unmapped | 502 |

326 non-human sites map to >1 human site (human isoforms in the same site group). 40 SGID-mapped rows change residue
(S↔T) between species. **Validation of the sequence fallback** (`psp_sequence_mapping_validation.tsv`): applied to the
3,751 SGID-mapped sites as positive control it returns a hit for 90.0% and agrees with SGID in **99.5%** of hits
(mouse 99.8%, rat 99.1%). Fallback-only mappings: 702 at 15/15 identity, 260 at 14, 140 at 13, 97 at 12, 60 at 11.

### 1.3 Gain from pooling species (`psp_pooling_gain.tsv`)

| scenario (all evidence) | unique human sites | KS pairs | kinases | kinases ≥5 | ≥20 | substrate proteins |
|---|---|---|---|---|---|---|
| human kinase × human substrate | 10,212 | 14,681 | 420 | 275 | 126 | 2,974 |
| any kinase × human substrate | 10,431 | 15,100 | 426 | 275 | 126 | 3,008 |
| HMR, SGID→KS-human sites only | 10,529 | 16,428 | 432 | 288 | 132 | 3,004 |
| HMR, SGID (KS + site dataset) | 12,136 | 18,452 | 437 | 301 | 143 | 3,355 |
| **HMR, SGID + sequence fallback** | **13,206 (+29%)** | **19,719 (+34%)** | **442** | **304** | **147** | 3,542 |
| all species, all mappings | 13,606 (+33%) | 20,473 (+39%) | 445 | 308 | 149 | 3,613 |

In-vivo-supported subsets: H×H 6,582 sites / 8,827 pairs; HMR 8,636 / 11,998. Using only SGIDs whose human member is in
the KS file itself adds almost no new *sites* (by construction) — the Phosphorylation_site_dataset SGIDs and the
sequence fallback are what add sites. 1,580 HMR pairs are supported in >1 species. Among HMR pairs, 4,712 come only
from non-human rows (both 1,375 / in vitro only 1,958 / in vivo only 1,379).

### 1.4 Isoforms
Substrate isoform accessions: 1,238 rows (5.2%), 265 distinct (human 1,013 rows / 193 accs; mouse 118; rat 104).
Most frequent suffixes -2 (496 rows), -8 (321), -3 (146), -1 (98; "-1" = explicitly the canonical in some entries).
Kinase isoform accessions: 239 rows (e.g. P05480-1 mouse SRC 96, P05771-2 PRKCB-II 72). The loader keeps
`sub_isoform`/`kin_isoform` separately; positions on isoform accessions are isoform coordinates.

### 1.5 Kinase distribution and concentration (`psp_kinase_site_counts.tsv`, `psp_kinase_thresholds.tsv`)

| pairs counted on | kinases ≥1 | ≥5 | ≥10 | ≥20 | ≥50 | ≥100 | Gini | top-10 share | top-20 | top-50 |
|---|---|---|---|---|---|---|---|---|---|---|
| human×human | 420 | 275 | 201 | 126 | 63 | 28 | 0.76 | 34.8% | 49.8% | 69.6% |
| HMR pooled | 442 | 304 | 217 | 147 | 76 | 44 | 0.78 | 35.8% | 51.2% | 70.5% |
| HMR, in vivo | 410 | 246 | 166 | 112 | 53 | 22 | 0.76 | 34.9% | 51.3% | 70.5% |

Top 10 (HMR): PRKACA 1,054, PRKCA 804, SRC 800, CDK1 796, MAPK1 778, CSNK2A1 775, MAPK3 567, GSK3B 522, CDK2 497, AKT1 472.
Pooling helps the tail (≥5 pairs: 275 → 304) but does not reduce concentration.

**Kinase groups** (UniProt pkinfam, `psp_kinase_group_counts.tsv`; kinases with ≥1 HMR pair / pkinfam total, pairs):
CMGC 55/62, 5,238 · AGC 52/57, 4,323 · Tyr 79/90, 2,998 · Other 48/81, 2,410 · CAMK 57/74, 1,973 · PIKK (ATM/ATR/MTOR/PRKDC/SMG1 + PI3Ks) 10/10, 869 ·
STE 50/56, 776 · CK1 12/12, 454 · TKL 23/34, 315 · NEK 10/11, 102 · Alpha 5/6, 81 · PDK/BCKDK 5/5, 29 · RIO 2/3 · ADCK 1/5 · RGC 0/5 · FAST 0/1.
**411 of 513 UniProt human protein kinase genes** (pkinfam ∪ KW-0418 with a "Protein kinase" domain) have ≥1 HMR pair.

**Non-protein-kinase "kinases"** (`psp_non_protein_kinase_entries.tsv`): 31 gene entries, 143 HMR pairs (0.7%):
fusions BCR/ABL (21), NPM/ALK (6); regulatory subunits CSNK2B (16), PRKAB1 (12), PRKAG2 (9), PIK3R1 (2); metabolic
kinases PKM (15), PGK1 (7), CKB, PCK1, KHK, PFKP, PFKFB3, PMVK, AK2, PIP4K2A; histidine kinase NME1; others JMJD6,
GTF2F1, DCAF1, BLVRA, BCR, TGM2, HSPA5, BRD4, BAZ1B, ENPP3, CERT1, PHKA1, STK19 (now WHR1 in UniProt — shown not to be a
kinase). FAM20C (14 pairs) is a genuine secretory-pathway Ser kinase but absent from pkinfam → keep it with a flag.
Lipid kinases in pkinfam's PI3/PI4 group with PSP protein substrates: PIK3CA 6, PIK3CG 5, PIK3C2A 5, PIK3CD 3, PIK3CB 1.
Recommendation: drop fusions/subunits/metabolic enzymes; keep PIKKs; decide explicitly on PI3Ks and FAM20C.

### 1.6 Multi-site / multi-label structure (HMR)
13,206 sites on 3,542 proteins; 2,217 proteins have >1 site, 288 have ≥10 (MAPT 86, LRRK2 53, TRPM7 47, TTK 43, YAP1 43);
the top 1% of proteins hold 8.3% of sites → protein-grouped splits are essential. 3,410 sites (25.8%) have >1 kinase
(max 24: AKT1 S473, H3 S10; GSK3B S9 23) → multi-label, and "negatives" for one kinase are frequently positives for another.

### 1.7 Autophosphorylation
1,283 rows have kinase gene = substrate gene (same organism); HMR: 1,104 pairs (5.6%) across 241 kinases
(LRRK2 52, TRPM7 47, CAMK2B 33, TTK 32, SRC 25, MELK 21, JAK2 19, BUB1 19). These may be cis or trans; motif-based
models should treat them with care (flag column `autophos`).

### 1.8 Windows: co-phosphorylated residues
71% of rows have ≥1 additional lower-case (modified) residue in the ±7 window (1: 6,216; 2: 4,360; ≥3: 6,370 rows)
— the lower-case marks reflect *any* known phosphosite, not co-occurrence on the same molecule; usable as a
"priming-site known" feature but beware leakage of site-annotation density.

### 1.9 15-mers vs UniProt 2026_03 (`psp_window_vs_uniprot.tsv`)
KS human sites (10,433 unique): **99.17% exact** window match at the stated position (canonical 9,810/9,892; isoform 536/541).
66 centre mismatches, all relocatable (exact 15-mer found elsewhere): shift +1 for 56 (histones H2A/H2B/H3/H4 — PSP numbers
without the initiator Met — and SLC3A2), −101 for 5, other drifts single. 20 sites lack a sequence: PSP uses VAR_ IDs
(O60504_VAR_055019, Q12809_VAR_014388), RefSeq/GenBank (NP_001184222, AAC50053 PAX3/FOXO1) or accessions no longer
in UP000005640 as given (P62988 ubiquitin, Q15300, Q17RH5, Q17R91, P18433-2). 1 position out of range.
Phosphorylation_site_dataset human sites: 99.23% exact (909 centre mismatches, 898 relocatable; 793 no sequence).

---------------------------------------------------------------------------------------------------------
## 2. PhosphoSitePlus Phosphorylation_site_dataset and Regulatory_sites (release 051925)

**Phosphorylation_site_dataset**: 379,676 rows (human 240,731, mouse 105,303, rat 31,745); all MOD_RSD "-p"; 50 non-STY.
Human S/T/Y: **240,698 sites** on 17,643 proteins (S 142,408 / T 58,879 / Y 39,411); 8,842 on isoform accessions;
11,923 flagged ambiguous. Every human site has ≥1 evidence count. `psp_sitedb_human_universe.tsv`:

| evidence definition | sites | S | T | Y | with KS (HMR) | dark fraction |
|---|---|---|---|---|---|---|
| all | 240,698 | 142,408 | 58,879 | 39,411 | 12,132 | 95.0% |
| LT_LIT ≥1 | 15,432 | 9,483 | 3,167 | 2,782 | 10,242 | 33.6% |
| MS_LIT+MS_CST ≥2 | 106,464 | 68,315 | 20,690 | 17,459 | 6,446 | 93.9% |
| **LT ≥1 or MS ≥2 (default `site_evidence_mask`)** | **114,275** | 72,822 | 22,588 | 18,865 | 11,532 | **89.9%** |
| LT ≥1 or MS_LIT ≥2 (no CST counts) | 81,005 | 59,716 | 14,153 | 7,136 | 11,293 | 86.1% |
| LT ≥1 or MS ≥5 | 62,176 | 39,582 | 10,707 | 11,887 | 11,102 | 82.1% |
| single MS record, no LT | 126,423 | 69,586 | 36,291 | 20,546 | 600 | 99.5% |

So **~90% of reasonably supported human phosphosites (~103k) have no kinase annotation** — the coverage target.
52.5% of PSP human sites rest on a single MS record. MS_CST counts strongly inflate Y sites (CST pY antibody enrichment):
LT≥1 or MS_LIT≥2 keeps 7.1k Y vs 18.9k with CST counts. Of 13,206 HMR KS sites, 12,132 are in the site dataset
(the rest: isoform coordinates, sequence-fallback positions, non-phospho-catalogued sites).

**Regulatory_sites**: 21,011 rows, columns GENE, PROTEIN, PROT_TYPE, ACC_ID, GENE_ID, HU_CHR_LOC, ORGANISM, MOD_RSD,
SITE_GRP_ID, SITE_+/-7_AA, DOMAIN, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT, PMIDs, LT_LIT, MS_LIT,
MS_CST, NOTES. All PTM types (p 16,321, ac 1,493, ub 1,425, sm 935, me/m1-3 645, gl 115, …). Human phospho:
**10,580 sites**, ON_FUNCTION 9,364, ON_PROCESS 6,706, ON_PROT_INTERACT 4,751, ON_OTHER_INTERACT 695; top functions:
molecular association (5,107), intracellular localization (2,705), enzymatic activity induced (1,444).
**75.1% of human regulatory phosphosites have a KS annotation** (vs 5% of all sites) — a strong study-bias signal:
"functional" and "kinase-annotated" are nearly the same sites. Use as a stratification/feature, not as labels.

---------------------------------------------------------------------------------------------------------
## 3. Perturbation phosphoproteomics compendium (data/01–03)

**01 conditions**: 1,937 conditions (1,329 distinct Condition strings), **178 PMIDs**, 204 samples/cell types (A549 165,
ARH-77/Ramos/SU-DHL-4 102 each, MCF7 87, RPE1 66). Source: "This study" 1,437, qPTM 461, Ochoa2016 39.
One study dominates: **PMID 36926954 = 790 conditions (40.8%)** (drug dose–response series, TMT10/11, Rituximab/Dasatinib/
MK2206/…); then Hijazi 31959955 (183), 32034161 (66), 34663829 (65). Labeling: LFQ 544, TMT10 477, TMT11 313, SILAC 292,
TMT 113, label-free 82, iTRAQ 46, other/NA 64. Perturbation type (keyword heuristic on Condition text, approximate —
`classify_condition` in `scripts/profile_compendium.py`): drug/inhibitor 1,217, ligand/cytokine/GF 272, baseline
sample contrast 111, genetic 85, physical/genotoxic stress 73, infection 43, other/unclear 136.

**02 sites**: 189,962 peptide_UIDs, all of the form `ACC_POS_M1`. **The `_M1` suffix is constant (all 189,962 = M1)**;
each (ACC, POS) occurs once and the windows carry no lower-case marks, so the data cannot tell us whether M = phospho
multiplicity (MaxQuant "___1") or a running index. Most likely it is multiplicity 1 / a de-duplicated site key; either way
it carries no information and sites can be keyed as (ACC, POS). 19,233 accessions, no isoform suffixes; S 133,604 /
T 42,354 / Y 14,004; 4,098 windows padded with `_`. Versus UniProt 2026_03: 98.81% exact window match; 1,204 centre
mismatches (sequence-version drift), 868 sites on 193 accessions absent from UP000005640 (obsolete/merged);
**2,887 accessions are TrEMBL** (16,153 Swiss-Prot).

**03 quantification**: 17,032,346 log2FC values, no NaN, no duplicate (site, condition) pairs; all IDs resolve to 01/02.
Wide matrix 189,962 × 1,937 float32, **density 4.63%**. Build (`scripts/build_compendium_matrix.py`, pyarrow csv,
16 threads): **11.6 s wall, 5.3 GB peak RSS**; `cache/compendium_wide.parquet` 104 MB (+ `compendium_long.parquet`).
Per-site coverage: median 11 conditions (IQR 4–52; 57,579 sites seen in <5 conditions); 102,467 sites ≥10, 48,689 ≥50,
32,758 ≥100, 9,743 ≥500. Per-condition: median 8,560 sites (min 1,981, max 29,577). Values: median 0.00, 1%/99%
quantiles −2.44/+2.52, 12% |log2FC|>1. Per-condition SD median 0.63 (5–95%: 0.23–1.42) → per-study scaling needed.

**PSP overlap** (`compendium_psp_overlap.tsv`):

| PSP set | sites | in compendium | ≥10 cond. | ≥100 cond. | median cond. | pairs with site in comp. | kinases ≥1 / ≥5 such pairs |
|---|---|---|---|---|---|---|---|
| human × human KS | 10,212 | 4,866 (47.6%) | 3,958 | 2,378 | 92 | 7,500 / 14,681 | 364 / 201 |
| HMR pooled KS | 13,206 | 6,079 (46.0%) | 4,913 | 2,908 | 86 | 9,842 / 19,719 | 384 / 229 |

PSP 15-mer = compendium 15-mer for 4,537/4,538 matched human sites (same coordinate system). Site-dataset overlap
(`compendium_psp_sitedb_overlap.tsv`): 104,160 compendium sites (54.8%) are PSP sites; 66.9% of PSP "LT≥1 or MS≥2"
sites are quantified (76,407), incl. 70,591 dark ones. Only 46% of LT-only-supported (LT≥1) sites are in the compendium.
85,802 compendium sites are not PSP sites by coordinate (10,741 of them do match a PSP 15-mer → coordinate/isoform
differences); the rest are sparsely measured (median 5 conditions vs 30 for PSP-known) — likely lower-confidence
localisations.

**Hijazi kinase-inhibitor matrix** (`kinhib_*.tsv`): 420 kinases × 60 inhibitors, binary, 762 links; **307 kinases have
≥1 inhibitor** (113 none; 109 exactly one); 286 of the 512 pkinfam human kinases are targeted. Inhibitor promiscuity:
PF3758309 98 targets, LY2835219 83, AZD1480 52; 15 inhibitors hit 1 kinase; U73122 (a PLC inhibitor) has 0.
Row names use legacy symbols (ADRBK1/2, CDC2L2, GSG2, ICK, KIAA0999, MST4, PAK7, ZAK, CABC1, SGK110, MGC42105, Q6ZSR9) and
include lipid kinases (PI4KB, PIK3C2B/G, PIK3C3, PIP4K2A–C, PIP5K1A/C, PIKFYVE) and BCR/ABL → a symbol-synonym map is
needed before joining with PSP. **Hijazi conditions** = PMID 31959955: 183 conditions = 61 treatments × 3 cell lines
(HL60, MCF7, NTERA2), all 60 min vs cell-line DMSO; Condition = `<cell>_<inhibitor>`; 60 treatments match kinhib columns
exactly, Edelfosine (lipid analogue) is the only treatment without a column (`kinhib_hijazi_conditions.tsv`).
17 kinhib inhibitors also occur in other studies (e.g. Dasatinib, MK2206, AZD8055 dose series in 36926954; KN62/KN93,
Torin1, PF4708671 in 26060313; CHIR99021, SP600125, JNK-IN-8 in 32034161) → `kinhib_condition_matches.tsv`.

---------------------------------------------------------------------------------------------------------
## 4. Human proteome reference / prediction universe (UniProt 2026_03, 02-Sep-2026)

Files, URLs, SHA-256: `DATA.md` and `external/uniprot/MANIFEST.json`. `proteome_sty_universe.tsv`:

| set | sequences | S/T/Y residues |
|---|---|---|
| Swiss-Prot canonical | 20,416 | **1,870,491** (S 953,038 / T 613,855 / Y 303,598) |
| Swiss-Prot canonical + isoforms | 42,547 | 4,046,421 (1,873,184 distinct 15-mers) |
| UP000005640 all canonical (incl. 127,104 TrEMBL) | 147,520 | 11,196,211 |
| UP000005640 incl. isoforms | 169,651 | 13,372,141 |

Known phosphosites as a share of the Swiss-Prot canonical S/T/Y universe: compendium 178,381 (9.5%); PSP human 240,698 (~12.9%);
KS-annotated HMR sites 13,206 (0.71%). Human protein kinase list: pkinfam 512 entries; KW-0418 reviewed 625
(459 with a "Protein kinase" domain); union by gene 513.

---------------------------------------------------------------------------------------------------------
## 5. Caveats / decisions for the PI
1. **Use the HMR pooled set with SGID → site-dataset mapping + sequence fallback** as the label universe candidate
   (+29% sites, +34% pairs vs human-only); keep `map_method` and organism as columns so benchmarks can be stratified
   (test sets should probably be human-row-supported pairs only).
2. In-vitro-only annotations are 44% of rows; decide whether test labels require in vivo, or report both.
3. Drop the 17 non-STY rows and the 31 non-protein-kinase entries (or flag); decide on PI3Ks / FAM20C explicitly.
4. Concentration is extreme (Gini 0.78; top 20 kinases = 51% of pairs) — per-kinase macro metrics are mandatory.
5. 26% of sites are multi-kinase; 63% of substrate proteins carry >1 site — split by protein (and homology), never by pair.
6. Regulatory and KS-annotated sites are the same well-studied ~10k sites (75% overlap) — strong ascertainment bias;
   "negatives" sampled from unannotated sites are mostly unlabelled, not negative.
7. Compendium covers only ~46% of KS sites (median ~86 conditions when present) → phosphoproteomics features will be
   missing for half of the training sites; models must handle missingness without it becoming a label proxy
   (annotated sites are much better covered: median 30 vs 5 conditions for PSP-known vs unknown compendium sites).
8. One study (36926954) is 41% of conditions; per-study normalisation/weighting needed. `_M1` carries no information.
9. Kinhib uses legacy gene symbols — map through UniProt synonyms before use.
10. UniProt REST streams are "current release" snapshots (2026_03); checksums pin what we used. 0.8% of PSP and 1.2%
    of compendium windows don't match UniProt at the stated position (Met-numbering, version drift, VAR_/RefSeq IDs).
