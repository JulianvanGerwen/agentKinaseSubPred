# Canonical IDs, kinase list, site universe and KS labels (2026-10-26)

Scope: build the identifiers, kinase list, observed-site universe (U_obs) and kinase–substrate (KS) label tables.
DESIGN §2, §3 and §5 are implemented with parameters, because the design is still under review. All numbers below
come from the scripts and were re-generated on 2026-10-26 (UniProt 2026_03, PSP 051925).

```
/workspace/.venv/bin/python -I scripts/download_uniprot.py      # adds ptm_tsv/ + kinase_detail/ (new, see DATA.md)
/workspace/.venv/bin/python -I scripts/build_kinase_list.py     # ~5 s  -> results/labels/kinase_list*.tsv, kinase_synonyms.tsv
/workspace/.venv/bin/python -I scripts/build_site_universe.py   # ~95 s (first run: Ochoa xlsx 40 s, cached) -> cache/sites/, results/sites/
/workspace/.venv/bin/python -I scripts/build_labels.py          # ~18 s -> cache/labels/, results/labels/label_*.tsv
/workspace/.venv/bin/python tests/test_ids.py; tests/test_labels_universe.py; tests/test_io_psp.py
```
Code: `src/ksp/ids.py` (store, windows, mapping, site_id), `src/ksp/kinases.py` (kinase list, synonyms,
`KinaseResolver`), `src/ksp/sites.py` (source loaders, `default_store()`), `src/ksp/labels.py` (PSP/OmniPath/UniProt
pairs). Manual decisions are in `config/kinase_curation.tsv`, which is data and not code. Each row there has a rationale.
Each run uses ≤16 cores, and all steps are single-threaded except pandas I/O.
Other changes: `openpyxl` was added to the venv and requirements.txt so the Ochoa/Cantley xlsx files can be read.
The HGNC complete set was already recorded in DATA.md (`external/hgnc/`).

## 1. Identifiers (`src/ksp/ids.py`)
- **Store:** Swiss-Prot canonical human sequences (20,416, from `UP000005640_reviewed_isoforms.fasta.gz`, isoform '').
  - Gene index: FASTA `GN=`, plus UniProt gene synonyms and HGNC previous symbols (+23,490 links). These never override a primary symbol.
  - Accession→gene index: all UP000005640 entries, including TrEMBL.
- **Windows:** `window(seq,pos,flank)`, padded with `_` (±7 gives a 15-mer, ±15 gives a 31-mer).
- **`site_id = ACC_RESPOS`:** `make_site_id` and `parse_site_id` validate the canonical accession (no isoform suffix) and that the residue is S/T/Y.
- **`map_window(acc,pos,15-mer,res,gene)`** maps a source site to the canonical:
  1. **Candidates:** the base accession if it is Swiss-Prot canonical. Otherwise the Swiss-Prot canonical(s) of the source gene, taken from the explicit gene or from the UniProt gene of a TrEMBL/obsolete accession.
  2. `exact_pos`: the window matches at the stated position.
  3. `shifted`: the window matches at exactly one other position.
  4. Otherwise unmapped, with reason `bad_window` / `no_candidate` / `no_match` / `multiple_matches`.
  - `_`/`X` in the source window are wildcards. This lets a window match when the source omits the initiator Met.
  - The canonical centre residue must equal the source residue.
- **`map_position(acc,pos,res,expected_length)`** is for sources without a window (Ochoa, OmniPath, UniProt MOD_RES). It checks the residue at the stated position on the canonical and, optionally, that the canonical length equals the source length.

## 2. Kinase list (`results/labels/kinase_list.tsv`, 568 rows; 503 `include=True`)
**Candidates** come from four sources:
- UniProt pkinfam (512);
- reviewed KW-0418 entries with a "Protein kinase" domain;
- KinHub (536; gives the Manning group/family/subfamily);
- curation `add` (FAM20A/B/C).

Every PSP kinase entry is also added as a row, so each one has a recorded decision.

**Exclusions** (`exclusion_category`, rationale per gene in `config/kinase_curation.tsv`):

| category | n | examples |
|---|---|---|
| disputed activity | 20 | BRD2/3/4/T, TAF1/1L, BAZ1A/1B, TRIM24/28/33, BCR, ABR, BLVRA, CERT1, GTF2F1, JMJD6, DCAF1, TGM2, HSPA5 |
| metabolic enzyme | 12 | PKM, PGK1, CKB, CKM, PCK1, KHK, PFKP, PFKFB3, PMVK, AK2, ALDOB, CAD |
| lipid kinase | 8 | PIK3CA/B/D/G, PIK3C2A, PIK3C3, PIKFYVE, PIP4K2A (in PSP / pkinfam) |
| regulatory subunit | 6 | CSNK2B, PRKAB1, PRKAG2, PIK3R1, PHKA1, PHKB |
| guanylate-cyclase pseudokinase | 6 | NPR1/2, GUCY2C/D/F, GUCY2E (rodent) |
| not a protein kinase | 5 | WHR1 (=STK19), NME1, ENPP1/3, NOTCH2 |
| glycan kinase | 2 | FAM20B (xylose kinase), POMK (mannose kinase) |
| fusion | 2 | BCR/ABL, NPM/ALK |
| pseudogene | 1 | PDPK2P |
| not in protein-kinase sources | 3 | viral UL97, PYV0001/PYV0094 |

**Included atypical kinases:**
- PIKKs: ATM, ATR, MTOR, PRKDC, SMG1, and TRRAP, which is a pseudokinase.
- FAM20C, and FAM20A, which is a pseudokinase.
- PDK1–4, BCKDK, ALPK1–3, EEF2K, TRPM6/7, RIOK1–3.
- ADCK1/2/5, COQ8A/B and FASTK, all flagged `activity_uncertain` (6).

Pseudokinases and uncertain kinases stay included by default: `--exclude-pseudokinases` and `--exclude-uncertain` turn them off.

- **`is_pseudokinase`** (38 included) is set by any of three rules:
  - a curation flag (FAM20A, TRRAP, ERBB3, KSR2);
  - the entry lacks UniProt KW-0418 (UniProt does not give "Kinase" to dead family members, e.g. STRADA/B, TRIB1–3, ILK, PEAK1, ULK4, RYK, PTK7);
  - UniProt CAUTION text says so (KSR1). Function text is not used because it mentions kinase-dead mutants.

  The flag is evaluated only for protein-kinase-family members.
- **`residue_class`** (DESIGN §5; the rule used is in `residue_class_rule`):
  - curation override: dual-specificity MAP2K1–7, DYRK1A/1B/2/3/4, TTK, CLK1–4, TESK1/2, PKMYT1 → STY (20); WEE1/WEE2 → Y;
  - otherwise Manning group TK (or pkinfam 'Tyr') → Y (90);
  - otherwise ST.

  Included: ST 391, Y 92, STY 20.
- **Kinase domain:** taken from UniProt DOMAIN features (Protein kinase, Alpha-type, PI3K/PI4K catalytic, Histidine kinase).
  - Primary domain = the one containing ACT_SITE; otherwise the last one. For JAKs this gives JH1, not JH2.
  - 498 of 503 included kinases have a domain. FAM20A/C, FASTK and others have no domain feature.
- **Manning group:** from KinHub (497 included). Otherwise from curation (FAM20*, CSNK2A3, CDK11A, PLK5, PEAK3), otherwise from pkinfam.
- **Synonyms** (`results/labels/kinase_synonyms.tsv`, 2,503 rows): own symbol; HGNC prev/alias; UniProt gene names; KinHub xName/Manning name; UniProt entry name; Cantley matrix names; PSP kinase names.
  - `use_for_join` = own symbol, OR the synonym points to exactly one kinase and is not the approved HGNC symbol of another gene.
  - Example: `PDK1` stays PDK1 and is never PDPK1; `PKD1` and `MST1` are not joinable.
  - `KinaseResolver` checks in this order: accession, own symbol, curated family names (never a gene), unambiguous synonym.
- **Join checks** (`kinase_list_summary.tsv`):
  - Hijazi kinhib rows: 411/420 resolved, including ADRBK1→GRK2, ZAK→MAP3K20, MST4→STK26, KIAA0999→SIK3, CDC2L2, GSG2, ICK, PAK7, CABC1, MGC42105, SgK110.
  - The 9 unresolved kinhib rows are lipid kinases outside the list (PI4KB, PIK3C2B/G, PIP4K2B/C, PIP5K1A/C), PRKY (not in the reviewed proteome) and Q6ZSR9 (an uncharacterised ORF).
  - Cantley matrices: 399/408 resolved by accession; the 9 unresolved are non-human (worm/yeast).

## 3. Site universe (`cache/sites/universe.parquet`; summary `results/sites/universe_summary.tsv`, `universe_overlap.tsv`)

| source | source sites | exact_pos | shifted | unmapped (main reasons) | unique canonical sites |
|---|---|---|---|---|---|
| PSP Phosphorylation_site_dataset (human, S/T/Y, -p) | 240,699 | 96.2% | 1.7% | 2.0% (no_match 4,516, of which 4,059 are isoform-only windows; no_candidate 322) | 232,789 |
| PSP KS human substrate sites (flag only) | 10,433 | 95.1% | 4.4% | 0.6% (61) | 10,201 |
| Compendium (02_phosphosite_information) | 189,962 | 94.9% | 3.7% | 1.5% (no_match 1,703; no_candidate 1,024, mostly obsolete or TrEMBL accessions without a Swiss-Prot gene) | 179,238 |
| CPTAC BCM v1.2 (10 cohorts, tumour+normal; ENSG→HGNC→UniProt, flank) | 126,946 | 97.0% | 2.8% | 0.3% (349) | 125,754 |
| Ochoa 2020 (no window: position+residue, length must match) | 116,258 | 95.6% (pos_residue) | – | 4.4% (length_changed 3,803; no_candidate 1,200) | 111,141 |

- **Universe size:** U_obs = 350,736 canonical sites on 18,702 proteins (S 217,233 / T 87,608 / Y 45,895).
  - 4 additional sites come only from PSP KS rows (`in_U_obs=False`).
  - 152,396 sites are seen in ≥2 sources.
- **Sites unique to one source:** PSP 103,773, compendium 50,225, CPTAC 28,002, Ochoa 12,052.
- **Largest source combinations:** PSP only 103,773; compendium only 50,225; all four 48,440 (+3,101 with KS).
- **Columns:**
  - `site_id, acc, pos, res, gene, window15, window31`;
  - flags: `src_{psp,psp_ks,compendium,cptac,ochoa}`, `n_sources`, `in_U_obs`;
  - per-source mapping method;
  - metadata: PSP `psp_lt_lit/ms_lit/ms_cst` (max over source rows), `psp_site_grp_id`, `psp_ambiguous`; `comp_n_conditions`; `cptac_cohorts`/`n_cohorts`; Ochoa `functional_score`/`biological_samples`/`best_loc_prob`.
  - **The metadata are popularity and ascertainment proxies. They must NOT be used as features** (DESIGN §5).
- **Mapping log:** `cache/sites/mapping_log.parquet` holds every source site with its method and failure reason.

## 4. Labels (`cache/labels/ks_pairs.parquet`, one row per (kinase, site_id, provenance); `ks_unmapped.parquet`)

**PSP** has 23,904 rows. The 17 non-STY rows are dropped.

**Kinase resolution** (all rows resolved):
- human rows go through the accession (19,045);
- non-human rows go through symbol/synonym. This assumes the same symbol means the human ortholog.
- 10 rows with an empty GENE field are resolved through the PSP kinase name.

**Substrate resolution, per row:**

| method | rows |
|---|---|
| human 15-mer, exact at stated position | 14,502 |
| human 15-mer, shifted (incl. the histone Met offset and isoform coordinates) | 928 |
| SITE_GRP_ID → human sites of the group (PSP site dataset + KS human rows, mapped by 15-mer) | 6,098 |
| ortholog fallback: same symbol, ungapped 15-mer, ≥11/15 identical, same residue (`--min-identity`) | 1,601 |
| unmapped | 758 (680 non-human without a human SGID member or ortholog hit; 75 human windows not in the canonical) |

After canonicalisation only 1 row maps to >1 human site.

**Flags** (aggregated over all PSP rows of a pair):
- `native_human`: ≥1 row where kinase and substrate organism are both human;
- `hmr`; `sub_native_human`;
- `species_origin` (e.g. `human>human;rat>rat`);
- `in_vivo`/`in_vitro` (any row), plus `_native` and `_hmr` versions, and `evidence_native`/`evidence_hmr` (both / in_vivo_only / in_vitro_only);
- `autophos`: kinase accession == substrate accession. PSP's same-gene/same-organism flag is kept as `psp_autophos_flag`;
- `kinase_included`, `kinase_residue_class`, `residue_in_class`, `in_universe`;
- `res_changed`: the SGID transfer changes S↔T;
- `sgid_multi`, `sub_map_methods`, `n_rows`.

**Variants** (each also requires `kinase_included & residue_in_class`):
- L1 = native rows;
- L2 = kinase and substrate organism ∈ {human, mouse, rat};
- L3 = all organisms;
- L2_vitro = HMR rows with in vitro evidence;
- L2_noauto = L2 minus autophosphorylation;
- EXT = external rows (OmniPath/UniProt);
- EXT_not_in_psp = an EXT pair that is absent from all PSP pairs.

The PSP pairs dropped by these requirements:
- 275 pairs have an excluded kinase.
- 214 pairs have a residue outside the kinase's class: 120 Y sites on ST kinases, 94 S/T sites on Y kinases (e.g. SRC 19, ABL1 16, AURKB 13 Y).

The ST/Y class rule removes all labels of 5 kinases: LMTK1, LMTK2, STYK1, ERN1, MAP3K21. **LMTK1–3 are in the Manning TK group but are reported as Ser/Thr kinases**, so the PI should review them.

| variant | pairs | sites | proteins | kinases | ≥5 | ≥20 | ≥50 | top-10 share |
|---|---|---|---|---|---|---|---|---|
| L1 | 14,197 | 9,850 | 2,931 | 378 | 265 | 120 | 62 | 35.2% |
| **L2 (default)** | **18,638** | **12,563** | 3,507 | **401** | **290** | **143** | 73 | 36.1% |
| L3 | 19,301 | 12,934 | 3,577 | 404 | 292 | 143 | 75 | 36.7% |
| L2_vitro | 14,182 | 10,033 | 2,961 | 379 | 255 | 118 | 62 | 37.4% |
| L2_noauto | 17,629 | 11,776 | 3,447 | 386 | 263 | 130 | 68 | 38.0% |
| L1 in vivo only / in vitro only / both | 3,488 / 5,641 / 5,068 | | | 273 / 298 / 305 | | | | |
| L2 − L1 (rodent-only pairs) | 4,441 | 3,444 | 1,385 | 259 | 130 | 48 | 19 | |
| EXT OmniPath (curated non-PSP) | 8,654 | 5,072 | 1,805 | 331 | 202 | 89 | 37 | |
| EXT UniProt (ECO:0000269 "by X") | 2,926 | 2,412 | 1,005 | 285 | 141 | 37 | 9 | |
| EXT any / of which not in PSP | 10,387 / 3,788 | 6,368 / 2,572 | | 371 / 291 | 243 / 149 | 107 / 58 | | |
| L2 + EXT | 22,574 | 14,023 | 3,741 | 417 | 310 | 166 | 92 | 33.8% |

Of the included kinases, 401 of 503 have ≥1 L2 label. Per-kinase counts are in `results/labels/kinase_label_counts.tsv`, and the full funnel is in `results/labels/label_funnel.tsv`.

**OmniPath enzsub** (`external/omnipath`, 39,037 phosphorylation rows):
- **Kept:** rows with ≥1 curated non-PSP source: SIGNOR, HPRD, phosphoELM, and NCI-PID / Reactome / SIGNOR via ProtMapper.
- **Dropped** (27,478): rows whose sources are only PhosphoSite*, ProtMapper, KEA, MIMP, text mining (Sparser, REACH, RLIMS-P, BEL) and/or PhosphoNetworks.
- dbPTM and Li2012 are aggregators that are off by default (`--omnipath-include-aggregators`).
- Columns: `ext_sources`, `ext_refs` (PMIDs of the curated sources only), `omnipath_lists_psp`, `omnipath_all_sources`, and `in_psp_ks` (our own PSP match).
- **Sites:** OmniPath has no windows, so positions are checked against the canonical residue. 548 residue mismatches are dropped.
- **Unresolved enzymes:** 2,038 rows have enzymes that are not kinases in our list, mostly complex members such as cyclins, AMPK β/γ, CAK/TFIIH subunits and phosphatases.
- 5,687 of 8,654 EXT pairs are also in PSP.

**UniProt** (`ptm_tsv`, MOD_RES "Phospho…; by X"):
- Only rows with ECO:0000269 evidence are used: 2,953 of 3,864 "by" annotations.
- Name resolution: whole name, then a normalised name (GSK3-beta→GSK3B), then "/" parts that all resolve to the same gene (PKB/AKT1→AKT1, STK4/MST1→STK4).
- "autocatalysis" maps to the entry's own gene when that gene is in the list.
- **Family names** (PKA, PKC, CK1, CK2, AMPK, MAPK, GSK3, CaMK2, PKB, SRC-type…; 56 names curated as `family_name`) are dropped: 585 annotations.
- Positions are current-release canonical, so they are mapped directly with a residue check.
- 1,999 of 2,926 EXT pairs are also in PSP.

## 5. Caveats and open decisions
1. **Orthology by symbol** for kinases (as in `io_psp`); the substrate ortholog fallback is ungapped (validated earlier at 99.5% agreement with SGID).
2. **Isoform-only sites are dropped.** In PSP this is 4,059 site rows and 54 KS human rows whose windows are not in the canonical. Isoform-coordinate sites whose window is in the canonical are kept and relocated.
3. **Ochoa 2020 has no windows.** It is mapped by position + residue on the 2026_03 canonical, and only when the length equals Ochoa's `prot_length`. This is a proxy for an unchanged sequence version; 3.3% are lost to length changes.
4. **CPTAC:** only the BCM-harmonised v1.2 files carry flanking 15-mers. The PDC Broad (RefSeq IDs, no flank) and UMich (peptide-level) tables are not used, so the CPTAC source is the BCM site list.
5. **OmniPath "not in PSP" ≠ independent literature.** SIGNOR and HPRD curate the same papers. Using EXT as an "independent" test needs its own leakage check by PMID.
6. **The residue-class rule is a pure group rule.** TK-group Ser/Thr kinases (LMTK1–3, STYK1?) lose all labels. Y labels on ST kinases (e.g. AURKB H3/Y, CSNK2A1 Y) are dropped; this is mostly autophosphorylation or noise. In-vivo-only rows are only 21% of the out-of-class rows.
7. **Rodent-transferred sites outside U_obs:** 487 L2 sites have no human observation. They are kept and flagged with `in_universe`.
8. **Restricted data:** row-level PSP-derived tables live only in `cache/` (CC BY-NC-SA). `results/` holds aggregates only.
9. **UniProt REST "current" snapshot** (2026_03) is pinned by SHA-256 in `external/uniprot/MANIFEST.json`.
