# F0 site-intrinsic features: AlphaFold structure and ESM-2 embeddings

Owner: structure/ESM agent, 2026-10-26. DESIGN §7 F0 (site) and F7 (kinase descriptors). No feature here is derived from labels.

## Scope
Every S/T/Y residue of every Swiss-Prot canonical human protein in UniProt 2026_03 (`ksp.uniprot`; 20,416 proteins, 11.41M residues, **1,870,491 S/T/Y**). Site id is `f"{acc}_{res}{pos}"` (1-based, canonical). The work does not depend on the site universe. The aligned exports cover the canonical index `cache/index/sites.parquet` (350,736 U_obs sites).

## Outputs
| path | content |
|---|---|
| `cache/structure/site_structure.parquet` | 1,870,491 rows: site_id, acc, pos, res plus the structure features below, with `af_status` (ok / not_in_afdb / seq_mismatch) |
| `cache/structure/residue_structure.parquet` | the same features for all 11.35M modelled residues (for window features) |
| `cache/structure/site_structure.meta.json` | counts, PAE calibration, mismatch and missing accession lists |
| `cache/esm2_650m/sty_<w>_<k>.npy` + `.site_ids.npy` | float16 (n, 1280) last-layer embeddings of S/T/Y residues. Each shard has an explicit site_id array, and `sty_index.parquet` maps site_id to (file, row). Loader: `ksp.features.esm_embed.load_sty` |
| `cache/esm2_650m/protein_mean.npy` (+ `_index.parquet`, `.accs.npy`) | 20,416 × 1280 float32 mean over all residues |
| `cache/esm2_650m/kinases.parquet` | one row per (kinase, kinase domain) from `results/labels/kinase_list.tsv` (562 accessions, 577 rows; `include` and `is_primary` flags). Columns: `emb_full`, `emb_domain_ctx` (domain mean taken from the full-protein run) and `emb_domain_iso` (the domain sequence embedded alone) |
| `cache/esm2_3b/...` | the same set for ESM-2 3B (D = 2560) |
| Feature store (`ksp.store`, ID-checked) | `structure/site_structure` (site × 18); `esm2_{650m,3b}/site_pca64` (site × 64, PCA fit on U_ms sites only, with `.components.npy` / `.mean.npy`); `esm2_{650m,3b}/protein_mean` (protein axis); `esm2_{650m,3b}/kinase_domain` (primary domain; 5/503 kinases have no domain, so NaN) and `kinase_full` (kinase axis) |
| `results/features/*.tsv/json` | validation |

Scripts:
- `src/ksp/features/structure.py`
- `src/ksp/features/esm_embed.py` (`run` / `merge`)
- `src/ksp/features/esm_kinases.py`
- `scripts/download_afdb_pae.py`
- `scripts/export_structure_esm_features.py`
- `scripts/validate_structure_esm_features.py`
- tests: `tests/test_features_structure_esm.py`

## Structure (AlphaFold DB v6)
**Input.** AFDB v6 human proteome tar (based on UniProt 2025_03; 23,586 models; SHA-256 in DATA.md). The tar is streamed and never extracted. mmCIF is parsed by our own minimal loop parser, and no downloaded code is run.

**Sequence check.** Every residue's `pdbx_sifts_xref_db_num/res` is compared with the 2026_03 canonical sequence:
- 20,249 proteins are ok;
- **41 have a sequence mismatch** (the UniProt sequence changed after 2025_03), so they are NA;
- **126 are not in AFDB**, so they are NA.

S/T/Y coverage is 1,860,126 / 1,870,491 = 99.45%.

**Proteins > 2700 aa.** 214 of 216 are assembled from AFDB fragments (1400 aa, step 200; UniProt numbering comes from SIFTS). Each residue takes its features from the fragment where it is most central. Contacts across fragments are not seen, so exposure is less reliable there (`af_fragmented` flag).

**Features** (a re-implementation of structuremap 0.0.10 / Bludau et al. 2022; the loop code is checked against a literal re-statement in the tests):
- `plddt`, `plddt_s5`, `plddt_s10`. The smoothed versions are window means over i±5 and i±10, truncated at the ends.
- `pse_{12,24}_{70,180}`: half-sphere exposure (θ = 70°) and full-sphere exposure (θ = 180°) at 12 Å and 24 Å. A residue counts as a neighbour when its Cα is within r Å and inside the angle θ of Cα→Cβ; glycine uses structuremap's pseudo-Cβ.
- **Caveat: no PAE.** These counts are computed without the PAE adjustment (pPSE with PAE = 0). AFDB serves PAE only per protein over HTTP and returned HTTP 429 even at 2–8 connections, so a full download (about 24 GB) was abandoned.
  - PAE was fetched for a random calibration set of 399/400 proteins (seed 0). True pPSE is stored for that set (`ppse_*`, `idr_sm`, `high_acc_5`; NA elsewhere).
  - Spearman correlation between PAE-free and PAE pPSE is 0.88–0.90.
- `idr` = `pse_24_180_s10` ≤ 71.25. The threshold was calibrated on the PAE set to reproduce structuremap's `ppse_24_180_s10` ≤ 34.27: accuracy 0.924, κ = 0.84. 38.5% of proteome residues are flagged as IDR.
- `high_acc` = `pse_12_70` ≤ 11, calibrated to reproduce `ppse_12_70` ≤ 5: accuracy 0.90, κ = 0.69.
- `disorder_plddt` = `plddt_s10` < 50.
- `sasa` and `rsa`:
  - SASA is our own vectorised Shrake–Rupley: 100 points, 1.4 Å probe, Bondi radii, k = 6 nearest-centre occlusion. Against Biopython ShrakeRupley it gives r = 0.997, with a mean absolute difference of about 2.7 Å².
  - RSA is SASA divided by the Tien 2013 maximum. It is computed on the isolated chain.
- `ss8` and `ss3`: DSSP assignments **as shipped in the AFDB mmCIF `_struct_conf`** (DSSP is not run locally). ss8 classes: H/G/I/P(PPII)/E/T/S/C. Proteome ss3 fractions: C 0.53, H 0.34, E 0.13.
- Skipped: distance to Pfam / kinase domain.

**Runtime.** 15.7 min on 15 CPU processes with single-threaded BLAS (14.6 min of that is parsing and geometry).

## ESM-2
**Models.** HF `facebook/esm2_t33_650M_UR50D` and `esm2_t36_3B_UR50D` (MIT; revisions in DATA.md). Settings:
- `EsmModel.last_hidden_state`, which is after the final layer norm and equals fair-esm `repr_layers=[last]`;
- BOS/EOS removed;
- fp16 with SDPA attention;
- batches of up to 16k tokens (8k for 3B).

**Long proteins.** 2,303 proteins are longer than 1022 aa. They use 1022-aa windows with stride 511, and the last window is right-aligned. Overlapping windows are averaged with equal weight.

**Precision check.** fp16 batched vs fp32 single-sequence: per-residue cosine ≥ 0.99999 for both models.

**Runtimes** (RTX 4090, GPUs 2+3, 2 workers):

| model | wall time | GPU-h | throughput | peak memory |
|---|---|---|---|---|
| 650M | 2.8 min | **0.093** | 34k residues/s per GPU (~61 proteins/s) | 2.1 GB |
| 3B | 7.9 min | **0.26** | 12.1k residues/s per GPU | 6.1 GB |

- The 200-protein benchmark measured 11.8k residues/s for 3B vs 31k for 650M, which projected 0.27 GPU-h for 3B. That is well under the 2 GPU-h limit, so the full 3B run was done.
- Kinase embeddings took 20 s (650M) and 55 s (3B).

**Disk.** 650M takes 4.8 GB and 3B takes 9.4 GB.

## Validation (`results/features/structure_esm_validation.json`; non-official sanity checks)
**Known disordered proteins vs folded controls** (`idp_panel_structure.tsv`):

| metric | 10 IDPs | folded controls |
|---|---|---|
| mean pLDDT | 61.7 | 90.6 |
| IDR fraction | 1.00 | 0.17 |
| RSA | 0.68 | 0.33 |

- `disorder_plddt` catches only 18% of IDP residues. AF often gives disordered proteins mid-range pLDDT with spurious helices (SNCA has pLDDT 75 and is 64% "helix"). The exposure-based `idr` is the better disorder flag.
- Expected failures:
  - INS (preproinsulin; the monomer prediction is low confidence) is flagged as IDR;
  - CALM1 (an elongated dumbbell) is 46% IDR.

**PSP Phosphorylation_site_dataset** (229,868 human canonical sites matched by exact 15-mer at the stated position) vs the other S/T/Y of the same proteins:
- IDR fraction is 0.59 vs 0.39 (OR 2.28; S 2.53, T 2.05, Y 1.86; low-throughput sites 2.56). This matches the enrichment of phosphosites in disordered regions reported by Bludau 2022.
- AUCs: −pLDDT 0.60, RSA 0.61, −`pse_24_180_s10` 0.63.
- Within ordered regions, phosphosites are more exposed (RSA 0.31 vs 0.27).

**ESM-2 650M linear probe** (protein-grouped 5-fold CV, 120k sites; sanity only):
- phosphosite vs not: AUC 0.76, against 0.63 for the structure features alone;
- the structural IDR flag predicted from ESM: AUC 0.97.

**Consistency.** Kinase full means agree with the proteome means (cos ≥ 0.999999). The stored ↔ store-aligned values were spot-checked.

The table broken down by residue × class (non-phospho / phospho / MS≥5 / PSP KS-annotated) is in `phospho_vs_nonphospho_structure.tsv`. The KS-annotated class uses all PSP labels and is descriptive only.

## Caveats
- **Version skew.** The AFDB v6 sequences are from 2025_03. 41 proteins with changed sequences are NA rather than remapped.
- **Isolated chains.** Structures are monomer predictions, so exposure ignores complexes, membranes and partners.
- **PAE-free exposure.** pPSE uses no PAE (see above). The thresholds were calibrated on 396 proteins.
- **Fragments.** Exposure for fragmented proteins (> 2700 aa) misses contacts across fragments.
- **ESM windows.** Embeddings in the window-overlap zones of proteins > 1022 aa are averages of embeddings with different context.
- **Protein means are unsplit.** They are computed on the whole protein and are not split-aware. They contain no label information, but homologous proteins have similar embeddings, so splits must stay homology-grouped as planned.
- **Kinase domains.** Domains come from `kinase_list.tsv` (`kinase_domains`). Some non-primary entries are not catalytic (e.g. CASK guanylate kinase-like), so use `is_primary`.
- **PCA-64** keeps 58% of the variance for 650M and 45% for 3B.
