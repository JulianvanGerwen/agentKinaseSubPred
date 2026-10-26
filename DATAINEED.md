# DATAINEED — datasets/materials requested from the user

Single place for all data requests. Each item: what, why, what we tried, current workaround, priority. Status: OPEN / RESOLVED / WITHDRAWN. The project never waits on these; workarounds are in use.

Last updated: 2026-10-26 15:00

## OPEN

### 2. decryptM (Zecha et al. 2023 Science, doi 10.1126/science.ade3925) supplementary drug/target table — priority MEDIUM
- **Need:** designated targets of each decryptM drug and experiment metadata to interpret the 790 decryptM contrasts in the compendium.
- **Tried:** science.org → 403; decryptM 2.0 Zenodo record is restricted.
- **Workaround:** infer drug identities from compendium condition names + kinobeads maps.

### 3. SIGNOR PhosphoSIGNOR / full human export — priority MEDIUM
- **Need:** curated KS relations with SIGNOR's "direct" flag (for additional labels and the independent external test set).
- **Tried:** signor.uniroma2.it unreachable from this machine (connection fails / proxy 503) across multiple attempts; no mirror found.
- **Workaround:** OmniPath's SIGNOR-derived rows (lack "direct" flag).
- **Ideal file:** SIGNOR "download all human data" TSV and/or the PTM (phosphorylation) export.

### 4. Bayer … Kuster, bioRxiv 2025 — kinase-inhibitor potency-coherence phosphoproteomics (17M dose–response curves, 133 inhibitors) — priority MEDIUM (nice to have)
- **Need:** site × inhibitor EC50 data with known affinities; claimed thousands of new KS pairs for 96 kinases (useful as features, not labels).
- **Tried:** not found in PRIDE; probably only in ProteomicsDB.
- **Workaround:** none (decryptM subset already in compendium).

### 5. Provenance of the provided phosphoproteomics compendium (data/01–03) — priority LOW (documentation)
- **Need:** citation and processing pipeline; in particular whether non-significant qPTM values were kept and how fold changes were normalised. Needed for DATA.md and to reason about study effects.
- **Workaround:** treated as user-provided; described from its metadata.

### 6. Human Cell Map (humancellmap.org) extra files — priority LOW
- **Need:** some files return HTTP 403 to scripted downloads (crapome matrix, supplementary prediction zips). Core preys/baits/SAINT files were obtained.

### 7. Hein et al. 2025 (Cell) final localisation table — priority LOW
- **Tried:** Cell supplement → 403. Workaround: Zenodo analysis repo (per-protein enrichment tables) obtained.

### 8. decryptM curve fits (PRIDE PXD037285, `*_Curves.zip`) — priority LOW-MEDIUM (feat-phospho, 2026-10-26)
- **Need:** per-site EC50 / curve class for the decryptM kinase-inhibitor dose series, to build a "potency coherence" F6b variant (|log EC50_site − log Kd(k, drug)|), which is closer to direct-substrate evidence than the per-dose log2FC we use now.
- **Tried:** not attempted (open CC0 data, ~3–5 GB; out of the time box for this round).
- **Workaround:** dose series used as per-dose log2FC with occupancy = dose/(dose+Kd) in the ridge (F6b).

### 9. CPTAC TMT plex / batch and ischemia time per tumour sample — priority LOW (feat-phospho, 2026-10-26)
- **Need:** confounders for the CPTAC residualisation (DESIGN KB-9). The LinkedOmics `*_meta.txt` and PDC `Clinical_meta_data_v1` have only ESTIMATE stromal/immune scores and WES purity (used).
- **Tried:** both metadata files inspected; no plex/batch/ischemia columns.
- **Workaround:** residualise on host protein + stromal + immune + purity only. Ideal file: per-cohort sample → TMT plex table (PDC study biospecimen manifests) and ischemia time.

## Licence decisions (for your information; no action unless you disagree)
- Phospho.ELM requires an academic click-through licence: **not accepted, skipped** (adds almost nothing beyond PSP).
- PSP, Kinase Library (Cantley) data and iPTMnet are CC BY-NC-SA / non-commercial: used for research; outputs derived from them inherit non-commercial terms.
- Phosformer code/weights are CC BY-NC-ND: we reimplement the architecture rather than modify their code; released weights used inference-only as a baseline.

## RESOLVED
- Klaeger et al. 2017 Kinobeads Tables S1–S3 — provided by user in data/Klaeger/ (2026-10-26 14:26). Thank you.
- PhosphoSitePlus `Phosphorylation_site_dataset` and `Regulatory_sites` (2025-05-19) — provided by user in data/ (2026-10-26 13:35).
- CPTAC phosphoproteomics — obtained publicly (LinkedOmics S3 v1.2 + PDC pan-cancer), no login needed.
