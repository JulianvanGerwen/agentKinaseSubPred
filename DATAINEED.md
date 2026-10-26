# DATAINEED — datasets/materials requested from the user

Single place for all data requests. Each item: what, why, what we tried, current workaround, priority. Status: OPEN / RESOLVED / WITHDRAWN. The project never waits on these; workarounds are in use.

Last updated: 2026-10-26 14:35

## OPEN

### 1. Klaeger et al. 2017 (Science, doi 10.1126/science.aan4368) Kinobeads supplementary tables — priority HIGH
- **Need:** inhibitor → kinase apparent Kd (Table S1/S2 "Kd_app", target lists for 243 drugs) to map inhibitor-perturbation conditions (decryptM, other drug studies in the compendium) to the kinases they inhibit (feature F6b).
- **Tried:** science.org supplementary → HTTP 403; PMC supplements → reCAPTCHA page; ProteomicsDB API (assay id not found within time box).
- **Workaround:** Reinecke 2024 (Nat Chem Biol) tables + CurveCurator reprocessed kinobeads (Zenodo 8399823) + Davis 2011 Kd + provided Hijazi matrix.
- **Ideal file:** the Science supplementary Excel files (aan4368_Tables_S1–S?.xlsx), dropped into `data/`.

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

## Licence decisions (for your information; no action unless you disagree)
- Phospho.ELM requires an academic click-through licence: **not accepted, skipped** (adds almost nothing beyond PSP).
- PSP, Kinase Library (Cantley) data and iPTMnet are CC BY-NC-SA / non-commercial: used for research; outputs derived from them inherit non-commercial terms.
- Phosformer code/weights are CC BY-NC-ND: we reimplement the architecture rather than modify their code; released weights used inference-only as a baseline.

## RESOLVED
- PhosphoSitePlus `Phosphorylation_site_dataset` and `Regulatory_sites` (2025-05-19) — provided by user in data/ (2026-10-26 13:35).
- CPTAC phosphoproteomics — obtained publicly (LinkedOmics S3 v1.2 + PDC pan-cancer), no login needed.
