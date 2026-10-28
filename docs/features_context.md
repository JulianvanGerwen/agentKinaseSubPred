# Context features: PPI (F3), localisation (F4), co-expression (F5)

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Owner: feat-context agent. Date: 2026-10-26. DESIGN v1 §7 F3/F4/F5.
Code: `src/ksp/features/{ppi_idmap,ppi,ppi_label,loc,coexpr}.py`; build scripts `scripts/features/context_{ppi,ppi_label,loc,coexpr,sanity}.py`.
Matrices are written only through `ksp.store.save_matrix` (store names `context/<module>/<name>`, files in `cache/features/context/<module>/`, each with a `.json` holding the index hashes). Each module also has an aggregate `descriptor.json` (copied to `results/features/context/<module>_descriptor.json`).
Shapes: kinase × protein = 503 × 18,702 (cache/index order); site vectors = 350,736 (sites.parquet order).

Rebuild everything (about 2 min, ≤10 threads):
```
export OMP_NUM_THREADS=10 OPENBLAS_NUM_THREADS=10 MKL_NUM_THREADS=10
.venv/bin/python -I scripts/features/context_ppi.py      # F3: about 40 s (dense 16k x 16k RWR inverse)
.venv/bin/python -I scripts/features/context_loc.py      # F4 + site topology: about 10 s
.venv/bin/python -I scripts/features/context_coexpr.py   # F5: about 10 s
.venv/bin/python -I scripts/features/context_ppi_label.py  # F3-label tests, timing, cross-fitted sanity
.venv/bin/python -I scripts/features/context_sanity.py   # enrichment / AUROC sanity, all modules
```

## 1. ID mapping (`ppi_idmap.py`)
Everything is mapped to the canonical Swiss-Prot accession (UniProt 2026_03). The mapping tables are HGNC complete set, the UniProt reviewed table with Entrez/RefSeq/topology (`external/uniprot_topology/`, downloaded for this module and recorded in DATA.md), and STRING v12 aliases. The STRING aliases are used only for ENSP→UniProt in COMPARTMENTS; no STRING interaction data is used.
- UniProt IDs have the isoform suffix stripped and must be in the reviewed set.
- Gene-level IDs (Entrez, ENSG, symbol) map to every reviewed accession of the gene.
- Symbol lookup order: HGNC approved → UniProt primary → HGNC previous (unique only) → HGNC alias (unique only).
- Every source logs its mapping rate to `results/features/context/id_mapping_rates.json`. A rate below 0.80 raises an error; nothing is dropped silently.
- Observed rates: BioPlex 293T 99.0 % / HCT116 99.5 % of rows; HuRI 98.8 %; OpenCell 99.8 %; Buljan 99.8 %; HCM 98.8 % of BFDR ≤ 0.01 rows (baits 191/192); HPA IF 98.8 % of genes; OpenCell localisation 100 %; HCM NMF 99.6 %; Hein 98.4 % of classified groups; COMPARTMENTS experiments 99.95 % / knowledge 99.1 %; HPA RNA 97.1 %; DepMap 99.99 %.
- Kinases are taken by accession from the index. No kinase-name resolution happens here; screen nodes are mapped like any other protein.

## 2. F3 PPI from systematic screens only (`ppi.py`, `context_ppi.py`)
| source | file | edges (mapped, undirected) | nodes | kinases (of 503) |
|---|---|---|---|---|
| bioplex_293t | BioPlex 3.0 293T network | 117,240 | 13,845 | 428 |
| bioplex_hct116 | BioPlex 3.0 HCT116 network | 70,777 | 10,072 | 325 |
| huri | HuRI.tsv (**not** Lit-BM or HI-union) | 51,286 | 8,143 | 182 |
| opencell | opencell-protein-interactions.csv | 29,231 | 5,217 | 207 |
| buljan | Buljan 2020 Table S2 (mmc3), the final high-confidence interactome | 7,254 | 2,376 | 351 |
| hcm | Human Cell Map SAINT bait–prey pairs with BFDR ≤ 0.01 | 34,251 | 4,439 | 164 |
| **combined** | union, unweighted | 279,677 | 16,178 | **478 (95 %)** |

The combined graph covers 15,490 of the 18,702 index proteins. It has 5 components; the giant component has 16,168 nodes.

Excluded: STRING (every channel), BioGRID, IntAct, HIPPIE, OmniPath, SIGNOR, BioPlex `interactionList_v4a`, and HuRI Lit-BM / HI-union.

Caveat: HCM is proximity labelling (BioID with compartment-marker baits), not physical binding. It is kept as its own flag so it can be ablated.

Features (kinase × protein):
- `ppi_edge_<src>` (int8): direct edge in that source.
- `ppi_edge_any` (int8): direct edge in any source.
- `ppi_n_studies` (0–5): number of studies with an edge; BioPlex counts once.
- `ppi_rwr`: random walk with restart, restart probability 0.5, computed exactly as R = r(I − (1−r)W)⁻¹ with W = AD⁻¹ (column-stochastic). Each column sums to 1.
- `ppi_rwr_dc`: degree-corrected version, log2((R+ε)/(bg+ε)), where bg[q] = mean RWR of q over all seeds (PageRank with uniform teleport) and ε = 0.01/N.
- `ppi_rwr_dc_pct`: percentile of `ppi_rwr_dc` within the kinase's row.
- `ppi_dist`: shortest path in hops, capped at 5. Lower means closer (descriptor `sign: -1`).

Missing values:
- RWR and distance features are NaN when the kinase or the protein is not in the graph.
- Edge flags are 0 when there is no edge. Use `ppi_kinase_in_source` (K × 6) and `ppi_protein_in_source` (P × 6) to tell "not screened" apart from "screened, no edge".
- `rwr_full.npy` (N × N float32, 1 GB), `rwr_background.npy` and `rwr_nodes.txt` are kept for F3-label.

## 3. F3-label: cross-fitted substrate neighbourhood (`ppi_label.py`)
`substrate_neighbourhood(train_pairs, query_pairs=None, exclude_same_cluster=cluster_map)` is a pure function of `train_pairs`; the graph itself is label-free.
- For each (k, q) it computes over k's training substrate proteins s: `nb_max`, `nb_mean` and `nb_top3` of `DC[q,s]` (the degree-corrected RWR from s), plus `nb_n_adj` (direct neighbours), `nb_frac_adj` and `nb_n_used`.
- Excluded substrates: s = q, s in the same cluster as q, and s = k (the autophosphorylation partner).
- Values are NaN if q is not in the graph or no training substrate is left.
- With `query_pairs=None` it returns dense K × P arrays. All 503 kinases × 18,702 proteins take **6–8 s** after a 1–2 s graph load.
- `PPINeighbourhood` (registered as `f3_ppi_neighbourhood`) implements the `ksp.crossfit` fit/transform contract. Its columns are prefixed `f3nb_`, and it uses `ctx.clusters` as a second guard.
- Unit tests in `context_ppi_label.py` cover self, same-cluster and autophos exclusion, the value of max/mean, and agreement between dense and pair modes. They pass.

## 4. F4 localisation (`loc.py`, `context_loc.py`)
Coarse vocabulary (15 terms): nucleoplasm, nucleolus, chromatin, nuclear_envelope, cytosol, plasma_membrane, cell_junction, mitochondrion, er, golgi, endolysosome_vesicle, peroxisome, centrosome_mt, actin_if, extracellular. The full mapping tables (`HPA_MAP`, `OPENCELL_MAP`, `HCM_MAP`, `HEIN_MAP`, `GO_ANCHORS`) in `src/ksp/features/loc.py` are the record. Summary:
- **HPA IF v25**: weight = main 1.0 / additional 0.5, multiplied by the per-location reliability (Enhanced 1, Supported 0.8, Approved 0.5, Uncertain 0.2).
  - Nuclear bodies and speckles → nucleoplasm.
  - Mitotic chromosome and kinetochore → chromatin.
  - Cytoplasmic bodies, aggresome and rods & rings → cytosol.
  - Vesicles, endosomes, lysosomes and lipid droplets → endolysosome_vesicle.
  - Centrosome, centriolar satellite, microtubules, spindle, cytokinetic bridge, midbody and cilia → centrosome_mt.
  - Actin and intermediate filaments → actin_if.
  - "Predicted to be secreted" → extracellular (0.5).
  - Sperm-specific compartments are dropped. Any location without a mapping raises an error.
- **OpenCell**: grades 3/2/1 → 1/0.6/0.3. `cytoskeleton` → actin_if 0.5 + centrosome_mt 0.5.
- **Human Cell Map**: the 20-dim NMF rank profile (BioID2021_NMF.xlsx, a copy of Go 2021 found in the Hein repo) is projected through rank → MMF-localisation name (`preys-latest.txt`). Mixed names are split 0.5/0.5.
- **Hein 2025**: one-hot of the consensus cluster annotation; `unclassified` counts as missing. 14-3-3 scaffold, translation, proteasome, stress granule and P-body → cytosol; ERGIC → er/golgi 0.5 each.
- **COMPARTMENTS** (score/5, from the ancestor-propagated full files): anchor GO terms per category.
  - Extracellular uses GO:0005615 (extracellular space), not "extracellular region" and not "vesicle", which would pull in the very common exosome (HDA) annotations.
  - Plain "nucleus" (GO:0005634) counts as nucleoplasm only when no nuclear sub-compartment scores as high.
  - The experiments channel has HPA as its only source, so it is redundant with `hpa`. It is kept as `loc_cos_comp_exp` but left out of the combined vector.
  - The knowledge channel is ablation-only (`ablation_only: true` in its metadata).

Features (kinase × protein, float16, NaN when either side is missing):
- `loc_cos_{hpa,opencell,hcm,hein,comp_exp,comp_know}`: cosine similarity per source.
- `loc_cos_combined`: cosine of the mean unit vector over HPA, OpenCell, HCM and Hein.
- `loc_cos_combined_know`: ablation version that adds the knowledge channel.
- `loc_n_sources` (int8): missingness summary, 0–4 primary sources with both sides annotated. Per-source missingness is the NaN pattern; `vec_<src>_{kinase,protein}` hold the raw vectors.
- `loc_hcm_nmf_cos`: cosine in the native NMF space.
- `loc_hein_profile_corr`: Pearson correlation of the z-scored enrichment profiles over 81 organelle-IP samples.
- `loc_hpa_min_reliability`.

Coverage (fraction of kinases / index proteins with a vector):

| source | kinases | proteins |
|---|---|---|
| HPA | 0.85 | 0.69 |
| OpenCell | 0.24 | 0.07 |
| HCM | 0.33 | 0.23 |
| Hein | 0.41 | 0.37 (profile: 0.52 / 0.43) |
| COMPARTMENTS experiments | 0.53 | 0.42 |
| COMPARTMENTS knowledge | 0.93 | 0.85 |

**Site topology** (`site_topo_code`, int8 per site, codes in the json):
- Sources: UniProt TOPO_DOM / TRANSMEM / INTRAMEM / SIGNAL / TRANSIT.
- Rule: a site in a non-TM protein that lies after the signal peptide counts as secreted/lumenal (code 5).
- `site_lumenal` (float16): 1 = extracellular/lumenal (codes 2, 5); 0 = cytoplasmic/nuclear/mitochondrial/peroxisomal side; NaN = TM segment, signal peptide, or a TM protein with no TOPO_DOM covering the site.
- Overall, 7.4 % of U_obs sites are flagged lumenal, against 0.7 % of DEV L2 KS sites. **FAM20C: 87.5 % lumenal (7/8 DEV sites)**.

## 5. F5 co-expression / co-dependency (`coexpr.py`, `context_coexpr.py`)
- `coexpr_hpa_tissue_spearman`: 51 HPA consensus tissues, log2(nTPM+1), ranked.
- `coexpr_depmap_pearson`: DepMap 26Q1 Chronos, 1,207 lines, NaN values mean-imputed.
- `coexpr_depmap_pearson_pc5`: the same after removing the top 5 PCs.
- Coverage: kinases 1.00 / 0.996; proteins 0.96 / 0.94.
- CPTAC protein co-abundance is not built.

## 6. Sanity: enrichment of KS pairs vs random pairs (not an official metric)
`results/features/context/sanity_all.tsv`, `sanity_ppi_label.tsv`.
- Positives: `ksp.labels.load_dev('L2')` minus autophosphorylation, protein level. That gives 6,436 pairs, 355 kinases and 2,322 proteins (DEV only; TEST is sealed).
- RAND background: same kinase × random DEV protein that is not a known pair in `load_dev('all')`.
- SHUF background: same kinase × a substrate of another kinase, which controls for study bias.
- Seeds 0–4. AUROC "all" imputes NaN as lowest; "obs" uses only pairs that are not missing.

| feature | KS pairs with value > 0 | RAND | SHUF | KS/SHUF | AUROC vs RAND (all/obs) | AUROC vs SHUF (all/obs) |
|---|---|---|---|---|---|---|
| edge BioPlex 293T | 1.6 % | 0.1 % | 0.2 % | 7.1 | 0.507 | 0.507 |
| edge BioPlex HCT116 | 0.9 % | 0.1 % | 0.1 % | 14 | 0.504 | 0.504 |
| edge HuRI | 0.3 % | 0.0 % | 0.1 % | 5.5 | 0.502 | 0.501 |
| edge OpenCell | 1.1 % | 0.1 % | 0.2 % | 5.6 | 0.505 | 0.504 |
| edge Buljan | 1.5 % | 0.1 % | 0.3 % | 5.3 | 0.507 | 0.506 |
| edge HCM | 0.4 % | 0.0 % | 0.1 % | 3.2 | 0.502 | 0.501 |
| edge any screen | 3.6 % | 0.3 % | 0.8 % | 4.5 | 0.516 | 0.514 |
| ppi_rwr_dc | – | – | – | – | 0.668 / 0.639 | 0.544 / 0.547 |
| ppi_dist (−) | – | – | – | – | 0.638 / 0.600 | 0.537 / 0.539 |
| F3-label nb_mean (cross-fitted, official CV folds + substrate clusters) | – | – | – | – | 0.681 / 0.662 | 0.547 / 0.556 |
| F3-label nb_max (cross-fitted) | – | – | – | – | 0.656 / 0.632 | 0.552 / 0.562 |
| F3-label nb_max (IN-SAMPLE, for contrast) | – | – | – | – | 0.683 | 0.583 |
| loc_cos_hpa | – | – | – | – | 0.646 / 0.592 | 0.556 / 0.563 |
| loc_cos_combined | – | – | – | – | 0.681 / 0.638 | 0.579 / 0.590 |
| loc_cos_comp_know (ablation) | – | – | – | – | **0.741** / 0.707 | 0.618 / 0.620 |
| loc_hein_profile_corr | – | – | – | – | 0.642 / 0.576 | 0.529 / 0.544 |
| coexpr HPA tissue | – | – | – | – | 0.596 / 0.581 | 0.560 / 0.560 |
| coexpr DepMap (pc5) | – | – | – | – | 0.545 / 0.518 | 0.515 / 0.514 |

Interpretation:
- Direct-edge enrichment vs SHUF is 3–14× (any screen: 4.5×). This matches the leakage audit's ≈ 5–8× for the screens (`results/leakage/`). Edges cover only 3.6 % of KS pairs, so they are high precision and low coverage.
- Much of the AUROC vs RAND comes from coverage, i.e. study bias. Well-studied substrates are in more screens and more localisation resources, and `loc_n_sources` alone gives 0.656. Against SHUF, which controls for study bias, every feature sits at 0.50–0.62.
- **No feature exceeds 0.8.** The strongest is the COMPARTMENTS knowledge channel at 0.74 vs RAND. That is expected: it imports GO/UniProt curation, often from the same papers as the KS labels. It stays ablation-only.
- The in-sample F3-label variant scores higher than the cross-fitted one (SHUF 0.583 vs 0.552). The gap is the leak that cross-fitting removes.
- The real substrate clusters change F3-label very little compared with self-only exclusion.
