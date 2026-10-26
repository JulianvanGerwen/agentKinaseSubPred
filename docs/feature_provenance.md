# Feature provenance register

One row per feature family: source, label-derived?, leakage risk, coverage. Each agent appends its own section.

## Context features (feat-context; details in docs/features_context.md)
| feature (store `context/...`) | source(s) | label-derived? | leakage risk | coverage |
|---|---|---|---|---|
| `ppi/ppi_edge_{bioplex_293t,bioplex_hct116,huri,opencell,buljan,hcm}`, `ppi_edge_any`, `ppi_n_studies` | BioPlex 3.0 (Huttlin 2021), HuRI (Luck 2020; not Lit-BM), OpenCell (Cho 2022), Buljan 2020 Table S2, Human Cell Map (Go 2021; SAINT BFDR ≤ 0.01) | no | low: systematic screens only. Buljan score cut-offs were tuned on IntAct literature PPIs. HCM is proximity, not binding. Enrichment vs SHUF is 3–14×, as in the leakage audit | edges on 3.6 % of DEV L2 KS protein pairs; kinases in graph 478/503 |
| `ppi/ppi_rwr`, `ppi_rwr_dc`, `ppi_rwr_dc_pct`, `ppi_dist` | combined screens graph (16,178 nodes, 279,677 edges) | no | low; study bias (degree) is partly removed by `_dc` | 95 % kinases, 83 % index proteins |
| `f3_ppi_neighbourhood` (crossfit family; `ksp.features.ppi_label`) | screens graph + training labels | **yes**: must be cross-fitted (ksp.crossfit) | high if not cross-fitted. In-sample vs cross-fitted AUROC vs SHUF is 0.583 vs 0.552 | about 94 % of KS pairs |
| `loc/loc_cos_hpa`, `loc_hpa_min_reliability`, `vec_hpa_*` | HPA subcellular IF v25 | no | very low (antibody IF); coverage is study-biased | 85 % kinases / 69 % proteins |
| `loc/loc_cos_opencell` | OpenCell live-cell localisation | no | very low | 24 % / 7 % |
| `loc/loc_cos_hcm`, `loc_hcm_nmf_cos` | Human Cell Map NMF (Go 2021) | no | low | 33 % / 23 % |
| `loc/loc_cos_hein`, `loc_hein_profile_corr` | Hein 2025 organelle IP-MS (Zenodo 14000085 repo tables) | no | low | 41–52 % / 37–43 % |
| `loc/loc_cos_comp_exp` | COMPARTMENTS experiments channel (HPA-derived only) | no | very low; redundant with HPA, not in the combined vector | 53 % / 42 % |
| `loc/loc_cos_comp_know`, `loc_cos_combined_know` | COMPARTMENTS knowledge channel (GO/UniProt curation) | no | **moderate**: curated, partly from the same papers as the KS labels. **Ablation only.** Highest single AUROC (0.74 vs RAND) | 93 % / 85 % |
| `loc/loc_cos_combined`, `loc_n_sources` | HPA + OpenCell + HCM + Hein | no | low | 92 % of KS pairs |
| `loc/site_topo_code`, `loc/site_lumenal` (site vectors) | UniProt 2026_03 topology (TOPO_DOM/TRANSMEM/SIGNAL/TRANSIT) | no | low to moderate: UniProt topology is curated, but it is not KS information | all 350,736 sites; 7.4 % lumenal; FAM20C KS sites 87.5 % lumenal |
| `coexpr/coexpr_hpa_tissue_spearman` | HPA RNA tissue consensus v25 | no | none | 100 % / 96 % |
| `coexpr/coexpr_depmap_pearson`, `_pc5` | DepMap 26Q1 Chronos (figshare) | no | none | 99.6 % / 94 % |

## Phosphoproteomics features F6 (feat-phospho; details in docs/features_phospho.md)
| feature (store `phospho/...`) | source(s) | label-derived? | leakage risk | coverage (non-NaN) |
|---|---|---|---|---|
| `f6b_ridge_beta`, `f6b_ridge_z` | compendium log2FC (692 inhibitor conditions ≤ 2 h: Hijazi 2020 + decryptM + 40 other studies) × inhibition matrix (Hijazi kinhib; Klaeger 2017 / Reinecke 2024 / CurveCurator / Davis 2011 Kd; curated aliases) | no | low–moderate: PSP in-vivo labels may come from the same inhibitor papers (DESIGN §13); study-removal variants `_clean_bench`, `_f6_clean`; ascertainment: labelled sites 2.2× more often non-missing | 384 kinases × 75k sites (10.0 M pairs); range: β scaled log2FC per unit occupancy, z ~ N(0,1)-like |
| `f6b_contrast` | same | no | same | 378 kinases × 63k sites |
| `f6b_kinase_maxcos` (kinase vector) | inhibition matrix | no | none (kinase-level) | 384 kinases; [−1, 1] |
| `f6c_prot` | CPTAC BCM v1.2 tumours (10 cohorts), site residualised on host protein + ESTIMATE stromal/immune + WES purity; kinase protein residualised on confounders | no | low; tumour benchmark (benchmarKIN) overlap → excluded in clean variant | 427 kinases × 69k sites; mean Fisher z |
| `f6c_aloop` | same + A-loop sites (DFG…APE in UniProt kinase domain) | no | low; kinase's own sites set NaN | 99 kinases × 69k sites |
| `f6a_coreg` family (crossfit; `ksp.features.phospho.CoregF6a`): `f6a_{comp,cptac}_{top,med}` | compendium study blocks / CPTAC residualised cohorts + training KS pairs | **yes**: cross-fit via ksp.crossfit; same-protein, same-cluster and self substrates excluded inside the function as well | high if not cross-fitted; ascertainment as above | ~300 kinases × 60k sites per fold; Fisher z |
| `sanity_f6a_*` | temporary protein-grouped 5-fold cross-fit | yes | **sanity only, never a model input** | — |
