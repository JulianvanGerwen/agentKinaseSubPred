> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

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

## Tyrosine-kinase features F9 (feat-tyr; details in docs/features_tyr.md)
| feature (store `tyr/...`) | source(s) | label-derived? | leakage risk | coverage |
|---|---|---|---|---|
| `sh2_dock_max`, `sh2_dock_local`, `sh2_dock_near`, `sh2_dock_delta`, `sh2_dock_spec` | Yaron-Barir 2024 Nature Suppl. Table 5 (MOESM7) `sh2_norm_scaled_matrices`: 76 SH2 domains x 6 positions (-2,-1,+1..+4) x 20 aa; scored on the **other** U_ms pY sites of the substrate protein, query site always excluded | no | low. Uses the U_ms observation status of *other* sites only, the same rule as the accepted F1 `primed_soft` variant (DESIGN §7 F1); the query site's own status is never read (A6/DL2-7). No KS label enters, so no crossfit family. **But** `sh2_dock_max/local/near` correlate with pY-cluster density, a site-popularity proxy: `sh2_dock_spec` (across-domain z-score) is the popularity-free version and must be the one reported if only one is kept | 24 of the 28 SH2-bearing kinases x 45,895 Y sites (998,496 non-NaN cells). JAK1/2/3 and TYK2 have no matrix. Range: log2 units, roughly [-6, +6]; `_spec` ~ N(0,1) |
| `rtk_pm` | HPA IF v25 + Human Cell Map protein vectors (existing F4 stores), plasma_membrane + cell_junction components; exposed only on transmembrane-receptor kinase rows | no | very low (same sources as F4, already audited) | 54 RTK rows x Y sites |
| `itam_gated`, `site_itam` (`itam_half`, `itam_full`, `itim`, `itam_order`) | UniProt 2026_03 canonical sequences; Yxx[LIV], tandem Yxx[LIV] at 6-12 aa, [ILVS]xYxx[LV] | no | none (sequence only) | all 45,895 Y sites; `itam_gated` on 26 non-receptor TK rows |
| `site_sh2` (`n_other_py_protein`, `n_other_py_50aa`, `adaptor_sh2_max`, `adaptor_sh2_local`) | U_ms site index; the 51 non-kinase SH2 matrices of MOESM7 | no | **moderate: site-popularity proxy.** B0-pop alone reaches 0.583 E-hard AUROC on Y kinases. Kept in a separate column block; **ablation only**, and they add nothing at n_pos>=10 (docs/features_tyr.md §3.2) | all Y sites |
| `kinase_desc` (`n_sh2`, `n_sh3`, `n_ptb`, `is_tm`, `is_receptor`, `has_sh2_matrix`) | UniProt 2026_03 `Domain [FT]` / `Keywords` | no | low; curated annotation but carries no KS information | 503 kinases. `n_ptb` is identically 0 (no human protein kinase has a UniProt PID/PTB domain) and should be dropped |

## Similarity-to-training-substrates features (feat-knn; details in docs/features_knn.md)
| feature (crossfit family `knn_sim`, `ksp.features.knn.KnnSim`) | source(s) | label-derived? | leakage risk | coverage |
|---|---|---|---|---|
| `knn_{esmc,esmcat,blosum}_{max,top3,kde,lr}`, `knn_log_ntrain` (default spaces; `esmw` optional) | site representations: ESM-2 650M per-residue embeddings in full-protein context (R3 window cache `cache/r3/win_emb.npy`) → centre residue (`esmc`), ±7 concatenated per-position PCA-128 (`esmcat`), ±7 mean-pooled (`esmw`); BLOSUM62 15-mer centre-weighted (`blosum`). Plus training KS pairs | **yes**: must be cross-fitted (ksp.crossfit). The query's own protein and substrate cluster are also excluded inside `transform` | **very high if not cross-fitted**: on fold-0 training rows a fit on all labels without exclusion gives macro AUROC 0.93–1.00, against 0.59–0.76 under the driver. With in-transform exclusion alone (no driver) it gives 0.58–0.76, close to the driver | U_ms sites (347,354) × kinases with ≥1 training site of the query's residue type; NaN otherwise. Range: cosine [−1, 1]; blosum mean weighted BLOSUM62 per position ≈ [−2, 8]; kde/lr in log units |
| `knn_{sp}_kde_bg` | the same representations; 2,000 random U_ms sites per residue type (seed 0) | no (kinase-free; reported only where the kinase has training sites) | low. It is a "typical-site" score, so it can act as a site-popularity proxy (fold-0 AUROC 0.50–0.54; 0.51–0.61 for Y kinases) | as above |
