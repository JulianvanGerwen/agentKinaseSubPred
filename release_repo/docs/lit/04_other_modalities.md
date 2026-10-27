> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# 04 — Non-phosphoproteomics modalities for kinase–substrate (KS) prediction

Author: literature/data agent (modalities review). Date: 2026-10-26.
Scope: data that can tell us *which kinase* phosphorylates *which S/T/Y site*, other than phosphoproteomics: PPIs, localisation, co-expression/co-dependency, site-level sequence/structure, kinase-side descriptors, kinase activity context, and structure-based (AlphaFold-type) KS modelling. For each source: content, size, URL, version, license, identifier type, **leakage risk** relative to PSP-derived labels, expected value.

Conventions: "Leakage" = the feature source itself contains (or was built from) curated literature KS knowledge such as PhosphoSitePlus (PSP), SIGNOR, HPRD, Reactome phosphorylation reactions, BioGRID "Biochemical Activity", or text-mined "A phosphorylates B" statements. Leakage inflates test metrics because the test labels were curated from the same papers. "Study bias" = a weaker problem: well-studied proteins have more edges and more labels, so a feature correlates with label density without being causal.

Things verified in this session are marked **[verified]**; numbers quoted from memory or secondary sources are marked *(verify)*.

---

## 0. Key empirical finding: a PPI leakage audit (new, reproducible)

Script: `scripts/audit_ppi_ks_leakage.py` (run: `python -I scripts/audit_ppi_ks_leakage.py --seed S`, ~40 s). Results: `results/leakage/ppi_ks_leakage_seed{0..4}.tsv`. Inputs (already in `external/`, downloaded by another agent): STRING v12.0 human `protein.links.detailed` and `protein.physical.links.detailed`, BioGRID 5.0.262 human tab3, BioPlex 3.0 (293T, HCT116), HuRI, OpenCell interactions. Labels: PSP Kinase_Substrate_Dataset (2025-05-19), human kinase × human substrate, protein-level, autophosphorylation removed → **7,711 KS protein pairs, 406 kinases, 2,915 substrate proteins** (7,653 after requiring STRING mapping).

For each source we compute the fraction of KS pairs with evidence, vs. two backgrounds of equal size per kinase: RAND (kinase × random proteome protein) and SHUF (kinase × PSP substrate of *another* kinase; controls for substrate study bias/degree). Mean over seeds 0–4 (SD of ratios ≤ 10 % except the sparse screens):

| Feature (evidence > 0 unless stated) | frac KS | frac RAND | frac SHUF | KS/SHUF |
|---|---|---|---|---|
| STRING full: experiments channel | 0.557 | 0.093 | 0.196 | 2.8 |
| STRING full: database channel | 0.174 | 0.012 | 0.021 | 8.1 |
| STRING full: textmining channel | **0.744** | 0.151 | 0.321 | 2.3 |
| STRING full: coexpression channel | 0.576 | 0.117 | 0.255 | 2.3 |
| STRING full combined ≥ 400 | 0.443 | 0.033 | 0.082 | 5.4 |
| STRING physical: experiments | **0.336** | 0.017 | 0.036 | 9.3 |
| STRING physical: database | 0.061 | 0.003 | 0.006 | 10.6 |
| STRING physical: textmining | 0.168 | 0.008 | 0.024 | 6.9 |
| STRING physical combined ≥ 400 | 0.196 | 0.008 | 0.017 | 11.4 |
| BioGRID any physical | 0.346 | 0.005 | 0.054 | 6.4 |
| **BioGRID "Biochemical Activity"** (in vitro kinase assays etc.) | **0.203** | 0.0004 | 0.008 | **25.8** |
| BioGRID physical excl. Biochemical Activity | 0.269 | 0.004 | 0.049 | 5.5 |
| BioGRID high-throughput only, excl. Biochem. Act. | 0.119 | 0.003 | 0.036 | 3.3 |
| BioPlex 3.0 293T | 0.015 | 0.0003 | 0.002 | 8.3 |
| BioPlex 3.0 HCT116 | 0.008 | 0.0001 | 0.001 | 8.5 |
| HuRI (Y2H) | 0.003 | 0.0002 | 0.0004 | (>8) |
| OpenCell (endogenous-tag IP-MS) | 0.011 | 0.0002 | 0.002 | 7.0 |

Conditional check (KS pairs only, seed 0): fraction with STRING *physical experiments* score > 0 is **0.90** for KS pairs whose *only* BioGRID evidence is "Biochemical Activity" (n = 592), 0.88 for pairs with other BioGRID physical evidence (n = 2,056), and only **0.048** for pairs with no BioGRID record (n = 5,005).

Interpretation:
1. **STRING's "experiments" channel is not leakage-free.** It imports curated BioGRID/IntAct/MINT/DIP/HPRD records, including in vitro kinase assays ("Biochemical Activity" in BioGRID, "phosphorylation reaction"/"enzymatic study" in IMEx). The usual advice "use only STRING experiments + physical mode" does **not** remove KS leakage: 90 % of BioGRID-biochem-only KS pairs carry a STRING experiments score. One third of all PSP KS pairs carry a STRING physical-experiments score, against 1–1.5 % in any single unbiased screen.
2. STRING textmining covers 74 % of KS pairs. STRING 12.5 adds a regulatory network with a "Catalysis of Phosphorylation" subtype, built by an LLM-based literature parser plus SIGNOR, KEGG and Reactome, and benchmarked on SIGNOR phosphorylation **[verified, Szklarczyk et al. NAR 2025 53:D730]**. That network is direct KS label leakage.
3. Unbiased systematic screens (BioPlex, HuRI, OpenCell, Buljan kinase AP-MS, Human Cell Map) show **genuine** enrichment (KS/SHUF ≈ 7–8×) at **low coverage** (≈ 1–2 % of KS pairs each). Buljan et al. 2020 report the same thing: 75 known KS pairs among 7,316 interactions, p < 1e-53. This is real signal: kinases do bind substrates and adaptors. It is sparse, so PPI is a precision feature, not a coverage feature.
4. The KS/SHUF ratio by itself cannot separate signal from leakage. Provenance decides that. Rule: **a PPI feature is admissible only if it comes from a systematic, hypothesis-free screen, or from curated records with all enzymatic/biochemical-activity evidence and all low-throughput literature records removed.** Ablation sanity check: if a curated-PPI feature raises held-out AUC much more than the screens-only feature does, treat the gain as leakage.

---

## 1. Protein–protein interactions

### 1.1 Resources

| Source | Content / size | URL / version | License | IDs | Leakage risk | Expected value |
|---|---|---|---|---|---|---|
| **STRING** | Functional (7 channels: neighborhood, fusion, co-occurrence, coexpression, experiments, database, textmining), physical, and (12.5) regulatory networks. Human: ~19.6k proteins, ~13.7M functional links. Also precomputed ProtT5 and network embeddings | https://string-db.org/cgi/download; web v12.5, files v12.0 (`external/string/`) **[verified]** | CC BY 4.0 **[verified]** | ENSP (aliases → UniProt) | **High**: textmining, database (KEGG, Reactome, GO complexes; regulatory parsing of SIGNOR/KEGG/Reactome), experiments (imports BioGRID/IntAct incl. kinase assays; see §0), regulatory network, and the **precomputed network embeddings**, which mix all channels | Usable features: coexpression channel (GEO, ProteomeHD; leakage-free but study-biased), possibly neighborhood/co-occurrence (weak in human). Use the rest only in a deliberately "leaky upper-bound" ablation |
| **BioGRID** | Curated physical/genetic interactions; experimental system and throughput per record; separate **BIOGRID-PTMS** file (PTM sites with enzymes) | https://downloads.thebiogrid.org ; release 5.0.262 (`external/biogrid/`) **[verified]** | MIT | Entrez, UniProt | **High**: "Biochemical Activity" = "phosphorylation of a substrate by a kinase in vitro" **[verified definition]**; 20 % of PSP KS pairs. PTMS file is a KS/site annotation set (candidate *label* source to check against PSP, never a feature) | After removing Biochemical Activity and keeping High-Throughput only: moderate, mostly duplicating the screens below |
| **IntAct / IMEx** | Curated molecular interactions (PSI-MI), interaction type and detection method, participant roles (enzyme / enzyme target) | https://www.ebi.ac.uk/intact/download ; release 2026-01-14 (`external/intact/`) **[verified date]** | CC BY 4.0 | UniProt | **High** unless filtered: drop interaction type MI:0217 *phosphorylation reaction* (and children/other enzymatic reactions), detection MI:0415 *enzymatic study*, participants with biological role *enzyme / enzyme target*. Also drop low-throughput publications (< N interactions per PMID, e.g. N = 100) | Low after filtering (overlaps BioPlex/HCM/Buljan, which IntAct also hosts) |
| **HIPPIE** | Integrated scored human PPIs (BioGRID, IntAct, MINT, HPRD, DIP, BIND, MIPS...) | https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/ ; v2.4, 2026-04-09 **[verified]** | free academic *(verify)* | Entrez/UniProt | **High**: inherits HPRD/BioGRID biochemical records; scoring rewards literature support | Do not use as a feature (no benefit over the raw screens) |
| **BioPlex 3.0** | AP-MS: 293T 118,162 interactions / 14,586 proteins (10,128 baits); HCT116 70,966 interactions (5,522 IPs) **[verified]** | https://bioplex.hms.harvard.edu ; Huttlin et al. Cell 2021, doi:10.1016/j.cell.2021.04.011 (`external/bioplex/`) | open download (cite) *(verify license)* | Entrez, UniProt | Low (systematic). Caveat: interaction probability (pInt) model trained on CORUM complexes, not KS | Medium precision, low coverage (1.5 % of KS pairs). Kinase–substrate binding in AP-MS is often stable adaptor/scaffold binding (e.g. CDK–cyclin–substrate, mTORC1–RPTOR) |
| **HuRI** | Binary Y2H: ~53k PPIs among ~8.3k proteins (HuRI); HI-union 64k **[verified file sizes]** | http://www.interactome-atlas.org ; Luck et al. Nature 2020, doi:10.1038/s41586-020-2188-x (`external/huri/`) | open *(verify)* | Ensembl gene | Very low for HuRI. **Do not use Lit-BM** (literature binary set) | Low: 0.3 % of KS pairs; Y2H misses transient enzyme–substrate contacts |
| **Buljan et al. 2020 kinase interactome** | AP-MS (Strep-HA, Flp-In T-REx HEK293, dox-inducible) of >300 kinase baits; **7,316 high-confidence interactions, 2,379 proteins**; >5,000 novel. BioID only for validation of selected kinases (PIM3, DYRK2), so it is *not* a kinome BioID map **[verified from full text]** | Mol Cell 2020 79:504, doi:10.1016/j.molcel.2020.07.001, PMC7427327; data in Table S2 | Elsevier OA (CC BY) *(verify)* | gene/UniProt | Low. Score cut-offs were tuned using IntAct literature PPIs (≥ 2 publications) as true positives (mild bias toward known complexes). 20 baits not expressed in HEK293 | **High value per edge**: kinase-centric, covers understudied kinases (90 baits with ≤ 30 citations). The paper raised the share of sites with a kinase assignment from 15 % to 26 % by intersecting with motif predictions. Top PPI source |
| **Human Cell Map (BioID)** | BioID of 192 compartment markers in HEK293; localises 4,145 proteins (NMF compartments) **[verified]** | https://humancellmap.org (download page returned 403 to curl; try browser/API) ; Go et al. Nature 2021, doi:10.1038/s41586-021-03592-2 | free academic *(verify)* | gene/UniProt | Low | Mainly as a **localisation** feature (NMF profiles; prey–prey co-proximity). Few kinase baits |
| **OpenCell** | 1,310 endogenously split-mNG-tagged proteins in HEK293T: live imaging localisation (graded annotations) + IP-MS, 29,922 interactions; abundance for 19k proteins **[verified counts]** | https://opencell.czbiohub.org/download ; Cho et al. Science 2022, doi:10.1126/science.abi6983 (`external/opencell/`) | CC BY *(verify)* | ENSG, UniProt | Low | Medium: endogenous expression, stoichiometry; 1.1 % of KS pairs |
| **Phospho-dependent binding: 14-3-3** | 14-3-3 interactome compendium (ANIA; Tinti et al., Database 2014) and 14-3-3-Pred (Madeira et al., Bioinformatics 2015); ELM LIG_14-3-3_CanoR/CterR regexes **[verified ELM classes]** | https://www.compbio.dundee.ac.uk/1433pred | academic | UniProt + site | Low–moderate (curation sometimes names the upstream kinase; use only binding, never the kinase) | Site-level prior: 14-3-3 sites are mostly basophilic (AKT/PKA/PKC/MARK/CAMK/PKD/RSK). Helps separate basophilic from acidophilic and proline-directed kinases |
| **Phospho-dependent binding: SH2/PTB** | SH2 specificity (peptide chips/arrays; e.g. Tinti et al. Cell Rep 2013 *(verify)*), ELM LIG_SH2_* / LIG_PTB_* regexes **[verified]** | http://elm.eu.org | ELM: free for academic use *(verify)* | site | Low | For Y kinases: substrates whose pY sites match the *kinase's own* SH2 specificity (SFK, ABL, SYK/ZAP70, TEC) are processively phosphorylated. Usable as a kinase-domain-architecture × substrate-motif interaction feature |
| **DirectContacts2** (Drew lab, Nat Commun 2026, doi:10.1038/s41467-026-75863-3) | Integrated map of *direct* physical human PPIs (ML over many datasets + structure) | paper | open *(verify)* | UniProt | Unknown: check whether training uses literature PPIs incl. enzymatic records | Optional. Audit before use |

### 1.2 Leakage filter recipe (recommended)
1. **Primary PPI graph ("screens-only")**: BioPlex 3.0 (293T ∪ HCT116), HuRI (not Lit-BM), OpenCell, Buljan 2020 kinase AP-MS, Human Cell Map prey–bait/prey–prey, plus optionally other proteome-scale screens in IntAct that are defined by *dataset/PMID with ≥ 1,000 interactions*. Record the PMIDs used.
2. **Never use** STRING textmining/database/regulatory channels or STRING network embeddings, HIPPIE, OmniPath, SIGNOR, Reactome/KEGG pathway co-membership, BioGRID "Biochemical Activity" or BIOGRID-PTMS, IntAct enzymatic interactions, or GO MF "kinase binding" (IPI) annotations as features. These may only appear in a clearly labelled "leaky upper bound" ablation.
3. **STRING experiments channel**: do not use the precomputed score. If an aggregate is needed, rebuild it from the screens in (1).
4. **Fold-internal label-derived graph features** (e.g. "fraction of the substrate's neighbours that are known substrates of kinase K", guilt-by-association): compute **inside each CV fold from training labels only**. Otherwise test labels leak through the neighbours.
5. Report a **leakage ablation**: model with screens-only PPI vs. model + curated PPI. Treat a large jump as leakage, not as a result.

### 1.3 Network-propagation features (on the screens-only graph)
- Kinase → protein: random walk with restart (RWR; restart 0.3–0.5) from the kinase node; personalised PageRank; shortest-path length (capped at 4); shared-neighbour count and Adamic–Adar; diffusion-state distance. Precompute one RWR vector per kinase: ~400 kinases × 20k nodes is a sparse matrix solve, minutes on CPU.
- Kinase-complex expansion: propagate from the kinase *and* its regulatory subunits and scaffolds (cyclins, RPTOR/RICTOR, 14-3-3, AKAPs) taken from the screens. This captures "kinase reaches substrate via scaffold".
- For leave-kinase-out, the features stay defined for unseen kinases as long as the kinase is a node in a screen (Buljan covers >300 kinases).
- Graph embeddings (node2vec/GNN) only on the screens-only graph, retrained per release. STRING's precomputed embeddings are excluded (leaky).
- Expected value: in published integrators (e.g. LinkPhinder, NetworKIN, SELPHI, KSMoFinder), network context improves protein-level ranking. At site level it adds little beyond motif, because it cannot tell apart sites on the same protein. Expect gains mainly for **kinase choice within a protein** in combination with motif, and for **understudied kinases with Buljan data**.

---

## 2. Subcellular localisation

| Source | Content / size | URL / version | License | IDs | Leakage risk | Value |
|---|---|---|---|---|---|---|
| **HPA subcellular (IF)** | Antibody IF in cell lines; 13,603 genes with main/additional location, reliability (Enhanced 1,230 / Supported 4,582 / Approved 6,999 / Uncertain 792), single-cell variation, cell-cycle dependence **[verified in file]** | https://www.proteinatlas.org/about/download ; HPA v25.1 (2026-05-25) **[verified]** (`external/hpa/subcellular_location.tsv.zip`) | CC BY-SA 3.0 *(verify)* | ENSG, gene | Very low (experimental) | **Top localisation source**: ~35 compartments, wide coverage. Use the reliability grade as a feature weight |
| **OpenCell localisation** | Live-cell endogenous localisation, 1,310 proteins, graded (3/2/1) annotations **[verified]** | as above | CC BY *(verify)* | ENSG | Very low | High quality, low coverage |
| **Human Cell Map** | BioID NMF compartment profiles, 4,145 proteins (HEK293) | humancellmap.org | academic *(verify)* | gene | Low | Continuous compartment vectors. Good for membrane-less and organelle-contact compartments |
| **Hein et al. 2025 organelle IP-MS** | Global organelle profiling, HEK293T, >7,600 proteins localised; remodelling on infection | Cell 2025 188:1137, doi:10.1016/j.cell.2024.11.028 **[verified via secondary listing]** | Elsevier OA *(verify)* | UniProt | Low | High: proteome-scale, quantitative, 2025 |
| **hyperLOPIT / LOPIT-DC (human)** | Density-gradient spatial proteomics: U2OS hyperLOPIT ~5k proteins, 12+ compartments (Geladaki et al. Nat Commun 2019 *(verify)*); SubCellBarCode 5 cell lines, ~12k proteins (Orre et al. Mol Cell 2019 *(verify)*); in Bioconductor `pRolocdata` | https://bioconductor.org/packages/pRolocdata | Artistic-2.0 (package) | UniProt | Low | Medium. Complements IF. Gives cytosol vs nucleus vs membrane fractions directly |
| **Nuclear/cytoplasmic fractionation** | e.g. SubCellBarCode N/C, CPTAC-style fractionations; Hein 2025 includes nuclear/cytosol markers | | | | Low | Simple N vs C feature. The most useful single axis for kinases (e.g. nuclear CDK/DYRK/HIPK vs cytosolic PKC/AKT) |
| **COMPARTMENTS (Jensen lab)** | Integrated localisation confidence (0–5) from channels: knowledge (UniProt/GO curated), experiments (HPA), text mining, predictions | https://compartments.jensenlab.org ; weekly builds (`external/compartments/` has experiments/knowledge/integrated) **[verified file list]** | CC BY 4.0 *(verify)* | ENSP, gene | Text mining: **no direct KS leakage** (protein–compartment co-mentions, not KS pairs), but strong **study bias**; the integrated score includes text mining. Knowledge channel = GO/UniProt curation (low KS leakage) | Use **knowledge + experiments** channels (as already downloaded), not *integrated* or *textmining*. Value: broad coverage, ontology-normalised (GO CC IDs) |
| **GO CC (GOA)** | Curated + IEA cellular-component annotations (280,334 human CC lines in GAF) **[verified count]** | http://current.geneontology.org (`external/go/`) | CC BY 4.0 | UniProt | Low for CC. (GO BP/MF are a different matter: "protein phosphorylation", "kinase binding" IPI annotations are KS-adjacent leakage; exclude BP/MF) | Use with evidence-code filtering (exclude IEA optionally). Propagate up the ontology to ~30–50 slim terms |

**Feature design.** For each (kinase K, substrate S) compute co-localisation: Jaccard/cosine of compartment vectors (HPA main+additional; HCM NMF; Hein 2025; COMPARTMENTS knowledge+experiments as a 50-term GO-slim vector), plus "shares main location" flags (nucleus, cytosol, plasma membrane, mitochondria, centrosome/spindle, ER/Golgi, cell junctions). Add the kinase's own compartment one-hot as a kinase-side covariate. Missing ≠ zero: add explicit missingness indicators, since coverage correlates with study bias.

Expected value: moderate and broadly useful. Localisation is a necessary condition, and it is the classical ingredient of NetworKIN/KinomeXplorer context. It cannot distinguish sites within a protein. Its main benefit is penalising motif-compatible but spatially implausible kinases (e.g. mitochondrial-matrix sites for cytosolic kinases).

---

## 3. Co-expression, co-abundance, co-dependency

| Source | Content / size | URL / version | License | IDs | Leakage | Value |
|---|---|---|---|---|---|---|
| **GTEx** | RNA-seq, ~54 tissues, ~1k donors | https://gtexportal.org/home/downloadsPage ; v10 *(verify current)* | open access (dbGaP for individual-level; summary data open) | ENSG | None | Tissue co-expression of K and S (Spearman across tissues), and the kinase-expression context (§6) |
| **HPA RNA tissue/cell line consensus** | Tissue consensus nTPM; 1,000+ cell lines (`rna_celline.tsv`, `rna_tissue_consensus.tsv` in `external/hpa/`) | HPA v25.1 | CC BY-SA 3.0 *(verify)* | ENSG | None | Same as GTEx, plus cell-line expression matching the perturbation-phosphoproteomics cell lines |
| **ProteomeHD** | 10,323 proteins × 294 SILAC perturbation ratios; treeClust co-regulation scores **[verified]** | https://www.proteomehd.net ; Kustatscher et al. Nat Biotech 2019, doi:10.1038/s41587-019-0298-5 | CC BY *(verify)* | UniProt | None (but it is also inside STRING coexpression) | Medium: protein co-regulation captures complex membership; the coregulation score is ready-made |
| **CPTAC protein co-abundance** | ~1,000+ tumours/NATs across 10 cancers (TMT); pan-cancer harmonised (Li et al. Cancer Cell 2023 *(verify)*) | https://pdc.cancer.gov ; LinkedOmics; Python `cptac` package | open (CPTAC data use) | gene/RefSeq | None | Medium: tumour co-abundance of K and S. Phospho-side use is covered in the phosphoproteomics review |
| **Tissue co-abundance atlas** | 7,873 proteomic samples, 11 tissues, 114M scored associations (Laman Trip, ... Beltrao, bioRxiv 2024, doi:10.1101/2024.05.15.594301) **[verified]** | preprint data | *(verify)* | UniProt | None | Medium: tissue-specific association probabilities, ready-made |
| **CCLE proteomics (Nusinow 2020)** | 375 cell lines, ~12.7k proteins (TMT) | Cell 2020, doi:10.1016/j.cell.2019.12.023 **[verified]**; DepMap portal | CC BY 4.0 (DepMap) | gene/UniProt | None | Medium: cell-line co-abundance; also kinase abundance per cell line |
| **ProCan-DepMapSanger** | 949 cell lines, ~8.5k proteins (DIA) | Gonçalves et al. Cancer Cell 2022, doi:10.1016/j.ccell.2022.06.010 **[verified]** | CC BY *(verify)* | UniProt | None | As above, more lines |
| **DepMap CRISPR (Chronos)** | ~1,100+ cell lines, genome-wide gene effect; co-essentiality (Pearson/GLS) links pathway partners (e.g. Wainberg et al. Nat Genet 2021 *(verify)*) | https://depmap.org/portal/data_page/ ; latest confirmed release **26Q1** (26Q2 thread exists) **[verified]** | CC BY 4.0 | Entrez/gene | None (CRISPR screens) | Medium for **pathway-level** KS (kinase and substrate co-essential when the substrate mediates the kinase's fitness effect: MTOR–RPTOR–RPS6KB, CDK–cyclins, PLK1 targets). Weak at site level. Only ~10–20 % of kinases have a strong fitness phenotype, so many pairs are flat (use GLS-corrected correlations) |

Feature design: per (K, S) correlation features computed once per resource: RNA tissue Spearman, CPTAC co-abundance (per cohort, then median), ProteomeHD coregulation, CCLE/ProCan co-abundance, DepMap GLS co-essentiality. All are protein-level and leakage-free, but they are protein-level, so they must be used with protein-grouped splits, which are already mandated.

Expected value: low–moderate individually. Co-expression is the weakest classical KS signal: KS pairs are often *not* co-expressed, because kinases are regulated post-translationally. Co-essentiality and co-abundance mostly identify complex/pathway partners. Use them as a few low-weight context features and test via ablation.

---

## 4. Substrate-site sequence and structure

### 4.1 Protein language models (PLMs)

HF availability checked live (2026-10-26) **[verified]**:

| Model | Params / dim | HF repo | Status here | License | Notes |
|---|---|---|---|---|---|
| ESM-2 8M / 35M | 320 / 480 | facebook/esm2_t6_8M_UR50D, esm2_t12_35M_UR50D | **cached** in `$HF_HOME/hub` (31 MB, 131 MB) | MIT | quick tests only |
| ESM-2 650M | 1280 | facebook/esm2_t33_650M_UR50D | downloadable, not gated; not cached | MIT | default workhorse; native `transformers` |
| ESM-2 3B | 2560 | facebook/esm2_t36_3B_UR50D | downloadable (2 shards, ~11 GB fp32); not cached | MIT | fits a 4090 in fp16 |
| ESM-2 15B | 5120 | facebook/esm2_t48_15B_UR50D | downloadable (7 shards, ~60 GB) | MIT | does **not** fit one 24 GB GPU in fp16. Needs sharding/8-bit. Transfer-learning gains over 650M–3B are small (Vieira et al. Sci Rep 2025, doi:10.1038/s41598-025-05674-x **[verified]**). Not recommended |
| **ESM-C 300M / 600M / 6B** | 960 / 1152 / 2560 *(6B verify)*; context 2048 | **biohub/ESMC-300M, biohub/ESMC-600M, biohub/ESMC-6B** (EvolutionaryScale repos now redirect to `biohub/*`) | downloadable, not gated. 600M = 2.3 GB fp32; 6B ≈ 29 GB fp32 (≈ 12 GB bf16) | **Now MIT** (Chan Zuckerberg Biohub, 2026; README tags `mit`+`other`, LICENSE.md is MIT) **[verified]**. *Earlier: Cambrian Open (300M/600M), Cambrian Non-Commercial (6B)* | Repos ship `modeling_esmc_remote.py` (custom code, loaded via `trust_remote_code`). **Rule: do not run downloaded code.** `transformers` ≥ 5.x has a **native ESMC port** used when `trust_remote_code` is off: `<plm-env>` (transformers 5.17.0) supports `model_type=esmc` **[verified]**; `.venv` (4.57.1) does not. So use `<plm-env>` without `trust_remote_code`, or `pip install esm` (PyPI) |
| ESM-3 small open | 1.4B | biohub/esm3-sm-open-v1 | downloadable | MIT (card) **[verified]** | generative multimodal; as an embedder no better than ESM-C. Skip |
| ProtT5-XL-U50 (encoder, half) | 1024 | Rostlab/prot_t5_xl_half_uniref50-enc | downloadable | MIT (ProtTrans GitHub) **[verified]** | strong per-residue embeddings, 1.2B encoder. STRING ships precomputed ProtT5 *protein-level* embeddings (sequence-only → leakage-free, unlike STRING network embeddings) |
| SaProt 650M (AF2) | 1280 | westlake-repl/SaProt_650M_AF2 | downloadable | MIT **[verified]** | needs Foldseek 3Di tokens from AFDB structures (Foldseek not installed; it is a binary from GitHub: request/approve before installing) |

PLM leakage: none from labels (UniRef pre-training only). Caveat for benchmarks: Phosformer-like models fine-tuned on PSP/other KS data must be retrained on our training split only.

**How to use for sites.** (a) Per-residue embedding at the phosphoacceptor (and mean of ±7 window) from **full-length** context. Proteins > 1,022 aa (ESM-2) or > 2,046 aa (ESM-C) need overlapping windows; take the window in which the site is most central. (b) Window-only embeddings (31-mer) are cheaper but lose context. Phosformer used kinase-domain + 15-mer peptide pairs. (c) Fine-tuning (LoRA) a PLM jointly over (kinase domain, site window) pairs is the Phosformer-style alternative, with higher compute.

**Compute estimates (4090, fp16/bf16, estimates to be timed):** human reference proteome ~20.4k canonical proteins, ~11.4M residues; S/T/Y ≈ 1.8M residues (KinoPlex count).
- ESM-2 650M full proteome: ~15–30 GPU-min. ESM-2 3B: ~1–2 GPU-h. ESM-C 600M: ~20–40 GPU-min. ESM-C 6B: ~2–4 GPU-h (bf16, 12 GB weights; long sequences need care at 2,048 tokens).
- Storage for S/T/Y residues only: 1.8M × 1280 × 2 B ≈ 4.6 GB (650M); ×2 for 2,560-d. Store per-site fp16 arrays in `cache/`, keyed by UniProt accession + position + UniProt release.
- Kinase domains (~437–500 kinases × ~300 aa): seconds.

### 4.2 AlphaFold-derived structural features

| Source | Content | URL / version | License | Leakage | Value |
|---|---|---|---|---|---|
| **AlphaFold DB human proteome** | ~23k predicted models (proteins > 2,700 aa split into fragments), pLDDT, PAE | https://alphafold.ebi.ac.uk/download (human tar ~5 GB; v4, check for newer v6 synced with recent UniProt *(verify)*) | CC BY 4.0 | None | **High for site-level "phospho-competence"** |
| **structuremap** (Bludau et al. PLoS Biol 2022, doi:10.1371/journal.pbio.3001636) | Python library: pPSE (prediction-aware part-sphere exposure; full sphere 24 Å & half sphere 12 Å), smoothed IDR flag (`nAA_24_180_pae` → IDR), secondary structure via DSSP, 3D PTM enrichment | https://github.com/MannLabs/structuremap | Apache-2.0 *(verify)* | None | Ready-made pipeline. Bludau showed phosphosites are enriched in IDRs/high-exposure regions and that pLDDT/pPSE separate regulatory sites |
| **KinoPlex** (Vanderwall DR, Huttlin EL, Mintseris J, Yaron-Barir TM, Johnson JL, Dong KD, et al. — Nat Biotech 2026, doi:10.1038/s41587-026-03239-5) **[verified abstract]** | ~20k AF models + PU transfer learning → ~567k "structurally phospho-competent" S/T/Y; intersected with Cantley PSSMs → ~250k high-confidence KS candidates; validated by deep K562 phosphoproteomics; "sequence-structure selective coupling" | paper (paywalled here; check supplementary data) | *(verify)* | Possibly uses known phosphosites (not KS labels) as PU positives. Check whether PSP KS pairs were used anywhere. Structural-competence scores are site-level and label-light | **Important**: both a feature (phospho-competence) and a **benchmark** to compare against |

Recommended structural site features (cheap, CPU): pLDDT at site and mean ±7; pPSE (12 Å half-sphere, 24 Å full-sphere); IDR flag; DSSP secondary structure (H/E/C) and relative solvent accessibility; distance to nearest folded domain boundary; Pfam/InterPro domain context (in kinase domain? activation loop? which domain family?); and residue-level conservation. Conservation sources: ortholog alignments from UniProt/Ensembl Compara, or ESM-derived pseudo-conservation (masked-LM entropy at the site, essentially free from the PLM pass). Activation-loop sites of kinases are a special class: autophosphorylation and upstream-kinase (PDK1, MAP2K, LKB1, CAMKK) sites, which can be flagged with KinCore activation-loop annotations.

Expected value: high for **site plausibility** (which S/T/Y can be phosphorylated at all, a large part of negatives), moderate for *kinase choice*. KinoPlex reports that kinase selectivity partly emerges from the structural accessibility of each kinase's preferred motif, which supports kinase × structure interaction features (e.g. PSSM score × pPSE).

### 4.3 Docking motifs and priming

| Mechanism | Kinases | Representation | Source |
|---|---|---|---|
| MAPK D-motif (ψ1-3–X3-5–φ-X-φ; DOC_MAPK_*; reverse D; FxFP DEF site DOC_MAPK_FxFP_2) | ERK1/2, p38, JNK, ERK5; also RSK/MK2 docking to MAPKs | Count/best-match of DOC_MAPK regexes on the *substrate protein* within ±N residues of the site and protein-wide; distance to site. Kinase side: "MAPK-group" flag → interaction feature | ELM 2024 (356 classes, 4,283 instances; Kumar et al. NAR 2024 52:D442, doi:10.1093/nar/gkad1058) **[verified]**; http://elm.eu.org/elms/elms_index.tsv |
| Cyclin RxL (DOC_CYCLIN_RxL_1), LxF, PxF, cyclin-D helix, CKS1 pT-P (DOC_CKS1_1) | CDK1/2/4/6 + cyclins | RxL match count and distance to site (Cy motif 10–40 aa from the site enhances CDK2-cycA); presence of an *upstream already-phosphorylated* TP for Cks1 (multisite processivity) | ELM; Örd/Loog yeast work |
| PLK1 polo-box: S-pS/pT-P (often primed by CDK1) | PLK1–3 | Substrate contains S-[ST]-P where the [ST] is observed phosphorylated (phospho-site catalogue) → "PBD docking available" flag + distance | literature; no dedicated ELM class found (MOD_Plk_* are phosphosite motifs only) **[verified]** |
| GSK3 priming: pS/pT at +4 (MOD_GSK3_1 `...([ST])...[ST]`) | GSK3A/B | (i) sequence flag S/T at +4; (ii) **observed** phosphorylation at +4 (from PSP *Phosphorylation_site_dataset* or the perturbation compendium; site catalogue, not KS labels); (iii) predicted priming-kinase score at +4 (max Cantley score over CDK/MAPK/CK2/DYRK) | Johnson 2023 PSPA matrices include phospho-primed positions |
| CK1 priming: pS/pT at −3 (MOD_CK1_1 `S..([ST])`), or acidic stretch; non-canonical SLS motif | CK1α/δ/ε/γ | Same three-level scheme at −3; acidic −3/−4 count | same |
| CK2 hierarchical (pS at +1/+3 acting as acidic determinant) | CK2A1/A2 | phosphorylated residue counted as acidic in window features | Kinase Library |
| AGC PIF pocket (DOC_AGCK_PIF_*) | PDK1 → AGC kinases (S6K, SGK, RSK) | on substrate kinase: hydrophobic motif present → PDK1-substrate flag | ELM |
| SRC-family/ABL processivity via SH2/SH3 | SFK, ABL, SYK, ZAP70, TEC | substrate has pY matching kinase SH2 specificity (LIG_SH2_SFK) or PxxP (LIG_SH3) → interaction feature | ELM |
| 14-3-3 binding | basophilic kinases | LIG_14-3-3_CanoR_1 match at the site → site prior for basophilic group | ELM, 14-3-3-Pred |
| PP1/PP2A/calcineurin docking (DOC_PP*) | phosphatases (negative/context) | optional context | ELM |

**How to represent priming in features (recommendation):**
1. *Sequence-only flag* (always available): S/T at +4 (GSK3), −3 (CK1), etc.
2. *Observed-phosphorylation flag*: whether the priming position is a detected phosphosite in a large site catalogue (PSP Phosphorylation_site_dataset or the provided compendium's site list; independent of KS labels; record counts of supporting studies). This is the most informative feature but carries study bias, so pair it with the sequence flag.
3. *Priming-kinase score*: max motif score (Cantley PSSMs) of plausible priming kinases at the priming position, giving a soft "primeable" probability.
4. *Phospho-aware sequence scoring*: when scoring with Cantley/Kinase Library PSSMs, pass primed residues as phosphorylated (the Kinase Library accepts lowercase s/t/y as phospho-priming residues *(verify API)*). For PLMs (no phospho token), add an explicit per-position channel. Alternatively, compute a second embedding with the primed S/T replaced by D/E (phosphomimetic) and use the difference.
5. *Kinase-side gating*: the model should learn kinase × priming interactions. In a one-model-for-all design, use kinase-family embeddings × priming features (or a per-kinase-group head). In GBDTs, add explicit product features for GSK3/CK1/PLK/CK2 groups.
6. *Evaluation*: report a dedicated "primed-kinase" slice (GSK3A/B, CK1s, PLK1–3, CK2, CDK multisite) as a benchmark subset.

**ELM leakage note:** ELM *instances* (and MOD_* instances) are literature-curated phosphosites, many with the kinase named. Use **class regexes only**; never use ELM instance membership as a feature.

---

## 5. Kinase-side features (sharing across kinases; leave-kinase-out / zero-shot)

| Source | Content | URL / version | License | Leakage | Value |
|---|---|---|---|---|---|
| **Kinase domain sequences / boundaries** | UniProt kinase-domain features; **KinCore** (Dunbrack): domain alignment, conformational classes, activation-loop boundaries, 437 catalytic typical kinases; **active-state AF2 models of all 437** (Gizzio, Faezov, Xu, Dunbrack, bioRxiv 2026, doi:10.64898/2026.02.19.706771) **[verified]** | https://dunbrack.fccc.edu/kincore/ (and /kincore/activemodels) **[verified URL]** | academic *(verify)* | None | **High**: the input for kinase PLM embeddings and specificity-determinant residues |
| **KinBase / Manning tree** | 518 (later ~540) human kinases, group/family/subfamily (AGC, CAMK, CK1, CMGC, STE, TK, TKL, Other, atypical) | http://kinase.com/web/current/kinbase/ (502 at time of check; mirror via Coral / UniProt "pkinfam") | academic | None | High: hierarchical one-hot / tree distances; enables partial pooling |
| **Kinase PLM embeddings** | ESM-2/ESM-C mean-pooled over kinase domain (and full length), ~400–540 kinases | computed locally | — | None | High for zero-shot: the Phosformer/DeepKinZero idea (Deznabi et al. Bioinformatics 2020, zero-shot via kinase embeddings *(verify)*) |
| **Specificity-determining residues (SDRs)** | Creixell et al. Cell 2015 *(verify)* SDR positions in the kinase domain alignment; one-hot residues at ~20–30 SDR positions | from KinCore alignment | — | None | Medium–high: interpretable kinase descriptors that generalise across families |
| **Cantley PSSMs as kinase descriptors** | Johnson et al. Nature 2023 (303 S/T kinases, doi:10.1038/s41586-022-05575-3); Yaron-Barir et al. Nature 2024 (93 Tyr kinases, 78 canonical + 15 non-canonical *(verify)*); Kinase Library | https://kinase-library.phosphosite.org **[verified URL]** | academic use; check terms | **Not PSP-derived** (positional scanning peptide arrays), so no label leakage. Johnson et al. used PSP only to validate | **Very high**: (i) site feature (score/percentile per kinase), (ii) kinase descriptor (flattened PSSM, ~10 positions × 23 residues), (iii) kinase–kinase similarity (PSSM correlation) for borrowing strength. Coverage ≈ 396 kinases → nearly all PSP kinases |
| **Kinase–kinase similarity kernels** | domain sequence identity, PSSM correlation, Manning-tree distance, PLM cosine | derived | — | None | For multi-task / kernel models: understudied kinases borrow from neighbours. Recommend multi-kernel (sequence + PSSM) |
| **Kinase "darkness"** | IDG dark-kinase list (Pharos/TCRD), PubMed counts | https://pharos.nih.gov | CC BY | None | For stratified evaluation (not as a feature) |

Leave-kinase-out design: the kinase must be represented by features available without its KS labels: PSSM (if profiled), domain PLM embedding, SDR one-hot, family/group, screens-only PPI neighbourhood, localisation and expression. Johnson/Yaron-Barir PSSMs are the strongest zero-shot descriptor. For non-profiled kinases (atypical, pseudokinases with residual activity), fall back to embedding/family.

---

## 6. Kinase activity / expression context

- **Which kinases are expressed where**: HPA tissue & cell-line RNA (v25.1), GTEx (tissues), CCLE/DepMap expression (cell lines; 26Q1), CPTAC/ProCan protein abundance. Uses: (a) a context feature for interpreting the provided perturbation phosphoproteomics (a site regulated in cell line X cannot be explained by a kinase not expressed in X); (b) a weak prior on KS plausibility (K and S both expressed in ≥ 1 common tissue/cell type); (c) filtering decoy negatives (do not treat "not annotated" as negative when the kinase is never co-expressed with the substrate in any assay context).
- **DepMap gene effect / co-dependency**: see §3. Kinase essentiality profiles also act as kinase descriptors (similar fitness profiles → shared pathway).
- **Kinase activity atlases** (inferred from phosphoproteomics, e.g. CPTAC kinase activities) belong to the phosphoproteomics review. Avoid circularity: activity inferred *from PSP substrates* leaks if then used to predict those substrates.
- Expected value: low as a direct KS feature, **high for the downstream use case** (cell-type-specific predictions and kinase activity estimation). Recommend offering an optional "context mask" at inference time rather than baking expression into the core model.

---

## 7. Structure-based KS prediction with AlphaFold-Multimer / AF3 (2024–2026)

What exists:
- **Kim, Huang, Johnson, Yaron-Barir, ..., Cantley, Hu, Perrimon (bioRxiv 2025, doi:10.1101/2025.10.10.681672)** **[verified]**: AlphaFold-Multimer v2.3.1 (LocalColabFold 1.5.5; 5 models × 5 recycles) on **>240,000 human pairs (207 S/T kinases × 1,162 TFs)** plus ~212k Drosophila pairs. Scored with iLIS, a PAE-local interface score: AUC 0.872 vs ipTM 0.822 on 774 positives / 6,004 controls, and 0.890 vs 0.776 in disordered complexes; high-confidence threshold iLIS ≥ 0.223 (10 % FDR). Recovers PSPA-like preferences and predicts a homeodomain docking motif for basophilic kinases. Run on a large heterogeneous HPC cluster (V100…H100); GPU-hours not reported. No site-level KS benchmark against Cantley motifs was reported. Data: Zenodo 10.5281/zenodo.17314335, code https://github.com/flyark/AFM-LIS **[verified URL]**.
- **KinoPlex (Nat Biotech 2026)**: structure (single-chain AF models) + PSSMs. This is the scalable structural route (§4.2), not complex modelling.
- **Active-state kinase AF2 models (Dunbrack lab 2026; AFEX, bioRxiv 2025.11.03.686188)**: active-conformation models of 436–437 kinases. AF-Multimer *can* model substrate-peptide-bound complexes with good templates and MSAs, but activation loops are not always active. Useful templates for peptide docking.
- **PTK specificity preprint (bioRxiv 2026, doi:10.1101/2026.05.10.724103)** *(details verify)*: PSPA-style tyrosine kinase profiling with decoys plus AF3 kinase–substrate complex models; notes AF3 training-set bias.
- **Cyclin–CDK docking with AF-Multimer** (fission yeast Cdc13-Cdc2 substrates; Willich, Kapadia, Nurse, bioRxiv *(verify)*): focused, not proteome-wide.
- **Subtimizer** (J Chem Inf Model 2026, doi:10.1021/acs.jcim.5c02430): AF-Multimer + ProteinMPNN substrate *design*; ipTM > 0.75 filter. Not a KS predictor.
- AF3 on phosphoproteins (bioRxiv 2025.04.14.648669): explicit phospho-residues gave only modest improvement over AF2/AF3 without phosphorylation. PepPCBench (bioRxiv 2025.04.08.647699): AF3-class peptide-complex accuracy is moderate, and confidence metrics are imperfect discriminators.
- I found no published kinome-wide AF-Multimer/AF3 screen that **benchmarks site-level KS prediction against PSP** with a positive result. Structure-based complex scoring therefore remains unproven for this task. The most convincing gains come from (i) docking-interface discovery (protein-level) and (ii) single-structure site accessibility (KinoPlex, structuremap).

Feasibility on 4× RTX 4090 (estimates; precomputed MSAs; MSAs for ~437 kinase domains computed once; substrate peptides run single-sequence):

| Setting | Pairs | Time/model | GPU-hours | Wall time (4 GPUs) | Verdict |
|---|---|---|---|---|---|
| Full-length K × full-length S (AF-M/Boltz-2), 1 model | 400 × 5,000 = 2M | 1–3 min | 33k–100k | 1–3 years | infeasible |
| Kinase domain × 15–25-mer site peptide, 1 model, 3 recycles | 400 × 10,000 sites = 4M | ~8–15 s | 9k–17k | 3–6 months | infeasible |
| Same, subsampled: 60 kinases (balanced across groups) × (100 PSP positives + 200 matched negatives) | 18k | ~10 s | ~50 | ~13 h | **feasible pilot** |
| Kinase domain × 200-aa substrate window (to capture docking motifs), 400 pre-ranked candidates per kinase from the motif model | 160k | ~30–60 s | 1.3k–2.7k | 2–4 weeks | possible only as a re-ranker for top candidates, if the pilot shows value |

Tools: AF3 weights require a request to Google DeepMind and have non-commercial terms (do not wait on it). Open alternatives with MIT/Apache licences: **Boltz-2**, **Protenix**, **Chai-1**; all support modified residues (SEP/TPO/PTR CCD codes) for priming experiments. ColabFold/AF-Multimer: Apache-2.0 code, CC BY 4.0 weights. Installing any of these means downloading code (repository/pip); this needs PI approval under the project rules.

Recommendation: run a **pilot** (≤ 1 GPU-day) measuring whether iLIS/ipSAE/ipTM of kinase-domain × peptide models adds AUROC/AUPRC *beyond* Cantley PSSM + structure features on held-out proteins. Include primed vs unprimed peptides for GSK3/CK1 (phosphoserine at +4/−3), the case where structure should help most. Abandon if the gain is < ~0.01 AUPRC. Watch for memorisation: kinases with peptide co-crystals in the PDB (PKA, CDK2, AKT, ...) will look better than dark kinases.

---

## 8. Prioritised recommendations

**Tier 1 (implement first)**
1. **Cantley PSSMs** (Johnson 2023 + Yaron-Barir 2024) as site scores and as kinase descriptors/similarities (no label leakage; ~396 kinases).
2. **Site structure features** from AFDB + structuremap (pLDDT, pPSE, IDR, DSSP) plus PLM pseudo-conservation. Compare against KinoPlex scores if their data are obtainable.
3. **PLM embeddings**: ESM-2 650M per-site + kinase-domain embeddings (native transformers, MIT). Then test ESM-C 600M (MIT now, native port in `<plm-env>`) and ESM-2 3B; choose by validation. Skip 15B and ESM-3.
4. **Screens-only PPI graph** (BioPlex 3.0 ∪ HuRI ∪ OpenCell ∪ Buljan 2020 ∪ Human Cell Map) with RWR/shared-neighbour features; fold-internal label-propagation features only.
5. **Localisation co-occurrence**: HPA v25.1 IF + OpenCell + Hein 2025 + COMPARTMENTS knowledge/experiments (not textmining/integrated), with missingness flags.
6. **Priming/docking features**: three-level priming representation; ELM DOC/LIG class regexes (never instances).

**Tier 2**: co-abundance/co-expression (ProteomeHD, CPTAC, CCLE/ProCan, tissue atlas), DepMap co-essentiality, kinase SDR one-hots, kinase expression context masks.

**Tier 3 / pilot only**: AF-Multimer/Boltz-2 kinase-domain × peptide pilot (~13 GPU-h), docking-window re-ranker if the pilot is positive.

**Excluded as features (leaky)**: STRING textmining/database/regulatory/experiments channels and STRING network embeddings; HIPPIE; BioGRID Biochemical Activity and BIOGRID-PTMS; IntAct enzymatic interactions; HuRI Lit-BM; Reactome/KEGG/SIGNOR/OmniPath; GO BP/MF; ELM instances. These belong in a clearly labelled leaky-upper-bound ablation, or are candidate *label* sources to be vetted separately.

## 9. Open items / requests
- Human Cell Map download returned HTTP 403 to curl. Try a browser or the API, or ask the user.
- Licenses marked *(verify)* should be confirmed when each dataset is recorded in DATA.md. The files in `external/` (STRING, BioGRID, IntAct, BioPlex, HuRI, OpenCell, HPA, COMPARTMENTS, GO, OmniPath) were downloaded by another agent and still need DATA.md entries.
- Installing Foldseek (for SaProt), Boltz-2/ColabFold and `pip install esm` needs PI approval (downloaded code).
- KinoPlex supplementary data (structural phospho-competence per site) would be a valuable feature and benchmark. Check availability and whether PSP KS pairs were used in its training.
