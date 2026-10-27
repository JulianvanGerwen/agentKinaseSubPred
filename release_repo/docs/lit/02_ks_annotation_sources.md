> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# 02 — Sources of experimentally verified kinase–substrate (site-level) annotations, and their quality

Author: literature agent (KS-sources), 2026-10-26. Scope: human (+ mouse/rat orthologous) **site-level** kinase→phosphosite (KS) annotations usable as labels, test sets or auxiliary features for the project (see `docs/AGENT_BRIEF.md`).
Project rule recap: *labels (train/val/test) must be literature-derived, experimentally verified KS annotations.* High-throughput (HT) screens, predictions and text-mining output are treated below as candidate **features** or **external sanity checks**, never as gold labels, unless explicitly argued otherwise.

All URL checks were made from this sandbox on 2026-10-26 with `curl` (HTTP code and `Content-Length` noted; nothing large downloaded). "Blocked here" means the sandbox proxy refused the connection; it is not evidence the resource is offline.
Numbers that I computed locally come from `scripts/lit/ks_source_overlap_stats.py` (reproducible; inputs listed in §9).

---

## 0. TL;DR verdict table

| Resource | Evidence type | Human site-level KS (approx.) | Net new vs PSP | License | Access (2026-10-26) | Train labels | Test labels | Features only |
|---|---|---|---|---|---|---|---|---|
| **PhosphoSitePlus KS dataset** (provided, release 2025-05-19) | LT literature, manual; in vivo + in vitro flags | 14,724 rows human→human; 23,904 all species | — | CC BY-NC-SA 3.0 | provided; downloads need free login (403 here) | **YES (core)** | **YES (core)** | — |
| **SIGNOR 4.0 / PhosphoSIGNOR** | LT literature, manual, causal; residue given in 76% of phospho reactions | 15,718 (de)phospho reactions; 6,903 sites, >400 kinases | ~3.1k human rows not in PSP (via OmniPath) | CC BY 4.0 | blocked here (proxy); OmniPath mirror works | **YES (after filtering)** | YES (as "external-curator" test) | — |
| **UniProtKB/Swiss-Prot MOD_RES "by KINASE"** | LT literature, manual (ECO:0000269), some by-similarity | 3,872 human sites with a "by" | ~960 exp.-evidence kinase–site pairs not in PSP (incl. 185 FAM20C) | CC BY 4.0 | REST API works | **YES (ECO:0000269 only, gene-level kinase names only)** | YES | — |
| **Phospho.ELM 9.0** | LT literature (+HTP sites without kinase) | ~2–3k KS (frozen 2010) | tiny (5 single-source rows in OmniPath) | academic, no redistribution, click-through | dump reachable (6.4 MB) | optional (redundant) | no | — |
| **HPRD 9** | LT literature (frozen 2010) | ~9.7k OmniPath rows touch HPRD | 17 single-source rows | academic | **offline (503)** | optional (via OmniPath) | no | — |
| **iPTMnet 6.2** | aggregator (PSP, UniProt, SIGNOR, PRO, …) + RLIMS-P/eFIP text mining | 24k+ enzyme–substrate–site (all PTM, all species) | text-mined part only | CC BY-NC-SA 4.0 (readme also says CC BY 4.0 — inconsistent) | ptm.txt 44 MB, score.txt 26 MB, OK | curated sub-sources only (PRO curated) | no | text-mined part as weak feature |
| **OmniPath enzyme–substrate** | aggregator of ~16 resources | 39,037 human phospho rows, 1,649 "enzymes" | 5,792 rows from LT curated DBs not in PSP; rest text-mining/HT/KEA | per-source (PSP part NC-SA) | API works (4.9 MB TSV) | **YES as a convenient source of SIGNOR/HPRD/PELM/NCI-PID/Reactome rows** (with filters) | partially | text-mining/HT rows as features |
| EPSD 2.0 / dbPTM 2025 / qPhos | aggregators of sites; KS mostly re-distributed PSP/UniProt/PELM | — | ~0 new curated KS | academic / unspecified | EPSD/GPS pages return empty to curl; dbPTM download page OK; qPhos needs registration | no | no | phosphosite universe / evidence counts (feature) |
| RegPhos 2.0 | aggregator (2014) | — | ~0 | academic | **offline** | no | no | no |
| **PhosphoNetworks (Newman 2013)** | HT in vitro protein microarray; **sites are motif-predicted** | 4,417 kinase→site (predicted site) | 3,283 single-source OmniPath rows | free academic/non-profit | OK (rawKSI 293 kB, refKSI 44 kB, comKSI 74 kB, highRes 154 kB) | **no** | no | protein-level in vitro KSR as feature |
| **Sugiyama 2019 / iKiP-DB** | HT in vitro, dephosphorylated HeLa lysate + 385 kinases, MS sites | 175,574 kinase–site pairs, 20,669 sites | ~96% not in PSP | CC BY 4.0 (Sci Rep SI); iKiP-DB in PTMsigDB v2 | SI + jPOST PXD011366 | **no** | no (at most "orthogonal support" analysis) | **YES (strong feature)** |
| **Kinase Library (Johnson 2023; Yaron-Barir 2024)** | HT in vitro peptide-library PSSMs (no proteins) | 303 S/T + ~93 Y kinases; scores for all sites | n/a (motifs) | CC BY-NC-SA 3.0 (python pkg) | site + pip `kinase-library` 1.8.0 | no | no | **YES (feature; also baseline to beat)** |
| GPS 5.0/6.0 "gold" set | aggregator (PSP+PELM+UniProt+…) | 30,043 KS sites (GPS 6.0, all species) | small; Müller-Dott: 544 unique interactions | academic | GPS download page lists only software | no (unclear provenance, contains PSP) | no | — |
| KinaseNET / PhosphoNET (Kinexus) | proprietary predictions + some curation | ~650k predicted sites | — | proprietary | kinasenet.ca blocked here; phosphonet.ca OK | no | no | no (closed, prediction) |
| RLIMS-P / eFIP / REACH / Sparser (text mining) | literature, automatic | thousands | ~5.8k text-mining-only OmniPath rows | varies | via OmniPath/iPTMnet | **no** | no | weak feature / curation queue |
| Chemical genetics (analog-sensitive, ATPγS) | direct, in lysate/cells, MS sites | tens–hundreds per kinase, ~10 human kinases | partly already in PSP | paper SI | paper SI | **yes if curated into PSP/SIGNOR; otherwise held-out "direct" test** | **YES (high-value direct test)** | — |
| iKinSub ("KinSub") | HT in vitro, rat IMCD lysate, few kinases | small | — | public | — | no | no | marginal |
| PhosPhAt | plant (Arabidopsis) | — | irrelevant | — | blocked here | no | no | no |

---

## 1. PhosphoSitePlus (PSP) — the provided `data/Kinase_Substrate_Dataset`

**Reference.** Hornbeck et al. 2012 NAR doi:10.1093/nar/gkr1122; 2015 NAR doi:10.1093/nar/gku1267; 2019 NAR doi:10.1093/nar/gky1159. Provided file: release header "051925" (2025-05-19), SHA-256 in DATA.md. License CC BY-NC-SA 3.0 (non-commercial; share-alike affects any redistributed derivative, e.g. a released label file or model trained on it should carry NC-SA terms — flag for REPORT.md). Static downloads (`Phosphorylation_site_dataset.gz`, `Regulatory_sites.gz`, etc.) require free registration; `https://www.phosphosite.org/staticDownloads.action` returns 403 to this sandbox.

**Content (computed, `ks_source_overlap_stats.py`).**
- 23,904 rows; human kinase→human substrate 14,724 rows (14,657 unique kinase-acc/substrate-acc/site keys), 421 human kinases, ~10k human sites, 3,103 substrate accessions.
- Residues: S 15,428, T 4,933, Y 3,526 plus 17 non-STY rows (H/K/R/D/I — histidine kinases etc.; **exclude**).
- Skew: top-10 kinases hold 34.7% and top-20 hold 49.7% of human→human rows; 126 kinases have ≥20 rows, 275 have ≥5, 146 have <5. Most-annotated substrate is tau `P10636-8` (241 rows), then TP53 (117), CTNNB1 (72).
- Species: human→human 14,724; human kinase→mouse substrate 2,426, →rat 1,282; mouse→mouse 2,334; rat→rat 983; plus cow/chicken/rabbit/viral.

**Fields.**
- `IN_VIVO_RXN` / `IN_VITRO_RXN` ("X" if any curated paper supports that context). PSP does not publish a formal definition, but in practice: *in vitro* = purified/recombinant kinase phosphorylates substrate/peptide (direct but potentially promiscuous); *in vivo* = the site's phosphorylation in cells/tissue depends on the kinase (overexpression, knock-down/out, inhibitor, activation), often combined with site mutants — **causal but not necessarily direct**. Counts (all rows): vitro-only 10,463, vivo-only 6,689, both 6,752. Human→human: vitro-only 5,885 (40%), vivo-only 3,668 (25%), both 5,171 (35%). PSP itself (2015) reported "in vitro outscore in vivo 11,396 to 6,231".
- `KIN_ORGANISM` / `SUB_ORGANISM`: species of the kinase protein used and of the substrate. Cross-species experiments exist (e.g., human kinase on rat substrate, 1,282 rows) — these are mostly in vitro with recombinant human kinase.
- `SITE_GRP_ID`: PSP "Modification Site Group" — a single ID shared by orthologous/isoform-equivalent sites across species and isoforms (PSP 2015: "Residue numbers cannot be used because they can vary between isoforms or orthologs"). This is PSP's own curated cross-species site map and is the cleanest way to pool mouse/rat annotations onto human sites.
- `SUB_ACC_ID` with isoform suffix (e.g. `P18433-2`, `P10636-8`): PSP stores each protein under its "parent" isoform, usually the **longest** isoform, which is often *not* UniProt canonical (tau: PSP uses 441-aa 2N4R `P10636-8`, UniProt canonical P10636 is the 758-aa PNS-tau). 1,239 rows (893 human→human) use isoform accessions; 239 kinase accessions also carry isoform suffixes; 8 kinase accessions are RefSeq/GenBank (`NP_…`, `AAA58698`, `YP_068414` viral). **Action:** map every site to UniProt canonical by locating the `SITE_+/-7_AA` 15-mer (lower-case = phospho) in the canonical sequence (exact, then gapped/fuzzy), never by residue number; drop/flag sites absent from canonical.
- `GENE` vs `KINASE`: `KINASE` is a PSP display name (e.g. "ERK1", "P38A", "CK2A1"); `GENE` is the official symbol but is blank for 7 kinase names — map via `KIN_ACC_ID`.

**Quality concerns specific to PSP.**
1. *Indirect in vivo annotations.* Of 65 S/T rows attributed to kinases whose annotations are ≥80% Tyr (43 such kinases, e.g. EGFR has 5 S/T of 101 rows), **45 (69%) are in-vivo-only**, versus a 25% in-vivo-only baseline — i.e. biochemically implausible annotations are ~2.8× enriched in the in-vivo-only class, consistent with downstream (indirect) effects being curated as KS. Petsalaki lab (Maier et al. 2025 MCP doi:10.1016/j.mcpro.2025.100994) makes the same point for EGFR serine sites.
2. *In vitro promiscuity.* Y sites attributed to S/T-dominant kinases (46 rows; CSNK2A1, CHEK1, GSK3B, LRRK2…) are ~50% in-vitro-only; some are real (GSK3B Y216 autophosphorylation), others plausibly artefacts.
3. *Family-level ambiguity.* Papers often use one isoform as representative (PRKACA for "PKA", PRKCA for "PKC"); the true in-vivo kinase may be a paralog. Treat paralog confusion as label noise; evaluation should report family-level as well as kinase-level metrics.
4. *Study bias* (see §6).

**Verdict.** Core train **and** test labels. Keep both evidence classes but carry `evidence ∈ {vivo_only, vitro_only, both}` as a label attribute: (i) primary metric on all; (ii) report metrics on the "both" and "vitro" (direct) subsets; (iii) consider down-weighting vivo-only or excluding residue-type-inconsistent pairs (an explicit, pre-registered filter, not post hoc).

---

## 2. Other manually curated literature resources (candidate label expansions)

### 2.1 SIGNOR 3.0 / 4.0 and PhosphoSIGNOR
- Refs: SIGNOR 3.0, Lo Surdo et al. 2023 NAR doi:10.1093/nar/gkac883; **SIGNOR 4.0**, Lo Surdo et al. 2026 NAR 54:D682 doi:10.1093/nar/gkaf1237 (PMID 41243950, online 2025-11-17).
- Content (4.0, Aug 2025): 40,940 causal relations, 15,718 (de)phosphorylation reactions; residue specified for 76% of phosphorylation reactions; PhosphoSIGNOR: >400 kinases, 90 phosphatases, 2,716 substrates, 6,903 unique phosphosites/peptides. Human-centred (also some mouse/rat curated with `TAX_ID`). Curation adds INDRA text-mining-assisted triage but every relation is curator-verified.
- Evidence: low-throughput literature; each relation has `MECHANISM=phosphorylation`, `RESIDUE`, `SEQUENCE` (±7), `TAX_ID`, `CELL_DATA`/`TISSUE_DATA`, `DIRECT` (YES/NO, the curator's judgement whether the interaction is direct), `SENTENCE`, `PMID`, `SCORE`. **Keep only `DIRECT=YES`, mechanism phosphorylation, residue non-empty, entity type protein (drop complexes/families or expand them carefully).**
- Overlap: From OmniPath (which re-distributes SIGNOR), 3,102 human SIGNOR-supported phospho rows are not in the provided PSP file (628 are SIGNOR-only across all OmniPath sources). Many SIGNOR entries are PSP-derived curation of the same papers; the new ones are a real but modest (~+15–20%) expansion.
- License: CC BY 4.0 (SIGNOR site; OmniPath records "CC BY 4.0, commercial OK") — the most permissive curated KS source.
- Access: `https://signor.uniroma2.it/downloads.php`, PhosphoSIGNOR TSV/JSON at `https://signor.uniroma2.it/PhosphoSIGNOR/downloads/` and `/apis/`. **Blocked from this sandbox** (proxy refused). Fallback: OmniPath enzsub (works) or ask user to fetch the TSV.
- Verdict: **train labels (yes, after filters)**, and good **independent-curator test** (sites only in SIGNOR, post-dating the PSP release, can serve as a temporal hold-out). Watch for complexes ("CDK1/CCNB1", "AMPK" subunits PRKAG1/PRKAB2 show up as enzymes in OmniPath) and phosphatase rows mislabelled as phosphorylation (PTPN1, PPP2CA appear among "phosphorylation" enzymes in OmniPath).

### 2.2 UniProtKB/Swiss-Prot — MOD_RES "Phosphoserine; by KINASE"
- Query used: `rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606) AND (reviewed:true) AND (ft_mod_res:by)&fields=accession,gene_primary,ft_mod_res&format=tsv` (1.6 MB).
- Content: 1,618 human reviewed entries, 3,872 phospho sites with a "by" kinase; 630 "autocatalysis"; many **family names** ("PKA" 243, "CK2" 184, "PKC" 152, "AMPK", "MAPK", "CK1") that cannot be mapped to a single gene. Evidence codes: ECO:0000269 (experimental, with PMID) on 2,959 sites; ECO:0000250 (by similarity, i.e. propagated from an ortholog) 621; ECO:0007744 (large-scale MS, mostly the site not the kinase) 691; ECO:0000255/0000305 (predicted/curator inference) 285.
- Overlap with PSP: of 2,567 rows whose kinase string equals a PSP human kinase gene symbol, 1,292 are already in PSP; **1,275 are not, 960 of them with ECO:0000269** (top: FAM20C 185 — largely Tagliabracci et al. 2015 secretome work; CDK1 67, MAPK1 47, PDPK1 30). ~297 of the new pairs hit sites PSP already links to other kinases. Some "new" pairs may be numbering mismatches (isoform), so re-map by sequence.
- License CC BY 4.0. Format TSV/JSON/flat file.
- Verdict: **train/test labels for ECO:0000269 rows with an unambiguous gene-level kinase**; drop by-similarity (orthology-propagated — would double-count rodent data), autocatalysis (handle separately as an autophosphorylation task), and family names (or keep as family-level labels only). "Catalytic activity" (Rhea) annotations on kinase entries are generic reactions (`L-seryl-[protein] + ATP = …`) and contain **no** substrate sites — not useful.

### 2.3 Phospho.ELM 9.0
- Diella et al. 2008 NAR doi:10.1093/nar/gkm772; Dinkel et al. 2011 NAR doi:10.1093/nar/gkq1104. Last release 9.0 / dump Sept 2010 (8,718 substrates; 42,573 instances incl. HTP sites; kinase only for a minority, LTP-only).
- Access: `http://phospho.elm.eu.org/dumps/phosphoELM_all_latest.dump.tgz` (6.4 MB, Last-Modified 2018-06-29) is reachable, but requires accepting an **academic license (no commercial use, no redistribution to third parties)** via form.
- Overlap: almost entirely subsumed by PSP — in OmniPath only 5 human phospho rows are Phospho.ELM-only. Sequences referenced to UniProt 12.3 (2007) → numbering drift.
- Verdict: redundant; not worth the license step. If used, only as extra provenance (PMIDs).

### 2.4 HPRD (Release 9, 2010)
- Keshava Prasad et al. 2009 NAR doi:10.1093/nar/gkn892. `http://www.hprd.org/` → **503** (effectively offline for years). OmniPath has 9,691 human rows touching HPRD but only 17 HPRD-only rows; 1,050 HPRD-supported rows absent from PSP (often coarse/indirect, e.g. kinase "family" assignments, old numbering).
- Verdict: use only via OmniPath, only rows corroborated by a second curated source; not as test labels.

### 2.5 iPTMnet (release 6.2) and PRO
- Huang et al. 2018 NAR doi:10.1093/nar/gkx1104. Integrates curated DBs (UniProt, PSP, PRO, SIGNOR, HPRD, Phospho.ELM, PhosphoGRID, P3DB…) and PIR text mining (RLIMS-P, eFIP). Reported scale ~654k PTM sites, >1,200 enzymes, >24k enzyme–substrate–site relations (all PTMs, all species).
- Files (`https://research.bioinformatics.udel.edu/iptmnet_data/files/current/`): `ptm.txt` 44.1 MB (Last-Modified 2025-10-14; columns ptm_type, source, substrate AC/gene, organism, site, enzyme AC/gene, note, pmid), `score.txt` 26.1 MB (2024-01-31), `protein.txt`. A sampled 18 MB slice of the phosphorylation block showed sources `sign` (SIGNOR), `unip` (UniProt), `rlim+` (RLIMS-P text-mined), `pro` (PRO curated proteoform; enzymes sometimes as PR: IDs = complexes/proteoforms).
- License: readme states **CC BY-NC-SA 4.0** at top and CC BY 4.0 at bottom — inconsistent; assume NC-SA.
- Verdict: no new curated content beyond its sources except PRO (small). Text-mined `rlim+` rows → features/curation queue only.

### 2.6 OmniPath enzyme–substrate (Türei et al. 2016 Nat Methods doi:10.1038/nmeth.4077; 2021 Mol Syst Biol doi:10.15252/msb.20209923)
Downloaded human phosphorylation rows (`https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references,curation_effort&types=phosphorylation&organisms=9606`, 4.9 MB TSV, SHA-256 in §9): **39,037 rows, 1,649 "enzymes"** (includes complex subunits, cyclins, phosphatases).
- Rows per primary source: ProtMapper 20,863; PhosphoSite 18,677; MIMP 12,539; KEA 10,999; HPRD 9,691; phosphoELM 9,195; SIGNOR 8,670; Sparser 5,427; REACH 5,099; PhosphoNetworks 4,212; RLIMS-P 2,145; NCI-PID 1,671; dbPTM 1,271; Reactome 989; BEL 905; Li2012 322.
- **It is mostly re-distribution:** 12,447 rows exactly match the provided PSP human→human keys (84.9% of PSP keys are present). Of the 20,590 rows not in our PSP file: 6,311 are PSP-labelled (older PSP snapshot, rodent-to-human homology translation, isoform numbering — not independent evidence); 5,841 text-mining only (ProtMapper over REACH/Sparser/RLIMS-P/BEL); 5,171 KEA-mixed (KEA aggregates incl. HT/predicted sets; "unspecified NC-SA"); 3,317 PhosphoNetworks-only (motif-predicted sites); **only 5,792 rows (482 enzymes, 2,650 sites) are backed by a low-throughput curated DB other than PSP** (SIGNOR 3,102, NCI-PID 1,189, HPRD 1,050, Reactome 661, Phospho.ELM 595, dbPTM 303), and of these 3,838 have an enzyme that is a PSP human kinase (the rest are complexes/phosphatases/cyclins).
- Note: MIMP and ProtMapper are not evidence sources — MIMP re-maps PSP/PELM/HPRD; ProtMapper (Bachman et al. 2022/2023, INDRA) normalises sites from DBs and text mining.
- Verdict: **useful as a single API to retrieve SIGNOR/HPRD/PELM/NCI-PID/Reactome curated rows**; filter to `sources ∩ {SIGNOR, HPRD, phosphoELM, NCI-PID, Reactome, dbPTM}` non-empty, enzyme ∈ human protein kinases, re-map sites by sequence. Never take PhosphoNetworks/KEA/text-mining rows as labels. Müller-Dott et al. 2025 (below) similarly dropped ProtMapper/KEA-only rows.

### 2.7 EPSD 2.0, dbPTM 2025, qPhos, RegPhos 2.0 (aggregators)
- **EPSD 2.0** (Chen et al. 2025 GPB doi:10.1093/gpbjnl/qzaf057): 2.77 M p-sites, 223 species; "regulator" annotations are re-collected from PSP/UniProt/PELM/etc. Download page returns an empty body to curl (JS); 36.2 GB total. No new curated KS.
- **dbPTM 2025** (NAR 53:D377 doi:10.1093/nar/gkae1005; `https://biomics.lab.nycu.edu.tw/dbPTM/download.php` OK): experimental phosphorylation file (1.6 M sites) and kinase-group "benchmark" sets (positive + negative sites; negatives = non-annotated S/T on the same proteins). KS content is aggregated from PSP/UniProt/HPRD/etc. No license stated.
- **qPhos** (Yu et al. 2019 NAR doi:10.1093/nar/gky1052): 3.5 M quantification events, 199k sites, 484 conditions — quantitative phosphoproteomics, *not* KS; "potential upstream kinases" are predictions. Download requires registration.
- **RegPhos 2.0** (Huang et al. 2014 Database doi:10.1093/database/bau034): KS integrated from PSP/PELM/HPRD/UniProt; all known URLs time out (offline).
- Verdict: none for labels. dbPTM/EPSD are useful as a **phosphosite universe** (which S/T/Y are observed phosphorylated — candidate set and negative sampling frame) and site-level evidence counts (feature). GPS/dbPTM-style negatives are a cautionary example (see §7).

---

## 3. High-throughput in vitro and chemical-genetic KS data (features, not labels)

### 3.1 Sugiyama, Imamura & Ishihama 2019 Sci Rep 9:10503 (doi:10.1038/s41598-019-46385-4)
- Design: HeLa cytoplasmic protein lysate, phosphatase-treated (thermosensitive AP, >95% dephosphorylation for most sites, heat-inactivated), then 385 recombinant kinases (354 WT protein kinases, 21 cancer mutants, 10 lipid kinases), 0.5 µg kinase/100 µg protein, 1 mM ATP, 3 h; dimethyl-label vs. no-kinase control, label-swap; TiO2 enrichment, LTQ-Orbitrap; substrate = kinase/control ratio >2.
- Output: 175,574 kinase–site pairs on 20,669 sites (9,998 pS, 4,698 pT, 5,973 pY); 1,427 motifs for 289 kinases. Data: SI Tables S1–S5 (CC BY 4.0 article) and jPOST/ProteomeXchange PXD011366.
- Quality: peptide FDR 0–0.3% (identification only — no KS-level FDR); authors state hits are not truly "specific" because of stochastic DDA sampling; caveats: heat-denatured proteome, no localisation/scaffolds/complexes, high kinase dose, residual phosphosites, pY in S/T-kinase reactions possibly from contaminating kinases in commercial preparations, HeLa proteome only. Recovers only ~4% of PSP relationships; only 4,913 sites overlap the deep HeLa in-vivo phosphoproteome; near-zero overlap with Newman 2013.
- Re-use: **iKiP-DB** (Mari et al. 2022 J Proteome Res doi:10.1021/acs.jproteome.2c00198) re-processed/filtered it (26,786 interactions per Müller-Dott 2025; included in PTMsigDB v2). Invergo 2022 PLoS Comput Biol doi:10.1371/journal.pcbi.1010110 (IV-KAPhE) used it to build per-kinase Naive-Bayes specificity models (349 kinases), then used these **as features** in a random forest trained on **PSP labels** — exactly the role we propose. Niinae, Sugiyama & Ishihama 2023 (bioRxiv doi:10.1101/2023.03.20.533483; later FEBS Open Bio, PMC11447832) showed lysate-derived PWMs predict PSP substrates about as well as peptide-library PWMs for CK2/ERK1/PKA, but between-cell-line correlations of intensities are only 0.64–0.79.
- Benchmark evidence: Müller-Dott et al. 2025 Nat Commun 16:4771 (doi:10.1038/s41467-025-59779-y) found adding iKiP-DB targets to PSP **decreased** kinase-activity-inference performance (tumour AUROC range 0.53–0.54 vs ≈0.65–0.68 for curated); Jaccard with any other resource ≤0.03.
- Verdict: **feature (per-kinase in vitro hit / ratio / motif score), not labels, not test**. Rationale: (i) rule requires literature-derived, experimentally verified annotations — a single HT screen with ratio>2 and no KS-level error control is not curated verification; (ii) in vitro promiscuity in denatured lysate; (iii) using it as test would reward motif-only models. It is acceptable for an *orthogonal sanity analysis* ("do our novel high-confidence predictions for dark kinases enrich among Sugiyama hits?"), clearly labelled as such. Leakage note: PSP curates some HT papers, so check whether Sugiyama's PMID (31324866) appears among PSP references (needs PSP reference download). If it does, those KS pairs are both label and feature: drop them from test sets or report test metrics stratified by 'is a Sugiyama hit' vs not.

### 3.2 PhosphoNetworks — Newman et al. 2013 Mol Syst Biol 9:655 (doi:10.1038/msb.2013.12)
- Protein microarrays: 4,191 full-length human proteins × 289 kinases → rawKSR 24,046 (protein-level), refKSR 3,656 (255 kinases, 742 substrates; Bayesian score), comKSR = refKSR + ~719–741 literature KSRs. **Site-level network (4,417 kinase→site, 230 kinases, 2,591 in vivo MS sites) is assigned by motif (M3) prediction** — sites are *not* measured. Validation in HeLa: ~76% of testable refKSRs showed kinase-dependent change (but ~85% of first-round changes were protein-stability effects); false-negative rate ~95% acknowledged.
- Files (`https://phosphonetworks.org/download.html`): rawKSI.csv 293 kB, refKSI.csv 44 kB, comKSI.csv 74 kB, highResolutionNetwork.csv 154 kB, motifs. Free for academic/non-profit.
- In OmniPath, 3,283 human rows are PhosphoNetworks-only — these are the motif-predicted sites. **Do not use as labels** (circular: site assignment *is* a motif prediction). rawKSR/refKSR (protein-level, in vitro) → **feature** "kinase phosphorylates this protein on an array".

### 3.3 Kinase Library — Johnson et al. 2023 Nature 613:759 (doi:10.1038/s41586-022-05575-3); Yaron-Barir et al. 2024 Nature 629:1174 (doi:10.1038/s41586-024-07407-y)
- Positional-scanning peptide libraries: 303 S/T kinases (>84% of predicted-active human S/T kinases) and the Tyr kinome (~78 canonical + non-canonical Tyr kinases; I could not re-verify the exact count in this session). Output = PSSMs + per-site percentile scores across the phosphoproteome; no protein substrates. Includes kinase-enrichment ("Kinase Library enrichment") for phosphoproteomics.
- Access: `https://kinase-library.phosphosite.org/` (OK), python `kinase-library` 1.8.0 (PyPI; **CC BY-NC-SA 3.0**; GitHub TheKinaseLibrary/kinase-library); PSSMs in Nature SI.
- Verdict: **feature and baseline** (the PI must retrain/compare fairly; note the motifs were derived without PSP so they are leakage-free w.r.t. our labels — but Johnson et al. *validated* against PSP, so do not tune thresholds on our test). Not labels.

### 3.4 GPS "gold standard" sets
- GPS 5.0 (Wang et al. 2020 GPB 18:72 doi:10.1016/j.gpb.2020.01.001) and GPS 6.0 (Chen et al. 2023 NAR 51:W243 doi:10.1093/nar/gkad383; 30,043 site-specific KS relations in 7,041 proteins, all species; 577 kinase/group predictors). The gold set is compiled from PSP, Phospho.ELM, UniProt, HPRD, SIGNOR etc. Müller-Dott 2025 extracted the GPS 5.0 set from SI: 352 kinases, only **544 interactions unique** vs other libraries. GPS download page offers only software/models (academic-only; 5 GB model zip).
- Verdict: not needed (redundant, unclear provenance/licensing, mixes species and family-level labels). Possibly useful for tracing a few extra PMIDs.

### 3.5 KinaseNET / PhosphoNET (Kinexus)
- PhosphoNET: ~650k human sites (≈14% experimentally observed), kinase predictions from proprietary matrices for ~500 kinases (Safaei et al. 2011 Proteome Sci). KinaseNET is Kinexus' kinase knowledge-base portal (kinasenet.ca blocked here; `www.kinexus.ca/kinasenet` 404). Proprietary, predictions → **exclude**.

### 3.6 Chemical-genetic (analog-sensitive kinase, ATPγS-analog) direct-substrate studies
- Principle (Shokat lab): gatekeeper-mutant kinase accepts N6-substituted ATPγS; thiophosphorylated peptides are covalently captured and sequenced → **direct** substrates with site localisation, in lysate or permeabilised cells. Examples: CDK1–cyclin B (Blethrow et al. 2008 PNAS doi:10.1073/pnas.0708966105), ERK2 (Carlson et al. 2011 Sci Signal doi:10.1126/scisignal.2002010), AMPK (Banko et al. 2011 Mol Cell doi:10.1016/j.molcel.2011.11.005), yeast Cdk1 (Holt et al. 2009 Science doi:10.1126/science.1172867; yeast — irrelevant except for orthology), CDK7/CDK9, PLK1 and others. KALIP (Xue et al. 2012 PNAS doi:10.1073/pnas.1119418109) combines in vitro lysate reactions with in vivo perturbation phosphoproteomics; OBIKA (on-bead in vitro kinase assays; STAR Protocols 2026) similarly intersects in vitro with in vivo hits. Kinase-catalysed biotinylation with ATP-biotin (K-BMAPS/K-BILDS, Pflum lab) is lower-throughput and protein-level.
- Quality: higher directness than PSP "in vivo"; still lysate context and analog-sensitive kinases can have altered specificity. Many of these papers are curated in PSP (as in vitro HTP). Per-study coverage: tens–hundreds of sites for ~10 human kinases.
- Verdict: if curated into PSP/SIGNOR → part of labels. A **hand-assembled "direct-evidence" external test set** from these papers (site-level, kinase-specific, with in-vivo corroboration à la KALIP) would be the most informative test of directness, but it is curation work; recommended as a later, small, high-value test (tagged by PMID, filtered against training sites).

### 3.7 Proximity labelling of kinase substrates
- TurboID–kinase fusions + phosphoproteomics: BIN2 (Arabidopsis; Kim et al. 2023 Plant Cell doi:10.1093/plcell/koad013 — plant), BioID kinase interactomes + kinase perturbation + motif for CK2/PKA in human cells (Niinae et al. 2021, PMC8325102: 24 and 35 putative substrates). SubMAPP (photo-activated TurboID phospho-mapping) is a method, no resource. I found **no 2023–2026 kinome-scale human proximity-labelling KS resource**. Use: PPI/proximity → feature (covered by the "other modalities" review).

### 3.8 Literature-mined resources: RLIMS-P, eFIP, REACH/Sparser (INDRA), BEL
- RLIMS-P (Torii et al. 2014 Database doi:10.1093/database/bau081): F-scores 0.76 kinase, 0.88 substrate, 0.92 site on BioCreative IV; eFIP full text: precision 92.4%, recall 76.5% (Tudor et al. 2015 Database doi:10.1093/database/bav020). Precision of the *kinase-substrate-site triple* is materially lower than the per-entity numbers and text mining cannot distinguish direct from indirect regulation or in vitro from in vivo.
- Verdict: **no labels**. Could be a feature ("a paper mentions kinase X with site Y") but it is literature-derived and highly correlated with the labels → strong leakage risk (it effectively encodes the answer for test sites). If used at all, only as a curation queue for humans, not as a model feature.

### 3.9 "KSD"/"KinSub"
- No current resource by these exact names found. Closest: **iKinSub** (NHLBI ESBL; in vitro kinase reactions on rat IMCD lysate for a few kinases, e.g. PKA, SGK1, CAMK2D; Bradford et al. 2014 Am J Physiol Cell Physiol) and the kinase-substrate logo database / KinasePredictor (Cell Commun Signal 2024 doi:10.1186/s12964-023-01436-2) — HT in vitro, rat, few kinases → marginal feature at most. "Yeast KID" (Sharifpoor et al. 2011, PMC3218865) is yeast. If the PI meant a specific resource, a citation is needed.

### 3.10 PhosPhAt (Arabidopsis; Heazlewood/Schulze) — plant; irrelevant. Host blocked here.

### 3.11 Other 2023–2026 developments worth knowing (not label sources)
- **SIGNOR 4.0 / PhosphoSIGNOR** (above) is the main new curated KS resource.
- Dabas et al. 2025 Science 390 (in vitro screen of 427 kinases on RNA Pol II CTD, 117 hits) — illustrates the promiscuity debate: an arXiv comment (Li & Zhou 2025, arXiv:2511.19444) argues short-peptide in vitro hits overstate in-cell targeting. Example of why HT in vitro ≠ labels.
- Large in vivo kinase-activity atlases (e.g. 33-cell-line phosphoproteome atlas, Nat Struct Mol Biol 2026 doi:10.1038/s41594-026-01877-6; CPTAC co-regulation maps such as CoPheeMap, bioRxiv 2024 doi:10.1101/2024.03.19.585786) and inhibitor-response data (Hijazi et al. 2020; decryptM) are **feature/benchmark** sources covered in the phosphoproteomics review.
- Benchmarks of KS libraries: Müller-Dott et al. 2025 (benchmarKIN) — PSP best; only 21.7% of interactions shared by ≥2 resources; median 8.5–18 targets per kinase for curated libraries vs 64–69 for NetworKIN/iKiP-DB.

---

## 4. Pooling human, mouse and rat

- Biology: kinase specificity is highly conserved (Yaron-Barir 2024: Tyr kinase specificity "fundamentally unchanged from worms to humans"); mouse/rat substrate sites with conserved ±7 windows are very likely phosphorylated by the human ortholog. PSP itself corroborates ~9% of kinases/8% of substrates across species (2015).
- Gain (computed): mouse/rat-substrate rows among H/M/R kinases: 7,234 → **4,669 new (kinase gene, site-group) pairs** not present with a human substrate. But **3,628 of these site groups have no human row in the KS file**, so the human residue cannot be obtained from `Kinase_Substrate_Dataset` alone. Options: (a) PSP `Phosphorylation_site_dataset.gz` (contains SITE_GRP_ID for human sites — **request from user**; requires PSP login); (b) align the rodent substrate to its human ortholog (UniProt/Ensembl Compara 1:1 orthologs; global alignment) and transfer the site only if the residue is S/T/Y-identical and the ±7 window identity ≥ ~60–70% (a pre-registered threshold), flagging the label as `orthology_transferred`.
- Kinase side: map rodent kinases to human by 1:1 ortholog (upper-cased gene symbol works for most; verify with accession-based orthology; watch paralog expansions, e.g. rodent-specific kinases).
- Leakage: orthologous sites of the same protein must fall in the same split as the human protein (they share SITE_GRP_ID → group by human protein after mapping). Report results with and without orthology-transferred labels; test sets should be human-only native annotations (transfers in training only), so that pooling is evaluated rather than assumed.
- 797 human-substrate rows used non-human kinases (mouse 328, rat 292, cow 115, chicken 12, viral 14) — keep mammalian ones with ortholog mapping; drop viral.

## 5. In vivo vs in vitro evidence

- In vitro (purified kinase) establishes the *capacity* for direct phosphorylation but is prone to promiscuity: loss of docking/scaffold/localisation constraints, high enzyme concentration, long incubations, peptide substrates with KM orders of magnitude above protein substrates (Johnson & Hunter 2005 Nat Methods doi:10.1038/nmeth731; Ubersax & Ferrell 2007 Nat Rev Mol Cell Biol doi:10.1038/nrm2203). Conversely, in vitro assays miss substrates requiring priming, adaptors or cofactors (Newman 2013 FN ~95%).
- In vivo (cellular perturbation) establishes physiological relevance but not directness: downstream kinases, phosphatase effects, protein-stability effects (Newman 2013: ~85% of kinase-dependent changes were abundance changes).
- Evidence from the provided data: in-vivo-only labels are 2.8× enriched among residue-type-implausible annotations (§1). Petsalaki/Maier 2025 note the same. Invergo 2022 trained on PSP (both classes) and used in vitro data as features.
- Recommendation: treat **"both" as the gold-standard core** (5,171 human rows), "in vitro only" as direct-but-possibly-promiscuous, "in vivo only" as possibly indirect. Pre-register: main metric on all PSP human labels; secondary metrics on "both" and on "in vitro"; an ablation training on (both + vitro) vs all. Do not choose the subset by test performance.

## 6. Study bias and the dark kinome

- PSP (provided): top-20 kinases = ~50% of human annotations; 146/421 kinases have <5 annotations; 421 kinases have ≥1 human→human annotation, so roughly 100–130 of the ~520–540 human protein kinases have none (exact number needs a curated kinome list). Literature: ~95% of known phosphosites lack a kinase; ~80% of annotated sites belong to ~20% of kinases; ~30% of kinases (~150) have no known substrate; PSP KS pairs have median 685 citations per kinase vs 194 for SELPHI predictions (Maier et al. 2025); Needham et al. 2019 Sci Signal doi:10.1126/scisignal.aau8645; Invergo & Beltrao 2018 Essays Biochem doi:10.1042/EBC20180019 (citation bias: well-studied kinases are tested on well-studied substrates). IDG "dark kinases" (~160) — Berginski et al. 2021 NAR doi:10.1093/nar/gkaa853 (Dark Kinase Knowledgebase).
- Consequences: (i) pooled AUC is dominated by CDK1/PKA/SRC/CK2/PKC/ERK → report **macro-averaged per-kinase** metrics with kinase-count strata; (ii) "hub" substrates (tau 241 rows, p53, β-catenin) leak through protein-level features unless split by protein; (iii) annotated sites are biased to well-studied proteins with abundant MS evidence → negatives sampled from all phosphosites are "easier" (protein-level features separate them) — sample negatives matched on protein study-level (e.g. from the same substrate proteins and from proteins with similar citation counts); (iv) sites annotated to only one kinase are not negatives for others (missing labels).
- Dark-kinase evaluation: leave-kinase-out (and leave-family-out) is the only honest proxy; Kinase Library motifs and Sugiyama in vitro data are the main label-independent information for dark kinases — another reason to keep them as features, not labels.

## 7. Negatives

- No database of true negative KS pairs exists. Practices in the literature:
  - Random kinase–site pairs from the phosphoproteome (Invergo 2022, residue-type matched; Maier et al. 2025 SELPHI2.0, 10–50× positives; results insensitive to ratio).
  - Non-annotated S/T/Y on the same substrate proteins (dbPTM/GPS/KinasePhos benchmarks; problematic: many are unobserved/unphosphorylatable, which inflates performance; and true unannotated sites).
  - Sites annotated to *other* kinases (cross-kinase negatives; common in Phosformer-style and kinase-classification setups) — biased toward kinases with opposite motifs, and contaminated by shared substrates/paralogs.
  - In vitro tested-but-not-phosphorylated sites (Sugiyama data provide quasi-negatives for in vitro specificity; Invergo used them for the in vitro model only).
  - Positive-unlabelled (PU) learning (e.g. Yang et al. 2016 Bioinformatics doi:10.1093/bioinformatics/btv550).
- Recommendation for this project: frame as PU; build negatives from **observed phosphosites** (from the phosphoproteomics compendium / PSP site dataset / dbPTM) not annotated to the kinase, residue-type matched (S/T vs Y; dual-specificity explicitly), with at least two pre-registered negative sets for evaluation: (A) random observed sites (easy), (B) sites on the *same proteins* as positives plus sites annotated to other kinases of the same group (hard). Metrics: per-kinase AUPRC with fixed positive:negative ratio and ranking metrics (top-k recovery among all observed sites of a held-out protein). Never treat unannotated pairs as confirmed negatives in error analysis.

## 8. Recommendations

**Training labels (in priority order).**
1. PSP human→human (provided), sites re-mapped to UniProt canonical by 15-mer; non-STY dropped; evidence class retained.
2. PSP mouse/rat (and human kinase on rodent substrate) transferred to human via SITE_GRP_ID / ortholog alignment (flagged; training only, ablated).
3. SIGNOR 4.0 (`DIRECT=YES`, residue given, protein kinase enzymes; family/complex entries expanded only if unambiguous) — via PhosphoSIGNOR TSV (user to fetch; blocked here) or OmniPath enzsub.
4. UniProt Swiss-Prot ECO:0000269 "by <GENE>" sites (gene-level names only; autocatalysis separate).
5. Optional: HPRD/Phospho.ELM/NCI-PID/Reactome rows from OmniPath only if corroborated by ≥1 other curated source and enzyme is a single protein kinase.

**Test labels.** Protein-grouped held-out PSP human sites (native, not transferred) as the main test; leave-kinase-out/leave-family-out; secondary: (i) "both"-evidence subset; (ii) SIGNOR/UniProt-only pairs (independent curators; subject to the same protein grouping — exclude proteins in training); (iii) later, a temporal hold-out from SIGNOR 4.0 entries added after the PSP release; (iv) optional curated chemical-genetics "direct" set.

**Features only.** Kinase Library PSSM scores; Sugiyama/iKiP-DB in vitro hits & motifs; PhosphoNetworks protein-level array hits; (later modalities: PPI, co-regulation, inhibitor perturbations). **Excluded:** text-mining KS (leakage), PhosphoNetworks site-level network (motif-predicted), KEA/GPS/NetworKIN/PhosphoNET predictions, aggregators' redistributed PSP rows (double counting).

**Requests to user.** (1) PSP `Phosphorylation_site_dataset.gz` (+ `Regulatory_sites.gz`) from the same 2025-05 release (needed for SITE_GRP_ID→human site mapping and the phosphosite universe); (2) PhosphoSIGNOR/SIGNOR download (proxy-blocked here) or permission to rely on OmniPath's SIGNOR copy; (3) confirmation that NC-SA licensing (PSP, Kinase Library, iPTMnet) is acceptable for model/label release; (4) whether to accept the Phospho.ELM academic click-through license (low value — recommend skip).

## 9. Reproducibility notes (transient query inputs, not project data)

Computed with `scripts/lit/ks_source_overlap_stats.py` on 2026-10-26. Inputs were fetched to a scratch directory and are not stored under `external/` (no project data derived from them yet); a data agent should re-acquire them into `external/<source>/` and record them in DATA.md before any modelling use.
- OmniPath enzsub (URL in §2.6), 4,868,762 bytes, SHA-256 `db90dc41315d73fef5b8276b9cdeb91a6da51745e98750e25203496774c9e8e4`.
- UniProt human reviewed MOD_RES "by" (URL in §2.2), 1,627,793 bytes, SHA-256 `16f15905cb48d1cdcbe2e49e96e57f68fbb1233e4a1d4fa78ca163886396407b`.
- iPTMnet `ptm.txt`: only byte-range samples inspected (source codes).
Caveats: kinase identity matched by gene symbol; isoform suffixes stripped for OmniPath/UniProt comparisons, so a minority of "new" pairs may be numbering artefacts; counts of resource contents quoted from papers were not independently re-derived.
