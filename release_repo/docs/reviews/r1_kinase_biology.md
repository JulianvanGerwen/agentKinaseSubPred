> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Review r1 — kinase signalling biology & phosphoproteomics

Reviewer role: kinase signalling biologist / phosphoproteomics. Reviewed: DESIGN.md v0, LITERATURE.md, docs/lit/03–04, docs/data_profile.md, docs/data_acquisition.md.

**Checks behind the numbers below.** I ran these read-only with `python` and ≤4 cores. The scripts are in my scratchpad and are not committed: `prep.py`, `a1.py`–`a4.py` in `<agent-scratchpad>/`. The PI should port the useful ones into `scripts/`.
- Label set: HMR PSP pairs (human, mouse and rat kinase × substrate), restricted to S/T/Y sites mapped to a human Swiss-Prot canonical sequence with the window verified. This gives 18,132 (kinase, site) pairs and 437 kinases.
- Paralogs: all-vs-all mmseqs2 on full-length kinase sequences, keeping alignments ≥250 aa long.

Overall verdict: the design is careful on leakage and statistics. Its biological weak spots are:
- the definition of "hard" negatives between close paralogs;
- the downstream benchmark, which as written is partly circular;
- inhibitor polypharmacology;
- docking and priming features as specified, which carry almost no information;
- residue-class and pseudokinase rules.

---

## Proposed changes (priority order)

### KB-1 (critical) — Mask paralogs in hard negatives (§5 "Evaluation negatives" E-hard, §5 "Training negatives", §6 primary endpoint, §6 site-centric)

**Problem.** E-hard defines a negative for *k* as any site annotated to another kinase. For close paralogs, those "negatives" are mostly unlabelled positives. PSP often assigns a generic kinase name to one gene:
- "PKA" → PRKACA;
- "CK2" → CSNK2A1;
- "AMPK" → PRKAA1;
- "PKC" → PRKCA.

**Measured.** Share of each kinase's E-hard negatives that are labelled to a paralog (≥50% identity, aln ≥250):

| kinase(s) | share of E-hard negatives labelled to a paralog |
|---|---|
| BLK, FGR, YES1, HCK, LCK, LYN (SRC family) | 0.51–0.54 |
| CDK3, CDK5, CDK2 | 0.13, 0.09, 0.09 |
| PRKCG, PRKCZ, PRKCI | ≈0.08 |
| MAPK7 | 0.08 |
| CSNK2A2 (vs CSNK2A1, ≥80% identity) | 0.07 |
| median over 123 kinases that have a paralog | 0.006 |

So the problem is concentrated in exactly the families where users most need specificity. Without masking, the macro AUROC for SFK members would reward reproducing PSP's naming habits, not biology.

**Change.**
1. Define P(k), the paralog set of *k*: kinases with ≥60% full-length identity over ≥250 aa (mmseqs2), plus any kinase in the same Manning subfamily. Write it to `results/labels/kinase_paralogs.tsv`.
2. **E-hard-pm**: E-hard with every site labelled to a kinase in P(k) removed from *k*'s negatives. Make the **primary endpoint macro AUROC under E-hard-pm**. This is acceptable now because no results exist yet. Keep unmasked E-hard as a secondary.
3. Training hard negatives also exclude P(k)-labelled sites, or carry a soft label of 0.5.
4. Site-centric Hit@k: also report **family-level Hit@k**, where a prediction counts as correct if a kinase in {k} ∪ P(k) is ranked, beside the gene-level value.

### KB-2 (critical) — Make the downstream kinase-activity benchmark non-circular (§6 "Downstream (secondary)")

**Problem 1: CPTAC.** As I understand Müller-Dott 2025, the benchmarKIN tumour benchmark (`benchmarking_gold_standard_top5per.Rds`) defines active/inactive kinases per tumour from the kinase's own protein abundance and activating-site phosphorylation in CPTAC. Verify this in its Methods. Two of our features are then circular with that ground truth:
- F6c correlates *s* with *k*'s own abundance and activation-loop phosphosites in CPTAC;
- F5 uses CPTAC co-abundance.

Seven of our ten CPTAC cohorts are the benchmark cohorts (BRCA, CCRCC, GBM, HNSCC, LSCC, LUAD, UCEC).

**Problem 2: perturbation studies.** 21 of the 24 Hernandez-Armenta PMIDs, and all of Hijazi, are inside the provided compendium. Excluding them by PMID is necessary but not sufficient: 17 Hijazi inhibitors also appear in other compendium studies (`kinhib_condition_matches.tsv`). The F6b "responds to inhibitors of K" signal therefore transfers by drug identity.

**Problem 3: indirect targets.** About 50 of the 173 Hernandez target rows are upstream kinases whose direct substrates are few: MAP2K1 14, MAP2K2 14, BRAF 14, RAF1 8. An indirect, pathway-level library wins on these, which pulls against our goal of predicting *direct* substrates.

**Change.**
1. Build a benchmark model variant "BK-clean":
   - for the tumour benchmark, it drops **all** CPTAC-derived features (F5-CPTAC and F6c);
   - for the perturbation benchmarks, it drops the benchmark studies **and** every compendium condition that uses a drug whose target set includes the benchmark target kinase (leave-drug/target-out).
2. Compare these libraries on the same kinase set, with N per kinase fixed a priori (N ∈ {50, 100, 250}):
   - (a) PSP-only;
   - (b) Cantley top-N;
   - (c) ours, predicted-only, with PSP pairs removed (this measures added value);
   - (d) PSP ∪ ours;
   - (e) a negative control: our library with sites shuffled within kinase in a degree-preserving way.
3. Report direct-substrate-poor upstream targets (MAP2K/RAF/MAP3K) separately, and state in advance that the library will not be tuned to them.
4. Libraries come from the DEV-trained evaluated model, not the deployment model.

### KB-3 (major) — Residue class and the kinase list: evidence-based rules, pseudokinases, ID hazards (§2 Kinases, §5 "Kinase residue class")

**Measured residue patterns.**
- **DYRK / CLK / HIPK / GSK3.** In PSP their Y labels are almost all cis activation-loop autophosphorylation:
  - DYRK1A/1B/3/4: 6 autophos Y, 0 trans Y;
  - CLK3/4: 1 trans Y each;
  - GSK3B: 1 trans Y.
  - Classing them as S/T/Y would add about 40k Y candidates per kinase with no trans biology behind them.
- **MAP2K1–7.** 25 trans Y labels (TEY/TGY/TPY loops); MAP2K1 alone has 10.
- **Cantley's 15 non-canonical Y kinases.** The Y atlas profiles 15 non-canonical tyrosine kinases. Those that are not already in DESIGN's list are LIMK1/2, BMPR2, PINK1, NEK10, PDHK1/3/4 and TNNI3K. These rest on in vitro peptide activity, and PSP has ≤1 trans Y label for each.
- **TK kinases on S/T.** 83 pairs, 66% in vivo-only. DESIGN's decision to drop them is correct.

**Change.**
1. **Residue classes:**
   - **Y only:** WEE1, WEE2.
   - **S/T/Y:** MAP2K1–7, PKMYT1 (CDK1 T14/Y15), TTK, TESK1/2, NEK10 (p53 Y327), PKDCC/VLK (a secreted Tyr kinase; missing from DESIGN).
   - **S/T, with Y labels dropped as cis autophosphorylation and the count reported:** DYRK1A/1B/2/3/4, CLK1–4, HIPK1–4, GSK3A/B.
   - **S/T by default, with Y as a sensitivity analysis:** the Cantley non-canonical set (LIMK1/2, BMPR2, PINK1, PDHK1/3/4, TNNI3K).
2. **Pseudokinases.** Add an explicit list (Manning 2002 / Kwon 2019) to `kinase_list.tsv` with the flag "no predictions, labels dropped":
   - TRIB1–3, STRADA/B, SCYL1–3, PEAK1/PRAG1, ULK4, NRBP1/2, TBCK, PXK, RYK, PTK7, EPHA10, MLKL, STK40, KSR1;
   - ILK: 9 PSP labels (GSK3B S9, AKT S473 …), legacy and now doubted;
   - POMK and FAM20B: sugar kinases.

   Keep CASK and FAM20C. Also add PRPF4B: like FAM20C it is missing from pkinfam, and it fell into the "NA" group in my check.
3. **Regulatory-subunit labels.** Map holoenzyme rows to the catalytic subunit instead of dropping them, and flag them: PRKAB1/PRKAG2 (21 pairs) → PRKAA1, CSNK2B (14) → CSNK2A1.
4. **Joining on IDs.** Join Cantley matrices and the kinhib table by **UniProt accession, never by name**. In `kinome_information.tsv`, 174 of 404 matrix names differ from the gene symbol. The worst case: matrix "PDK1" is PDPK1, while gene *PDK1* is PDHK1. PSP has 71 "PDK1"-named rows that are PDPK1, and 14 rows that are PDHK1.

### KB-4 (major) — Docking-motif features as written carry almost no information (§9 "Docking-motif kinases", F-family in §7)

**Measured on human Swiss-Prot canonical proteins.**

| motif (raw regex) | proteome | all KS substrate proteins | kinase-specific substrates |
|---|---|---|---|
| ELM D-motif (DOC_MAPK_gen_1) | 76% | 84% | MAPK1: 86% |
| RxL | 77% | — | CDK2: 88% |
| S-[ST]-P | 32% | 45% | PLK1: 53% |

The small excess is mostly a protein-length and "studied protein" effect.

**Change.** Replace "ELM class regexes, distance to site" with features that require context:
1. A motif instance counts only if it lies in an IDR (AF pLDDT < 70 or predicted disorder) and is conserved (ortholog or PLM pseudo-conservation). Use a D-motif PSSM (e.g. Zeke et al. 2015) rather than a regex.
2. Encode the signed distance to the site. Orientation matters: MAPK D-motifs usually lie N-terminal to the site, and the DEF motif F-x-F-P lies C-terminal.
3. **Phospho-dependent docking.** Require the docking phosphosite to be an *observed* phosphosite:
   - CDK–Cks1: pT-P 10–40 aa N-terminal of the site;
   - PLK1 polo-box: S-[pS/pT]-P;
   - TK SH2 processivity: the best score of the kinase's own SH2 domain (Yaron-Barir MOESM7 SH2 matrices, already downloaded) over observed pY sites on the substrate protein, for the SRC, ABL, SYK/ZAP70, TEC and CSK families;
   - PDPK1 PIF: the substrate carries an AGC hydrophobic motif F-x-x-F-[ST]-[FY].
4. Pre-register per-family ablations: does each docking feature improve AUROC for its own family under E-hard-pm?

### KB-5 (major) — Inhibitor-target mapping: polypharmacology, cascades, paradoxical activation (§7 F6b)

**Measured in the Hijazi matrix.**
- 177 of the 307 inhibited kinases have an inhibitor profile identical to at least one other kinase, so they cannot be separated by this matrix.
- Only 54 kinases are hit by an inhibitor with ≤5 targets.
- The matrix is binary, with no potency.

**Known biology the mean/min contrast misreads.**
- **Cascades:** MEK inhibition lowers ERK substrates; mTORC1 inhibition lowers S6K substrates.
- **Paradoxical or feedback increases:**
  - RAF inhibitors activate ERK in RAS-mutant cells;
  - ATP-competitive AKT inhibitors hyperphosphorylate AKT T308/S473;
  - rapalogs and AZD8055 relieve feedback onto AKT and ERK.
- **Time:** 24 h contrasts mostly reflect changes in protein levels.

**Change.**
1. Compute a per-site deconvolution across inhibitor contrasts: non-negative regression of −log2FC(s,c) on the kinase-inhibition matrix I(k,c).
   - Where Kd_app is available, I(k,c) is occupancy from dose and Kd_app: decryptM drugs via Kinobeads/Reinecke, or LINCS/Davis via KB-12.
   - Otherwise I(k,c) is the binary value.
2. Collapse kinases whose I-vectors are collinear (|r| > 0.95) into a "co-inhibited group". Each member gets the group's value plus an identifiability flag.
3. Use only contrasts ≤ 2 h, or down-weight longer ones by a fixed rule.
4. Keep "decrease" and "increase" as separate features.
5. Add decryptM potency coherence (|log EC50_site − log Kd_app|) once the curve zips are downloaded.
6. Keep a curated drug→target table (`results/labels/drug_targets.tsv`, with source, version and Kd) covering decryptM's 31 drugs and the drug contrasts in qPTM, not only Hijazi.

### KB-6 (major) — Add a site-within-protein discrimination metric (§6, key secondary)

**Rationale.** Users interpreting phosphoproteomics need to know *which site* on a protein a kinase targets. Several feature families cannot answer that, because they are constant within a protein:
- PPI (F3);
- localisation (F4);
- co-expression (F5);
- docking-motif presence (KB-4), except its distance term.

E-all gains can come entirely from those protein-level features.

**Change.**
1. For each kinase *k* and each test protein with ≥1 positive for *k*, score the positives against all other residue-matched U_obs sites on that protein. Report the macro over kinases of the pooled within-protein AUROC.
2. Pre-register a flag: any model or ablation whose E-all gain over B1/B2 is not matched by a gain on this metric is labelled "protein-level gain".
3. Report the canonical-pair panel (KB-7a) in the same way.

### KB-7 (major) — Biologically meaningful positive and negative controls (§11)

All of these are report-only and pre-registered now. None of them is used to tune the model.

**a. Canonical-pair panel.** Freeze `tests/canonical_ks_pairs.tsv` now (~40 pairs). Examples:
- MAP2K1 → MAPK1 T185/Y187; PDPK1 → AKT1 T308; MTOR → RPS6KB1 T389 and EIF4EBP1 T37/T46;
- AKT1 → GSK3B S9; GSK3B → CTNNB1 S33/S37/T41; CSNK1A1 → CTNNB1 S45;
- PRKACA → CREB1 S133; PRKAA1 → ACACA S79; STK11 → PRKAA1 T172;
- ATM → TP53 S15; ATR → CHEK1 S345; CHEK1 → CDC25A S123;
- WEE1 → CDK1 Y15; PKMYT1 → CDK1 T14; CDK7 → CDK1 T161;
- AURKB → H3C1 S10; IKBKB → NFKBIA S32/S36; TBK1 → IRF3 S386; RIPK3 → MLKL T357/S358;
- PINK1 → PRKN S65; LRRK2 → RAB10 T73; EEF2K → EEF2 T56;
- PDK1–4 → PDHA1 S293/S300/S232; BCKDK → BCKDHA S337; ALPK1 → TIFA T9;
- ZAP70 → LAT Y220; LCK → ZAP70 Y493.

Report, from out-of-fold predictions (CV on DEV; TEST at final), each pair's within-kinase U_obs percentile and its within-protein rank.

**b. Motif-fraction checks.** Take the top-100 *unlabelled* predictions per kinase from a model variant **without F1/Cantley**, so the check is not circular. Expected values, measured on PSP labels:

| kinase(s) | feature | label fraction |
|---|---|---|
| CDK1/2 | P at +1 | 0.90 |
| MAPK1/3 | P at +1 | 0.89 |
| PRKACA / AKT1 | R/K at −3 or R at −2 | 0.79 / 0.86 |
| ATM | Q at +1 | 0.88 |
| CSNK2A1 | acidic at +1..+3 | 0.86 |
| FAM20C | E at +2 (S-x-E) | 0.93 (n = 14) |
| GSK3B | S/T at +4 | 0.56 (0.41 known pS) |

Background across all S/T labels is P+1 0.30, basic −3/−2 0.36, Q+1 0.06.

Flag any kinase whose prediction fraction falls more than 0.15 below its label fraction. Note that MTOR is not an SQ kinase: P+1 0.50, Q+1 low. Its predictions must not be dominated by S/T-Q.

**c. PWM concordance.** For atlas kinases, compute Pearson r between the log-odds PWM of top-N predictions and the Cantley PSSM, again from the no-F1 variant. The reference ceiling is r(PWM of training labels, Cantley).

**d. Biological negative controls.**
- Predictions on lumenal or extracellular sites (KB-8) for cytosolic kinases should sit near the bottom.
- For FAM20C, ≥80% of its top-100 predictions should be on signal-peptide or TM proteins.
- Mitochondrial-matrix sites should not be predicted for non-mitochondrial kinases.
- No predictions are emitted for pseudokinases (KB-3).

**e. Cascade confusion.** For upstream kinases (MAP2K1/2, MAP3Ks, BRAF/RAF1, PDPK1, STK11, CAMKK2), report the fraction of their top-50 predictions that are PSP substrates of their downstream kinases (ERK, AGC, AMPK-related). A high fraction means co-regulation and inhibitor features have learned the pathway, not direct catalysis.

**f. Popularity of predictions.** Report the fraction of top-N predictions that fall on the 5% most-annotated proteins, compared with the background fraction.

### KB-8 (major) — Add site-level topology and compartment compatibility (§7 F4; acquisition)

**Problem.** Compartment vectors at protein level miss the most decisive constraint: *which side of the membrane the residue is on*.
- Cytosolic kinases cannot reach lumenal or extracellular segments.
- FAM20C, and PKDCC for Y, act only in the secretory lumen.
- PDK1–4 and BCKDK act only in the mitochondrial matrix.

**Change.**
1. Download UniProt per-residue features: Topological domain, Transmembrane, Intramembrane, Signal peptide, Transit peptide, Subcellular location. This is one REST stream, minutes of work.
2. Derive a per-site class: cytosolic/nuclear-facing, lumenal/extracellular, matrix, or unknown.
3. Add a kinase × site compatibility feature: the kinase's compartment class (cytosolic/nuclear, secretory lumen, matrix) crossed with the site class.
4. Use it for reporting unorthodox kinases and for KB-7d.
5. For a FAM20C external test, add the FAM20C-dependent secretome phosphosites of Tagliabracci et al. 2015 (Cell 161:1619). These are KO/KD phosphoproteomics, so they are external and evidence-flagged, not training labels.

### KB-9 (major) — CPTAC and compendium co-phosphorylation confounders (§7 F5, F6a, F6c)

**Problem.** Beyond protein abundance, tumour co-phosphorylation is driven by several confounders:
- tumour purity and immune/stromal content: immune-kinase substrates co-vary with infiltration;
- TMT plex batch;
- cold-ischemia time: stress kinase sites rise with delayed fixation (Mertins 2014);
- shared or multi-site peptides: BCM site IDs carry the peptide, and sites on the same peptide share quantification.

**Change.**
1. Within each cohort, take residuals after regressing out:
   - the site's own protein abundance;
   - purity, and ESTIMATE immune and stromal scores;
   - TMT plex;
   - ischemia time where available;
   - the top 10 PCs, a number fixed now.
2. Use Spearman correlation on tumour samples only.
3. Collapse sites that come from the same peptide before computing correlations.
4. In F6a, for both compendium and CPTAC, **exclude training substrates on the same protein, and those within ±30 aa**, from the aggregation. Otherwise deployment predictions cluster around known sites (a "halo"): same-peptide or same-protein co-variation is not evidence of the same kinase.

### KB-10 (minor) — Priming features should be soft and tied to specific priming partners (§7 F1, §9)

**Measured.**
- GSK3B labels: S/T at +4 in 56%, and a known pS at +4 in 41%.
- CK1α/δ/ε labels: known pS/pT at −3 in only 24%, and D/E at −3 in 20%. The rest follow non-canonical SLS or acidic-cluster patterns, e.g. CTNNB1 S45.
- Ascertainment of the "known neighbour" status is small: 1.41 vs 1.22 known sites within ±7 for KS-labelled vs other PSP sites.
- Compendium coverage of GSK3B sites is not depleted (52% vs 49% for all labelled sites). CK1 (32–39%) and SRC (37%) are less well covered.

**Change.**
1. Keep both the primed and unprimed Cantley scores, as in the design.
2. Set the "observed phosphosite at offset" flag from MS evidence only (MS_LIT + MS_CST ≥ 2), not from LT-only.
3. Score "predicted priming" only for the canonical partner pairs:
   - GSK3 +4 ← CDK, MAPK, DYRK1A, CK2, PKA;
   - CK1 −3 ← CK2, PKA, CK1 itself;
   - CK2 hierarchical ← CK2;
   - PLK1 polo-box ← CDK1 (KB-4).
4. Report the priming stratum separately for GSK3A/B, CK1α/δ/ε and PLK1–4.

### KB-11 (minor) — Be realistic about what can be evaluated for unorthodox kinases (§6 stratification, §9)

**Measured** (human-native sites, with expected TEST positives at 20% in brackets):

| kinase(s) | human-native sites |
|---|---|
| FAM20C | 12 (≈2) |
| PDK1–4 | 4–7 each |
| BCKDK | 2 |
| ALPK1–3 | ≤5 |
| RIOK1–3 | ≤1 |
| SMG1 | 6 |
| EEF2K | 13 |
| TRPM7 | 57, but alpha-kinase labels overall are 73% autophosphorylation |

Only about 158 kinases have ≥15 human-native sites (≈ ≥3 test positives). The "Atypical" stratum is therefore really the PIKKs ATM, ATR, PRKDC and MTOR.

**Change.**
1. Label the Atypical stratum honestly: report "PIKK-SQ" (ATM, ATR, PRKDC, SMG1) separately from MTOR, which is proline/hydrophobic-directed.
2. Exclude autophosphorylation from alpha-kinase test positives.
3. For kinases with fewer than 3 test positives, report pre-registered case studies from the KB-7a panel instead of metrics.
4. For atlas-absent atypicals, use the nearest-atlas fallback only within the same fold class (ALPK and TRPM → EEF2K; PDK2/3 → PDHK1/4). Set it to NA for RIOKs.

### KB-12 (minor) — Data worth acquiring now (value per effort)

| # | Item | Why | Effort / access |
|---|---|---|---|
| 1 | UniProt per-residue topology, signal and transit peptides | KB-8; decisive for FAM20C, PDHK, BCKDK, PKDCC; a negative control | minutes (REST) |
| 2 | **Buljan 2020** Mol Cell Table S2 (7,316 kinase-bait interactions, >300 kinases) | Best per-edge PPI, covers understudied kinases, enables LKO PPI features | small; open access (verify CC BY) |
| 3 | Inhibitor target affinities: **HMS LINCS KINOMEscan** (open), **Davis 2011** supplement, **Drug Target Commons** API, **ChEMBL** (a connector is available here), Klaeger 2017 / Reinecke 2024 via the ProteomicsDB API or Reinecke MassIVE MSV000092248 | KB-5 occupancy model; maps decryptM and qPTM drugs | small–moderate |
| 4 | decryptM `*_Curves.zip` (PXD037285, CC0, ~3–5 GB) | EC50 for potency coherence | moderate download |
| 5 | **KinCore** domain boundaries and activation-segment positions — `external/kincore/` is **empty** | Needed by §4.4 LKO clustering and F6c activation-loop sites | small |
| 6 | CPTAC sample annotations: purity, ESTIMATE, TMT plex, ischemia time (LinkedOmics meta / PDC clinical) | KB-9 | small |
| 7 | Analog-sensitive / ATPγS direct-substrate supplements: CDK1-as (Blethrow 2008), CDK2-as (Chi 2008), ERK2 (Carlson 2011), CDK9 (Sansó 2016), AMPK (Schaffer 2015); plus Tagliabracci 2015 for FAM20C | High-directness **external test** set, stratified as "chemical-genetic" | moderate (SI tables); check PSP overlap first |
| 8 | KinoPlex supplementary data (per-site structural phospho-competence and KS candidates) | Comparator and possible feature; check whether PSP KS pairs entered its PU training | small if released |
| — | Already in hand but unused: Yaron-Barir MOESM7 SH2 matrices | KB-4 SH2 processivity | none |

Lower value: LINCS P100 (about 96 sites), ProteomeHD (covered by CPTAC co-abundance), and STRING in any form (leaks).

---

## Risks I could not resolve
- **Directness of PSP labels.** In vivo-only labels (22% of S/T pairs) include indirect relations. Evaluating every evidence class equally may favour models that learn pathway context. Suggestion: pre-register in vitro-supported test positives as a co-reported stratum, without changing the primary endpoint.
- **benchmarKIN tumour ground truth.** I described its construction from memory of Müller-Dott 2025. Confirm in the paper's Methods before finalising KB-2.
- **Bayer/Kuster 2025** (17M dose curves, 96 kinases) remains the largest potential source of quasi-direct evidence. If it is obtained, its KSR calls use Kinase Library motif plausibility, so they are not independent of F1.
- **Specificity of analog-sensitive kinases** and the in vitro peptide context of the Cantley atlases may differ from cellular specificity. Docking-dominated kinases (LRRK2 on Rab switch-II, PDPK1, MAP2K) will look poor on motif-only baselines whatever we do.
- **Statistical power.** With about 158 kinases at ≥3 test positives, per-kinase AUROC for the many kinases near the threshold is very noisy. That is for the statistics reviewer, but it limits any biological conclusion about lowly studied kinases.
