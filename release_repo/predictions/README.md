# Released predictions

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

**Licence: CC BY-NC-SA 3.0, derived from PhosphoSitePlus labels.** Non-commercial use only, share-alike, with
attribution to PhosphoSitePlus. See `LICENSE-PREDICTIONS.md`.

## Files

| File | Rows | Size | SHA-256 |
|---|---|---|---|
| `kinase_top200.tsv.gz` | 100,600 (503 kinases × top 200 sites) | 1,952,136 B | `98f348ed330784c8dbe52f9d327a14226975a74ab07b1e7612ebb6cf58d64567` |
| `site_top5.tsv.gz` | 1,753,680 (350,736 sites × top 5 kinases) | 34,750,123 B | `61e20f51b4eb043365ef492a0f20c31168f30bc6bacd637d40154fd0a8c3729f` |
| `coverage.tsv` | 3 methods | 659 B | — |
| `RELEASE_MANIFEST.json` | — | 1,606 B | — |

`RELEASE_MANIFEST.json` is a byte-identical copy of the file that the release run wrote. It contained no absolute
paths or user names, so nothing had to be removed. Its `path` fields are relative to the original project
directory: `results/release/<file>` there is `predictions/<file>` here. The `release_scores` entry pins the full
score table (131,144,410 kinase–site pairs, 4.1 GB parquet). That table is not included because of its size.
`python examples/read_predictions.py` checks both SHA-256 values.

## Version and provenance

- Release run `r2r3_release__20261027-012616__s0`, written on 2026-10-27. Inputs: R2 deployment run
  `r2_deploy__20261027-010348__s0` and the R3 deployment predictions.
- Model: the frozen final configuration (FREEZE 721dc197), retrained on DEV+TEST labels (L2 = human plus mouse and
  rat PSP pairs; 15,202 training pairs). R2 seeds 0–4 and R3 seeds 0–4, fused with the same weights as in CV.
- Labels: PhosphoSitePlus Kinase_Substrate_Dataset, release header "051925" (2025-05-19). Sequences and site
  numbering: UniProt human reference proteome UP000005640, release 2026_03.
- Scored universe: the 350,736 sites of U_obs (observed phosphosites) × the 503 kinases whose residue class matches
  the site.
- Fusion weights (0.6 for `score`, 0.5 for `score_site`), per-kinase calibrators, the softmax temperature and the
  ST/Y isotonic maps all come from the DEV cross-validation run `r2_fusion_r3__20261026-191901__s0`.
- Code: `scripts/models/r2_deploy.py` and `scripts/models/r2_release.py`; see `docs/model_r2.md` §7.

## Identifiers

- `site_id` = `<UniProt accession>_<residue><position>`, for example `P04637_S15`. Positions are 1-based on the
  UniProt 2026_03 canonical sequence and include the initiator methionine.
- `kinase_acc` = UniProt accession of the kinase; `gene` = its HGNC symbol as used in the project kinase list.

## Columns

`kinase_top200.tsv.gz` (sorted by kinase, then rank):

| Column | Meaning |
|---|---|
| `kinase_acc`, `gene` | kinase |
| `rank` | 1–200 within the kinase, by `score` |
| `site_id` | site |
| `score` | fused score: 0.6 · z(R2 logit) + 0.4 · z(R3 logit). Each model's logits are z-scored over all of its 131 M deployment rows. Unitless; use it for ranking within a kinase. |
| `prob` | per-kinase calibrated probability, logistic(a_k + b_k · score), with a_k and b_k fitted on DEV out-of-fold predictions. It estimates P(curated PSP-type positive \| score) under PU labels, so it is a lower bound. It is not comparable across kinases with different label counts. |
| `pct` | percentile of `score` among the kinase's candidate U_obs sites (100 = best) |
| `assignment_prob` | site-assignment probability for this kinase at this site (see below) |

`site_top5.tsv.gz` (sorted by site, then rank):

| Column | Meaning |
|---|---|
| `site_id` | site |
| `rank` | 1–5 within the site, by `assignment_prob` (descending) |
| `kinase_acc`, `gene` | kinase |
| `assignment_prob` | see below |
| `prob`, `pct` | as above, for this kinase–site pair |
| `score_site` | fused score with equal weights, 0.5 · z(R2) + 0.5 · z(R3). These weights were chosen on CV for the site-centric endpoint (M2 Hit@5). |

Only 316 of the 503 kinases appear in `site_top5`; the others are never among the top five at any site.

## What `assignment_prob` means, and its caveats

The method (DESIGN amendments A17 and A31b):
1. At each site, a softmax with temperature T is taken over the calibrated logits of all modelled kinases whose
   residue class matches the site.
2. The top kinase's softmax value is mapped through an isotonic curve. There are separate curves for S/T and for Y
   sites, fitted on labelled DEV sites.
3. The other kinases share the remaining mass in proportion to their softmax values.

Caveats:
- **It is conditional.** It assumes that the site is a direct substrate of one of the 503 modelled kinases. It is
  not the probability that the site is phosphorylated, or that any kinase targets it.
- **There is no ≥0.8 tier.** No site reaches 0.8; the S/T curve tops out at about 0.756.
- **On TEST, the ≥0.5 tier for S/T sites was calibrated**: 0.615 claimed vs 0.603 observed, n = 345. That check
  used the frozen DEV-trained model. The deployment model's maps were not re-tested, because no out-of-fold
  predictions exist for DEV+TEST.
- **Y sites are not trusted.** Tyr-kinase specificity was not demonstrated on TEST. The Y curve tops out at about
  0.445 for the top kinase.
- **Rank 1 is not always the model's top kinase.** After step 2 lowers the top kinase's mass, step 3 can give
  another kinase more than the top kinase. In this file, at 25,532 of the 350,736 sites (7.3%; S 17,826, T 7,663,
  Y 43), the rank-1 kinase by `assignment_prob` is not the kinase with the highest `prob` among the five listed. For
  the model's own ordering, re-sort the five rows by `prob`, or by `score_site` as in the M2 evaluation.
- **This also affects the "no Y site ≥ 0.5" statement.** `RELEASE_MANIFEST.json` and the report say that no Y-site
  assignment reaches 0.5. That holds for the top kinase (the quantity that `coverage.tsv` counts). But 19 Y sites in
  this file have a non-top kinase with `assignment_prob` 0.55, through step 3. All 19 are INSR/IGF1R paralog pairs
  (INSR 13, IGF1R 6). This is a property of the assignment formula; the released numbers were not changed.
- For the same reason, the count of rank-1 rows with `assignment_prob` ≥ 0.5 in this file (29,905) is 19 higher
  than Coverage@0.5 in `coverage.tsv` (29,886). The coverage table counts the top kinase's mapped value only.

## Coverage (`coverage.tsv`)

| Method | Coverage@0.5 (sites with no PSP kinase) | Coverage@0.8 | Sites in some kinase's top 1% |
|---|---|---|---|
| Final (deployment) | 29,886 (26,762) | 0 | 162,783 |
| Cantley only (B1, its own DEV-OOF calibration) | 30,258 (29,220) | 1,828 | 170,065 |
| Random scores | 6 (6) | 0 | 329,349 |

Coverage counts are not precision, and the denominator is all 350,736 U_obs sites. Cantley's ≥0.8 tier comes from
the tail of an isotonic map that was fitted on very few sites, so it is not evidence of higher precision
(`docs/REPORT_internal.md` §9).

## Seeds

Splits and annealing: seed 0. R2 and R3 final models: seeds 0–4. Random-score baseline in `coverage.tsv`: seed 0.
Seeds are recorded per run in the run manifests (`results/runs/<run_id>/manifest.json` in the project directory;
not included here).

## Spot checks (from `docs/REPORT_internal.md` §9)

- TP53 S15: ATM ranks first, assignment probability 0.70 (the known kinase).
- AKT1 S473: PDPK1, then AKT1, then MTOR third. MTOR (mTORC2) is the known kinase.
- EIF2S1 S52: PKC family; a miss. The known kinases are EIF2AK1–4.
