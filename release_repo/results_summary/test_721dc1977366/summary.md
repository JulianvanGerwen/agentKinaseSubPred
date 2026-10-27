> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

**Standalone metrics under missing='worst' (unscored pairs ranked last); NOT head-to-head - see the matched comparisons below.**

| config | role | M1 | M2 Hit@5 (score col) | M3 top1% | M4 | T3 M1 | T3 missing kinases | missing/required | train s | infer s |
|---|---|---|---|---|---|---|---|---|---|---|
| final | final | 0.7940 | 0.3627 (score_site) | 0.3632 | 0.7987 | 0.7518 | 0 | 0/13230299 | 1136.6 | 316.9 |
| R2 | ablation | 0.7783 | 0.2637 (score) | 0.3118 | 0.7813 | 0.7436 | 0 | 0/13230299 | 417.9 | 298.9 |
| R3 | ablation | 0.7748 | 0.3223 (score) | 0.2920 | 0.7691 | 0.7030 | 0 | 0/13230299 | 718.7 | 3.1 |
| B1_i | baseline | 0.7137 | 0.1628 (score) | 0.1607 | 0.6989 |  |  | 1199434/13230299 | 0.0 | 4.36 |
| B1_ii | baseline | 0.7137 | 0.1888 (score) | 0.1607 | 0.6989 |  |  | 1199434/13230299 | 0.0 | 0.97 |
| B2 | baseline | 0.6848 | 0.1549 (score) | 0.1383 | 0.6640 | 0.6157 | 0 | 0/13230299 | 3.33 | 39.74 |
| B7 | baseline | 0.7249 | 0.2432 (score_site) | 0.1998 | 0.7166 | 0.7083 | 0 | 0/13230299 | 89.12 | 287.03 |
| B0_prior | control | 0.5000 | 0.1088 (score_site) | 0.0100 | 0.5000 |  |  | 0/13230299 | 0.0 | 0.75 |
| B0_pop | control | 0.5785 | 0.0203 (score) | 0.1273 | 0.5878 |  |  | 0/13230299 | 0.0 | 0.41 |
| B0_miss | control | 0.5654 | 0.0367 (score) | 0.0659 | 0.5795 |  |  | 0/13230299 | 19.02 | 30.42 |
| B3 | baseline | 0.7079 | 0.1634 (score) | 0.1040 | 0.6878 |  |  | 0/13230299 | 3152.41 | 6478.84 |
| B4 | baseline | 0.7047 | 0.2327 (score) | 0.3044 | 0.7008 |  |  | 0/13230299 | 737.0 | 2043.4 |

T3: each kinase scored by the model trained with its LKO family removed; T3 M2 mixes LKO models.

**Primary superiority tests (m = 9; Holm over these; B = 10000)**

| comparison | endpoint | status | diff | 95% CI | p_boot | p_holm | beats | matched kinases | frac A kept |
|---|---|---|---|---|---|---|---|---|---|
| final vs B1_i | M1 | tested | 0.0745 | [0.0465, 0.0991] | 0.0002 | 0.0018 | True | 386 | 0.7413 |
| final vs B1_i | M3 | tested | 0.2051 | [0.1650, 0.2346] | 0.0002 | 0.0018 | True | 386 | 0.7413 |
| final vs B1_ii | M2 | tested | 0.1787 | [0.1469, 0.2130] | 0.0002 | 0.0018 | True | 386 | 0.7413 |
| final vs B3 | M1 | tested | 0.0861 | [0.0594, 0.1051] | 0.0002 | 0.0018 | True | 503 | 1.0000 |
| final vs B3 | M2 | tested | 0.1993 | [0.1507, 0.2290] | 0.0002 | 0.0018 | True | 503 | 1.0000 |
| final vs B3 | M3 | tested | 0.2592 | [0.2131, 0.2908] | 0.0002 | 0.0018 | True | 503 | 1.0000 |
| final vs B4 | M1 | tested | 0.0894 | [0.0630, 0.1172] | 0.0002 | 0.0018 | True | 503 | 1.0000 |
| final vs B4 | M2 | tested | 0.1300 | [0.0908, 0.1565] | 0.0002 | 0.0018 | True | 503 | 1.0000 |
| final vs B4 | M3 | tested | 0.0588 | [0.0096, 0.1198] | 0.0178 | 0.0178 | True | 503 | 1.0000 |

**Secondary comparisons (Holm separately)**

| comparison | endpoint | status | diff | 95% CI | p_boot | p_holm | beats | matched kinases | frac A kept |
|---|---|---|---|---|---|---|---|---|---|
| final vs R2 | M1 | tested | 0.0158 | [0.0090, 0.0224] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs R2 | M2 | tested | 0.0990 | [0.0731, 0.1277] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs R2 | M3 | tested | 0.0514 | [0.0335, 0.0714] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs R2 | M4 | tested | 0.0174 | [0.0068, 0.0246] | 0.0010 | 0.0034 | True | 503 | 1.0000 |
| final vs R3 | M1 | tested | 0.0192 | [0.0076, 0.0302] | 0.0010 | 0.0034 | True | 503 | 1.0000 |
| final vs R3 | M2 | tested | 0.0405 | [0.0074, 0.0595] | 0.0130 | 0.0130 | True | 503 | 1.0000 |
| final vs R3 | M3 | tested | 0.0712 | [0.0406, 0.0970] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs R3 | M4 | tested | 0.0296 | [0.0168, 0.0404] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B7 | M1 | tested | 0.0692 | [0.0438, 0.0881] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B7 | M2 | tested | 0.1195 | [0.0833, 0.1459] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B7 | M3 | tested | 0.1634 | [0.1230, 0.1975] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B7 | M4 | tested | 0.0821 | [0.0506, 0.0978] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B2 | M1 | tested | 0.1093 | [0.0761, 0.1395] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B2 | M2 | tested | 0.2078 | [0.1683, 0.2385] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B2 | M3 | tested | 0.2250 | [0.1706, 0.2626] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B2 | M4 | tested | 0.1346 | [0.0984, 0.1573] | 0.0002 | 0.0034 | True | 503 | 1.0000 |
| final vs B1_i | M4 | tested | 0.0865 | [0.0573, 0.1034] | 0.0002 | 0.0034 | True | 386 | 0.7413 |

Footnote: 95% CIs are unadjusted paired percentile intervals; 'beats' requires Holm-adjusted p < 0.05 AND CI lower bound > 0. Comparisons use the matched universe (PRIOR-1; matched-kinase M2, A27).
