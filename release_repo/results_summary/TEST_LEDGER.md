> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*


## STARTED TEST batch 2026-10-27 00:10:01 CET (attempt 1)

<!-- KSP-LEDGER {"attempt": 1, "batch_id": "721dc1977366-1793056201", "code_identity": "264d94313576608c77678f337f55d7450323c64daf549355340ae819f2030991", "code_sha256": {"__init__.py": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "canary.py": "511c41b156d40a1fe34020a85fb3784963645f9f39cf628f13e359c0d41a78af", "crossfit.py": "58fec1579a64cc7c42143e0b9d078345fe394875a9786b26fa8d2717cf40afe6", "deploy.py": "a4d5ede76e7144c24f8dcec3d6da80ef8fa07e4cfdd551107743a4b3d2a9d2c0", "eval.py": "84457e83dec9123cfe7cf0d54eb6ba05fb608f384ba5331a0cb9c835c6f4de16", "eval_test.py": "738a8c416aabf2247cbdbeaa13cb614c5dfe7449980b5d316ea6663152cdea17", "ids.py": "fc57f5b553ecb76950d5966815f5457d89ec5d91762edc6a137d7dca3034b10c", "index.py": "90806a63d3dc8947953e93b4dd603097784b8b723a45fcfd52c95371ef758a4f", "io_psp.py": "27978501a31ca757f74b7b90e30d5dc7fdae6204eb04a6a411220a125eeb08b0", "kinases.py": "e538a7e75c2638d5f33bc19641ada5fcd2299cec1025e8d3d523d006bfe2c4b4", "labels.py": "11dc2660e7f7a9ea182d2bf93f6ce592b0d30708223167337a9c85364d552e87", "manifest.py": "f63e634710f2a504cb82c71e4fadc50c5cf2709403e10bd7cd9e6bb713613c3d", "sites.py": "97408ec6bb2ba0820e64f2d20e41495aa433cb0f7724ed8db7d208449c073c11", "splits.py": "4e391b3e23841b62b24ce3b95bdab925b42772e914f8f357c6f6ad8dc42e27c2", "store.py": "a713e9103746e1c5660f649f1831aea4d337e7da1954762bbf86042f429637c5", "strata.py": "a683e85d23486a94a661403d644b41939b86f8fbab3db7b5fa9c7b8426fff625", "uniprot.py": "9f183f51f33d8675f41e8d7ad2803490aa256c594d69620fbe6582a619d724ac"}, "event": "STARTED", "freeze_sha": "721dc1977366ee7bfb9af187345f72180c537d025b78cac4a7cbfb3e6f7b8138", "git_head": "c3e2a1448eb78e4a030d4498dadbd82db827676e", "inputs": {"cache/runs/b3_retrained_test__20261026-232229__s0/test_predictions.parquet": "a468cdbc8d57cfabcef58b7880dcb985f9ab2fa9ee54d44f61458dee3bc464a3", "cache/runs/b4_selphi2_final__20261026-212619__s0/test_predictions.parquet": "89365e66920f602558c0370f2ff4bb705d12b8aec14e64220a63b93ea5b4cbfb", "cache/runs/r2_final_t3__20261026-221039__s0/test_t3_predictions.parquet": "344433cbb4928ce7f91ffd19c705c186fa996fec7a3c1bdff813cefbb37f053f", "cache/runs/r2_final_test__20261026-221039__s0/test_predictions.parquet": "2aa381c9fbccaba1b2da022af453276b0ca960e878a72e3be2017ff67ee9f17b", "cache/runs/r2r3_final__20261026-230138__s0/test_predictions.parquet": "1d95b271672b1db20ae63c2918a32a911e6800dc969b53b6514a3dceff7ce9c4", "cache/runs/r2r3_final__20261026-230138__s0/test_t3_predictions.parquet": "d62195b8f3a66df64cb0f91f5261469d1c80ed1204922bd016470aaafc8ad1b1", "cache/runs/r3_final_distill_test__20261026-220424__s0/test_predictions.parquet": "6223c72c3a93d8082c957753c541dc0737a8e9b5921e90ce19e29af71828d384", "cache/runs/r3_final_distill_test__20261026-220424__s0/test_t3_predictions.parquet": "9c6c403ea27690d8f5c337b47525a4cb84452ad0cc9be1cf31de78a4b78b8874", "cache/runs/simple_test_b0_miss__20261026-212754__s0/test_predictions.parquet": "07c434d093827d27d413e4563bd15c77dba097128db85ee03f7f9f6b53f21798", "cache/runs/simple_test_b0_pop__20261026-212754__s0/test_predictions.parquet": "234ee7e5c075aef7e82b628ecf8d7e469b9c1a8eed981b0fa642933789b82a31", "cache/runs/simple_test_b0_prior__20261026-212754__s0/test_predictions.parquet": "0bf51f64c36273dca115522292571dc4eccd647effd9a9f50ec953545f12fdf3", "cache/runs/simple_test_b1_arm_i__20261026-214815__s0/test_predictions.parquet": "609988b9e84a48dc9857403071f6f56a1641172d249ea7608e68eda80d493655", "cache/runs/simple_test_b1_arm_ii__20261026-214850__s0/test_predictions.parquet": "f3213ed39439e5b995f13c34817e9447870912b2c087c9271f816750cae776f4", "cache/runs/simple_test_b2_pssm__20261026-212755__s0/test_predictions.parquet": "bc78031604405008ae6de8b7cfd4161460166faeab234e4367980a051ea52879", "cache/runs/simple_test_b2_pssm__20261026-212755__s0/test_t3_predictions.parquet": "25731be1ec6a214e66ce4ef1063e6e3a58ac78a66edc10d804f3ace3743e175a", "cache/runs/simple_test_b7_lr__20261026-212809__s0/test_predictions.parquet": "f9977e688f1956311e92b1b1209bfd4ec7cc7fff623015c5c74123418c14da6a", "cache/runs/simple_test_b7_lr__20261026-212809__s0/test_t3_predictions.parquet": "2058cc30380d7f791c07d86c485cc06906763043bf0482a53158793308919eeb"}, "tag": "TEST"} -->

FREEZE.json sha256 721dc1977366ee7bfb9af187345f72180c537d025b78cac4a7cbfb3e6f7b8138; git HEAD c3e2a1448eb78e4a030d4498dadbd82db827676e; harness 974d520; 12 configs; pre-flight {"final": 25922452, "final:T3": 25922452, "R2": 25922452, "R2:T3": 25922452, "R3": 25922452, "R3:T3": 25922452, "B1_i": 19215814, "B1_ii": 19215814, "B2": 25922452, "B2:T3": 25922452, "B7": 25922452, "B7:T3": 25922452, "B0_prior": 25922452, "B0_pop": 25922452, "B0_miss": 25922452, "B3": 25922452, "B4": 25922452}

## COMPLETED TEST batch 2026-10-27 00:57:29 CET

<!-- KSP-LEDGER {"attempt": 1, "batch_id": "721dc1977366-1793056201", "event": "COMPLETED", "freeze_sha": "721dc1977366ee7bfb9af187345f72180c537d025b78cac4a7cbfb3e6f7b8138", "outputs": "results/test/721dc1977366", "tag": "TEST"} -->

Configs 12; B=10000; runtime 2838.1s; outputs results/test/721dc1977366

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
