## 10. Limitations

1. **No demonstrated Tyr-kinase specificity.**
   - On TEST, Y kinases score M1 0.595 against the site-popularity floor B0-pop at 0.585. The final model beats that floor on only 10 of 23 Y kinases.
   - The review rates this **NOT TRUSTED** as a skill claim. Cantley is at chance (0.496).
   - T3 on Y kinases is 0.531.
   - Y-kinase M1 fell from CV to TEST: 0.667 → 0.595 (paired Δ −0.072 ± 0.033). SYK, IGF1R and LYN collapse.
   - No Y site can receive an assignment probability above 0.445.
   - On CV the Y gap was already large (0.63 vs 0.81). F9 did not help. Part of the gap is endpoint composition (hubs, indirect in-vivo labels).
   - Sources: `docs/audit/audit_test.md` §0 #10, §4; `results/test/721dc1977366/final_t3_per_kinase.tsv`; `docs/strata_cv.md`.
2. **No high-confidence assignment tier.**
   - After A31 there is no ≥0.8 tier. TEST coverage@0.8 is 0, and coverage@0.5 is 5,773 sites, all S/T.
   - The S/T ≥0.5 tier is honest on labelled TEST sites (0.615 claimed vs 0.603 observed, n = 345).
   - Y probabilities are mildly overconfident (0.41 vs 0.34 in [0.2, 0.5)).
   - These are precisions on labelled sites only (PU); they do not transfer to unlabelled sites.
   - Sources: `final_reliability.tsv`; `audit_test.md` §6.
3. **M2 favours well-studied kinases.**
   - CV Hit@5 is 0.142 at 5–19 labels against Cantley's 0.162. A32(b) debiasing failed.
   - TEST cannot refute this: it holds only 13 such kinases, the best-labelled of the bin.
   - TEST M2 (0.363) is K_eval composition, not improvement; CV on the same 117 kinases is 0.348 (`audit_test.md` §2).
4. **The margin over Cantley is partly a weak-baseline draw.**
   - The TEST M1 margin (+0.075 [0.046, 0.099]) exceeds the CV margin (+0.063 [0.048, 0.074]).
   - That is because Cantley fell by 0.025 ± 0.013 on the same kinases, mostly on Y kinases. Quote both margins (`audit_test.md` §3).
5. **Priming.**
   - The TEST primed-site advantage (0.897 vs 0.704; 6 kinases, 65 pairs) is only against non-priming-aware Cantley.
   - On unprimed positives the final model and Cantley are near-tied (0.740 vs 0.725).
   - On CV, priming-aware Cantley beat the model on primed sites (0.945 vs 0.893), and it was not a TEST comparator.
6. **Baselines are our reimplementations.**
   - B3: the fidelity check passed (0.888 vs 0.885), but it has 1 seed and a PI-fixed 1:16 negative ratio.
   - B4: 10 of 100 bags, a substituted PSSM/coreg, and 66.6% of TEST pairs on NA-default site features.
   - The B4 M3 win is modest: Δ +0.059 [0.010, 0.120], Holm p 0.018.
   - The superiority claims are against "our faithful retrained reimplementations under this benchmark", not against the published systems as distributed.
7. **TEST strata carry no CIs, so they are descriptive.**
   - Small strata are anecdotal: PIKK (4 kinases), FAM20C (1 kinase), STY (2 in K_eval), TKL (3) and the 5–19 bin (13).
   - Hub-site and external-tier-B gains are mostly popularity (B0-pop 0.788 and 0.724).
8. **Items the single batch could not provide, and which cannot be added without a second sealed access:**
   - hub-excluded Y;
   - A16 calibration on the temporal and external strata;
   - TEST coverage baselines;
   - a paired T3 test;
   - stratum CIs.

   See §7.9.
9. **Phosphoproteomics adds almost nothing to CV M1** (F6 ≈ +0.002; the no-phospho variant scores 0.771). Downstream, these features may help on unseen drugs, but the evidence is indirect.
10. **Downstream parity, not superiority, with PSP**, plus an unseen-drug gap of 0.04–0.06 (`docs/downstream.md` §7).
11. **The canary gate passes only pooled.**
    - Two of five CV folds exceed 0.57 (0.575 and 0.676).
    - Every model sits above 0.5, including label-free controls. The sampler explanation is plausible but unproven.
12. **LKO generalisation is weaker on novel sites.**
    - On CV T2, sites with no other-kinase label score 0.703 against a 0.757 control.
    - T2 and T3 are "unseen labels", not kinase-naive: label-free Cantley percentiles enter R3 pretraining.
    - TEST T3 (0.752) has no paired test against a baseline.
13. **Optimism and post-freeze changes.**
    - About 90 configurations were tried on CV. TEST confirms no M1 optimism: paired Δ on common kinases +0.001 ± 0.012.
    - A24 was a post-freeze metric correction applied to all methods; it widened learned-vs-Cantley CV margins by about 0.01.
14. **Label noise and PU.** In-vivo-only labels are enriched for indirect events. All precision and probability statements are lower bounds.
15. **Unbuilt items.**
    - F8 (Sugiyama), docking motifs, F6d and U_all (unobserved S/T/Y).
    - The deployment model is still running (§9).
    - `scripts/run_all.sh`, promised in DESIGN §11.8, does not exist.
16. **Data gaps** (DATAINEED): the SIGNOR direct export, decryptM curve fits, CPTAC plex/ischemia confounders and PTB matrices were unavailable; workarounds are in use.
