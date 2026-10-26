# B3 Phosformer: recipe, deviations, fairness rules

**Source.** Zhou et al. 2023, *Bioinformatics* 39(2):btad046. This is `data/paper2.pdf`, not paper1, which is SELPHI2.0.

**Reference material.** I read the released artefacts in `external/phosformer/` (CC BY-NC-ND 4.0) for reference only. None of their code was copied or modified. Our implementation is in `src/ksp/baselines/phosformer/` and `scripts/baselines/phosformer_*.py`. It is built on the `transformers` library classes `RobertaModel` and `RobertaForMaskedLM`.

**How to read the Source column.**
- *Paper* means the paper states the value.
- *Released* means we inferred it from their released files.
- *Ours* means the paper is silent and we made the choice.

## 1. Architecture and input

| Item | Value | Source |
|---|---|---|
| Encoder | RoBERTa-style: 6 layers, 12 heads, d=768, FFN 3072, GELU, dropout 0.1, absolute positions, max_pos 4000, vocab 95 (character level) | Paper §2.3.1 (6/12/768); released config.json (the rest) |
| Parameters | 46.3 M (MLM model) | Measured |
| Tokens | One token per residue, plus `<s>` and `</s>`. The vocab layout equals the released vocab: `<s>`0 `<pad>`1 `</s>`2 `<unk>`3 `<mask>`4, then ASCII 33–122 → 5–94 | Paper; released vocab.json |
| Pair input | `<s>` 11-mer `</s></s>` kinase-domain `</s>`, i.e. the RoBERTa pair convention, so the central peptide residue is token 6. Our tokenizer matches `transformers.RobertaTokenizer` on 70/70 pairs, including padded ones | Paper §2.3.2 |
| Peptide | 11-mer (±5) centred on S/T/Y. If the site is within 5 residues of a terminus, the missing positions get the `<pad>` token. As in the released input pipeline, these keep attention_mask=1, and RoBERTa position ids skip them | Paper §2.1, §2.3.2 |
| Kinase input | Unaligned kinase-domain sequence | Paper |
| Head | Hidden state at token 6 → dropout → Linear(768,768) → tanh → dropout → Linear(768,2) → softmax | Paper §2.3.2 ("feed-forward layer … binary classification layer … softmax"); token index and head shape confirmed from the released state dict |

**Our kinase-domain inputs.**
- Source: the UniProt "Protein kinase" domain coordinates in `cache/index/kinases.parquet`, cut from the canonical UniProt 2026_03 sequence. Median length is 265 aa.
- **Deviation:** 5 of 503 kinases have no domain annotation (FASTK, PAN3, PEAK3, FAM20C, FAM20A). For these we use the full canonical sequence (473–887 aa). This is logged as a fallback in every manifest.
- Their own domains came from their HMM-based boundaries. The ±5 boundary-shift augmentation covers the difference.

## 2. MLM pretraining

| Item | Paper | Ours |
|---|---|---|
| Corpus | 295,320 NCBI-nr kinase domains (86,944,788 aa) + 5,338 substrate proteins (3,845,761 aa) | **Their released corpus file.** Counts are identical to the paper: 295,320 / 86,944,788 and 5,338 / 3,845,761. **It is filtered for TEST leakage (§2.1).** |
| Objective | MLM, 15% of residue tokens masked | Same. Replacement is 80% `<mask>` / 10% random amino acid / 10% unchanged. The paper says only "replace with mask tokens"; 80/10/10 is the RoBERTa/HF default. |
| Optimiser | AdamW, lr 5e-5 | AdamW lr 5e-5, betas (0.9, 0.999), eps 1e-8, weight decay 0. LR decays linearly to 0 with no warmup (HF Trainer defaults). Gradient clipping at 1.0. |
| Batch / epochs | 48/device × 8 A5000 = 384 sequences per step; 50 epochs | 384 sequences per optimiser step, 50 epochs (39,050 steps), on 1 RTX 4090. Batches are length-bucketed and split into micro-batches of ≤49k tokens. Loss = sum over masked tokens / number of masked tokens in the full batch. |
| Max length | Not stated. Released max_pos is 4000; the longest kinase domain is 4,984 aa | Sequences longer than 1022 aa are randomly cropped to 1022 each epoch. This affects 0.1% of kinase domains and 9% of substrates. |
| Validation | None | 1,000 sequences (0.33%) are held out to monitor MLM loss and accuracy. |
| Precision | Not stated | bf16 autocast plus torch.compile. |

**Time pilot** (60–100 steps on GPU 0):
- eager mode: 0.53 s/step, 220k tok/s, so 50 epochs ≈ 5.7 h;
- compiled, after warm-up: about 0.38 s/step, so ≈ 4.2 h.

The final run (filtered corpus, compiled, 49k-token micro-batches) measures **0.38 s/step**, so 38,900 steps finish in about 4.3 h (ETA 19:15). That is under the 8 h cap, so we keep the **faithful 50 epochs**. Note that DESIGN §8.1 and PRIOR-3 said "do not redo / ~6 GPU-h"; CMP-2 and DESIGN v1 asked us to reproduce it.

**Fallback if pretraining overruns.**
- Trigger: the measured ETA at step 2,000 exceeds 6 h, or the GPU-0 contention needed for Phase B makes it exceed G3 − 6 h.
- Action: restart with `--epochs E`, where E = floor(50 × available_h / ETA_h), and with the linear LR schedule rescaled to the new total. Never truncate a 50-epoch schedule mid-way: its LR would not have decayed.
- If a run is already past 60% of a 50-epoch schedule when time runs out, use the last 5-epoch checkpoint instead, and label it "B3a-partial (N/50 epochs)".
- Either way, the change is recorded here and in `docs/baseline_phosformer.md`.

### 2.1 TEST-leakage filtering of the corpus (PI AIFM2-8)
The 5,338 substrate proteins are "all substrates in their curated dataset". They therefore overlap our TEST proteins. The first launch on the unfiltered corpus was **stopped after about 1 minute and kept only as a timing pilot** (`cache/phosformer/pilot_unfiltered_aborted/`).

The final corpus is built by `scripts/baselines/phosformer_corpus_filter.py`, once `results/splits/test_proteins.tsv` exists. The TEST set T is the TEST proteins plus every protein in a substrate cluster that contains a TEST protein.

- **substrates.fa:** a sequence is removed if either
  - its UniProt accession is in T, or
  - mmseqs2 finds a hit to any T protein at ≥40% identity with ≥50% coverage of both query and target. This is the same edge rule as the substrate clusters, so non-human orthologs of TEST proteins are removed too.
- **kinases.fa** (domains only): a sequence is removed if it hits a T protein at ≥90% identity over ≥80% of the domain. This catches TEST substrate proteins that are themselves kinases (autophosphorylation and kinase-on-kinase substrates) and their near-identical orthologs.
  - The 40% cluster rule is deliberately **not** applied here. All ePK domains are about 25–50% identical to each other, so it would delete most of the corpus and change the method.
  - Counts of hits at 40% and 50% are reported in `corpus_filter.json` for transparency.
- **Residual exposure, stated rather than fixed.** DEV proteins of every CV fold stay in the corpus. This is sequence only, with no labels, so CV out-of-fold sequences were seen unlabelled during MLM. The situation is the same as for any pLM feature (ESM-2/UniRef). Per-fold pretraining (5 × 4 h) is out of budget.

## 3. Fine-tuning

| Item | Paper / released | Ours (B3a, on our folds) |
|---|---|---|
| Loss | Focal loss, γ=2, softmax over 2 classes, no α | Same (`model.focal_loss`) |
| Optimiser | AdamW lr 2e-5 | Same, with weight decay 0 and linear decay over the planned steps |
| Early stopping | "early stopping steps as 5" (unit unclear) | Patience 5 evaluations on the **inner** validation split (`ksp.splits.train_valid_split_inner`, grouped by cluster, 20%). Never on the outer fold. |
| Batch size / epochs | Not stated. The released train.csv is one pre-augmented pass of 2.41M rows | Batch 64. Each epoch is a freshly resampled and augmented pass. Validation runs every 1/4 epoch; at most 3 epochs |
| Negative ratio | 1:16 total; released train.csv has 141,972 pos : 2,271,552 neg, of which 45% hard / 55% easy, i.e. about **8 hard + 8 easy** | Choose between **1:16 (8 hard + 8 easy)** and **1:50 (25 + 25)** on CV fold 0 by M1-CV; the winner is used on all folds. Each kinase gets the same total budget |
| Hard negatives | Peptides phosphorylated by *another* kinase | Drawn from training-fold sites in `NegativeMasks.hard_candidates(k)` minus `excluded(k)`. Excluded sites are KP(k), Hamming-≤2 neighbours and paralog-labelled sites, all from `ksp.labels` |
| Easy negatives | S/T/Y 11-mers from substrate proteins with no evidence of phosphorylation by any kinase | U_ms sites of the matching residue class in training-fold proteins, minus `excluded(k)`. Our universe holds observed phosphosites, so "easy" here means **no kinase label** rather than "never phosphorylated". This is a deviation forced by the shared negative generator |
| Upsampling | Kinases with ≥50 positives upsampled to the maximum (paper: 1118; released: 1012 ≈ maximum in the training split); kinases with fewer than 50 upsampled to 50 | Same rule applied to our fold's training positives: ≥50 → max count in that fold; <50 → 50 |
| Boundary shift | Both domain ends independently ±5 | Same, resampled per example |
| Kinase masking | 5% of kinase tokens → `<mask>` (released: 4.996%) | Same, resampled per example |
| Residue class | Not enforced | Negatives match the kinase's residue class (shared generator) |

### Fairness contract (CMP-2 / PRIOR-3)
B3a and our own models R2/R3 share:
- the label set (`ksp.labels.load_dev('L2')` for training, `cv_eval_rows` for evaluation);
- the splits (`results/splits/cv_folds.tsv`);
- the residue-class policy;
- the negative masks (`NegativeMasks`: KP, Hamming-2 and paralog);
- the evaluation candidate sets (E-hard and E-all from `ksp.eval`);
- the seed list [0..4];
- the early-stopping rule (inner validation only).

Only architecture, inputs, pretraining corpus and the method's own augmentations differ.
- The augmentations (resampling, shift, masking, multilevel negatives) are part of the method. They are kept in B3a and ablated as one block if time allows.
- The negative ratio is re-tuned on our CV (F1). The window-matched arm (F2) is tracked in docs/baseline_phosformer.md as optional.

## 4. Arms
- **B3a (headline):** our MLM run on the filtered corpus, fine-tuned on our folds.
- **B3-pub (contaminated reference):** the released fine-tuned checkpoint, run zero-shot. Its training data overlaps our TEST labels, so it is never compared on equal footing.
- **B3b:** ESM-2 initialisation (ablation; optional, only if GPU time remains).
- **B3c:** random initialisation (optional).
