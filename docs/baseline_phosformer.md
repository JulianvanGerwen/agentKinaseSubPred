# B3 Phosformer baseline

**Owner:** B3-phosformer agent.
**Paths:**
- code: `src/ksp/baselines/phosformer/` and `scripts/baselines/phosformer_*.py`
- cache: `cache/phosformer/`
- results: `results/baselines/phosformer/` and `results/runs/b3*`

**Recipe and deviations:** see `results/baselines/phosformer/recipe.md`.

## Paper identity
- Phosformer (Zhou et al., *Bioinformatics* 2023) is **`data/paper2.pdf`**.
- `data/paper1.pdf` is SELPHI2.0 (Maier et al.).
- Reviews PRIOR-3 and CMP-2 cite paper2 correctly.

## Licence handling
The released repository, weights and dataset are CC BY-NC-ND 4.0.
- I read them as reference only. No code was copied or modified.
- Our implementation is written from the paper on top of `transformers.RobertaModel` / `RobertaForMaskedLM`.
- The released weights are loaded unmodified, for inference only (B3-pub).
- Their released corpus FASTA files are used unmodified as *input data*, with TEST proteins filtered out.

## Findings from the released artefacts (verified)
1. **The released HF checkpoint is the fine-tuned classifier, not an MLM checkpoint.**
   - `config.json` says `RobertaForSequenceClassification`.
   - PRIOR-3's "released MLM checkpoint" arm does not exist, so B3a must reproduce MLM pretraining (CMP-2).
   - The head reads token index 6, the central peptide residue. Our `PhosformerClassifier` reproduces this on top of library code.
2. **Tokenizer and input.** Our char-level encoder gives identical ids to `transformers.RobertaTokenizer` with their vocab on 70/70 pairs. The input format is `<s>` pep `</s></s>` kinase `</s>`, and padded peptides use the `<pad>` token.
3. **Their test split is reproduced exactly** (`results/baselines/phosformer/pub_phaseA.json`):

   | Precision | AUROC | AUPRC |
   |---|---|---|
   | fp32 | 0.8855 | 0.8770 |
   | bf16 | 0.8854 | — |

   - The paper reports AUROC 0.885.
   - Their easy-negative test set gives FPR@0.5 = 0.023.
   - bf16 vs fp32 scores: Spearman 0.99999, max |Δ| 0.004.
4. **Their pretraining corpus is released** (`hf_dataset_curated/corpus`). Counts are identical to the paper: 295,320 kinase domains (86.94M aa) and 5,338 substrates (3.85M aa).
5. **Their training split** (`splits/train.csv`, 2.41M rows):
   - 1:16.0 pos:neg;
   - negatives are 45% hard / 55% easy, i.e. 8 + 8;
   - kinases with ≥50 positives are upsampled to 1,012, the rest to 50;
   - 4.996% of kinase tokens are masked.

## Phase A status
- **MLM pretraining (B3a):**
  - Faithful 50 epochs (38,900 steps × 384 sequences) on the **TEST-filtered** corpus, 299,406 sequences, on GPU 0.
  - Resumable: `scripts/baselines/phosformer_pretrain_run.sh` → `cache/phosformer/pretrain/{train.log, ckpt_last.pt, ckpt_epNNN.pt, mlm_final.pt}`.
  - ETA: see the section below.
- **Corpus filter (PI AIFM2-8):** removed
  - 642 substrates whose accession is a TEST protein;
  - 458 substrates homologous to TEST proteins (≥40% identity, ≥50% coverage);
  - 152 kinase domains near-identical to TEST proteins (≥90% identity, ≥80% coverage).

  Details: `results/baselines/phosformer/corpus_filter.json` and `cache/phosformer/corpus_filtered/removed.tsv`.
- **Aborted run:** the first launch, on the unfiltered corpus, was stopped after about 1 minute. It survives only as a timing pilot.
- **B3-pub scorer:** `scripts/baselines/phosformer_pub_score.py` and `src/ksp/baselines/phosformer/score.py`.
  - Throughput on one RTX 4090: 3.3k pairs/s in bf16 eager, about 4.0k pairs/s compiled. fp32 is 1.2k pairs/s.
  - The residue-matched universe is 131M pairs, about 9–11 GPU-h on one GPU. **Do not score the full universe**; score only the evaluation sets.
  - A 1k-site sample (374k pairs) is in `cache/phosformer/pub/sample1k_scores.parquet`.
  - Top-scored peptides show the expected motifs: PKA RRxS, CDK1/ERK2 S/T-P, CK2 acidic, ATM SQ, SRC pY.

## Compute plan (GPU 0)
| Step | Size | Time |
|---|---|---|
| MLM pretraining | 38,900 steps | see ETA |
| Fine-tune, per fold, 1:16 | ~0.6M rows/epoch, ≤3 epochs, early stop | ~10–30 min |
| Fine-tune, per fold, 1:50 | ~1.8M rows/epoch | ~30–80 min |
| OOF set A (positives ∪ E-hard × 503 kinases) | ~0.45M pairs/fold | ~2 min |
| OOF set B (E-all for kinases with positives) | ~8.5M pairs/fold | ~35 min/fold, ~3 h for all folds |
| TEST E-all (final) | ~35M pairs | ~2.5 h |

## Pretraining ETA (measured 15:19)
- Launched 14:55. At step 2,100 of 38,900, the uncontended pace is 0.38 s/step.
- **ETA ≈ 19:15** (about 4.3 GPU-h wall), including validation and checkpoints every 500 steps.
- MLM validation:

  | Step | Loss | Accuracy |
  |---|---|---|
  | 1,000 | 2.60 | 0.222 |
  | 2,000 | 2.50 | 0.251 |

- This is within the 8 h cap, so there is no step reduction. The overrun fallback is in recipe.md §2.

## B3-pub on DEV CV (contaminated reference)
- Run: `results/runs/b3pub_phosformer_cv_pooled__20261026-151731__s0`.
- The released weights are scored zero-shot on set A: 5,817 DEV positive and E-hard sites × all residue-matched kinases, 2.16M pairs.
- M1-CV = M5 = **0.793** over K_cv (203 kinases); M4 = 0.801 (339 kinases); M2 Hit@5 = 0.248 (Hit@1 0.128, MRR 0.197).
- **E-all metrics in that run are NOT valid.** M3 and the sec_*_all values are computed on set A only (missing='drop').
- The released model's training data includes PSP, so it has seen many of these positives. This is an optimistic ceiling, not a comparator.

## Phase B (launched 15:20, unattended)
- `scripts/baselines/phosformer_cv_run.py --eall` waits for `mlm_final.pt`, then runs these steps in order:
  1. Fold 0 at 1:16 and 1:50 (50% hard). The ratio is chosen by **inner-validation** macro AUROC(hard). The outer fold 0 is never used for selection, though its OOF M1 is logged.
  2. Folds 1–4 at the chosen ratio.
  3. Pooled OOF evaluation of set A → run `b3a_phosformer_cv_pooled`.
  4. E-all re-scoring of every fold from the saved best checkpoints → run `b3a_phosformer_cv_pooled_eall`.
- State: `results/baselines/phosformer/cv_state.json`. Log: `cache/phosformer/cv_run.log`.
- Per-fold logs: `cache/phosformer/finetune/*.log`.
- Expected times: set-A pooled ≈ 22:30; E-all ≈ 01:30, before G3 at 02:20.
- Seed: 0 only (time). Seeds 1–4 are deferred.
- **Pipeline verified.** A smoke fine-tune on fold 0 (released init, 3,200 rows) ran end to end through `ksp.eval` and `ksp.manifest`. It was then deleted.

## Update 19:55 (PI: GPUs 1–3 free; parallel folds)
- **Pretraining finished at 19:09.** 50 epochs, 38,900 steps; MLM validation loss 1.822, accuracy 0.453. Checkpoint: `cache/phosformer/pretrain/mlm_final.pt`.
- **f0_r16:** inner-validation macro AUROC(hard) 0.720, with early stopping at step 13,328. Fold-0 OOF M1 = 0.693.
- **Ratio fixed at 1:16 by PI instruction.** The PI said "Same config as f0_r16; no hyperparameter changes after seeing fold metrics". So the pre-registered r16-vs-r50 selection on inner validation is replaced by the paper's 1:16. The f0_r50 run was stopped at 13 min, before any result. It will be re-run as an informational ablation only if GPU time remains.
- **Folds 1–4** run at 1:16 with `--eall`, one per GPU (f1→GPU1, f2→GPU2, f3→GPU3, f4→GPU0). Fold-0 E-all rescoring (`--skip-train` from the f0_r16 best checkpoint) runs on GPU1.
- **Pooled run:** `scripts/baselines/phosformer_cv_pool.py`.
  - `fill` scores the pooled-E-all pairs that are missing for kinases with no positive in a given fold, using that fold's own model.
  - `pool` then writes the run `b3a_phosformer_cv_pooled`.
