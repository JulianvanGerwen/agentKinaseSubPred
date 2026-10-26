# R5 pre-registered selection rule (written 2026-10-26 18:45 CET, before any fold-0 R5 metric was read)

Harness: frozen at commit 974d520 (DESIGN A28). Metrics only from `ksp.eval`. Labels only via `ksp.labels.load_dev()`.

## Configurations tried on the eval fold (fold 0) — the complete list, fixed in advance

| # | config | description |
|---|---|---|
| 0 | `r5li:li35_base` (SMOKE) | time pilot only: 200 steps, 1 seed. Recorded in REGISTRY with note `SMOKE`; **not eligible for selection** |
| 1 | `r5li:li35_base` | ESM-2 35M trunk, late interaction + R3 pair head, 3 seeds |
| 2 | `r5li:li35_cpt_aux` | 1 + Cantley-percentile distillation pretraining (1500 steps) and auxiliary loss 0.2 |
| 3 | `r5li:li150_cpt_aux` | 2 with the ESM-2 150M trunk (batch 256) |

No further eval-fold configurations will be run for selection. Diagnostic ablations (`li35_cpt_aux_noli`,
`li35_cpt_aux_lionly`) quantify the late-interaction term and are **not eligible for selection**, exactly as R3's
`noctx` ablation was not.

## Selection metric

**M1-CV = macro AUROC on E-hard over K_cv**, as reported by `ksp.eval` on fold 0. Ties (|ΔM1| < 0.002) broken by
M2 Hit@5. This is the pre-declared CV model-selection rule of DESIGN §6.8.

## Adoption gate (given by the PI brief)

R5 is carried to all 5 folds **only if**, on fold 0:

- **M1 ≥ R3 fold-0 M1 + 0.01 = 0.741 + 0.01 = 0.751**, **or**
- **M2 gain ≥ 0.02**, i.e. **M2 ≥ 0.221 + 0.02 = 0.241**.

R3's fold-0 reference is its final selected variant `reg1x5:w7_nce05_cpt_aux` (docs/model_r3.md §1, per-fold table):
fold 0 M1 **0.741**, M2 **0.221**, M3 0.267.

If no configuration clears the gate, R5 is reported as a negative result and nothing is carried forward.

## Canary gate (applies to whatever is adopted, on pooled CV)

`ksp.eval.canary_check` on the pooled 5-fold OOF predictions: mean **negative** canary AUROC on E-hard **≤ 0.57**
and **positive** canary AUROC **≥ 0.90**. Per-fold values reported alongside (DESIGN A27).

## Optimism note

Counting the time pilot, R5 will have used **4** eval-fold runs over **3** distinct selectable configurations; all
are logged in `results/cv/REGISTRY.tsv` (`model=R5`). R3 used 10 distinct configurations for comparison.
