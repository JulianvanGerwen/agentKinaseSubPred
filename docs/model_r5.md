# R5: fine-tuned ESM-2 late-interaction kinase–site model

Owner: R5 agent, 2026-10-26. Harness frozen at commit 974d520 (DESIGN A28); every metric on this page comes from
`ksp.eval`. Labels only through `ksp.labels.load_dev()`.

Paths: code `src/ksp/models/r5/{data,model,train}.py`, drivers `scripts/models/r5_{pilot,kinase_embed,cv}.py` and
`r5_queue.sh`, cache `cache/r5/`, results `results/models/r5/`, runs `results/runs/r5_*`.

## 1. What was asked, and what the compute allowed

The brief asked for a **cross-encoder** that jointly attends over the kinase domain and the substrate window, with a
fallback to **late interaction** if full cross-encoding could not score the ~21.3M candidate pairs of a CV fold
within about 1.5 GPU-h.

**Time pilot first** (`scripts/models/r5_pilot.py` → `results/models/r5/pilot_timing.json`, one RTX 4090, bf16,
inference only, joint sequence = 300 kinase-domain + 31 window + 3 special tokens = 334):

| Model | batch | pairs/s | GPU-h for 21.3M pairs (1 fold) |
|---|---|---|---|
| esm2_t12_35M | 64 | 1,637 | **3.6** |
| esm2_t12_35M | 256 | 1,383 | **4.3** |
| esm2_t30_150M | 64 | 447 | **13.2** |
| esm2_t30_150M | 256 | 414 | **14.3** |

A full cross-encoder is therefore **2.4×–9.5× over budget for inference alone**, before any training, and ×5 folds.

The two permitted fallbacks:
- **Two-stage (top-K rerank).** Rejected. `ksp.eval.evaluate` scores the *exhaustive* E-all universe and ranks
  unscored pairs last (`missing="worst"`); M3 (top-1% recall over E-all) and the E-all secondary metrics would be
  computed against a universe the model never saw, and the leaderboard already blanks M3 for runs with partial
  E-all coverage. The harness is frozen, so it cannot be adapted. A two-stage R5 would not be comparable to R2/R3.
- **Late interaction (ColBERT-style).** Adopted. It keeps a *token-level* kinase↔substrate interaction but
  factorises the two encoders, so each side is encoded once per fold and only the MaxSim head runs per pair.
  Measured cost of scoring a whole fold: **4.4 s** for 21.28M pairs (0.0012 GPU-h), i.e. three orders of magnitude
  inside budget.

## 2. Model

**Site tower (fine-tuned ESM-2).** ESM-2 (35M or 150M) over the ±15 window (`window31` of
`cache/index/sites.parquet`; 31 residues + CLS/EOS, terminus overhang padded and attention-masked). Three
**zero-initialised** additive inputs are injected at the embedding layer, so training starts from the unmodified
pretrained trunk:
1. the frozen ESM-2 650M per-residue embedding of the same residue **in full-protein context** (`cache/r3`, exactly
   R3's site input) — the fine-tuned trunk refines R3's representation rather than replacing it;
2. a phospho-channel embedding (a flanking residue is itself a U_ms phosphosite). **The centre site's own status is
   never used** (forced 0; DESIGN A6/DL2-7);
3. the protein-context vector (ESM-2 650M protein mean), broadcast over the window.

ESM-2's `token_dropout` scaling (×0.88) and its `position_ids` construction are reproduced explicitly, because the
`inputs_embeds` path in `transformers` skips both.

Outputs: per-token late-interaction vectors `S` (B, 33, 64) and a pooled site vector `z` (centre ‖ masked mean).

**Kinase tower.** Per-residue ESM-2 650M embeddings of the catalytic domain, cached once
(`scripts/models/r5_kinase_embed.py` → `cache/r5/kin_tok.npy`, 503 × 320 × 1280 fp16, built in 4 s). On top of them
a trainable Linear → 2-layer transformer (d=192) → `m=8` learned query tokens attend over the domain → `Q` (K, 8, 64).
The R3 descriptor blocks (kinase-domain mean embedding, Cantley score-profile PCA, Manning group/family, descriptor
dropout p=0.2) give a global kinase vector `u` (K, 128). Caching the kinase *tokens* is what makes the kinase tower
cost ~1% of a step instead of dominating it.

**Head.**

```
score(s,k) = gamma * sum_{j<8} max_t <Q[k,j], S[s,t]> / sqrt(64)     <- late interaction (new)
           + <W z_s, u_k>/sqrt(h) + MLP(z_s * u_k) + b(u_k)          <- R3's bilinear pair head
```

The second line makes R5's head a strict superset of R3's, so the late-interaction term is an addition that can be
ablated (`use_li` / `use_pair`).

**Losses (unchanged from R3, re-used from `ksp.models.r3.model`).** Balanced masked BCE (positives weighted
∝ n_k^-0.5) + 0.5 × filtered listwise InfoNCE over kinases per positive (denominator excludes the site's other true
kinases, masked entries and Par(k)) + optional Cantley-percentile distillation. Targets, masks (KP, Hamming ≤ 2,
Par(k)) and the training/validation protocol come from `ksp.models.r3.data.build_targets` and
`ksp.splits.train_valid_split_inner` verbatim, so R5 and R3 are directly comparable.

**Canaries.** Canaries have no kinase domain, so they get random-but-fixed pseudo-domain token blocks (16 tokens,
seeded, scaled to the real token s.d.) — the same role `ksp.canary.canary_descriptors` plays for R3. They carry no
sequence information.

**Efficiency note (identical metric, cheaper).** The inner-validation early-stopping metric only ever reads
label-carrying rows (positives, and negatives restricted to rows carrying some positive), so R5 scores only those
rows for it. The value is identical to R3's; the cost is ~50× lower.

## 3. Results

(filled in below once the runs complete)

## 4. Runtime

(filled in below)

## 5. Deviations and assumptions

1. **Late interaction instead of a full cross-encoder**, for the measured reason in §1. The kinase and the substrate
   still interact at token level, but not through joint self-attention.
2. **Kinase-domain embeddings are computed on the isolated domain sequence** (the `emb_domain_iso` convention of
   `src/ksp/features/esm_kinases.py`), not in full-protein context.
3. **Long domains** (35 of 503, e.g. SRPK1/2, MASTL, CDC7) keep the N- and C-terminal 160 residues and drop the
   middle insert. **5 of 503 kinases have no annotated catalytic domain** and fall back to the first 320 residues;
   both are flagged in `cache/r5/kin_meta.json`.
4. **GPUs 1–3** as assigned (GPU 0 is Phosformer's), ≤ 12 cores (3 jobs × 4 threads).
5. **3 seeds per run**, not R3's 5, for time. This is conservative: it can only understate R5 against R3's 5-seed
   numbers.
6. **No git** (agent rule), so manifests carry the hash the manifest module queries.
