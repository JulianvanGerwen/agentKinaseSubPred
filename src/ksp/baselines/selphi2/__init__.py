"""B4 SELPHI2.0 retrained on our splits (DESIGN §8.1 B4, PRIOR-4). Owner: B4-selphi2 agent.

Reference: Maier, Petursson, Lussana, Petsalaki, MCP 2025 24(8):100994 (data/paper1.pdf) and the Apache-2.0
code in external/selphi2/selphi_2.0 (read as reference; logic ported, nothing executed).

Modules
  features  : the 44 columns of SELPHI's released feature table (residue + 43), regenerated / substituted from our
              data; label-derived ones (PSP PSSM `score`, KSEA-activity co-regulation `coreg*`) are
              ksp.crossfit LabelDerivedFeature families.
  model     : negative sampling (their random kinase-site pairs, adapted to our masks), random-forest grid search
              and bagging, OOF prediction.
Docs: docs/baseline_selphi2.md.
"""
