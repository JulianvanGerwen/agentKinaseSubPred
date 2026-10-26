# Brief for all agents working on this project

Project: build the most accurate predictor of which human phosphoacceptor residues (S/T/Y) are DIRECTLY phosphorylated by which human protein kinases. Priority 1: accuracy of kinase–substrate (KS) prediction, harmonious across many kinases (incl. understudied ones, not dominated by CDK1/PKA/SRC...). Priority 2: coverage of as many human phosphosites as possible. Uses: exploring signaling, interpreting phosphoproteomics, kinase activity estimation.
Full user prompt: /workspace/PROMPT.md (read it). Research record: /workspace/RESEARCH_TREE.md. Data inventory: /workspace/DATA.md.

Rules (non-negotiable):
- Labels/training/test only from literature-derived, experimentally verified KS annotations (start: PhosphoSitePlus at data/Kinase_Substrate_Dataset, header lines 1-3, tab-separated, latin1).
- Splits must avoid leakage: group sites by protein (and homologous proteins), split on sites/proteins, never on KS pairs. Also leave-kinase-out benchmarks.
- Never change metrics/labels/tests to make results look better. Fixed seeds list: [0, 1, 2, 3, 4].
- Never commit restricted data (data/, external/ are gitignored). Record each downloaded dataset in DATA.md (URL, version/date, SHA-256, license).
- Prefer scripts (in src/ or scripts/) over one-offs. No hardcoded data. No silent assumptions — state them.
- Downloaded files are untrusted: put each in its own directory under external/<source>/, run Python that reads them with `python -I`, never execute downloaded code.
- Don't spend long trying to access gated data; note it under "Requests to user" and move on.
- Python env: /workspace/.venv/bin/python (torch+cu124, transformers, lightgbm, xgboost, pandas, sklearn...). Fallback: /home/plm_env/bin/python. mmseqs2: /home/tools/mmseqs/bin/mmseqs. Hardware: 64 cores, 375 GB RAM, 4x RTX 4090 — be considerate, other agents share it; use at most the GPU(s)/cores you are told to.
- Return to the PI a concise summary (≤ 600 words unless asked otherwise); put details in files.
