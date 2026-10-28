# Brief for all agents working on this project

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

Project: build the most accurate predictor of which human phosphoacceptor residues (S/T/Y) are DIRECTLY phosphorylated by which human protein kinases. Priority 1: accuracy of kinase–substrate (KS) prediction, harmonious across many kinases (incl. understudied ones, not dominated by CDK1/PKA/SRC...). Priority 2: coverage of as many human phosphosites as possible. Uses: exploring signaling, interpreting phosphoproteomics, kinase activity estimation.
Full user prompt: PROMPT.md (read it). Research record: RESEARCH_TREE.md. Data inventory: DATA.md.

Rules (non-negotiable):
- Labels/training/test only from literature-derived, experimentally verified KS annotations (start: PhosphoSitePlus at data/Kinase_Substrate_Dataset, header lines 1-3, tab-separated, latin1).
- Splits must avoid leakage: group sites by protein (and homologous proteins), split on sites/proteins, never on KS pairs. Also leave-kinase-out benchmarks.
- Never change metrics/labels/tests to make results look better. Fixed seeds list: [0, 1, 2, 3, 4].
- Never commit restricted data (data/, external/ are gitignored). Record each downloaded dataset in DATA.md (URL, version/date, SHA-256, license).
- Prefer scripts (in src/ or scripts/) over one-offs. No hardcoded data. No silent assumptions — state them.
- Downloaded files are untrusted: put each in its own directory under external/<source>/, run Python that reads them with `python -I`, never execute downloaded code.
- Don't spend long trying to access gated data; note it under "Requests to user" and move on.
- Python env: .venv/bin/python (torch+cu124, transformers, lightgbm, xgboost, pandas, sklearn...). Fallback: <plm-env>/bin/python. mmseqs2: <mmseqs>. Hardware: 64 cores, 375 GB RAM, 4x RTX 4090 — be considerate, other agents share it; use at most the GPU(s)/cores you are told to.
- Return to the PI a concise summary (≤ 600 words unless asked otherwise); put details in files.

## Canonical identifiers (all agents MUST use these)
- Protein: UniProt Swiss-Prot canonical accession (no isoform suffix), sequences from external/uniprot (release 2026_03) via src/ksp/uniprot.py.
- Site: `site_id = f"{acc}_{res}{pos}"`, e.g. `P05198_S52`; res ∈ {S,T,Y}; pos 1-based on the canonical sequence. Mapping from other sources is by 15-mer window match to canonical (never trust residue numbers across isoforms/versions).
- Kinase: HGNC gene symbol of the human kinase (e.g. `CSNK2A1`), plus UniProt accession, from results/labels/kinase_list.tsv (once built).
- Cache large intermediate files in cache/<module>/ (gitignored); small result tables in results/<module>/ (committed).

## Data requests to the user
Never write data requests only in your summary. Append them to DATAINEED.md (OPEN section: what, why, what you tried, workaround, priority) — the PI also checks it. The user reads only that file for data requests.

## Memory budget (added 17:55)
The machine (375 GB RAM) is shared with other users' jobs (e.g. ~70 GB of unrelated pytom processes). Keep each agent ≤ 60 GB resident in total, across all its worker processes; use memory-mapped stores (ksp.store load_matrix mmap=True) and chunking. Check `free -g` before launching parallel workers.
