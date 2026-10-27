> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.
>
> *Copied from the project record for this release. Absolute paths were replaced by relative paths or placeholders; paths refer to the original project directory, and most `cache/` and `results/` files they name are not in this repository.*

# Design review brief (round N)
You are an independent reviewer of DESIGN.md for a kinase–substrate prediction project.
Read: PROMPT.md (fixed goal & rules), DESIGN.md, LITERATURE.md, and as needed docs/lit/*.md, docs/data_profile.md, docs/data_acquisition.md.
Rules for reviewers:
- The project goal and its parts (PROMPT.md) are fixed. Propose only concrete changes within scope, or state that none are needed in your area.
- Each proposed change: ID (e.g., STAT-1), severity (critical / major / minor), the exact DESIGN.md section, what to change (concrete wording/parameter/procedure), and why (evidence or reasoning). Max ~12 changes; prioritise.
- Also list "risks I could not resolve" briefly.
- Be specific and quantitative where possible (e.g., thresholds, sample-size estimates, compute estimates). You may run quick read-only checks on data (python, ≤4 cores, no GPU) to support a point, e.g., counting how many kinases would have ≥3 test positives.
- Write your review to docs/reviews/<round>_<role>.md and return a ≤300-word summary of your top changes.
