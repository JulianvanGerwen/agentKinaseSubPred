# Design review brief (round N)
You are an independent reviewer of /workspace/DESIGN.md for a kinase–substrate prediction project.
Read: /workspace/PROMPT.md (fixed goal & rules), /workspace/DESIGN.md, /workspace/LITERATURE.md, and as needed docs/lit/*.md, docs/data_profile.md, docs/data_acquisition.md.
Rules for reviewers:
- The project goal and its parts (PROMPT.md) are fixed. Propose only concrete changes within scope, or state that none are needed in your area.
- Each proposed change: ID (e.g., STAT-1), severity (critical / major / minor), the exact DESIGN.md section, what to change (concrete wording/parameter/procedure), and why (evidence or reasoning). Max ~12 changes; prioritise.
- Also list "risks I could not resolve" briefly.
- Be specific and quantitative where possible (e.g., thresholds, sample-size estimates, compute estimates). You may run quick read-only checks on data (/workspace/.venv/bin/python, ≤4 cores, no GPU) to support a point, e.g., counting how many kinases would have ≥3 test positives.
- Write your review to /workspace/docs/reviews/<round>_<role>.md and return a ≤300-word summary of your top changes.
