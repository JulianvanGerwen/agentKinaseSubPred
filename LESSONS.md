# LESSONS — process rules for running this project efficiently
(One-line rule + evidence. Pruned at each milestone.)

- Let reviewers run quick read-only data checks: in r1, 6/7 reviewers backed changes with measured numbers (paralog-contaminated hard negatives 51–54% for SRC family; 68% same-protein co-labelled pairs → leave-site-out leaks) that pure reading would have missed.
- Run literature, data profiling and data acquisition in parallel from the start: all three finished within ~20 min wall-clock of each other (evidence: M0→M1 in <1 h).
- Benchmark compute before planning around it: the design guessed co-regulation was the bottleneck; measurement showed 3.6 s on GPU, while LightGBM on 400k×160 took 15 min (r1 algorithms review).
- Key everything by UniProt accession, not symbol: Cantley "PDK1"=PDPK1 vs gene PDK1=PDHK1; 174/404 Cantley names ≠ gene symbols (r1 reviews).
