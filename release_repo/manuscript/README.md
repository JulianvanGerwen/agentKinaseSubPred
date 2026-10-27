# Manuscript

> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

This directory holds the manuscript prepared for *PLOS Computational Biology*. It is written in Markdown.
Converting it to Word or LaTeX for submission is left to the publishing researcher. Authors are a placeholder.

| File | What it is |
|---|---|
| `manuscript.md` | main text, tables and figure legends |
| `supplement.md` | S1 Text: approaches that were tried and did not work, and other negative results |
| `references.md` | the reference list, with the source each entry was checked against; unverified details are marked `[verify]` |
| `NUMBER_MAP.tsv`, `SUPP_NUMBER_MAP.tsv` | every number in `manuscript.md` / `supplement.md`, with the project file it comes from |
| `figures/` | Fig1–Fig5 as PNG (300 dpi) and PDF, and `make_figures.py`, which builds them |
| `INDEPENDENT_REVIEW_OF_DRAFT.md` | an independent agent's review of an earlier draft; its line numbers refer to that draft. Kept for transparency |
| `CLARITY_CHANGES.md` | log of the final clarity edit of `manuscript.md` and `supplement.md` |

Paths in these files (`results/...`, `cache/...`, `REPORT.md`) refer to the original project directory. The metric
tables they cite are in `../results_summary/` (`results/test/721dc1977366/` → `results_summary/test_721dc1977366/`)
and the release files in `../predictions/`; `REPORT.md` is `../docs/REPORT_internal.md`. Label-level files are not
in this repository.

Rebuild the figures from the repository root (needs numpy, pandas and matplotlib):

```
python manuscript/figures/make_figures.py            # overwrites the figures here
python manuscript/figures/make_figures.py /some/dir  # writes them elsewhere
```

Where the manuscript and the internal report differ, the manuscript is the later, reviewed version; the numbers in both
come from REPORT.md v3.
