> This work was produced as an exercise in testing whether AI agents can run a bioinformatics research project. A lead AI agent (the 'PI', Claude Opus 5.5) directed sub-agents (Claude Opus 5.5 and Sonnet 5.5) under a human-written brief. The human researcher supplied the brief, the data and PhosphoSitePlus access; they did not write the code, run the analyses, or write the text. All results should be checked before use.

# Licence of the released predictions

The files in `predictions/` (`kinase_top200.tsv.gz`, `site_top5.tsv.gz`, `coverage.tsv`, `RELEASE_MANIFEST.json`)
and the metric tables in `results_summary/` were produced by models trained and evaluated on kinase–substrate labels
from **PhosphoSitePlus®** (PSP; Cell Signaling Technology, Inc.; www.phosphosite.org). PSP data are distributed under
the **Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported licence (CC BY-NC-SA 3.0)**,
https://creativecommons.org/licenses/by-nc-sa/3.0/. Because these files are derived from PSP, they are released under
the same licence:

- **Attribution.** Cite PhosphoSitePlus (www.phosphosite.org) and its publications. The project used Hornbeck et al.
  2012 *Nucleic Acids Res.* doi:10.1093/nar/gkr1122, Hornbeck et al. 2015 *Nucleic Acids Res.*
  doi:10.1093/nar/gku1267 and Hornbeck et al. 2019 *Nucleic Acids Res.* doi:10.1093/nar/gky1159 (references as listed in `docs/lit/02_ks_annotation_sources.md`). Also cite
  this work (`CITATION.cff`).
- **Non-commercial.** Do not use these files for commercial purposes. For commercial use of PSP-derived material,
  contact Cell Signaling Technology for a PSP licence.
- **Share-alike.** If you remix, transform or build on these files, distribute your contributions under the same
  licence.

Other non-commercial inputs also shaped these predictions: the Cantley kinase-library matrices (kinase-library 1.8.0,
CC BY-NC-SA 3.0) and iPTMnet (CC BY-NC-SA 4.0; used only to build a list of excluded studies). Their terms are
compatible with the licence above. See `NOTICE` and `DATA_SOURCES.md`.

No PSP data are included in this repository. The predictions are scores for kinase–site pairs. They do not mark
which pairs are PSP positives, and they do not contain PSP sequences, evidence fields or identifiers.

The code in this repository is under the MIT licence (`LICENSE`). The MIT licence does not apply to the files listed
above.

The predictions are provided "as is", without warranty of any kind. They were produced by AI agents and have not been
validated experimentally; see the limitations in `README.md` and `predictions/README.md`.
