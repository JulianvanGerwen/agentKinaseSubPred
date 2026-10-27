"""Downstream kinase-activity-inference evaluation (benchmarKIN; DESIGN §6.7, §16 A21). See docs/downstream.md.

Modules
  bench_data  benchmarKIN perturbation (Hernandez-Armenta, Hijazi) and CPTAC tumour benchmarks -> site_id / kinase acc
  libraries   kinase -> substrate libraries (PSP all / DEV, Cantley top-N, shuffled, predictions top-N, union / minus)
  inference   kinase activity methods (RoKAI z-score, KSEA, decoupler-style ULM and wmean)
  evaluate    benchmarKIN metrics (scaled rank, PHit(k), subsampled AUROC; tumour AUROC with 80% subsampling)
  api         run_benchmarks(libraries, ...) -> tidy result tables (the entry point for later model libraries)
"""
