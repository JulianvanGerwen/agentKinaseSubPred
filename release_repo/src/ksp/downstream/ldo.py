"""Leave-drug-out (A21) stratification of the perturbation benchmark.

No leave-drug-out *model* exists (it would need F6a/F6b rebuilt with `exclude_drugs` and an R2 retrain). The harness
instead stratifies the benchmark experiments by whether the clean model's phospho features could have seen the
same drug:
  * benchmark drug identities: `phospho_inhib.find_drugs` on the experiment description (Hernandez), drug name
    (Hijazi) or treatment (tyrosine), using the same union-find drug identity (Klaeger / Reinecke / Davis synonyms +
    results/features/phospho/drug_aliases.tsv) as the F6b feature builder. Class tokens ('MEKi', 'p38i') count only
    when no named compound is found, and then only by class name.
  * "seen": some compendium condition OUTSIDE the clean_bench PMIDs (i.e. kept in the clean variant's F6a/F6b) names
    the same drug identity, in Condition or Control text (f6b_conditions.tsv `drugs`). Time filters are ignored
    (F6a uses all conditions), so "seen" is conservative.
  * strata: drug_unseen (all named drugs unseen), drug_seen (>= 1 seen), no_drug (stimuli, genetic, cell cycle,
    DNA damage; leave-drug-out does not apply).
Evaluating on drug_unseen experiments approximates leave-drug-out: for these drugs the clean model had no
phosphoproteomic response data at all.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache

import pandas as pd

from .bench_data import ROOT

PH = ROOT / "results" / "features" / "phospho"


@lru_cache(maxsize=1)
def _drug_tools():
    from ksp.features import phospho_inhib as pi
    from ksp.kinases import load_resolver
    _, uf, _ = pi.build_affinity(load_resolver(), verbose=False)
    return pi, uf, pi.read_aliases()


def _ids_from_text(text: str) -> tuple[set, set]:
    pi, uf, al = _drug_tools()
    hits = pi.find_drugs(str(text), set(uf.p.keys()), uf, al)
    named = {d for d, _ in hits if not d.startswith("class:")}
    cls = {d for d, _ in hits if d.startswith("class:")}
    return named, cls


def compendium_drug_ids(exclude_pmids) -> tuple[set, set]:
    pi, uf, _ = _drug_tools()
    t = pd.read_csv(PH / "f6b_conditions.tsv", sep="\t", dtype=str, keep_default_na=False)
    t = t[~t.PMID.isin({str(p) for p in exclude_pmids})]
    named, cls = set(), set()
    for txt in t.drugs:
        for tok in re.findall(r"(?:cond|ctrl):([^@|]+)@", txt):
            k = pi.norm(tok)
            if k in uf.p:
                named.add(uf.find(k))
            else:
                cls.add(f"class:{k}")
    return named, cls


def _compendium_text(exclude_pmids) -> str:
    """All normalised Condition/Control strings of compendium conditions outside `exclude_pmids` (text fallback)."""
    pi, _, _ = _drug_tools()
    c = pd.read_csv(ROOT / "data" / "01_phosphoproteomics_conditions.csv", dtype=str, keep_default_na=False)
    c = c[~c.PMID.isin({str(p) for p in exclude_pmids})]
    return "|".join(pi.norm(a) + "|" + pi.norm(b) for a, b in zip(c.Condition, c.Control))


def drug_strata(experiments: pd.DataFrame) -> pd.DataFrame:
    """experiments: bench.experiments. Returns experiment, dataset, description, drugs, drug_classes, stratum.
    A drug also counts as seen if any synonym (>= 5 characters) occurs as a substring of a kept compendium
    Condition/Control text (fallback for conditions the F6b parser did not resolve)."""
    pi, uf, _ = _drug_tools()
    excl = json.loads((PH / "exclude_lists.json").read_text())["clean_bench"]["exclude_pmids"]
    seen_named, seen_cls = compendium_drug_ids(excl)
    text = _compendium_text(excl)
    groups: dict = {}
    for k in list(uf.p):
        groups.setdefault(uf.find(k), set()).add(k)
    rows = []
    for e, d, desc in zip(experiments.experiment, experiments.dataset, experiments.description):
        named, cls = _ids_from_text(desc)
        ids = named or cls
        seen = (named & seen_named) if named else (cls & seen_cls)
        for i in ids:                                       # text fallback over all synonyms
            syn = groups.get(i, set()) | ({i[6:]} if i.startswith("class:") else {i})
            if any(len(s) >= 5 and s in text for s in syn):
                seen = seen | {i}
        stratum = ("drug_seen" if seen else "drug_unseen") if ids else "no_drug"
        rows.append((e, d, desc, ";".join(sorted(named)), ";".join(sorted(cls)), ";".join(sorted(seen)), stratum))
    return pd.DataFrame(rows, columns=["experiment", "dataset", "description", "drug_ids", "drug_classes",
                                       "seen_in_clean_compendium", "stratum"])
