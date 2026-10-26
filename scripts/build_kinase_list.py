#!/usr/bin/env python
"""Build results/labels/kinase_list.tsv and results/labels/kinase_synonyms.tsv (+ kinase_list_summary.tsv).

Rules: src/ksp/kinases.py docstring and docs/labels_universe.md. Manual decisions: config/kinase_curation.tsv.
Usage: /workspace/.venv/bin/python -I scripts/build_kinase_list.py [--exclude-pseudokinases] [--exclude-uncertain]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ksp import io_psp, kinases  # noqa: E402


def psp_kinase_entries(hg: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PSP kinase genes (upper-case human symbol after io_psp.map_kinases_to_human) with a human accession:
    human rows' KIN_ACC_ID (most frequent), else HGNC approved/previous/alias symbol -> UniProt."""
    d = io_psp.filter_psp(io_psp.load_psp(), sty_only=True)
    d = io_psp.map_kinases_to_human(d)
    d = d[d.kin_human_mapped | d.kin_gene_uc.ne("")]
    d = d.assign(g=d.kin_human_gene.fillna(d.kin_gene_uc).str.upper())
    hum = d[d.kin_organism.eq("human")]
    acc = hum.groupby("g").kin_acc.agg(lambda s: s.value_counts().index[0])
    hsym = {}
    for col in ["symbol", "prev_symbol", "alias_symbol"]:
        for r in hg.itertuples(index=False):
            for s in str(getattr(r, col)).split("|"):
                if s and r.uniprot_ids:
                    hsym.setdefault(s.upper(), r.uniprot_ids.split("|")[0])
    out = []
    for g, n in d.g.value_counts().items():
        out.append((g, acc.get(g, hsym.get(g, "")), int(n)))
    names = hum[["kin_name", "g"]].drop_duplicates().rename(columns={"kin_name": "name", "g": "gene"})
    return pd.DataFrame(out, columns=["gene", "acc", "n_rows"]), names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclude-pseudokinases", action="store_true")
    ap.add_argument("--exclude-uncertain", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "results" / "labels"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    hg = kinases.load_hgnc()
    psp_k, psp_names = psp_kinase_entries(hg)
    kl = kinases.build_kinase_list(psp_k, a.exclude_pseudokinases, a.exclude_uncertain)
    # PSP row counts per kinase-list gene (documentation only)
    acc2n = psp_k[psp_k.acc.ne("")].groupby("acc").n_rows.sum()
    kl["psp_ks_rows"] = kl.uniprot_acc.map(lambda a: acc2n.get(a, 0) if isinstance(a, str) and a else 0).astype(int)
    psp_extra_n = dict(zip(psp_k.gene, psp_k.n_rows))
    m = kl.psp_ks_rows.eq(0) & kl.in_psp
    kl.loc[m, "psp_ks_rows"] = kl.loc[m, "gene"].str.upper().map(psp_extra_n).fillna(0).astype(int)
    kl.to_csv(out / "kinase_list.tsv", sep="\t", index=False)

    # psp kinase names: map their PSP gene to the kinase-list gene through accession
    acc_of = dict(zip(psp_k.gene, psp_k.acc))
    a2g = dict(zip(kl.uniprot_acc, kl.gene))
    psp_names = psp_names.assign(gene=psp_names.gene.map(lambda g: a2g.get(acc_of.get(g, ""), g)))
    syn = kinases.build_synonyms(kl, psp_names)
    syn.to_csv(out / "kinase_synonyms.tsv", sep="\t", index=False)

    # summary + verification of joins with other sources
    res = kinases.KinaseResolver(kl.assign(include=kl.include.astype(bool)), syn)
    rows = [("candidates", len(kl)), ("include", int(kl.include.sum())),
            ("include_pseudokinase", int((kl.include & kl.is_pseudokinase).sum())),
            ("include_activity_uncertain", int((kl.include & kl.activity_uncertain).sum()))]
    for c, n in kl[~kl.include].exclusion_category.value_counts().items():
        rows.append((f"excluded:{c}", int(n)))
    for c, n in kl[kl.include].residue_class.value_counts().items():
        rows.append((f"include_residue_class:{c}", int(n)))
    for c, n in kl[kl.include].manning_group.replace("", "(none)").value_counts().items():
        rows.append((f"include_group:{c}", int(n)))
    rows.append(("include_with_kinase_domain", int(kl[kl.include].kinase_domain_start.notna().sum())))
    rows.append(("synonyms_rows", len(syn)))
    rows.append(("synonyms_use_for_join", int(syn.use_for_join.sum())))
    # kinhib (Hijazi) legacy symbols
    kh = pd.read_csv(ROOT / "data" / "kinhib_inhib_comb__1.csv", usecols=[0]).iloc[:, 0].astype(str)
    r = [res.resolve(n) for n in kh]
    rows.append(("kinhib_rows", len(kh)))
    rows.append(("kinhib_resolved", sum(g is not None for g, _ in r)))
    rows.append(("kinhib_resolved_included", sum(g is not None and res.included.get(g, False) for g, _ in r)))
    unres = sorted(n for n, (g, _) in zip(kh, r) if g is None)
    rows.append(("kinhib_unresolved_names", " ".join(unres)))
    cn = kinases.load_cantley_names()
    rc = [res.resolve(m, acc=a)[0] for m, a in zip(cn.matrix_name, cn.acc)]
    rows.append(("cantley_matrices", len(cn)))
    rows.append(("cantley_resolved", sum(g is not None for g in rc)))
    rows.append(("cantley_unresolved", " ".join(sorted(m for m, g in zip(cn.matrix_name, rc) if g is None))))
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(out / "kinase_list_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows, columns=["metric", "value"]).to_string(index=False))


if __name__ == "__main__":
    main()
