#!/usr/bin/env python
"""Positive control for the sequence-based ortholog site mapping used as a fallback in profile_psp.py:
apply map_sites_by_sequence() to non-human KS sites that ARE mapped by SITE_GRP_ID (to a human site of
the Phosphorylation_site_dataset / KS dataset) and measure agreement with the SITE_GRP_ID mapping.
Also reports the identity distribution of fallback-only mappings.

Output: results/data_profile/psp_sequence_mapping_validation.tsv
Usage: python scripts/verify_psp_sequence_mapping.py  (needs cache/psp_human_mapped_rows.parquet + external/uniprot)
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp import io_psp, uniprot  # noqa: E402


def main(min_identity=11):
    m = pd.read_parquet("cache/psp_human_mapped_rows.parquet")
    fa = uniprot.read_fasta()
    ann = uniprot.read_annotations()
    acc2gene = dict(zip(ann.Entry, ann["Gene Names (primary)"].str.upper()))
    canon = fa[(fa.isoform == "") & fa.db.eq("sp")].copy()
    canon["gene_uc"] = canon.acc.map(acc2gene).fillna(canon.gene.str.upper())
    sg = m[m.map_method.str.startswith("site_grp_id") & m.h_sub_acc_full.notna()]
    # one SGID target per non-human site, canonical human targets only (sequence mapping uses canonical seqs)
    sg = sg[~sg.h_sub_acc_full.str.contains("-")]
    tgt = sg.groupby(["sub_acc_full", "position"]).apply(
        lambda g: set(zip(g.h_sub_acc_full, g.h_position.astype(int))), include_groups=False).rename("sgid_targets")
    sites = sg.drop_duplicates(["sub_acc_full", "position"])[["sub_organism", "sub_acc_full", "position", "residue",
                                                                "window", "sub_gene_uc"]]
    res = io_psp.map_sites_by_sequence(sites, canon[["gene_uc", "acc_full", "acc", "seq"]], min_identity)
    res = res.join(tgt, on=["sub_acc_full", "position"])
    res["hit"] = res.seq_h_sub_acc_full.notna()
    res["agree"] = [(a, int(p)) in t if pd.notna(a) else False
                    for a, p, t in zip(res.seq_h_sub_acc_full, res.seq_h_position.fillna(-1), res.sgid_targets)]
    rows = [dict(check="SGID-mapped non-human sites (canonical human target)", n=len(res),
                 n_seq_hit=int(res.hit.sum()), n_agree=int(res.agree.sum()),
                 precision=round(res.agree.sum() / max(res.hit.sum(), 1), 4),
                 recall=round(res.agree.mean(), 4))]
    for o in ["mouse", "rat"]:
        r = res[res.sub_organism == o]
        rows.append(dict(check=f"  {o}", n=len(r), n_seq_hit=int(r.hit.sum()), n_agree=int(r.agree.sum()),
                         precision=round(r.agree.sum() / max(r.hit.sum(), 1), 4), recall=round(r.agree.mean(), 4)))
    fb = m[m.map_method.eq("sequence")].drop_duplicates(["sub_acc_full", "position"])
    for k, v in fb.seq_identity.value_counts().sort_index().items():
        rows.append(dict(check=f"fallback-only sites with identity {int(k)}/15", n=int(v)))
    for k, v in res[res.hit].seq_identity.value_counts().sort_index().items():
        rows.append(dict(check=f"control sites with identity {int(k)}/15 (agree)", n=int(v),
                         n_agree=int(res[res.hit & (res.seq_identity == k)].agree.sum())))
    out = pd.DataFrame(rows)
    out.to_csv("results/data_profile/psp_sequence_mapping_validation.tsv", sep="\t", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
