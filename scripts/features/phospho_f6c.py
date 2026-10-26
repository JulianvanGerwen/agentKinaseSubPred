#!/usr/bin/env python
"""F6c CPTAC kinase-state correlation (not label-derived), DESIGN §7 F6c.

Outputs (kinase x site float16 via ksp.store):
  cache/features/phospho/f6c_prot{_v}.npy   mean Fisher z over cohorts: residualised site vs kinase protein abundance
  cache/features/phospho/f6c_aloop{_v}.npy  mean Fisher z over cohorts: residualised site vs kinase A-loop phospho
results/features/phospho/f6c_aloop_sites.tsv (A-loop range + observed sites per kinase), f6c_report{_v}.json

Usage: CUDA_VISIBLE_DEVICES=1 PYTHONPATH=src .venv/bin/python -I scripts/features/phospho_f6c.py
       [--variant v --exclude-cptac BRCA,LUAD | --exclude-cptac all] [--min-n 20]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import store  # noqa: E402
from ksp.ids import CanonicalStore  # noqa: E402
from ksp.features import phospho as ph  # noqa: E402
from ksp.features import phospho_cptac as pc  # noqa: E402


def aloop_table(kin: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    st = CanonicalStore.from_uniprot(all_fasta=False)
    by_acc = {a: g for a, g in sites.groupby("acc")}
    recs = []
    for i, r in enumerate(kin.itertuples()):
        seq = st.seqs.get(r.uniprot_acc, "")
        al = pc.activation_loop(seq, r.kinase_domain_start, r.kinase_domain_end) \
            if pd.notna(r.kinase_domain_start) else None
        rec = dict(kinase_idx=i, kinase=r.gene, uniprot_acc=r.uniprot_acc, aloop_start=None, aloop_end=None,
                   motif="", aloop_seq="", site_idx="", site_ids="")
        if al:
            s, e, how = al
            rec.update(aloop_start=s, aloop_end=e, motif=how, aloop_seq=seq[s - 1:e])
            g = by_acc.get(r.uniprot_acc)
            if g is not None:
                g = g[(g.pos >= s) & (g.pos <= e)]
                rec.update(site_idx=";".join(map(str, g.index)), site_ids=";".join(g.site_id))
        recs.append(rec)
    return pd.DataFrame(recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="")
    ap.add_argument("--exclude-cptac", default="", help="comma list of cohorts, or 'all'")
    ap.add_argument("--min-n", type=int, default=20)
    a = ap.parse_args()
    t0 = time.time()
    sfx = f"_{a.variant}" if a.variant else ""
    excl = True if a.exclude_cptac == "all" else [x for x in a.exclude_cptac.split(",") if x]
    kin, sites = ph.kinase_index(), ph.site_index()
    blocks = ph.load_blocks(kinds=("cptac",), exclude_cptac=excl)
    kp = np.load(ph.CACHE / "cptac_kinase_protein.npz")
    al = aloop_table(kin, sites)
    al.to_csv(ph.RESULTS / "f6c_aloop_sites.tsv", sep="\t", index=False)
    prof_prot, prof_al = {}, {}
    n_al = {}
    for b in blocks:
        prof_prot[b.cohort] = (kp[f"{b.cohort}_kidx"], kp[f"{b.cohort}_X"])
        rpos = {int(s): j for j, s in enumerate(b.rows)}
        ki, P = [], []
        for r in al.itertuples():
            if not r.site_idx:
                continue
            js = [rpos[int(s)] for s in r.site_idx.split(";") if int(s) in rpos]
            if not js:
                continue
            Z = pc.zrows(b.X[js])
            with np.errstate(all="ignore"):
                p = np.nanmean(Z, 0)
            if np.isfinite(p).sum() >= a.min_n:
                ki.append(r.kinase_idx)
                P.append(p.astype(np.float32))
        prof_al[b.cohort] = (np.array(ki, np.int64), np.stack(P) if P else np.zeros((0, b.n_cols), np.float32))
        n_al[b.cohort] = len(ki)
    own = {i: sites.index[sites.acc == a_].to_numpy() for i, a_ in enumerate(kin.uniprot_acc)}
    own = {k: v for k, v in own.items() if len(v)}
    rep = dict(variant=a.variant, exclude_cptac=a.exclude_cptac, cohorts=[b.cohort for b in blocks], min_n=a.min_n,
               n_kinases_with_aloop_range=int(al.aloop_start.notna().sum()),
               n_kinases_with_aloop_sites_in_universe=int((al.site_idx != "").sum()),
               n_kinases_with_aloop_profile_per_cohort=n_al,
               n_kinases_with_protein_per_cohort={c: int(len(v[0])) for c, v in prof_prot.items()})
    for name, prof in (("f6c_prot", prof_prot), ("f6c_aloop", prof_al)):
        Fz, nc = pc.f6c_correlations(blocks, prof, len(kin), len(sites), own, min_n=a.min_n)
        store.save_matrix(f"phospho/{name}{sfx}", Fz.astype(np.float16), "kinase", "site",
                          meta=dict(variant=a.variant, exclude_cptac=a.exclude_cptac, min_n=a.min_n,
                                    label_derived=False))
        rep[name] = dict(n_finite=int(np.isfinite(Fz).sum()), n_kinases=int(np.isfinite(Fz).any(1).sum()),
                         n_sites=int(np.isfinite(Fz).any(0).sum()),
                         mean_abs_z=float(np.nanmean(np.abs(Fz))))
        print(name, rep[name], f"{time.time()-t0:.0f}s", flush=True)
    rep["seconds"] = round(time.time() - t0, 1)
    (ph.RESULTS / f"f6c_report{sfx}.json").write_text(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()
