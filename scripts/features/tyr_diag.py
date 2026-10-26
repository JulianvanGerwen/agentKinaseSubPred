"""Diagnose the tyrosine-kinase E-hard gap (DEV CV only; labels via ksp.labels.load_dev).

Questions (owner: feat-tyr):
  D1  Is the Y macro-AUROC deficit driven by tiny-n kinases / few negatives, or is it uniform?
  D2  Is it label noise? Per-kinase E-hard AUROC of Cantley split by evidence class (in vitro vs in vivo-only).
  D3  Composition of the Y label set and of the Y E-hard pool: autophos-adjacent sites, kinase-protein
      substrates, activation-loop windows, near-duplicate 15-mers, sites shared by many kinases.
  D4  Does the Cantley Tyr atlas separate kinase substrate sets at all?  Cross-kinase AUROC matrix:
      score kinase k's matrix on (own eval positives) vs (eval positives of a reference kinase j).

Outputs -> results/features/tyr/
  diag_per_kinase.tsv, diag_evidence.tsv, diag_pool_composition.tsv, diag_crosskinase_auroc.tsv,
  diag_summary.json
Run: /workspace/.venv/bin/python -I scripts/features/tyr_diag.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import index, labels, store  # noqa: E402

OUT = ROOT / "results/features/tyr"
OUT.mkdir(parents=True, exist_ok=True)
SCORE = "motif/cantley_unprimed_pct_published"   # official B1 arm (A23)


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    return roc_auc_score(y, np.r_[pos, neg])


def main():
    ki, si = index.load_kinases(), index.load_sites()
    kpos = pd.Series(np.arange(len(ki)), index=ki.uniprot_acc)
    spos = pd.Series(np.arange(len(si)), index=si.site_id)
    ums = si.in_U_ms.to_numpy().astype(bool)
    resY = si.res.to_numpy() == "Y"
    M, _ = store.load_matrix(SCORE)

    allp = labels.load_dev("all")          # masks (A6: all rows, all sites)
    ev = labels.load_dev("eval")           # evaluation positives (U_ms)
    nm = labels.NegativeMasks(allp)
    par = labels.paralog_sets()
    kmeta = ki.set_index("uniprot_acc")

    # ---------------------------------------------------------------- D1 + D2: per kinase, by evidence
    rows = []
    evY = ev[ev.res.eq("Y")]
    pos_idx_of = {}
    for acc, g in evY.groupby("kinase_acc"):
        if acc not in kpos.index:
            continue
        kr = kpos[acc]
        hard = nm.hard_candidates(acc)
        hard = hard[ums[hard] & resY[hard]]                       # official universe: U_ms, class-matched
        neg = np.setdiff1d(hard, nm.excluded(acc))
        pi = spos.reindex(g.site_id).dropna().astype(int).to_numpy()
        neg = np.setdiff1d(neg, pi)
        pos_idx_of[acc] = pi
        r = M[kr]
        sn = r[neg]
        ec = g.set_index("site_id").evidence_native.reindex(si.site_id.to_numpy()[pi]).to_numpy()
        sp = r[pi]
        vit, viv = ec != "in_vivo_only", ec == "in_vivo_only"
        rows.append(dict(
            kinase_acc=acc, gene=kmeta.gene[acc], manning_family=kmeta.manning_family[acc],
            n_pos=len(pi), n_neg_hard=len(neg), n_paralogs=len(par.get(acc, ())),
            auroc_hard=auc(sp, sn),
            n_vitro=int(vit.sum()), auroc_vitro=auc(sp[vit], sn),
            n_vivo_only=int(viv.sum()), auroc_vivo_only=auc(sp[viv], sn)))
    per_k = pd.DataFrame(rows).sort_values("n_pos", ascending=False)
    per_k.to_csv(OUT / "diag_per_kinase.tsv", sep="\t", index=False)

    summ = {"score": SCORE, "n_Y_kinases_scored": int(len(per_k)),
            "macro_auroc_hard_all": float(per_k.auroc_hard.mean())}
    for mn in (1, 3, 5, 10, 20, 50):
        s = per_k[per_k.n_pos >= mn]
        summ[f"macro_auroc_hard_npos_ge_{mn}"] = float(s.auroc_hard.mean())
        summ[f"n_kinases_npos_ge_{mn}"] = int(len(s))
    # pooled (micro) AUROC over all Y (k, site) pairs, for contrast with the macro mean
    summ["median_n_neg_hard"] = float(per_k.n_neg_hard.median())
    summ["frac_kinases_npos_le_3"] = float((per_k.n_pos <= 3).mean())
    summ["frac_pos_in_kinases_npos_le_3"] = float(per_k.loc[per_k.n_pos <= 3, "n_pos"].sum() / per_k.n_pos.sum())

    pair = per_k[(per_k.n_vitro >= 5) & (per_k.n_vivo_only >= 5)]
    ev_tab = pair[["gene", "n_vitro", "auroc_vitro", "n_vivo_only", "auroc_vivo_only"]]
    ev_tab.to_csv(OUT / "diag_evidence.tsv", sep="\t", index=False)
    summ["evidence_paired_n_kinases"] = int(len(pair))
    summ["evidence_paired_macro_vitro"] = float(pair.auroc_vitro.mean())
    summ["evidence_paired_macro_vivo_only"] = float(pair.auroc_vivo_only.mean())
    summ["evidence_paired_delta"] = summ["evidence_paired_macro_vitro"] - summ["evidence_paired_macro_vivo_only"]

    # ---------------------------------------------------------------- D3: pool composition (Y vs S/T)
    kin_accs = set(ki.uniprot_acc)
    comp = []
    for res_name, mask in (("Y", resY), ("ST", ~resY)):
        evr = ev[ev.res.eq("Y")] if res_name == "Y" else ev[ev.res.isin(["S", "T"])]
        pool = np.unique(np.concatenate([nm.hard_idx]))
        pool = pool[ums[pool] & mask[pool]]
        sids = si.site_id.to_numpy()[pool]
        accs = si.acc.to_numpy()[pool]
        wins = si.window15.to_numpy()[pool]
        nlab = evr.groupby("site_id").kinase_acc.nunique().reindex(sids).fillna(0)
        wc = Counter(wins)
        comp.append(dict(
            residue=res_name, n_pool_sites=len(pool),
            frac_substrate_is_kinase=float(np.isin(accs, list(kin_accs)).mean()),
            mean_kinases_per_site=float(nlab.mean()), max_kinases_per_site=int(nlab.max()),
            frac_sites_ge4_kinases=float((nlab >= 4).mean()),
            frac_window_dup=float(np.mean([wc[w] > 1 for w in wins])),
            n_distinct_proteins=int(len(set(accs))),
            frac_top10_proteins=float(pd.Series(accs).value_counts().head(10).sum() / len(pool))))
    comp = pd.DataFrame(comp)
    comp.to_csv(OUT / "diag_pool_composition.tsv", sep="\t", index=False)

    # Y-specific: the SFK activation-loop window family, and sites labelled to many kinases
    aloop = evY[evY.site_id.map(si.set_index("site_id").window15).fillna("").str.contains("IEDNEYTAR")]
    summ["Y_aloop_IEDNEYTAR_pairs"] = int(len(aloop))
    summ["Y_aloop_IEDNEYTAR_sites"] = int(aloop.site_id.nunique())
    summ["Y_aloop_kinases_labelled"] = int(aloop.kinase_acc.nunique())
    shared = evY.groupby("site_id").kinase_acc.nunique().sort_values(ascending=False)
    summ["Y_top_shared_sites"] = {str(k): int(v) for k, v in shared.head(10).items()}
    summ["Y_frac_pairs_on_kinase_proteins"] = float(evY.acc.isin(kin_accs).mean())
    evST = ev[ev.res.isin(["S", "T"])]
    summ["ST_frac_pairs_on_kinase_proteins"] = float(evST.acc.isin(kin_accs).mean())

    # ---------------------------------------------------------------- D4: cross-kinase separation
    refs = [g for g in ("SRC", "FYN", "LCK", "LYN", "ABL1", "ABL2", "EGFR", "ERBB2", "JAK1", "JAK2",
                        "SYK", "ZAP70", "INSR", "IGF1R", "MET", "FGFR1", "PTK2", "BTK") if g in set(ki.gene)]
    g2a = dict(zip(ki.gene, ki.uniprot_acc))
    mat = pd.DataFrame(index=refs, columns=refs, dtype=float)
    for kg in refs:                      # scoring kinase (matrix)
        ka = g2a[kg]
        if ka not in pos_idx_of:
            continue
        r = M[kpos[ka]]
        for jg in refs:                  # substrate set
            ja = g2a[jg]
            if ja not in pos_idx_of:
                continue
            own, other = pos_idx_of[ka], pos_idx_of[ja]
            other = np.setdiff1d(other, own)     # sites shared by both are uninformative
            mat.loc[kg, jg] = auc(r[own], r[other])
    mat.to_csv(OUT / "diag_crosskinase_auroc.tsv", sep="\t")
    off = mat.to_numpy()[~np.eye(len(refs), dtype=bool)]
    summ["crosskinase_mean_offdiag_auroc"] = float(np.nanmean(off))
    summ["crosskinase_frac_offdiag_ge_0.6"] = float(np.nanmean(off >= 0.6))
    # within- vs between-family contrast
    famof = {g: kmeta.manning_family[g2a[g]] for g in refs}
    wi, be = [], []
    for i, a in enumerate(refs):
        for j, b in enumerate(refs):
            if i == j or np.isnan(mat.loc[a, b]):
                continue
            (wi if famof[a] == famof[b] else be).append(mat.loc[a, b])
    summ["crosskinase_mean_within_family"] = float(np.mean(wi)) if wi else None
    summ["crosskinase_mean_between_family"] = float(np.mean(be)) if be else None

    (OUT / "diag_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))
    print("\n--- pool composition ---\n", comp.round(3).to_string(index=False))
    print("\n--- cross-kinase AUROC (row = matrix, col = substrate set) ---\n", mat.round(2).to_string())


if __name__ == "__main__":
    main()
