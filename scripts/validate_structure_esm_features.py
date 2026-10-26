#!/usr/bin/env python
"""Sanity validation of the F0 site-intrinsic features (cache/structure/site_structure.parquet, cache/esm2_650m/).

1. Known intrinsically disordered proteins vs well-folded controls: pLDDT, IDR flag, RSA.
2. PSP Phosphorylation_site_dataset (human, mapped to canonical by exact 15-mer window match at the stated
   canonical position; isoform-specific rows skipped) vs non-phospho S/T/Y of the same proteins: IDR fraction,
   pLDDT, RSA, secondary structure; odds ratios.
3. ESM-2 linear probe: logistic regression on 650M S/T/Y embeddings, phosphosite vs not, protein-grouped 5-fold CV
   (sanity only; NOT a model result). Also ESM -> structural IDR flag probe.
4. Consistency: kinases.parquet emb_full vs protein_mean.npy cosine.
Writes results/features/structure_esm_validation.json and results/features/phospho_vs_nonphospho_structure.tsv.
Usage: python -I scripts/validate_structure_esm_features.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ksp.features.esm_embed import canonical_swissprot, load_sty  # noqa: E402

OUT = ROOT / "results/features"
IDP = {"P37840": "SNCA", "P10636": "MAPT", "P46527": "CDKN1B", "P38936": "CDKN1A", "P10451": "SPP1",
       "P17096": "HMGA1", "P02686": "MBP", "O60829": "PAGE4", "P06454": "PTMA", "P05814": "CSN2",
       "P49023": "PXN?", "Q16186": "ADRM1?"}
IDP = {k: v for k, v in IDP.items() if not v.endswith("?")}
FOLDED = {"P61626": "LYZ", "P02144": "MB", "P68871": "HBB", "P00918": "CA2", "P02766": "TTR", "P62937": "PPIA",
          "P00441": "SOD1", "P62988": "UBB?", "P69905": "HBA1", "P00698": "LYZ_chicken?", "P01308": "INS",
          "P62158": "CALM?", "P0DP23": "CALM1"}
FOLDED = {k: v for k, v in FOLDED.items() if not v.endswith("?")}


def psp_sites(seqs):
    s = pd.read_parquet(ROOT / "cache/psp_human_sites.parquet")
    s = s[(s.mod_type == "p") & s.is_sty & (s.isoform == "") & s.acc.isin(seqs.keys())].copy()

    def win_ok(r):
        q = seqs[r.acc]
        p = r.position - 1
        w = "".join(q[i] if 0 <= i < len(q) else "_" for i in range(p - 7, p + 8))
        return w == r.window
    s["ok"] = [win_ok(r) for r in s.itertuples()]
    n_all = len(s)
    s = s[s.ok]
    s["site_id"] = s.acc + "_" + s.residue + s.position.astype(str)
    return s.drop_duplicates("site_id"), n_all


def odds(a, b):
    """a, b boolean arrays (feature true in phospho / non-phospho)."""
    pa, pb = a.mean(), b.mean()
    return float(pa / (1 - pa) / (pb / (1 - pb)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seqs = dict(canonical_swissprot()[["acc", "seq"]].values)
    st = pd.read_parquet(ROOT / "cache/structure/site_structure.parquet")
    meta = json.loads((ROOT / "cache/structure/site_structure.meta.json").read_text())
    rep = {"structure_meta": {k: v for k, v in meta.items() if not k.endswith("accs") and k != "other_fail"}}
    rep["coverage"] = dict(n_sty=len(st), n_with_features=int(st.plddt.notna().sum()),
                           frac=float(st.plddt.notna().mean()),
                           af_status=st.af_status.value_counts().to_dict())

    # 1. IDP panel
    res = pd.read_parquet(ROOT / "cache/structure/residue_structure.parquet",
                          columns=["acc", "pos", "plddt", "rsa", "idr", "disorder_plddt", "ss3"])
    rows = []
    for grp, d in (("IDP", IDP), ("folded", FOLDED)):
        for acc, g in d.items():
            r = res[res.acc == acc]
            if len(r):
                rows.append(dict(group=grp, acc=acc, gene=g, L=len(r), plddt=r.plddt.mean(), idr=r.idr.mean(),
                                 disorder_plddt=r.disorder_plddt.mean(), rsa=r.rsa.mean(),
                                 helix=(r.ss3 == "H").mean(), strand=(r.ss3 == "E").mean()))
    panel = pd.DataFrame(rows)
    panel.to_csv(OUT / "idp_panel_structure.tsv", sep="\t", index=False, float_format="%.3f")
    rep["idp_panel"] = panel.groupby("group")[["plddt", "idr", "disorder_plddt", "rsa"]].mean().round(3).to_dict()
    rep["proteome_residue_means"] = dict(plddt=float(res.plddt.mean()), idr=float(res.idr.mean()),
                                         disorder_plddt=float(res.disorder_plddt.mean()))
    del res

    # 2. phospho vs non-phospho S/T/Y
    ps, n_psp = psp_sites(seqs)
    rep["psp"] = dict(n_canonical_rows=n_psp, n_window_matched=len(ps))
    st["phos"] = st.site_id.isin(set(ps.site_id))
    lt = set(ps.loc[ps.lt_lit > 0, "site_id"])
    hi = set(ps.loc[ps.ms_total >= 5, "site_id"])
    ks = set(ps.loc[ps.ks_any, "site_id"])
    d = st[st.acc.isin(set(ps.acc)) & st.plddt.notna()].copy()
    d["cls"] = np.where(d.site_id.isin(ks), "KS-annotated (PSP)",
                        np.where(d.site_id.isin(hi), "phospho MS>=5",
                                 np.where(d.phos, "phospho other", "non-phospho")))
    tab = []
    for (res_, cls), g in d.groupby(["res", "cls"]):
        tab.append(dict(res=res_, cls=cls, n=len(g), plddt=g.plddt.mean(), plddt_s10=g.plddt_s10.mean(),
                        rsa=g.rsa.mean(), idr=g.idr.mean(), disorder_plddt=g.disorder_plddt.mean(),
                        high_acc=g.high_acc.mean(), helix=(g.ss3 == "H").mean(), strand=(g.ss3 == "E").mean(),
                        pse_12_70=g.pse_12_70.mean()))
    for cls, g in d.groupby("cls"):
        tab.append(dict(res="STY", cls=cls, n=len(g), plddt=g.plddt.mean(), plddt_s10=g.plddt_s10.mean(),
                        rsa=g.rsa.mean(), idr=g.idr.mean(), disorder_plddt=g.disorder_plddt.mean(),
                        high_acc=g.high_acc.mean(), helix=(g.ss3 == "H").mean(), strand=(g.ss3 == "E").mean(),
                        pse_12_70=g.pse_12_70.mean()))
    tab = pd.DataFrame(tab)
    tab.to_csv(OUT / "phospho_vs_nonphospho_structure.tsv", sep="\t", index=False, float_format="%.4f")
    P, N = d[d.phos], d[~d.phos]
    from sklearn.metrics import roc_auc_score
    rep["phospho_vs_non"] = dict(
        n_phos=len(P), n_non=len(N), n_lt=int(d.site_id.isin(lt).sum()),
        or_idr=odds(P.idr == 1, N.idr == 1), or_disorder_plddt=odds(P.disorder_plddt == 1, N.disorder_plddt == 1),
        or_idr_by_res={r: odds(P[P.res == r].idr == 1, N[N.res == r].idr == 1) for r in "STY"},
        or_idr_LT=odds(d[d.site_id.isin(lt)].idr == 1, N.idr == 1),
        auc_neg_plddt=float(roc_auc_score(d.phos, -d.plddt)), auc_rsa=float(roc_auc_score(d.phos, d.rsa)),
        auc_neg_plddt_s10=float(roc_auc_score(d.phos, -d.plddt_s10)),
        auc_neg_pse_24_180_s10=float(roc_auc_score(d.phos, -d.pse_24_180_s10)),
        frac_idr_phos=float(P.idr.mean()), frac_idr_non=float(N.idr.mean()),
        # ordered regions only: phosphosites in structured regions are more exposed?
        rsa_in_ordered_phos=float(P[P.idr == 0].rsa.mean()), rsa_in_ordered_non=float(N[N.idr == 0].rsa.mean()),
        auc_rsa_in_ordered=float(roc_auc_score(d[d.idr == 0].phos, d[d.idr == 0].rsa)))

    # 3. ESM linear probe (sanity only)
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    rng = np.random.default_rng(0)
    e = d.sample(n=min(120000, len(d)), random_state=0)
    idx, M = load_sty(ROOT / "cache/esm2_650m", e.site_id.values)
    X = M.astype(np.float32)
    y = e.phos.values.astype(int)
    yi = e.idr.values.astype(int)
    groups = e.acc.values
    probe = {}
    for name, yy in (("phospho", y), ("idr", yi)):
        pred = np.zeros(len(yy))
        for tr, te in GroupKFold(5).split(X, yy, groups):
            clf = LogisticRegression(C=0.05, max_iter=300)
            clf.fit(X[tr], yy[tr])
            pred[te] = clf.decision_function(X[te])
        probe[name] = float(roc_auc_score(yy, pred))
    # structure-only logistic baseline for phospho, same folds
    F = e[["plddt", "plddt_s10", "rsa", "pse_12_70", "pse_24_180_s10", "idr"]].fillna(0).values
    F = (F - F.mean(0)) / (F.std(0) + 1e-9)
    pred = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(F, y, groups):
        pred[te] = LogisticRegression(max_iter=300).fit(F[tr], y[tr]).decision_function(F[te])
    probe["phospho_structure_only"] = float(roc_auc_score(y, pred))
    probe["n"] = len(e)
    probe["note"] = "protein-grouped 5-fold CV on a 120k random sample of S/T/Y of PSP-phosphoproteins; sanity only"
    rep["esm_probe"] = probe

    # 4. consistency
    km = pd.read_parquet(ROOT / "cache/esm2_650m/kinases.parquet").drop_duplicates("acc")
    pi = pd.read_parquet(ROOT / "cache/esm2_650m/protein_mean_index.parquet").set_index("acc")
    V = np.load(ROOT / "cache/esm2_650m/protein_mean.npy")
    a = np.stack(km.emb_full.values)
    b = V[pi.loc[km.acc, "row"].values]
    cos = (a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    rep["kinase_vs_proteome_mean_cos_min"] = float(cos.min())
    (OUT / "structure_esm_validation.json").write_text(json.dumps(rep, indent=1, default=float))
    print(json.dumps(rep, indent=1, default=float))


if __name__ == "__main__":
    main()
