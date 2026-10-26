"""Build the F9 tyrosine-kinase feature arrays (all via ksp.store; none label-derived).

Outputs (cache/features/tyr/, loaded with ksp.store.load_matrix("tyr/<stem>")):

  kinase x site (503 x 350,736), float32, NaN outside {Y site, kinase of Y/STY residue class}:
    sh2_dock_max      max SH2 log2 score over the OTHER observed (U_ms) pY sites of the same protein
    sh2_dock_local    log2 sum_j exp2(score_j) * exp(-|dpos_j| / TAU)  -- distance-weighted local
                      "effective pY-docking affinity" around the query site (TAU = 50 aa)
    sh2_dock_near     max SH2 score over other observed pY within +/-50 aa (NaN if none)
    sh2_dock_delta    sh2_dock_max minus the SH2 score of the query site's own window
                      (does the kinase prefer a *neighbour* over this site as a docking anchor?)
    sh2_dock_spec     z-score of sh2_dock_max across all 76 SH2 domains at the same site: is this
                      pY neighbourhood preferred by THIS kinase's SH2 specifically?  This is the
                      kinase-specific component, with the generic "good docking nearby" /
                      site-popularity component divided out.
    rtk_pm            substrate plasma-membrane/cell-junction score, exposed ONLY for
                      transmembrane-receptor kinases (NaN for non-RTKs) -- explicit RTK x membrane
                      interaction
    itam_gated        site ITAM order (0/1/2), exposed ONLY for Syk/Tec/Src-family kinases

  site vectors (350,736 x n):
    site_sh2 [4]      n_other_py_protein, n_other_py_50aa, adaptor_sh2_max, adaptor_sh2_local
                      (max / distance-weighted docking score over the 51 NON-kinase SH2 domains:
                       how strongly the neighbourhood recruits SH2 adaptors at all)
    site_itam [4]     itam_half, itam_full, itim, itam_order

  kinase vectors (503 x 6):
    kinase_desc       n_sh2, n_sh3, n_ptb, is_tm, is_receptor, has_sh2_matrix

Run: /workspace/.venv/bin/python -I scripts/features/tyr_build.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import index  # noqa: E402
from ksp.features import tyr  # noqa: E402
from ksp.store import load_matrix, save_matrix  # noqa: E402
from ksp.uniprot import DEFAULT_DIR, read_fasta  # noqa: E402

OUT = ROOT / "cache/features/tyr"
RES = ROOT / "results/features/tyr"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

UNIPROT_FASTA = DEFAULT_DIR / "proteome" / "UP000005640_reviewed_isoforms.fasta.gz"
TAU = 50.0          # aa decay length for the distance-weighted docking score
NEAR = 50           # aa window for the "near" variant
SYK_TEC_SRC = {"Syk", "Tec", "Src", "Csk", "Abl", "Fer", "Ack", "FAK"}   # Manning families, non-receptor TKs


def main():
    t0 = time.time()
    ki, si = index.load_kinases(), index.load_sites()
    nk, ns = len(ki), len(si)
    kpos = pd.Series(np.arange(nk), index=ki.uniprot_acc)
    y_cls = ki.residue_class.isin(["Y", "STY"]).to_numpy()

    # ---------------------------------------------------------------- SH2 matrices
    mats = tyr.load_sh2_matrices()
    names, tab = tyr.sh2_lookup(mats)
    name_pos = {n: i for i, n in enumerate(names)}
    k2doms = {}
    for sh2, gene in tyr.SH2_TO_KINASE.items():
        if sh2 in name_pos:
            k2doms.setdefault(gene, []).append(name_pos[sh2])
    adaptor_rows = [i for n, i in name_pos.items() if n not in tyr.SH2_TO_KINASE]
    print(f"SH2 matrices: {len(names)} domains, {len(k2doms)} kinases, {len(adaptor_rows)} adaptors")

    # ---------------------------------------------------------------- observed pY sites (U_ms)
    isY = (si.res.to_numpy() == "Y")
    ums = si.in_U_ms.to_numpy().astype(bool)
    obs = np.flatnonzero(isY & ums)                 # anchor candidates: observed pY only
    allY = np.flatnonzero(isY)                      # query rows: every Y site in the index
    print(f"Y sites: {len(allY)} total, {len(obs)} observed (U_ms)")

    # SH2 score of every observed pY window, for every domain
    wins = si.window15.to_numpy()
    S_obs = tyr.score_windows(tab, wins[obs])               # (n_dom, n_obs)
    S_all = tyr.score_windows(tab, wins[allY])              # (n_dom, n_allY)

    acc_obs, pos_obs = si.acc.to_numpy()[obs], si.pos.to_numpy()[obs].astype(np.int64)
    acc_all, pos_all = si.acc.to_numpy()[allY], si.pos.to_numpy()[allY].astype(np.int64)
    # group observed anchors by protein
    order = np.argsort(acc_obs, kind="stable")
    acc_sorted = acc_obs[order]
    bounds = {a: (lo, hi) for a, lo, hi in zip(*_runs(acc_sorted))}

    # per-kinase aggregates over query Y sites.  Collapse each kinase's tandem SH2 domains to their
    # max up front, so the per-site loop is fully vectorised over kinases.
    kin_genes = sorted(g for g in k2doms if g in set(ki.gene))
    g2a = dict(zip(ki.gene, ki.uniprot_acc))
    Skin_obs = np.stack([S_obs[k2doms[g]].max(axis=0) for g in kin_genes])     # (n_kin, n_obs)
    Skin_own = np.stack([S_all[k2doms[g]].max(axis=0) for g in kin_genes])     # (n_kin, n_allY)
    Sad_obs = S_obs[adaptor_rows].max(axis=0) if adaptor_rows else None
    E_obs = np.exp2(np.clip(Skin_obs, -20, 20))
    Ead_obs = np.exp2(np.clip(Sad_obs, -20, 20)) if adaptor_rows else None

    nq = len(allY)
    agg = {m: np.full((len(kin_genes), nq), np.nan, np.float32)
           for m in ("max", "local", "near", "delta", "spec")}
    n_other = np.zeros(nq, np.float32)
    n_other_near = np.zeros(nq, np.float32)
    ad_max = np.full(nq, np.nan, np.float32)
    ad_local = np.full(nq, np.nan, np.float32)

    for qi in range(nq):
        b = bounds.get(acc_all[qi])
        if b is None:
            continue
        sel = order[b[0]:b[1]]                         # positions inside `obs` / `S_obs` columns
        d = pos_obs[sel] - pos_all[qi]
        sel, d = sel[d != 0], d[d != 0]                # exclude the query site itself
        if len(sel) == 0:
            continue
        ad = np.abs(d)
        n_other[qi] = len(sel)
        near = ad <= NEAR
        n_other_near[qi] = int(near.sum())
        w = np.exp(-ad / TAU).astype(np.float32)
        sc = Skin_obs[:, sel]                          # (n_kin, m)
        agg["max"][:, qi] = sc.max(axis=1)
        agg["local"][:, qi] = np.log2(E_obs[:, sel] @ w + 1e-12)
        if near.any():
            agg["near"][:, qi] = sc[:, near].max(axis=1)
        agg["delta"][:, qi] = agg["max"][:, qi] - Skin_own[:, qi]
        # specificity contrast: is this pY neighbourhood preferred by THIS kinase's SH2 relative to
        # the whole SH2 panel?  Removes the generic "many/good docking sites nearby" component.
        dm = S_obs[:, sel].max(axis=1)                 # best anchor per SH2 domain (all 76)
        sd = dm.std()
        agg["spec"][:, qi] = (sc.max(axis=1) - dm.mean()) / (sd if sd > 1e-6 else np.nan)
        if adaptor_rows:
            ad_max[qi] = Sad_obs[sel].max()
            ad_local[qi] = np.log2(Ead_obs[sel] @ w + 1e-12)
        if qi % 20000 == 0:
            print(f"  {qi}/{nq}  {time.time()-t0:.0f}s", flush=True)

    # ---------------------------------------------------------------- scatter into kinase x site
    meta_common = dict(label_derived=False, source="Yaron-Barir 2024 Suppl. Table 5 SH2 matrices",
                       phospho_input="U_ms membership of OTHER sites only (never the query site)")
    for m, stem in (("max", "sh2_dock_max"), ("local", "sh2_dock_local"),
                    ("near", "sh2_dock_near"), ("delta", "sh2_dock_delta"),
                    ("spec", "sh2_dock_spec")):
        A = np.full((nk, ns), np.nan, np.float32)
        for ri, g in enumerate(kin_genes):
            kr = kpos[g2a[g]]
            if not y_cls[kr]:
                continue
            A[kr, allY] = agg[m][ri]
        save_matrix(f"tyr/{stem}", A, "kinase", "site",
                    meta=dict(description=f"SH2 docking ({m}) on other observed pY of the substrate",
                              tau_aa=TAU, near_aa=NEAR, n_kinases=int(np.isfinite(A).any(1).sum()),
                              **meta_common))
        print(f"saved tyr/{stem}  kinases={int(np.isfinite(A).any(1).sum())} "
              f"nonnan_cells={int(np.isfinite(A).sum())}")

    # ---------------------------------------------------------------- site vectors
    sv = np.full((ns, 4), np.nan, np.float32)
    sv[allY, 0] = n_other
    sv[allY, 1] = n_other_near
    sv[allY, 2] = ad_max
    sv[allY, 3] = ad_local
    save_matrix("tyr/site_sh2", sv, "site", "none",
                col_ids=["n_other_py_protein", "n_other_py_50aa", "adaptor_sh2_max", "adaptor_sh2_local"],
                meta=dict(description="pY-neighbourhood descriptors; adaptor_* use the 51 non-kinase SH2 domains",
                          caution="n_other_py_* is a site-popularity proxy; ablate it",
                          **meta_common))

    fa = read_fasta(UNIPROT_FASTA)
    can = fa[(fa.isoform == "") & (fa.db == "sp")]
    seqs = dict(zip(can.acc, can.seq))
    it = tyr.itam_features(seqs, si.acc.to_numpy(), si.pos.to_numpy().astype(np.int64))
    iv = np.array(it.to_numpy(dtype=np.float32), copy=True)
    iv[~isY] = np.nan
    save_matrix("tyr/site_itam", iv, "site", "none", col_ids=list(it.columns),
                meta=dict(description="ITAM/ITIM-like context (Yxx[LIV], tandem pair, [ILVS]xYxx[LV])",
                          label_derived=False, source="UniProt 2026_03 canonical sequences"))

    # ---------------------------------------------------------------- kinase descriptors
    kd = tyr.kinase_adaptor_domains(ki.uniprot_acc.tolist())
    kd["has_sh2_matrix"] = kd.acc.map({g2a[g]: 1 for g in kin_genes}).fillna(0).astype(int)
    kcols = ["n_sh2", "n_sh3", "n_ptb", "is_tm", "is_receptor", "has_sh2_matrix"]
    save_matrix("tyr/kinase_desc", kd[kcols].to_numpy(dtype=np.float32), "kinase", "none",
                col_ids=kcols,
                meta=dict(description="kinase adaptor-domain descriptors", label_derived=False,
                          source="UniProt 2026_03 Domain [FT] / Keywords"))
    kd.to_csv(RES / "kinase_adaptor_domains.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- RTK x plasma membrane
    vh, mh = load_matrix("context/loc/vec_hpa_protein")
    vc, mc = load_matrix("context/loc/vec_hcm_protein")
    voc = list(np.load(ROOT / "cache/features/context/loc/vec_hpa_protein.cols.npy", allow_pickle=True))
    pm = [voc.index("plasma_membrane"), voc.index("cell_junction")]
    prot = index.load_proteins()
    p_pm = np.nanmax(np.stack([np.asarray(vh)[:, pm].sum(1), np.asarray(vc)[:, pm].sum(1)]), axis=0)
    ppos = pd.Series(np.arange(len(prot)), index=prot.acc)
    site_prot = ppos.reindex(si.acc).to_numpy()
    site_pm = p_pm[site_prot].astype(np.float32)
    is_rtk = (kd.set_index("acc").is_tm.reindex(ki.uniprot_acc).to_numpy() > 0) & y_cls
    A = np.full((nk, ns), np.nan, np.float32)
    A[np.ix_(is_rtk, allY)] = site_pm[allY]
    save_matrix("tyr/rtk_pm", A, "kinase", "site",
                meta=dict(description="substrate plasma-membrane/cell-junction score, RTK rows only",
                          source="HPA IF v25 + Human Cell Map (F4 stores)", label_derived=False,
                          n_rtk=int(is_rtk.sum())))

    # ---------------------------------------------------------------- ITAM gated to Syk/Tec/Src etc.
    gate = (ki.manning_family.isin(SYK_TEC_SRC).to_numpy() & y_cls)
    A = np.full((nk, ns), np.nan, np.float32)
    A[np.ix_(gate, allY)] = iv[allY, 3]
    save_matrix("tyr/itam_gated", A, "kinase", "site",
                meta=dict(description="site ITAM order, non-receptor TK rows only",
                          gated_families=sorted(SYK_TEC_SRC), label_derived=False,
                          n_kinases=int(gate.sum())))

    json.dump(dict(n_sh2_domains=len(names), n_sh2_kinases=len(kin_genes),
                   sh2_kinases=sorted(kin_genes), n_adaptor_domains=len(adaptor_rows),
                   n_rtk_rows=int(is_rtk.sum()), n_itam_gated_rows=int(gate.sum()),
                   tau_aa=TAU, near_aa=NEAR, runtime_s=round(time.time() - t0, 1)),
              open(RES / "tyr_build_meta.json", "w"), indent=2)
    print(f"done in {time.time()-t0:.0f}s")


def _runs(sorted_keys):
    """(unique keys, start, stop) of equal runs in a sorted array."""
    u, st = np.unique(sorted_keys, return_index=True)
    sp = np.r_[st[1:], len(sorted_keys)]
    return u, st, sp


if __name__ == "__main__":
    main()
