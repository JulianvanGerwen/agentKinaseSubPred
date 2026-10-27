"""Post-TEST exploratory (DESIGN A34, pre-declared target (c) "Tyr specificity"). DEV CV ONLY.

Reproduces the frozen final R2+R3 late fusion on the pooled CV out-of-fold predictions (R2/R3 run
artefacts are read-only and never modified) and tests two cheap Tyr-specific interventions with
**nested** selection: every tunable is chosen on the four training folds of each outer fold.

  A. Y-specific fusion weight (brief (ii)+(iii))
       score = w_Y * z(R2) + (1 - w_Y) * z(R3) on Y-class kinase rows, w_Y picked per outer fold.
       ST rows keep the frozen global w = 0.6.

  B. Tyr motif repair (brief (i))
       score' = score_fused + gamma * z(C_k) on Y-class kinase rows, where C_k is kinase k's own
       Cantley Tyr percentile row if that row ranks k's positives above its E-hard negatives on the
       inner folds, and otherwise the row of k's nearest kinase-domain neighbour that owns an
       in-class Tyr matrix. Both the per-kinase own/neighbour decision and gamma are taken from the
       inner folds only, so the "which matrices are anti-correlated?" question is itself answered
       inside the folds and no label reaches the held-out fold.

Selection metric (inner folds): mean over the inner folds of the macro per-kinase within-fold E-hard
AUROC over Y kinases. Within one fold the A24 mid-rank transform is monotone, so a within-fold AUROC
on raw scores equals the official one exactly; only the aggregation differs from pooled-OOF, and it
matches the `M1_fold_mean` convention already used in the fusion sweep.

Reporting metric: the official `ksp.eval.evaluate` pooled-OOF path with A24, plus a paired
substrate-cluster bootstrap restricted to the Y stratum (the same `_boot_matrix` route that
docs/strata_cv.md uses for stratum CIs).

Outputs -> results/models/posttest_tyr/{explore.json, per_kinase_Y.tsv, selection_trace.tsv}
Run: /workspace/.venv/bin/python -I scripts/models/posttest_tyr/tyr_fusion_explore.py --threads 8
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from ksp import eval as kev  # noqa: E402
from ksp import index, labels, splits, store  # noqa: E402

OUT = ROOT / "results/models/posttest_tyr"
OUT.mkdir(parents=True, exist_ok=True)

R2_RUN = "r2_final_ens5__20261026-191213__s-1"            # R2 side of the frozen final fusion
R3_POOLED = "r3_w7_nce05_cpt_aux_cv_pooled__20261026-162205__s0"
W_FROZEN = 0.6                                             # frozen final fusion weight for `score`
WS = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
GAMMAS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
SEED = 0


def load_preds(run: str) -> pd.DataFrame:
    m = json.load(open(ROOT / "results" / "runs" / run / "manifest.json"))
    if "predictions" not in m["outputs"]:
        R3F = json.load(open(ROOT / "results/models/r3/FINAL_OOF.json"))
        return pd.concat([load_preds(x) for x in R3F["cv_fold_runs"]], ignore_index=True)
    p = pd.read_parquet(ROOT / m["outputs"]["predictions"]["path"],
                        columns=["kinase_idx", "site_idx", "score"])
    return p.astype({"kinase_idx": np.int32, "site_idx": np.int32, "score": np.float64})


def eval_items(data, k, kcls):
    """(positive site indices, E-hard negative site indices) for kinase k -- evaluate()'s own logic."""
    nS = len(data.sites)
    pos = data.positives.get(k, np.zeros(0, np.int64))
    cand = data.universe & kev._res_ok(data.sites.res.to_numpy(), kcls)
    excl = np.zeros(nS, bool); excl[data.masks.excluded(k)] = True
    posm = np.zeros(nS, bool); posm[pos] = True
    hard = np.zeros(nS, bool); hard[data.masks.hard_candidates(k)] = True
    return pos, np.flatnonzero(cand & ~excl & ~posm & hard)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    t0 = time.time()
    rep = {"design": "DESIGN A34 post-TEST exploratory; DEV CV only; released model unchanged",
           "r2_run": R2_RUN, "r3_run": R3_POOLED, "w_frozen": W_FROZEN, "boot_B": a.boot}

    ki, si = index.load_kinases(), index.load_sites()
    kacc = ki.uniprot_acc.to_numpy()
    kcls_of = dict(zip(kacc, ki.residue_class))
    y_cls = ki.residue_class.isin(["Y", "STY"]).to_numpy()
    y_kidx = set(np.flatnonzero(y_cls).tolist())

    print("loading R2 / R3 OOF predictions ...", flush=True)
    A, Bp = load_preds(R2_RUN), load_preds(R3_POOLED)
    M = A.merge(Bp, on=["kinase_idx", "site_idx"], how="inner", suffixes=("_r2", "_r3"))
    if len(M) != len(A):
        raise ValueError(f"R2/R3 rows differ: {len(A)} vs {len(M)}")
    del A, Bp
    print(f"  {len(M):,} pooled OOF rows ({time.time()-t0:.0f}s)", flush=True)

    data = kev.prepare("cv", fold=None)
    site_fold = data.sites.acc.map(splits.load_cv_folds()).to_numpy()
    kidx = M.kinase_idx.to_numpy(); sidx = M.site_idx.to_numpy()
    fold = site_fold[sidx]

    # per-fold z standardisation over ALL rows, exactly as scripts/models/r2_fusion.py
    for c in ("score_r2", "score_r3"):
        z = np.array(M[c].to_numpy(), dtype=np.float64, copy=True)
        for f in range(5):
            mk = fold == f
            z[mk] = (z[mk] - z[mk].mean()) / z[mk].std()
        M[c] = z
    z2, z3 = M.score_r2.to_numpy(), M.score_r3.to_numpy()
    base_score = W_FROZEN * z2 + (1 - W_FROZEN) * z3
    fold = np.nan_to_num(fold, nan=-1).astype(np.int8)

    def run_eval(score, tag):
        P = pd.DataFrame({"kinase_idx": kidx, "site_idx": sidx, "score": score.astype(np.float32)})
        res = kev.evaluate(P, data, missing="error", keep_items=True)
        print(f"  [{tag}] overall M1 = {res.endpoints['M1_auroc_hard']:.4f} ({time.time()-t0:.0f}s)",
              flush=True)
        return res

    res_base = run_eval(base_score, "baseline (frozen w=0.6)")
    pk0 = res_base.per_kinase.set_index("kinase_acc")
    # Headline Y stratum = residue_class 'Y', matching docs/strata_cv.md §4 and the 0.63 reference.
    # Y+STY is the set the interventions actually touch and is reported as a secondary row.
    pure_y = set(kacc[ki.residue_class.eq("Y").to_numpy()])
    yset = set(kacc[y_cls])
    y_accs = [k for k in pk0.index[pk0.in_K_primary] if k in yset]          # touched set (Y+STY)
    y_pure = [k for k in y_accs if k in pure_y]                              # headline stratum
    st_accs = [k for k in pk0.index[pk0.in_K_primary] if k not in yset]
    rep["n_Y_kinases_K_cv_pure"] = len(y_pure)
    rep["n_Y_kinases_K_cv_incl_STY"] = len(y_accs)
    rep["n_ST_kinases_K_cv"] = len(st_accs)

    def macro(res, accs):
        p = res.per_kinase.set_index("kinase_acc")
        return float(np.nanmean(p.loc[[k for k in accs if k in p.index], "auroc_hard"]
                                .to_numpy(dtype=float)))

    rep["baseline"] = {"M1_all": float(res_base.endpoints["M1_auroc_hard"]),
                       "M1_Y": macro(res_base, y_pure), "M1_Y_inclSTY": macro(res_base, y_accs),
                       "M1_ST": macro(res_base, st_accs)}
    print("  baseline:", json.dumps({k: round(v, 4) for k, v in rep["baseline"].items()}), flush=True)

    # ------------------------------------------------- map each kinase's eval rows into M
    print("mapping per-kinase evaluation rows ...", flush=True)
    key = kidx.astype(np.int64) * (len(si) + 1) + sidx
    order = np.argsort(key, kind="stable")
    ksorted = key[order]
    rows = {}
    for k in y_accs:
        kx = int(np.flatnonzero(kacc == k)[0])
        ip, ih = eval_items(data, k, kcls_of[k])
        for nm, arr in (("pos", ip), ("hard", ih)):
            q = kx * (len(si) + 1) + np.asarray(arr, dtype=np.int64)
            loc = np.searchsorted(ksorted, q)
            ok = (loc < len(ksorted)) & (ksorted[np.clip(loc, 0, len(ksorted) - 1)] == q)
            if not ok.all():
                raise ValueError(f"{k}: {(~ok).sum()}/{len(q)} {nm} rows absent from the OOF frame")
            rows[(k, nm)] = order[loc]
    n_items = sum(len(v) for v in rows.values())
    print(f"  {n_items:,} Y evaluation rows mapped ({time.time()-t0:.0f}s)", flush=True)
    ymask = np.isin(kidx, list(y_kidx))

    def sel_metric(score, folds_sel):
        """Mean over folds_sel of the macro within-fold Y E-hard AUROC (exact under A24)."""
        per_fold = []
        for f in folds_sel:
            v = []
            for k in y_accs:
                rp, rh = rows[(k, "pos")], rows[(k, "hard")]
                rp, rh = rp[fold[rp] == f], rh[fold[rh] == f]
                if len(rp) and len(rh):
                    x = kev.auroc(score[rp], score[rh])
                    if np.isfinite(x):
                        v.append(x)
            if v:
                per_fold.append(np.mean(v))
        return float(np.mean(per_fold)) if per_fold else np.nan

    trace, transfer = [], []

    # ------------------------------------------------- A. Y-specific fusion weight
    print("\n== A. Y-specific fusion weight (nested) ==", flush=True)
    scoreA = base_score.copy()
    chosen_w = {}
    for f in range(5):
        inner = [g for g in range(5) if g != f]
        sc = {w: sel_metric(w * z2 + (1 - w) * z3, inner) for w in WS}
        best = max(sc, key=sc.get)
        chosen_w[f] = best
        sel = ymask & (fold == f)
        scoreA[sel] = (best * z2 + (1 - best) * z3)[sel]
        trace += [{"variant": "A", "outer_fold": f, "param": "w_Y", "value": w,
                   "inner_metric": sc[w], "chosen": w == best} for w in WS]
        transfer.append({"variant": "A", "outer_fold": f, "chosen": best,
                         "inner_gain": sc[best] - sc[W_FROZEN],
                         "outer_base": sel_metric(base_score, [f]),
                         "outer_new": sel_metric(best * z2 + (1 - best) * z3, [f])})
        print(f"  fold {f}: w_Y = {best}  (inner {sc[best]:.4f}, w=0.6 gives {sc[0.6]:.4f}; "
              f"outer {transfer[-1]['outer_base']:.4f} -> {transfer[-1]['outer_new']:.4f})", flush=True)
    rep["A_chosen_w_per_fold"] = chosen_w
    res_A = run_eval(scoreA, "A: Y-specific w")
    rep["A"] = {"M1_all": float(res_A.endpoints["M1_auroc_hard"]),
                "M1_Y": macro(res_A, y_pure), "M1_Y_inclSTY": macro(res_A, y_accs),
                "M1_ST": macro(res_A, st_accs)}

    # ------------------------------------------------- B. Tyr motif repair
    print("\n== B. Tyr motif repair (nested) ==", flush=True)
    C, _ = store.load_matrix("motif/cantley_unprimed_pct_published")
    C = np.asarray(C)
    isY_site = si.res.to_numpy() == "Y"
    have_own = ~np.isnan(C[:, isY_site]).all(axis=1)
    ident = labels.kinase_domain_identity()
    kpos_of = {k: i for i, k in enumerate(kacc)}

    # kinase_domain_identity() is a LONG edge table (acc_a, acc_b, identity, method), not a matrix.
    # Build the symmetric neighbour ranking from it.
    ed = ident.dropna(subset=["identity"])
    ed = pd.concat([ed[["acc_a", "acc_b", "identity"]],
                    ed[["acc_b", "acc_a", "identity"]].rename(
                        columns={"acc_b": "acc_a", "acc_a": "acc_b"})], ignore_index=True)
    ed = ed[ed.acc_a != ed.acc_b].sort_values("identity", ascending=False)
    nbr_rank = {k: g.acc_b.tolist() for k, g in ed.groupby("acc_a", sort=False)}
    nbr_id = {(r.acc_a, r.acc_b): r.identity for r in ed.itertuples()}

    fam_of = dict(zip(kacc, ki.manning_family))

    def neighbour(acc, cross_family: bool):
        """Nearest Y-class kinase-domain neighbour that owns an in-class Tyr matrix.
        cross_family=True skips the kinase's own Manning family, because same-family matrices tend
        to share the same defect (FGFR1's nearest neighbour is FGFR2, AXL's is MERTK, ...)."""
        for j in nbr_rank.get(acc, []):
            jx = kpos_of.get(j)
            if jx is None or not (y_cls[jx] and have_own[jx]):
                continue
            if cross_family and fam_of.get(j) == fam_of.get(acc):
                continue
            return j
        return None

    gene_of = dict(zip(kacc, ki.gene))
    #: the three repair sources the inner folds choose between (plus "none")
    ARMS = ("own", "nbr", "nbr_xfam")
    src_of = {}
    for k in y_accs:
        src_of[(k, "own")] = k
        src_of[(k, "nbr")] = neighbour(k, False)
        src_of[(k, "nbr_xfam")] = neighbour(k, True)
    rep["B_neighbour_of"] = {
        gene_of[k]: {arm: (gene_of.get(src_of[(k, arm)]),
                           round(float(nbr_id.get((k, src_of[(k, arm)]), np.nan)), 3))
                     for arm in ("nbr", "nbr_xfam") if src_of[(k, arm)]}
        for k in y_accs}
    print(f"  repair sources: nbr for {sum(1 for k in y_accs if src_of[(k,'nbr')])}/{len(y_accs)}, "
          f"nbr_xfam for {sum(1 for k in y_accs if src_of[(k,'nbr_xfam')])}/{len(y_accs)}", flush=True)

    # Cantley channel per kinase over ALL of that kinase's OOF rows (not just pos/hard): A24 ranks
    # within (fold, kinase) over the whole candidate set, so the repair has to move the whole set.
    # z-scored over those rows so one gamma has a common scale across kinases.
    krows_all = {}
    o2 = np.argsort(kidx, kind="stable")
    kb = np.searchsorted(kidx[o2], np.arange(len(kacc) + 1))
    for k in y_accs:
        kx = int(np.flatnonzero(kacc == k)[0])
        krows_all[k] = o2[kb[kx]:kb[kx + 1]]
    chan = {}
    for k in y_accs:
        ra = krows_all[k]
        for arm in ARMS:
            src = src_of[(k, arm)]
            if src is None or src not in kpos_of or not have_own[kpos_of[src]]:
                continue
            v = C[kpos_of[src]][sidx[ra]].astype(np.float64)
            fin = np.isfinite(v)
            if fin.sum() < 2 or np.std(v[fin]) < 1e-9:
                continue
            mu, sd = v[fin].mean(), v[fin].std()
            chan[(k, arm)] = np.nan_to_num((v - mu) / sd, nan=0.0)
    for arm in ARMS:
        rep[f"B_n_kinases_with_{arm}"] = sum(1 for k in y_accs if (k, arm) in chan)

    scoreB = base_score.copy()
    gammas, picks = {}, {}
    for f in range(5):
        inner = [g for g in range(5) if g != f]
        # (1) per kinase: own vs neighbour matrix, judged on the inner folds only
        pick = {}
        for k in y_accs:
            ra = krows_all[k]
            pos_in_ra = np.searchsorted(ra, rows[(k, "pos")])
            hard_in_ra = np.searchsorted(ra, rows[(k, "hard")])
            best_u, best_s = None, -np.inf
            for arm in ARMS:
                if (k, arm) not in chan:
                    continue
                v = chan[(k, arm)]
                rp, rh = rows[(k, "pos")], rows[(k, "hard")]
                mp, mh = np.isin(fold[rp], inner), np.isin(fold[rh], inner)
                if mp.sum() == 0 or mh.sum() == 0:
                    continue
                s = kev.auroc(v[pos_in_ra][mp], v[hard_in_ra][mh])
                if np.isfinite(s) and s > best_s:
                    best_u, best_s = arm, s
            # only repair when the inner folds say the chosen channel carries signal at all
            pick[k] = best_u if (best_u is not None and best_s > 0.5) else None
        # (2) one shared gamma, also from the inner folds
        def apply(g, folds_sel):
            s = base_score.copy()
            if g == 0:
                return s
            for k in y_accs:
                if pick.get(k) is None:
                    continue
                ra = krows_all[k]
                m = np.isin(fold[ra], folds_sel)
                s[ra[m]] = base_score[ra[m]] + g * chan[(k, pick[k])][m]
            return s
        sc = {g: sel_metric(apply(g, inner), inner) for g in GAMMAS}
        gbest = max(sc, key=sc.get)
        gammas[f] = gbest
        picks[f] = {k: (pick[k] or "none") for k in y_accs}
        s_f = apply(gbest, [f])
        sel = ymask & (fold == f)
        scoreB[sel] = s_f[sel]
        trace += [{"variant": "B", "outer_fold": f, "param": "gamma", "value": g,
                   "inner_metric": sc[g], "chosen": g == gbest} for g in GAMMAS]
        transfer.append({"variant": "B", "outer_fold": f, "chosen": gbest,
                         "inner_gain": sc[gbest] - sc[0.0],
                         "outer_base": sel_metric(base_score, [f]),
                         "outer_new": sel_metric(s_f, [f])})
        cnt = {arm: sum(1 for k in y_accs if pick.get(k) == arm) for arm in ARMS}
        n_own, n_nb = cnt["own"], cnt["nbr"] + cnt["nbr_xfam"]
        print(f"  fold {f}: gamma = {gbest} (inner {sc[gbest]:.4f}, gamma=0 gives {sc[0.0]:.4f}; "
              f"outer {transfer[-1]['outer_base']:.4f} -> {transfer[-1]['outer_new']:.4f}); "
              f"own {n_own}, nbr {cnt['nbr']}, nbr_xfam {cnt['nbr_xfam']}, "
              f"dropped {len(y_accs)-n_own-n_nb}", flush=True)
    rep["B_gamma_per_fold"] = gammas
    rep["B_picks_per_fold"] = picks
    res_B = run_eval(scoreB, "B: motif repair")
    rep["B"] = {"M1_all": float(res_B.endpoints["M1_auroc_hard"]),
                "M1_Y": macro(res_B, y_pure), "M1_Y_inclSTY": macro(res_B, y_accs),
                "M1_ST": macro(res_B, st_accs)}

    for tag, r in (("A", res_A), ("B", res_B)):
        rep[tag]["ST_drift"] = float(abs(macro(r, st_accs) - rep["baseline"]["M1_ST"]))
        if rep[tag]["ST_drift"] > 1e-9:
            print(f"  WARNING: {tag} moved ST M1 by {rep[tag]['ST_drift']:.2e} (expected 0)", flush=True)

    # ------------------------------------------------- paired cluster bootstrap on the Y stratum
    print("\n== paired substrate-cluster bootstrap, Y stratum ==", flush=True)
    cl = kev._universe_clusters(res_base)
    mats = {t: kev._boot_matrix(r, cl, a.boot, SEED, y_pure)
            for t, r in (("base", res_base), ("A", res_A), ("B", res_B))}
    db = np.nanmean(mats["base"]["M1"], 1)
    rep["baseline"]["M1_Y_lo95"] = float(np.percentile(db, 2.5))
    rep["baseline"]["M1_Y_hi95"] = float(np.percentile(db, 97.5))
    for t in ("A", "B"):
        d = np.nanmean(mats[t]["M1"], 1) - db
        d = d[np.isfinite(d)]
        p = min(1.0, 2 * min(((d <= 0).sum() + 1) / (len(d) + 1), ((d >= 0).sum() + 1) / (len(d) + 1)))
        rep[t].update(Y_delta=rep[t]["M1_Y"] - rep["baseline"]["M1_Y"],
                      Y_delta_lo95=float(np.percentile(d, 2.5)),
                      Y_delta_hi95=float(np.percentile(d, 97.5)), Y_delta_p_boot=float(p))
        print(f"  {t}: Y M1 {rep['baseline']['M1_Y']:.4f} -> {rep[t]['M1_Y']:.4f}  "
              f"delta {rep[t]['Y_delta']:+.4f} [{rep[t]['Y_delta_lo95']:+.4f}, "
              f"{rep[t]['Y_delta_hi95']:+.4f}] p={p:.3f}", flush=True)

    # ------------------------------------------------- outputs
    tab = pd.DataFrame({"kinase_acc": y_accs})
    km = ki.set_index("uniprot_acc")
    tab["gene"] = tab.kinase_acc.map(km.gene)
    tab["manning_family"] = tab.kinase_acc.map(km.manning_family)
    tab["n_pos"] = tab.kinase_acc.map(pk0.n_pos)
    tab["n_neg_hard"] = tab.kinase_acc.map(pk0.n_neg_hard)
    for t, r in (("base", res_base), ("A", res_A), ("B", res_B)):
        tab[f"auroc_{t}"] = tab.kinase_acc.map(r.per_kinase.set_index("kinase_acc").auroc_hard)
    tab["delta_A"] = tab.auroc_A - tab.auroc_base
    tab["delta_B"] = tab.auroc_B - tab.auroc_base
    tab["B_pick_fold0"] = tab.kinase_acc.map(picks[0])
    tab.sort_values("n_pos", ascending=False).round(4).to_csv(OUT / "per_kinase_Y.tsv",
                                                              sep="\t", index=False)
    pd.DataFrame(trace).round(5).to_csv(OUT / "selection_trace.tsv", sep="\t", index=False)
    TR = pd.DataFrame(transfer)
    TR["outer_gain"] = TR.outer_new - TR.outer_base
    TR.round(5).to_csv(OUT / "inner_outer_transfer.tsv", sep="\t", index=False)
    rep["transfer"] = {v: {"mean_inner_gain": float(g.inner_gain.mean()),
                           "mean_outer_gain": float(g.outer_gain.mean()),
                           "folds_with_outer_gain": int((g.outer_gain > 0).sum())}
                       for v, g in TR.groupby("variant")}
    print("\n--- inner-fold gain vs held-out gain ---\n", TR.round(4).to_string(index=False))
    rep["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "explore.json").write_text(json.dumps(rep, indent=2, default=str))
    print("\n" + json.dumps({k: rep[k] for k in ("baseline", "A", "B")}, indent=2, default=str))
    print(f"\ndone in {rep['runtime_s']:.0f}s")


if __name__ == "__main__":
    main()
