#!/usr/bin/env python
"""CV strata report (DESIGN §6.4, §9, §17 A24/A26/A27): M1/M2/M3/M5 per stratum with cluster-bootstrap CIs for the
main pooled-OOF runs, using the frozen ksp.eval functions and the canonical stratum definitions in ksp.strata.

Per run (worker mode, `--tags`):
  1. predictions -> ksp.eval.map_predictions -> ksp.eval.fold_rank_normalise (A24, exactly the step evaluate() does
     for pooled CV; done once so the frame can be reused), then ksp.eval.evaluate(..., fold_normalise=False).
     The full-run M1/M2/M3 are checked against the run's metrics_a24.json (leaderboard).
  2. kinase-level strata: per-kinase metrics of that evaluation, macro over K_cv ∩ stratum (M1, M2, M3) and over
     K_any ∩ stratum (M5 = leaderboard 'M5bin' convention); CIs from ksp.eval._boot_matrix (the matrix
     ksp.eval.bootstrap averages; same weights; B=1000 substrate-cluster resamples, seed 0).
  3. site/pair-level strata (evidence, hub, primed): positive subsets. M1 = AUROC(stratum positives, the kinase's
     unchanged E-hard negatives) - identical to re-evaluating with subset labels because removed positives stay in
     KP(k) and are therefore masked, never negatives (verified against a subset-label evaluate() for R2). M2 = mean
     Hit@5 of the stratum's site pairs (independent of other positives). M3 = recall of the stratum's positives in
     the kinase's top-1% list of the full evaluation (fixed cut-off; differs from a subset re-evaluation only by
     the cut-off count, see sanity file).
  4. A26 hub-excluded diagnostic: re-evaluation with hub-site positives removed AND hub sites removed from E-hard.
Combine mode (`--combine`): results/strata/cv_strata.tsv, cv_strata_site_per_kinase.tsv, cv_strata_sanity.tsv.
Usage: python -I scripts/strata_cv_report.py --tags R2,B1ii   |   --combine
"""
import argparse
import dataclasses
import json
import os
import sys
import time
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_MAX_THREADS"):
    os.environ.setdefault(v, "8")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import labels  # noqa: E402
from ksp import strata as ST  # noqa: E402

RUNS = ROOT / "results" / "runs"
PRED = ROOT / "cache" / "runs"
WORK = ROOT / "cache" / "strata"
OUT = ROOT / "results" / "strata"
B_BOOT = 1000
SEED = 0
SMALL_K, SMALL_POS = 5, 30          # pre-declared: a stratum is 'small' if < 5 kinases (with defined M1) or < 30 pos.

R2_ID = "r2_final_knn_s0_cv5__20261026-181711__s0"
B1II_ID = "simple_b1_universe_reeval__20261026-155703__s0"
RUNSPEC = {   # tag: (run_id, description, kinase-restriction source run or None)
    "R2": (R2_ID, "R2 final (kNN, seed 0) pooled 5-fold, latest R2 final run on LEADERBOARD", None),
    "R3": ("r3_final_distill_cv_pooled__20261026-173414__s0", "R3 final (distilled) pooled CV", None),
    "B1i": ("simple_b1_published_reeval__20261026-155352__s0", "B1 Cantley arm (i) published percentile", None),
    "B1ii": (B1II_ID, "B1 Cantley arm (ii) universe-recalibrated", None),
    "B1prim": ("simple_b1_primed_hard_reeval__20261026-155750__s0", "B1 Cantley, hard priming (pp=True on U_ms)",
               None),
    "B2": ("simple_b2_pssm_reeval__20261026-160718__s0", "B2 PSSM", None),
    "B4": ("b4_selphi2_noleak__20261026-171245__s0", "B4 SELPHI2 (no-leak retrain; 33% coverage)", None),
    "B7": ("simple_b7_lr_reeval__20261026-160451__s0", "B7 simple-strong LR", None),
    "B0prior": ("simple_b0_prior_reeval__20261026-155334__s0", "B0-prior floor", None),
    "B0pop": ("simple_b0_pop_reeval__20261026-155346__s0", "B0-pop floor", None),
    "B0miss": ("simple_b0_miss_reeval__20261026-155334__s0", "B0-miss floor", None),
    # A27 matched-kinase view vs Cantley: both restricted to the kinases B1 arm (ii) scores
    "R2_kmB1": (R2_ID, "R2 restricted to kinases scored by B1 (A27 matched-kinase)", B1II_ID),
    "B1ii_kmB1": (B1II_ID, "B1 arm (ii) on its own scored kinases (A27 matched-kinase)", B1II_ID),
}


# ============================================================================================== helpers
class HardDrop:
    """NegativeMasks view with some sites removed from E-hard (A26 hub-excluded diagnostic); excluded() unchanged."""

    def __init__(self, m, drop_idx):
        self._m = m
        self._drop = np.asarray(sorted(drop_idx), np.int64)

    def excluded(self, k):
        return self._m.excluded(k)

    def hard_candidates(self, k):
        h = self._m.hard_candidates(k)
        return h[~np.isin(h, self._drop)]

    def __getattr__(self, n):
        return getattr(self._m, n)


def scored_kinases(run_id: str, data) -> set:
    import pyarrow.parquet as pq
    f = PRED / run_id / "predictions.parquet"
    names = pq.ParquetFile(f).schema_arrow.names
    p = pd.read_parquet(f, columns=["kinase_idx"] if "kinase_idx" in names else ["kinase_acc"])
    if "kinase_idx" in p:
        ki = np.unique(p.kinase_idx.values)
        return set(data.kinases.uniprot_acc.values[ki])
    return set(p.kinase_acc.unique())


def load_normalised(run_id: str, data, kinase_filter: set | None, chunk: int = 48) -> pd.DataFrame:
    """Memory-lean A24 normalisation: ksp.eval.map_predictions + ksp.eval.fold_rank_normalise applied per block of
    kinases (normalisation groups are (kinase, fold), so blocking by kinase gives identical values). After
    normalisation, rows that no endpoint reads are dropped: canary kinases, non-universe sites and sites outside the
    kinase's residue class (required pairs, E-hard, E-all and M2 candidates are all universe- and class-restricted).
    Returns kinase-sorted int-ID predictions with score (fold percentile) and score_site (raw, for M2)."""
    p = pd.read_parquet(PRED / run_id / "predictions.parquet")
    cols = [c for c in ["kinase_acc", "site_id", "kinase_idx", "site_idx", "score", "score_site"] if c in p]
    p = p[cols]
    kacc = data.kinases.uniprot_acc.values
    if "kinase_idx" not in p:      # string IDs: map once (validated again by map_predictions below)
        p = E.map_predictions(p, data)
        p = p.rename(columns={"k": "kinase_idx", "s": "site_idx"})
    kk = p.kinase_idx.values
    real = ~pd.Series(kacc).str.startswith("CANARY_").values
    allowed = real.copy()
    if kinase_filter is not None:
        allowed &= np.isin(kacc, list(kinase_filter))
    res_arr = data.sites.res.values
    cls_ok = {c: E._res_ok(res_arr, c) & data.universe for c in ["ST", "Y", "STY"]}
    kcls = data.kinases.residue_class.values
    parts = []
    for lo in range(0, len(kacc), chunk):
        ks_ = [i for i in range(lo, min(lo + chunk, len(kacc))) if allowed[i]]
        if not ks_:
            continue
        m = (kk >= lo) & (kk < lo + chunk)
        sub = p[m]
        sub = sub[np.isin(sub.kinase_idx.values, ks_)]
        if not len(sub):
            continue
        P = E.map_predictions(sub, data)
        del sub
        P = E.fold_rank_normalise(P, data)
        keep = np.zeros(len(P), bool)
        for i in np.unique(P.k.values):
            mi = P.k.values == i
            keep[mi] = cls_ok[kcls[i]][P.s.values[mi]]
        P = P[keep]
        parts.append(pd.DataFrame({"kinase_idx": P.k.values.astype(np.int16), "site_idx": P.s.values.astype(np.int32),
                                   "score": P.score.values, "score_site": P.score_site.values}))
        del P
    del p, kk
    out = pd.concat(parts, ignore_index=True)
    del parts
    return out.sort_values("kinase_idx", kind="stable", ignore_index=True)


def ci_row(boot_cols: dict, point: dict) -> dict:
    out = {}
    for m, v in point.items():
        out[m] = v
        d = boot_cols.get(m)
        if d is None or not len(d):
            out[f"{m}_lo"] = out[f"{m}_hi"] = np.nan
            continue
        d = d[~np.isnan(d)]
        out[f"{m}_lo"] = float(np.percentile(d, 2.5)) if len(d) else np.nan
        out[f"{m}_hi"] = float(np.percentile(d, 97.5)) if len(d) else np.nan
    return out


def macro(pk: pd.DataFrame, ks, col):
    ks = [k for k in ks if k in pk.index]
    return E._nanmean(pk.loc[ks, col]) if ks and col in pk else np.nan


def stratum_summary(pk, bm_p, bm_a, kp_order, ka_order, kset, label: dict) -> dict:
    """pk indexed by kinase_acc (per-kinase metrics); bm_*: _boot_matrix outputs over kp_order / ka_order."""
    kp = [k for k in kp_order if k in kset]
    ka = [k for k in ka_order if k in kset]
    jp = [i for i, k in enumerate(kp_order) if k in kset]
    ja = [i for i, k in enumerate(ka_order) if k in kset]
    point = {"M1": macro(pk, kp, "auroc_hard"), "M2": macro(pk, kp, "hit5"), "M3": macro(pk, kp, "top1_recall"),
             "M5": macro(pk, ka, "auroc_hard")}
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            bc = {"M1": np.nanmean(bm_p["M1"][:, jp], 1) if jp else None,
                  "M2": np.nanmean(bm_p["M2"][:, jp], 1) if jp else None,
                  "M3": np.nanmean(bm_p["M3"][:, jp], 1) if jp else None,
                  "M5": np.nanmean(bm_a["M1"][:, ja], 1) if ja else None}
    r = dict(label)
    r["n_kinases_Kcv"] = int(pk.loc[kp, "auroc_hard"].notna().sum()) if kp else 0
    r["n_kinases_Kany"] = int(pk.loc[ka, "auroc_hard"].notna().sum()) if ka else 0
    r["n_pos_Kcv"] = int(pk.loc[kp, "n_pos"].sum()) if kp else 0
    r["n_pos_Kany"] = int(pk.loc[ka, "n_pos"].sum()) if ka else 0
    r.update(ci_row(bc, point))
    r["small"] = bool(r["n_kinases_Kcv"] < SMALL_K or r["n_pos_Kcv"] < SMALL_POS)
    r["small_Kany"] = bool(r["n_kinases_Kany"] < SMALL_K or r["n_pos_Kany"] < SMALL_POS)
    r["genes_Kcv"] = ",".join(sorted(pk.loc[kp, "gene"].astype(str))) if len(kp) <= 15 else f"{len(kp)} kinases"
    return r


def subset_result(res, data, keep: dict):
    """keep: kinase_acc -> bool array over data.positives[k] (= order of res.items[k]['pos_s']). Returns an
    EvalResult whose items/site_level/per_kinase hold only the kept positives (negatives unchanged)."""
    items, rows = {}, []
    pk = res.per_kinase.set_index("kinase_acc")
    sl = res.site_level
    sl_keys = []
    for k, m in keep.items():
        it = res.items.get(k)
        if it is None or not m.any():
            continue
        assert len(it["pos_s"]) == len(m) == len(data.positives[k]), k
        items[k] = {"pos_s": it["pos_s"][m], "pos_c": it["pos_c"][m], "hard_s": it["hard_s"], "hard_c": it["hard_c"],
                    "top_hits": it["top_hits"][m], "all_idx_hash": it["all_idx_hash"]}
        sl_keys.append(pd.DataFrame({"kinase_acc": k, "site_idx": data.positives[k][m]}))
        rows.append({"kinase_acc": k, "n_pos": int(m.sum()),
                     "auroc_hard": E.auroc(it["pos_s"][m], it["hard_s"]),
                     "top1_recall": float(it["top_hits"][m].mean())})
    keys = pd.concat(sl_keys, ignore_index=True) if sl_keys else pd.DataFrame(columns=["kinase_acc", "site_idx"])
    sls = sl.merge(keys, on=["kinase_acc", "site_idx"], how="inner") if len(sl) else sl
    new = pd.DataFrame(rows)
    if len(new):
        h = sls.groupby("kinase_acc").hit5.mean() if len(sls) else pd.Series(dtype=float)
        new["hit5"] = new.kinase_acc.map(h)
        for c in ["gene", "in_K_primary", "in_K_any", "residue_class", "manning_group"]:
            new[c] = new.kinase_acc.map(pk[c])
    return dataclasses.replace(res, items=items, site_level=sls, per_kinase=new)


def boot_pair(res, cl, kp_order, ka_order):
    """_boot_matrix over K_any ∪ K_cv once; the K_cv matrix is a column subset (per-kinase columns are independent
    and the cluster weights depend only on (clusters, B, seed), so this equals a separate call over kp_order)."""
    allk = list(ka_order) + [k for k in kp_order if k not in set(ka_order)]
    bm = E._boot_matrix(res, cl, B_BOOT, SEED, allk)
    pos = {k: i for i, k in enumerate(allk)}
    jp = [pos[k] for k in kp_order]
    ja = [pos[k] for k in ka_order]
    return {m: v[:, jp] for m, v in bm.items()}, {m: v[:, ja] for m, v in bm.items()}


# ============================================================================================== one run
def run_one(tag: str, data, ktab: pd.DataFrame, site_attr: dict, hub_idx: np.ndarray, check_equiv: bool) -> None:
    run_id, desc, kfrom = RUNSPEC[tag]
    t0 = time.time()
    kfilter = scored_kinases(kfrom, data) if kfrom else None
    Pn = load_normalised(run_id, data, kfilter)
    t_load = time.time() - t0
    res = E.evaluate(Pn, data, missing="worst", fold_normalise=False, keep_items=True)
    t_eval = time.time() - t0 - t_load
    sanity = {"tag": tag, "run_id": run_id, "desc": desc, "n_pred_rows": int(len(Pn)),
              "M1_here": res.endpoints["M1_auroc_hard"], "M2_here": res.endpoints["M2_hit5"],
              "M3_here": res.endpoints["M3_top1_recall"], "M1_fold_mean_here": res.endpoints.get("M1_fold_mean")}
    mf = RUNS / run_id / "metrics_a24.json"
    if mf.exists() and kfilter is None:
        e = json.loads(mf.read_text())["endpoints"]
        sanity.update({"M1_leaderboard": e["M1_auroc_hard"], "M2_leaderboard": e["M2_hit5"],
                       "M3_leaderboard": e["M3_top1_recall"]})
    pk = res.per_kinase[~res.per_kinase.kinase_acc.str.startswith("CANARY_")].set_index("kinase_acc")
    kset_filter = kfilter if kfilter is not None else set(pk.index)
    kp_order = [k for k in data.k_primary if k in pk.index and k in kset_filter]
    ka_order = [k for k in data.k_any if k in pk.index and k in kset_filter]
    cl = E._universe_clusters(res)
    bm_p, bm_a = boot_pair(res, cl, kp_order, ka_order)
    # cross-check against ksp.eval.bootstrap on the unrestricted run
    if kfilter is None:
        bt = E.bootstrap(res, B=B_BOOT, seed=SEED).set_index("endpoint")
        sanity["boot_M1_lo_eval"] = float(bt.loc["M1", "lo95"])
        sanity["boot_M1_hi_eval"] = float(bt.loc["M1", "hi95"])
    rows = []
    base = {"tag": tag, "run_id": run_id}
    rows.append(stratum_summary(pk, bm_p, bm_a, kp_order, ka_order, set(kp_order) | set(ka_order),
                                {**base, "level": "all", "stratum": "all", "value": "all"}))
    for (s, v), ks in ST.kinase_sets(ktab).items():
        rows.append(stratum_summary(pk, bm_p, bm_a, kp_order, ka_order, ks,
                                    {**base, "level": "kinase", "stratum": s, "value": v}))
    t_kin = time.time() - t0 - t_load - t_eval
    # ---- site/pair-level strata (positive subsets)
    site_pk = []
    for (s, v), (kin_ok, attr) in site_attr.items():
        keep = {}
        for k in pk.index:
            if k not in kset_filter or k not in kin_ok or k not in data.positives:
                continue
            keep[k] = attr[k] == v
        sub = subset_result(res, data, keep)
        spk = sub.per_kinase.set_index("kinase_acc") if len(sub.per_kinase) else pd.DataFrame()
        if not len(spk):
            continue
        kpo = [k for k in kp_order if k in spk.index]
        kao = [k for k in ka_order if k in spk.index]
        sbp, sba = boot_pair(sub, cl, kpo, kao)
        rows.append(stratum_summary(spk, sbp, sba, kpo, kao, set(spk.index),
                                    {**base, "level": "site", "stratum": s, "value": v}))
        x = sub.per_kinase.assign(tag=tag, stratum=s, value=v)
        site_pk.append(x[["tag", "stratum", "value", "kinase_acc", "gene", "in_K_primary", "in_K_any", "n_pos",
                          "auroc_hard", "hit5", "top1_recall"]])
    t_site = time.time() - t0 - t_load - t_eval - t_kin
    # ---- A26 hub-excluded re-evaluation (hub positives removed and hub sites removed from E-hard)
    hub_set = set(hub_idx.tolist())
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
    pm = pm[~pm.site_idx.isin(hub_set)].reset_index(drop=True)
    posd = {k: np.sort(g.site_idx.values) for k, g in pm.groupby("kinase_acc")}
    dh = dataclasses.replace(data, positives=posd, pos_meta=pm, masks=HardDrop(data.masks, hub_idx))
    resh = E.evaluate(Pn, dh, missing="worst", fold_normalise=False, keep_items=True)
    pkh = resh.per_kinase.set_index("kinase_acc")
    kpo = [k for k in kp_order if k in pkh.index]
    kao = [k for k in ka_order if k in pkh.index]
    hbp, hba = boot_pair(resh, cl, kpo, kao)
    rows.append(stratum_summary(pkh, hbp, hba, kpo, kao, set(pkh.index),
                                {**base, "level": "hub_excluded", "stratum": "hub_excluded", "value": "all"}))
    for v in ["ST", "Y", "STY"]:
        ks = set(ktab.kinase_acc[ktab.residue_class.eq(v)])
        rows.append(stratum_summary(pkh, hbp, hba, kpo, kao, ks, {**base, "level": "hub_excluded",
                                                                   "stratum": "hub_excluded:residue_class",
                                                                   "value": v}))
    t_hub = time.time() - t0 - t_load - t_eval - t_kin - t_site
    # ---- equivalence check (R2 only): subset-label evaluate() vs post-hoc subset, primed (curated) + hub=yes
    if check_equiv:
        eq = []
        for (s, v) in [("primed", "primed"), ("hub", "hub")]:
            kin_ok, attr = site_attr[(s, v)]
            pmx = data.pos_meta[data.pos_meta.kinase_acc.isin(kin_ok)].copy()
            pmx = pmx[[bool(attr[k][np.searchsorted(data.positives[k], si)] == v)
                       for k, si in zip(pmx.kinase_acc, pmx.site_idx)]].reset_index(drop=True)
            posx = {k: np.sort(g.site_idx.values) for k, g in pmx.groupby("kinase_acc")}
            dx = dataclasses.replace(data, positives=posx, pos_meta=pmx)
            kx = np.where(np.isin(data.kinases.uniprot_acc.values, list(posx)))[0]
            Px = Pn[np.isin(Pn.kinase_idx.values, kx) | np.isin(Pn.site_idx.values, pmx.site_idx.values)]
            rx = E.evaluate(Px, dx, missing="worst", fold_normalise=False, keep_items=False).per_kinase
            rx = rx.set_index("kinase_acc")
            keep = {k: attr[k] == v for k in posx}
            sp = subset_result(res, data, keep).per_kinase.set_index("kinase_acc")
            j = rx.index.intersection(sp.index)
            eq.append({"stratum": f"{s}={v}", "n_kinases": len(j),
                       "max_abs_diff_auroc_hard": float(np.nanmax(np.abs(rx.loc[j, "auroc_hard"] - sp.loc[j, "auroc_hard"]))),
                       "max_abs_diff_hit5": float(np.nanmax(np.abs(rx.loc[j, "hit5"] - sp.loc[j, "hit5"]))),
                       "max_abs_diff_top1": float(np.nanmax(np.abs(rx.loc[j, "top1_recall"] - sp.loc[j, "top1_recall"]))),
                       "macro_diff_top1": float(np.nanmean(rx.loc[j, "top1_recall"] - sp.loc[j, "top1_recall"]))})
        sanity["equivalence"] = eq
    sanity["runtime_s"] = {"load_normalise": round(t_load, 1), "evaluate": round(t_eval, 1),
                           "kinase_strata_boot": round(t_kin, 1), "site_strata_boot": round(t_site, 1),
                           "hub_excluded": round(t_hub, 1), "total": round(time.time() - t0, 1)}
    WORK.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(WORK / f"{tag}.strata.tsv", sep="\t", index=False)
    if site_pk:
        pd.concat(site_pk).to_csv(WORK / f"{tag}.site_pk.tsv", sep="\t", index=False)
    # per-kinase table of the full evaluation (for partition checks / priming tables)
    res.per_kinase.drop(columns=[c for c in res.per_kinase if c.startswith("auroc_hard_f")]).to_csv(
        WORK / f"{tag}.per_kinase.tsv", sep="\t", index=False)
    (WORK / f"{tag}.sanity.json").write_text(json.dumps(sanity, indent=1, default=E._json_default))
    print(f"{tag}: M1 {sanity['M1_here']:.4f} (leaderboard {sanity.get('M1_leaderboard')}) "
          f"runtime {sanity['runtime_s']}", flush=True)


def site_attributes(data, ktab):
    """{(stratum, value): (kinase set, {kinase_acc: array of stratum labels over data.positives[k]})}."""
    S = data.sites
    sid = S.site_id.values
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
    ev = dict(zip(zip(pm.kinase_acc, pm.site_idx), ST.evidence_class(pm)))
    dev_all = labels.load_dev("all")
    hp = dev_all[labels.evaluable_positive(dev_all)]
    hubs = ST.hub_sites(hp, data.par)
    hub_idx = np.where(np.isin(sid, list(hubs)))[0]
    hubm = np.zeros(len(sid), bool)
    hubm[hub_idx] = True
    cur = ST.priming_curated(data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY_")])
    can = dict(zip(ktab.kinase_acc, ktab.priming_cantley_offsets))
    can = {k: tuple(int(o) for o in v.split(",")) for k, v in can.items() if v}
    st_nonpriming = set(ktab.kinase_acc[ktab.residue_class.isin(["ST", "STY"]) & ~ktab.priming & ~ktab.priming_cantley])
    ctrl = {k: (-3, 4) for k in st_nonpriming}
    observed = set(sid[S.in_U_ms.values.astype(bool)])
    real = [k for k in data.positives if not k.startswith("CANARY_")]
    out = {}
    evd, hubd, prd, prc, prx = {}, {}, {}, {}, {}
    for k in real:
        pos = data.positives[k]
        evd[k] = np.array([ev.get((k, s), "other") for s in pos], dtype=object)
        hubd[k] = np.where(hubm[pos], "hub", "non_hub").astype(object)
        q = pd.DataFrame({"kinase_acc": k, "site_id": sid[pos]})
        if k in cur:
            prd[k] = np.where(ST.primed_mask(q, cur, S, observed), "primed", "unprimed").astype(object)
        if k in can:
            prc[k] = np.where(ST.primed_mask(q, can, S, observed), "primed", "unprimed").astype(object)
        if k in ctrl:
            prx[k] = np.where(ST.primed_mask(q, ctrl, S, observed), "primed", "unprimed").astype(object)
    for v in ST.EVIDENCE_CLASSES:
        out[("evidence", v)] = (set(evd), evd)
    for v in ["hub", "non_hub"]:
        out[("hub", v)] = (set(hubd), hubd)
    for v in ["primed", "unprimed"]:
        out[("primed", v)] = (set(prd), prd)
        out[("primed_cantley", v)] = (set(prc), prc)
        out[("primed_ctrl_nonpriming_ST", v)] = (set(prx), prx)
    return out, hub_idx, {"n_hub_sites_dev": len(hubs), "n_hub_sites_universe": int(hubm[data.universe].sum())}


def worker(tags):
    t0 = time.time()
    data = E.prepare("cv", fold=None)
    # canaries are never reported here: drop their positives so evaluate() does not loop over them
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].reset_index(drop=True)
    data = dataclasses.replace(data, pos_meta=pm,
                               positives={k: v for k, v in data.positives.items() if not k.startswith("CANARY_")})
    ktab = ST.kinase_table()
    site_attr, hub_idx, hub_info = site_attributes(data, ktab)
    print(f"prepare + strata attributes {time.time() - t0:.0f}s; {hub_info}", flush=True)
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "hub_info.json").write_text(json.dumps(hub_info))
    for tag in tags:
        try:
            run_one(tag, data, ktab, site_attr, hub_idx, check_equiv=(tag == "R2"))
        except Exception as ex:   # reported, never silent
            import traceback
            traceback.print_exc()
            print(f"FAILED {tag}: {type(ex).__name__}: {ex}", flush=True)


# ============================================================================================== combine + sanity
PARTITIONS_KIN = ["manning_group", "residue_class", "train_bin", "priming", "priming_cantley", "pikk", "secretory"]


def combine():
    OUT.mkdir(parents=True, exist_ok=True)
    tabs, pks, sans = [], [], []
    for tag in RUNSPEC:
        f = WORK / f"{tag}.strata.tsv"
        if not f.exists():
            print("missing", tag)
            continue
        tabs.append(pd.read_csv(f, sep="\t"))
        if (WORK / f"{tag}.site_pk.tsv").exists():
            pks.append(pd.read_csv(WORK / f"{tag}.site_pk.tsv", sep="\t"))
        sans.append(json.loads((WORK / f"{tag}.sanity.json").read_text()))
    T = pd.concat(tabs, ignore_index=True)
    T.to_csv(OUT / "cv_strata.tsv", sep="\t", index=False, float_format="%.4f")
    SP = pd.concat(pks, ignore_index=True)
    SP.to_csv(OUT / "cv_strata_site_per_kinase.tsv", sep="\t", index=False, float_format="%.4f")
    chk = []
    for s in sans:
        tag = s["tag"]
        x = T[T.tag.eq(tag)]
        allr = x[x.stratum.eq("all")].iloc[0]
        r = {"tag": tag, "check": "full-run M1/M2/M3 = leaderboard metrics_a24.json"}
        if "M1_leaderboard" in s:
            r["max_abs_diff"] = max(abs(s["M1_here"] - s["M1_leaderboard"]), abs(s["M2_here"] - s["M2_leaderboard"]),
                                    abs((s["M3_here"] or 0) - (s["M3_leaderboard"] or 0)))
            chk.append(r)
        if "boot_M1_lo_eval" in s:
            chk.append({"tag": tag, "check": "stratum 'all' M1 CI = ksp.eval.bootstrap M1 CI",
                        "max_abs_diff": max(abs(allr.M1_lo - s["boot_M1_lo_eval"]),
                                            abs(allr.M1_hi - s["boot_M1_hi_eval"]))})
        # kinase partitions recombine: sum_g n_g M_g / sum n_g == M_all
        for part in PARTITIONS_KIN:
            g = x[x.stratum.eq(part)]
            for m, n in [("M1", "n_kinases_Kcv"), ("M5", "n_kinases_Kany")]:
                gg = g[g[n] > 0]
                rec = (gg[m] * gg[n]).sum() / gg[n].sum()
                chk.append({"tag": tag, "check": f"kinase partition {part}: n-weighted {m} recombines",
                            "max_abs_diff": abs(rec - allr[m])})
        # site partitions recombine per kinase (M1, M2 exact; M3 exact under the fixed cut-off)
        pkf = pd.read_csv(WORK / f"{tag}.per_kinase.tsv", sep="\t").set_index("kinase_acc")
        sp = SP[SP.tag.eq(tag)]
        for st_ in ["evidence", "hub", "primed", "primed_cantley"]:
            g = sp[sp.stratum.eq(st_)]
            if not len(g):
                continue
            for col in ["auroc_hard", "hit5", "top1_recall"]:
                w = g.dropna(subset=[col])
                if col == "hit5":
                    w = w.assign(nn=w.n_pos)
                rec = (w[col] * w.n_pos).groupby(w.kinase_acc).sum() / w.n_pos.groupby(w.kinase_acc).sum()
                full = pkf.loc[rec.index, col]
                # only kinases whose positives are fully covered by the partition's defined values
                ntot = g.groupby("kinase_acc").n_pos.sum()
                okk = ntot[ntot.eq(pkf.loc[ntot.index, "n_pos"])].index.intersection(rec.index)
                okk = [k for k in okk if g[g.kinase_acc.eq(k)][col].notna().all()]
                d = (rec.loc[okk] - full.loc[okk]).abs()
                chk.append({"tag": tag, "check": f"site partition {st_}: per-kinase n-weighted {col} recombines "
                            f"({len(okk)} kinases)", "max_abs_diff": float(d.max()) if len(d) else np.nan})
        for e in s.get("equivalence", []):
            for c in ["max_abs_diff_auroc_hard", "max_abs_diff_hit5", "max_abs_diff_top1", "macro_diff_top1"]:
                chk.append({"tag": tag, "check": f"post-hoc vs subset-label evaluate(), {e['stratum']}: {c} "
                            f"({e['n_kinases']} kinases)", "max_abs_diff": e[c]})
        if tag.startswith("B0"):
            ks = x[(x.level.isin(["all", "kinase"])) & (x.n_kinases_Kcv >= SMALL_K)]
            chk.append({"tag": tag, "check": f"B0 floor: max |M1-0.5| over non-small kinase strata ({len(ks)})",
                        "max_abs_diff": float((ks.M1 - 0.5).abs().max())})
            chk.append({"tag": tag, "check": "B0 floor: |M1(all)-0.5|", "max_abs_diff": float(abs(allr.M1 - 0.5))})
        chk.append({"tag": tag, "check": "runtime_s total", "max_abs_diff": s["runtime_s"]["total"]})
    C = pd.DataFrame(chk)
    C.to_csv(OUT / "cv_strata_sanity.tsv", sep="\t", index=False, float_format="%.6g")
    print(C.to_string())


# ============================================================================================== markdown tables
KEY_STRATA = [("all", "all", "all"), ("kinase", "residue_class", "ST"), ("kinase", "residue_class", "Y"),
              ("kinase", "residue_class", "STY"),
              ("kinase", "train_bin", "1-4"), ("kinase", "train_bin", "5-19"), ("kinase", "train_bin", "20-99"),
              ("kinase", "train_bin", ">=100"), ("kinase", "priming", "True"), ("kinase", "priming", "False"),
              ("kinase", "priming_cantley", "True"), ("site", "primed", "primed"), ("site", "primed", "unprimed"),
              ("site", "primed_cantley", "primed"), ("site", "primed_cantley", "unprimed"),
              ("site", "primed_ctrl_nonpriming_ST", "primed"), ("site", "primed_ctrl_nonpriming_ST", "unprimed"),
              ("kinase", "pikk", "True"), ("kinase", "pikk_sq", "True"), ("kinase", "fam20c", "True"),
              ("site", "evidence", "in_vitro_only"), ("site", "evidence", "in_vivo_only"), ("site", "evidence", "both"),
              ("site", "hub", "hub"), ("site", "hub", "non_hub"), ("hub_excluded", "hub_excluded", "all"),
              ("hub_excluded", "hub_excluded:residue_class", "Y"), ("hub_excluded", "hub_excluded:residue_class", "ST")]
TAB_TAGS = ["R2", "R3", "B1i", "B1ii", "B1prim", "B2", "B4", "B7", "B0prior", "B0pop", "B0miss"]


def tables():
    T = pd.read_csv(OUT / "cv_strata.tsv", sep="\t")
    SP = pd.read_csv(OUT / "cv_strata_site_per_kinase.tsv", sep="\t")
    f = lambda x: "" if pd.isna(x) else f"{x:.3f}"  # noqa: E731
    L = []
    for m, kn, pn in [("M1", "n_kinases_Kcv", "n_pos_Kcv"), ("M5", "n_kinases_Kany", "n_pos_Kany"),
                      ("M2", "n_kinases_Kcv", "n_pos_Kcv"), ("M3", "n_kinases_Kcv", "n_pos_Kcv")]:
        tags = [t for t in TAB_TAGS if t in set(T.tag)]
        L += [f"### {m} by stratum ({'K_any' if m == 'M5' else 'K_cv'}; point [95% cluster-bootstrap CI] for R2)", "",
              "| stratum | n kin | n pos | " + " | ".join(tags) + " |", "|---|---|---|" + "---|" * len(tags)]
        for lv, s, v in KEY_STRATA:
            x = T[T.level.eq(lv) & T.stratum.eq(s) & T.value.astype(str).eq(v)].set_index("tag")
            if "R2" not in x.index:
                continue
            r2 = x.loc["R2"]
            small = r2["small" if m != "M5" else "small_Kany"]
            cells = []
            for t in tags:
                if t not in x.index:
                    cells.append("")
                    continue
                c = f(x.loc[t, m])
                if t == "R2" and not pd.isna(x.loc[t, m]):
                    c += f" [{f(x.loc[t, m + '_lo'])}, {f(x.loc[t, m + '_hi'])}]"
                cells.append(c)
            name = f"{s}={v}" + (" (small)" if small else "")
            L.append(f"| {name} | {int(r2[kn])} | {int(r2[pn])} | " + " | ".join(cells) + " |")
        L.append("")
    # A27 matched-kinase view vs Cantley
    L += ["### A27 matched-kinase comparison (both restricted to kinases scored by B1 arm ii)", "",
          "| stratum | n kin (K_cv) | R2 M1 | B1ii M1 | R2 M2 | B1ii M2 | R2 M3 | B1ii M3 |", "|---|---|---|---|---|---|---|---|"]
    for lv, s, v in KEY_STRATA:
        x = T[T.level.eq(lv) & T.stratum.eq(s) & T.value.astype(str).eq(v)].set_index("tag")
        if not {"R2_kmB1", "B1ii_kmB1"} <= set(x.index):
            continue
        a, b = x.loc["R2_kmB1"], x.loc["B1ii_kmB1"]
        L.append(f"| {s}={v} | {int(a.n_kinases_Kcv)} | {f(a.M1)} | {f(b.M1)} | {f(a.M2)} | {f(b.M2)} | {f(a.M3)} | "
                 f"{f(b.M3)} |")
    L.append("")
    # paired primed - unprimed within priming kinases (>= 3 of each)
    L += ["### Within-kinase primed vs unprimed (kinases with >= 3 positives of each; E-hard AUROC)", "",
          "| run | stratum | n kinases | mean AUROC primed | mean AUROC unprimed | mean diff | kinases with diff < 0 |",
          "|---|---|---|---|---|---|---|"]
    for t in [x for x in TAB_TAGS if x in set(SP.tag)]:
        for s in ["primed", "primed_cantley", "primed_ctrl_nonpriming_ST"]:
            g = SP[SP.tag.eq(t) & SP.stratum.eq(s)].pivot_table(index="kinase_acc", columns="value",
                                                                 values=["auroc_hard", "n_pos"])
            if not len(g) or ("n_pos", "primed") not in g:
                continue
            ok = g[(g[("n_pos", "primed")] >= 3) & (g[("n_pos", "unprimed")] >= 3)]
            d = ok[("auroc_hard", "primed")] - ok[("auroc_hard", "unprimed")]
            L.append(f"| {t} | {s} | {len(ok)} | {f(ok[('auroc_hard', 'primed')].mean())} | "
                     f"{f(ok[('auroc_hard', 'unprimed')].mean())} | {f(d.mean())} | {int((d < 0).sum())} |")
    L.append("")
    # per-kinase table for curated priming kinases
    pk = {t: pd.read_csv(WORK / f"{t}.per_kinase.tsv", sep="\t").set_index("kinase_acc")
          for t in ["R2", "R3", "B1i", "B1prim", "B7"] if (WORK / f"{t}.per_kinase.tsv").exists()}
    sp = SP[SP.stratum.eq("primed")]
    L += ["### Curated priming kinases, per kinase (E-hard AUROC; primed / unprimed positives)", "",
          "| kinase | n pos (primed/unprimed) | " + " | ".join(f"{t} all / primed / unprimed" for t in pk) + " |",
          "|---|---|" + "---|" * len(pk)]
    for k in sorted(set(sp.kinase_acc), key=lambda a: -pk["R2"].loc[a, "n_pos"]):
        g = sp[sp.kinase_acc.eq(k)]
        n = {v: int(g[g.tag.eq("R2") & g.value.eq(v)].n_pos.sum()) for v in ["primed", "unprimed"]}
        cells = []
        for t, P in pk.items():
            a = {v: g[g.tag.eq(t) & g.value.eq(v)].auroc_hard for v in ["primed", "unprimed"]}
            cells.append(f"{f(P.loc[k, 'auroc_hard']) if k in P.index else ''} / "
                         f"{f(a['primed'].iloc[0]) if len(a['primed']) else '-'} / "
                         f"{f(a['unprimed'].iloc[0]) if len(a['unprimed']) else '-'}")
        L.append(f"| {pk['R2'].loc[k, 'gene']} | {n['primed']}/{n['unprimed']} | " + " | ".join(cells) + " |")
    L.append("")
    # PIKK / FAM20C per kinase
    L += ["### PIKK and FAM20C, per kinase (E-hard AUROC / Hit@5 / top-1% recall)", "",
          "| kinase | n pos | " + " | ".join(pk) + " |", "|---|---|" + "---|" * len(pk)]
    for g in ["ATM", "ATR", "PRKDC", "SMG1", "MTOR", "FAM20C"]:
        P0 = pk["R2"]
        k = P0.index[P0.gene.eq(g)]
        if not len(k):
            continue
        k = k[0]
        cells = [f"{f(P.loc[k, 'auroc_hard'])} / {f(P.loc[k, 'hit5'])} / {f(P.loc[k, 'top1_recall'])}"
                 if k in P.index else "" for P in pk.values()]
        L.append(f"| {g} | {int(P0.loc[k, 'n_pos'])} | " + " | ".join(cells) + " |")
    (OUT / "cv_strata_tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", default=None)
    ap.add_argument("--combine", action="store_true")
    a = ap.parse_args()
    if a.combine:
        combine()
        tables()
    else:
        worker([t for t in a.tags.split(",") if t])
