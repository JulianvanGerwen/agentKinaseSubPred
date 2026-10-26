"""ksp.eval - THE ONLY metric code (DESIGN §6, §11.1, §16 A18-A20). Owner: core-infra.

    data = prepare('cv')                       # DEV OOF evaluation (labels via ksp.labels.load_dev)
    res  = evaluate(preds, data)               # preds: DataFrame[kinase_acc, site_id, score(, score_site)]
    res.endpoints, res.per_kinase, res.strata(), res.missing, res.canary
    bootstrap(res, B=1000) / paired_bootstrap(resA, resB, B) / holm(pvals) / matched_universe(pA, pB)

TEST evaluation is only reachable through ksp.eval_test (it passes the sealed labels to prepare('test', ...)).

Definitions (frozen at G1)
  Universe      : U_ms sites on the evaluated split's proteins (A6). For split='cv', fold=None -> all DEV proteins
                  (pooled out-of-fold predictions); fold=f -> proteins of CV fold f.
  Positives(k)  : human-native PSP pairs of k, in class, not autophosphorylation, site in U_ms and in the universe.
  Masks(k)      : ksp.labels.NegativeMasks: KP(k) (all sources/species, identical 15-mers), Hamming<=2 15-mers of KP(k)
                  sites, sites labelled to Par(k). Mask labels: CV = all DEV rows (+canaries); TEST = DEV U TEST rows.
  E-all(k)      : universe sites of k's residue class, minus masks (and minus positives) = negatives.
  E-hard(k)     : E-all(k) sites carrying >=1 human-native non-autophos PSP label to a kinase outside {k} U Par(k).
  M1  : macro over K_eval of AUROC(E-hard)            M4 : macro over K_any of AUROC(E-hard)
  M2  : site-centric filtered Hit@5 (score_site if given, else score). Candidates for true pair (k,s) = real kinases
        whose class covers res(s), minus the other kinases labelled to s (any source) and minus Par(k) (A18).
        Expected value under random tie-breaking. Averaged within kinase, then macro over K_eval.
        Secondary: Hit@1, Hit@10, MRR (1/expected rank), family-level Hit@5 (Manning family; family score = max).
  M3  : macro over K_eval of recall of positives within the top 1% (ceil(0.01 * |E-all U P|)) of E-all U P;
        ties at the cut-off get fractional credit. Secondary: pAUROC(FPR<=0.01, McClish-standardised, as
        sklearn max_fpr), AUPRC (average precision) and its lift (AP / prevalence), AUROC on E-all.
  CV  : the same endpoints with K_cv in place of K_eval (M1-CV == M5) and K_any_cv (>=1 DEV positive) for M4.
  Missing predictions: policy 'worst' (default; missing scores = -inf, counted), 'drop' (only scored pairs; use for
  matched-universe comparisons), or 'error'. Unknown kinase/site IDs and duplicate pairs always raise.
  Canaries (kinase_acc CANARY_*) are evaluated in res.canary only, never in endpoints/strata.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import canary as _canary
from . import index, labels, splits

PRIMING = {"GSK3A", "GSK3B", "CSNK2A1", "CSNK2A2", "CSNK2A3", "PLK1"}   # + CK1 family (added in kinase_strata)
TRAIN_BINS = [-1, 0, 4, 19, 99, 1e12]
TRAIN_BIN_LABELS = ["0", "1-4", "5-19", "20-99", ">=100"]
NEG_CANARY_MAX = 0.57
POS_CANARY_MIN = 0.90


# ================================================================================================ metric primitives
def _groups_desc(scores: np.ndarray):
    """Sort descending; return order, unique-score group starts."""
    o = np.argsort(-scores, kind="mergesort")
    s = scores[o]
    starts = np.r_[0, np.where(np.diff(s) != 0)[0] + 1] if len(s) else np.zeros(0, int)
    return o, starts


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUROC, ties count 1/2. NaN if either side empty."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.r_[pos, neg]
    from scipy.stats import rankdata
    r = rankdata(allv)
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def roc_points(pos, neg):
    """ROC vertices (fpr, tpr) over unique-score thresholds (descending), starting at (0,0)."""
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[np.asarray(pos, float), np.asarray(neg, float)]
    o, st = _groups_desc(s)
    yy = y[o]
    ends = np.r_[st[1:], len(s)] - 1
    tp = np.cumsum(yy)[ends]
    fp = np.cumsum(1 - yy)[ends]
    return np.r_[0, fp / max(len(neg), 1)], np.r_[0, tp / max(len(pos), 1)]


def pauroc(pos, neg, max_fpr: float = 0.01, standardise: bool = True) -> float:
    """Partial AUROC up to max_fpr (linear interpolation at max_fpr), McClish-standardised like sklearn
    roc_auc_score(max_fpr=...): 0.5 = random, 1 = perfect."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    fpr, tpr = roc_points(pos, neg)
    stop = np.searchsorted(fpr, max_fpr, "right")
    x, y = fpr[:stop], tpr[:stop]
    if stop < len(fpr):
        x = np.r_[x, max_fpr]
        y = np.r_[y, np.interp(max_fpr, [fpr[stop - 1], fpr[stop]], [tpr[stop - 1], tpr[stop]])]
    area = float(np.trapezoid(y, x)) if hasattr(np, "trapezoid") else float(np.trapz(y, x))
    if not standardise:
        return area
    mn, mx = max_fpr ** 2 / 2, max_fpr
    return 0.5 * (1 + (area - mn) / (mx - mn))


def average_precision(pos, neg) -> float:
    """AP = sum_t (R_t - R_{t-1}) P_t over unique-score thresholds (sklearn definition)."""
    if len(pos) == 0:
        return float("nan")
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[np.asarray(pos, float), np.asarray(neg, float)]
    o, st = _groups_desc(s)
    yy = y[o]
    ends = np.r_[st[1:], len(s)] - 1
    tp = np.cumsum(yy)[ends]
    n = ends + 1
    prec = tp / n
    rec = tp / len(pos)
    return float(np.sum(np.diff(np.r_[0, rec]) * prec))


def topfrac_hits(pos, neg, frac: float = 0.01) -> np.ndarray:
    """Per-positive expected hit (0..1) of being within the top ceil(frac*N) of pos U neg (fractional tie credit)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    N = len(pos) + len(neg)
    if len(pos) == 0:
        return np.zeros(0)
    n = max(1, math.ceil(frac * N))
    allv = np.r_[pos, neg]
    t = np.sort(allv)[::-1][n - 1]
    above = (allv > t).sum()
    tied = (allv == t).sum()
    credit = (n - above) / tied
    return np.where(pos > t, 1.0, np.where(pos == t, credit, 0.0))


def topfrac_recall(pos, neg, frac: float = 0.01) -> float:
    h = topfrac_hits(pos, neg, frac)
    return float(h.mean()) if len(h) else float("nan")


def expected_hit_at(n_greater: np.ndarray, n_equal: np.ndarray, K: int) -> np.ndarray:
    """P(rank <= K) when the true item ties with n_equal others (uniform tie-breaking); rank = 1 + n_greater + U."""
    n_greater, n_equal = np.asarray(n_greater, float), np.asarray(n_equal, float)
    return np.clip((K - n_greater) / (n_equal + 1), 0, 1)


def weighted_auroc_matrix(pos_s, neg_s, Wp: np.ndarray, Wn: np.ndarray) -> np.ndarray:
    """AUROC for B weight vectors at once (bootstrap). Wp: B x n_pos, Wn: B x n_neg. Ties count 1/2."""
    s = np.r_[pos_s, neg_s]
    o, st = _groups_desc(-s)        # ascending order groups
    W = np.concatenate([Wp, np.zeros_like(Wn)], 1)[:, o]
    V = np.concatenate([np.zeros_like(Wp), Wn], 1)[:, o]
    gp = np.add.reduceat(W, st, axis=1)
    gn = np.add.reduceat(V, st, axis=1)
    below = np.cumsum(gn, 1) - gn
    num = (gp * (below + 0.5 * gn)).sum(1)
    den = Wp.sum(1) * Wn.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


# ================================================================================================ evaluation data
@dataclass
class EvalData:
    split: str
    fold: int | None
    sites: pd.DataFrame                       # canonical index sites
    universe: np.ndarray                      # bool over index sites
    site_cluster: np.ndarray                  # cluster id per index site (-1 = no cluster)
    kinases: pd.DataFrame                     # real kinases + canaries (index order, canaries appended)
    positives: dict                           # kinase_acc -> np.ndarray of index positions
    pos_meta: pd.DataFrame                    # one row per positive pair: kinase_acc, site_idx, strata flags
    ext_positives: dict                       # kinase_acc -> index positions of external-curated test pairs
    masks: labels.NegativeMasks
    par: dict
    site_kinases: dict                        # site_idx -> set of kinase_acc labelled (any source) to the site
    kinase_sets: pd.DataFrame                 # eval_kinase_sets.tsv
    k_primary: list                           # K_eval (test) / K_cv (cv)
    k_any: list                               # K_any (test) / kinases with >=1 positive in universe (cv)
    train_counts: pd.Series                   # kinase_acc -> DEV L2 count
    info: dict = field(default_factory=dict)


def _res_ok(res: np.ndarray, cls: str) -> np.ndarray:
    if cls == "STY":
        return np.ones(len(res), bool)
    if cls == "ST":
        return (res == "S") | (res == "T")
    if cls == "Y":
        return res == "Y"
    raise ValueError(f"unknown residue class {cls!r}")


def prepare(split: str = "cv", fold: int | None = None, labels_df: pd.DataFrame | None = None,
            with_canaries: bool = True) -> EvalData:
    """Build the evaluation data for split 'cv' (DEV; labels from ksp.labels.load_dev) or 'test' (labels_df = the
    sealed TEST label rows, passed only by ksp.eval_test)."""
    S = index.load_sites()
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    clmap = splits.load_clusters()
    s_cl = S.acc.map(clmap).fillna(-1).astype(np.int64).values
    ks = splits.load_eval_kinase_sets()
    dev_all = labels.load_dev("all")
    if with_canaries:
        dev_all = pd.concat([dev_all, _canary.canary_labels()], ignore_index=True)
    test_prot = splits.load_test_proteins()
    if split == "cv":
        folds = splits.load_cv_folds()
        prot = set(folds.index) if fold is None else set(folds.index[folds == fold])
        lab = dev_all
        mask_pairs = dev_all
        k_primary = ks.kinase_acc[ks.in_K_cv].tolist()
    elif split == "test":
        if labels_df is None:
            raise ValueError("split='test' requires labels_df (only ksp.eval_test may supply sealed labels)")
        prot = test_prot
        lab = labels_df if with_canaries else labels_df[~labels_df.get("is_canary", False).astype(bool)]
        mask_pairs = pd.concat([dev_all, lab], ignore_index=True)
        k_primary = ks.kinase_acc[ks.in_K_eval].tolist()
    else:
        raise ValueError(f"split must be 'cv' or 'test', got {split!r}")
    universe = S.acc.isin(prot).values & S.in_U_ms.values
    lab = lab[lab.acc.isin(prot)]
    ev = lab[labels.evaluable_positive(lab) & lab.site_id.isin(spos.index)].copy()
    ev["site_idx"] = spos[ev.site_id].values
    ev = ev[universe[ev.site_idx.values]]
    ev = ev.drop_duplicates(["kinase_acc", "site_id"])
    positives = {k: np.sort(g.site_idx.values) for k, g in ev.groupby("kinase_acc")}
    meta_cols = [c for c in ["kinase_acc", "site_idx", "evidence_native", "new_since_2021"] if c in ev]
    pos_meta = ev[meta_cols].reset_index(drop=True)
    ext = lab[lab.get("ext_test", pd.Series(False, index=lab.index)).astype(bool) & lab.site_id.isin(spos.index)]
    ext = ext.assign(site_idx=lambda d: spos[d.site_id].values)
    ext = ext[universe[ext.site_idx.values]]
    ext_pos = {k: np.unique(g.site_idx.values) for k, g in ext.groupby("kinase_acc")}
    par = labels.paralog_sets()
    masks = labels.NegativeMasks(mask_pairs, paralogs=par)
    mp = mask_pairs[mask_pairs.site_id.isin(spos.index)]
    sk = mp.assign(i=spos[mp.site_id].values).groupby("i").kinase_acc.agg(set).to_dict()
    K = _canary.kinases_with_canaries() if with_canaries else index.load_kinases()
    if split == "cv":
        k_any = sorted(k for k in positives if not k.startswith("CANARY_"))
    else:
        k_any = ks.kinase_acc[ks.in_K_any].tolist()
    tc = pd.Series(ks.n_dev_L2.values, index=ks.kinase_acc.values)
    info = {"n_universe_sites": int(universe.sum()), "n_positive_pairs": int(len(ev)),
            "n_positive_pairs_real": int((~ev.kinase_acc.str.startswith("CANARY_")).sum()),
            "n_ext_pairs": int(len(ext))}
    return EvalData(split, fold, S, universe, s_cl, K, positives, pos_meta, ext_pos, masks, par, sk, ks,
                    k_primary, k_any, tc, info)


# ================================================================================================ predictions
def map_predictions(preds: pd.DataFrame, data: EvalData, score_col: str = "score") -> pd.DataFrame:
    """Validate and map predictions to (k_i, s_i, score[, score_site]). Raises on unknown IDs, duplicates, NaN.
    Either string IDs (kinase_acc, site_id) or integer positions (kinase_idx into kinases_with_canaries() order =
    cache/index/kinases.parquet + canaries appended; site_idx into cache/index/sites.parquet) - the latter is the
    efficient form for pooled-CV E-all (~10^8 rows)."""
    if {"kinase_idx", "site_idx"} <= set(preds.columns):
        ki = pd.Series(preds.kinase_idx.values)
        si = pd.Series(preds.site_idx.values)
        if (ki < 0).any() or (ki >= len(data.kinases)).any() or not np.issubdtype(ki.dtype, np.integer):
            raise KeyError("kinase_idx out of range of the kinase index (+canaries)")
        if (si < 0).any() or (si >= len(data.sites)).any() or not np.issubdtype(si.dtype, np.integer):
            raise KeyError("site_idx out of range of the site index")
    else:
        need = {"kinase_acc", "site_id"}
        if not need <= set(preds.columns):
            raise ValueError(f"predictions need kinase_acc/site_id (or kinase_idx/site_idx); got {list(preds.columns)}")
        kpos = pd.Series(np.arange(len(data.kinases)), index=data.kinases.uniprot_acc.values)
        ki = preds.kinase_acc.map(kpos)
        if ki.isna().any():
            bad = preds.kinase_acc[ki.isna()].unique()[:10]
            raise KeyError(f"{ki.isna().sum()} prediction rows with unknown kinase_acc, e.g. {list(bad)}")
        spos = pd.Series(np.arange(len(data.sites)), index=data.sites.site_id.values)
        si = preds.site_id.map(spos)
        if si.isna().any():
            bad = preds.site_id[si.isna()].unique()[:10]
            raise KeyError(f"{si.isna().sum()} prediction rows with unknown site_id, e.g. {list(bad)}")
    if score_col not in preds.columns:
        raise ValueError(f"predictions need a {score_col!r} column")
    out = pd.DataFrame({"k": ki.values.astype(np.int64), "s": si.values.astype(np.int64),
                        "score": preds[score_col].values.astype(np.float64)})
    if "score_site" in preds.columns:
        out["score_site"] = preds.score_site.values.astype(np.float64)
    key = out.k.values * np.int64(len(data.sites)) + out.s.values
    if len(np.unique(key)) != len(key):
        raise ValueError(f"{len(key) - len(np.unique(key))} duplicate (kinase, site) prediction rows")
    if not np.isfinite(out.score).all():
        raise ValueError("non-finite scores in predictions (encode 'not scored' by omitting the row)")
    return out


# ================================================================================================ evaluation
@dataclass
class EvalResult:
    split: str
    fold: int | None
    endpoints: dict
    per_kinase: pd.DataFrame
    site_level: pd.DataFrame                  # one row per (kinase, positive site) for M2
    missing: dict
    canary: dict
    items: dict                               # kinase -> dict of arrays for bootstrap
    pair_strata: pd.DataFrame
    data_info: dict

    def to_json(self) -> dict:
        return {"split": self.split, "fold": self.fold, "endpoints": self.endpoints, "missing": self.missing,
                "canary": self.canary, "data_info": self.data_info,
                "strata": strata(self).to_dict(orient="records"),
                "pair_strata": self.pair_strata.to_dict(orient="records")}

    def save_json(self, path) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_json(), fh, indent=1, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _nanmean(x) -> float:
    x = np.asarray(list(x), float)
    x = x[~np.isnan(x)]
    return float(x.mean()) if len(x) else float("nan")


def evaluate(preds: pd.DataFrame, data: EvalData, missing: str = "worst", frac_top: float = 0.01,
             kinases: list | None = None, keep_items: bool = True) -> EvalResult:
    """All endpoints for one prediction set. `kinases` restricts the per-kinase loop (default: every kinase with
    >= 1 positive in the universe, incl. canaries)."""
    if missing not in ("worst", "drop", "error"):
        raise ValueError("missing must be 'worst', 'drop' or 'error'")
    P = map_predictions(preds, data)
    nS = len(data.sites)
    res_arr = data.sites.res.values
    kacc = data.kinases.uniprot_acc.values
    kcls = data.kinases.residue_class.values
    kpos = {a: i for i, a in enumerate(kacc)}
    P = P.sort_values("k", kind="stable")
    kb = np.searchsorted(P.k.values, np.arange(len(kacc) + 1))
    todo = kinases if kinases is not None else sorted(data.positives)
    buf = np.full(nS, np.nan)
    rows, items, miss = [], {}, {"pairs_required": 0, "pairs_missing": 0, "per_kinase": {}}
    pair_rows = []
    pm = data.pos_meta.copy()
    for k in todo:
        if k not in kpos:
            raise KeyError(f"kinase {k} not in kinase index")
        i = kpos[k]
        pos = data.positives.get(k, np.zeros(0, np.int64))
        sl = slice(kb[i], kb[i + 1])
        buf[:] = np.nan
        buf[P.s.values[sl]] = P.score.values[sl]
        cand = data.universe & _res_ok(res_arr, kcls[i])
        excl = np.zeros(nS, bool)
        excl[data.masks.excluded(k)] = True
        posm = np.zeros(nS, bool)
        posm[pos] = True
        neg_all = cand & ~excl & ~posm
        hard = np.zeros(nS, bool)
        hard[data.masks.hard_candidates(k)] = True
        neg_hard = neg_all & hard
        req = posm | neg_all
        nmiss = int(np.isnan(buf[req]).sum())
        miss["pairs_required"] += int(req.sum())
        miss["pairs_missing"] += nmiss
        miss["per_kinase"][k] = nmiss
        if nmiss and missing == "error":
            raise ValueError(f"kinase {k}: {nmiss} required (kinase, site) pairs have no prediction")
        sc = buf.copy()
        if missing == "worst":
            sc[np.isnan(sc)] = -np.inf
        ip = pos[~np.isnan(sc[pos])]
        ih = np.where(neg_hard & ~np.isnan(sc))[0]
        ia = np.where(neg_all & ~np.isnan(sc))[0]
        sp, sh, sa = sc[ip], sc[ih], sc[ia]
        # infinities -> finite floor below everything (keeps rank-based metrics exact)
        fl = np.nanmin(np.r_[sp[np.isfinite(sp)], sa[np.isfinite(sa)], 0.0]) - 1.0
        sp, sh, sa = (np.where(np.isfinite(x), x, fl) for x in (sp, sh, sa))
        hits = topfrac_hits(sp, sa, frac_top)
        r = {"kinase_acc": k, "n_pos": len(ip), "n_pos_required": len(pos), "n_neg_hard": len(ih),
             "n_neg_all": len(ia), "n_missing": nmiss,
             "auroc_hard": auroc(sp, sh), "auroc_all": auroc(sp, sa), "pauroc01_all": pauroc(sp, sa, 0.01),
             "ap_all": average_precision(sp, sa), "top1_recall": float(hits.mean()) if len(hits) else np.nan,
             "n_pos_clusters": len(np.unique(data.site_cluster[ip]))}
        prev = len(sp) / max(len(sp) + len(sa), 1)
        r["ap_lift_all"] = r["ap_all"] / prev if prev > 0 else np.nan
        rows.append(r)
        if keep_items:
            items[k] = {"pos_s": sp, "pos_c": data.site_cluster[ip], "hard_s": sh, "hard_c": data.site_cluster[ih],
                        "top_hits": hits, "all_idx_hash": hash(ia.tobytes()) ^ hash(ip.tobytes())}
        # pair-level strata (positives subsets vs the same E-hard negatives)
        sub = pm[pm.kinase_acc.eq(k)]
        sub_s = pd.Series(sc[sub.site_idx.values], index=sub.site_idx.values)
        strata_defs = []
        if "evidence_native" in sub:
            for ev in ["in_vivo_only", "in_vitro_only", "both"]:
                strata_defs.append((f"evidence={ev}", sub.site_idx.values[sub.evidence_native.eq(ev).values]))
        if "new_since_2021" in sub:
            nw = sub.new_since_2021.astype(bool).values
            strata_defs += [("temporal=new_since_2021", sub.site_idx.values[nw]),
                            ("temporal=pre_2021", sub.site_idx.values[~nw])]
        e = data.ext_positives.get(k)
        if e is not None:
            e = e[~posm[e]]
            strata_defs.append(("external_curated", e))
        for name, idx in strata_defs:
            if len(idx) == 0:
                continue
            v = sc[idx]
            v = v[~np.isnan(v)]
            v = np.where(np.isfinite(v), v, fl)
            pair_rows.append({"kinase_acc": k, "stratum": name, "n_pos": len(v), "auroc_hard": auroc(v, sh),
                              "auroc_all": auroc(v, sa), "top1_recall": topfrac_recall(v, sa, frac_top)})
        del sub_s
    pk = pd.DataFrame(rows)
    pk = _attach_kinase_attrs(pk, data)
    site_level = _site_centric(P, data, missing)
    if len(site_level):
        g = site_level.groupby("kinase_acc")
        agg = g.agg(hit1=("hit1", "mean"), hit5=("hit5", "mean"), hit10=("hit10", "mean"), mrr=("rr", "mean"),
                    fam_hit5=("fam_hit5", "mean"), n_site_pairs=("hit5", "size"))
        pk = pk.merge(agg, left_on="kinase_acc", right_index=True, how="left")
    is_can = pk.kinase_acc.str.startswith("CANARY_")
    real = pk[~is_can].set_index("kinase_acc")
    E = _endpoints(real, data)
    E["n_missing_pairs"] = miss["pairs_missing"]
    E["n_required_pairs"] = miss["pairs_required"]
    can = canary_check(pk[is_can]) if is_can.any() else {}
    ps = pd.DataFrame(pair_rows)
    ps_tab = _pair_strata_table(ps, data)
    miss["per_kinase"] = {k: v for k, v in miss["per_kinase"].items() if v}
    miss["site_centric_missing_scores"] = int(site_level.attrs.get("n_missing", 0)) if len(site_level) else 0
    return EvalResult(data.split, data.fold, E, pk, site_level, miss, can, items, ps_tab, data.info)


def _attach_kinase_attrs(pk: pd.DataFrame, data: EvalData) -> pd.DataFrame:
    K = data.kinases.set_index("uniprot_acc")
    pk = pk.copy()
    pk["gene"] = pk.kinase_acc.map(K.gene)
    pk["manning_group"] = pk.kinase_acc.map(K.manning_group).fillna("(none)")
    pk["manning_family"] = pk.kinase_acc.map(K.manning_family)
    pk["residue_class"] = pk.kinase_acc.map(K.residue_class)
    tc = pk.kinase_acc.map(data.train_counts).fillna(0)
    pk["n_train_dev_L2"] = tc.astype(int)
    pk["train_bin"] = pd.cut(tc, TRAIN_BINS, labels=TRAIN_BIN_LABELS).astype(str)
    pk["priming"] = pk.gene.isin(PRIMING) | pk.manning_family.eq("CK1")
    pk["pikk"] = pk.manning_family.eq("PIKK")
    pk["in_K_primary"] = pk.kinase_acc.isin(set(data.k_primary))
    pk["in_K_any"] = pk.kinase_acc.isin(set(data.k_any))
    return pk


def _endpoints(real: pd.DataFrame, data: EvalData) -> dict:
    kp = [k for k in data.k_primary if k in real.index]
    ka = [k for k in data.k_any if k in real.index]
    E = {"K_primary_name": "K_eval" if data.split == "test" else "K_cv",
         "n_K_primary": len(data.k_primary), "n_K_primary_evaluated": len(kp),
         "n_K_any": len(data.k_any), "n_K_any_evaluated": len(ka)}

    def m(col, ks_):
        if col not in real:
            return float("nan")
        return _nanmean(real.loc[ks_, col]) if ks_ else float("nan")
    E["M1_auroc_hard"] = m("auroc_hard", kp)
    E["M1_n_kinases_defined"] = int(real.loc[kp, "auroc_hard"].notna().sum()) if kp else 0
    E["M2_hit5"] = m("hit5", kp)
    E["M3_top1_recall"] = m("top1_recall", kp)
    E["M4_auroc_hard_Kany"] = m("auroc_hard", ka)
    if data.split == "cv":
        E["M5_auroc_hard_Kcv"] = E["M1_auroc_hard"]
    for c in ["hit1", "hit10", "mrr", "fam_hit5", "pauroc01_all", "ap_all", "ap_lift_all", "auroc_all"]:
        E[f"sec_{c}"] = m(c, kp)
    return E


def _site_centric(P: pd.DataFrame, data: EvalData, missing: str) -> pd.DataFrame:
    """Filtered site-centric ranks for every real positive pair (k, s) (A18, A20)."""
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
    if not len(pm):
        return pd.DataFrame()
    K = data.kinases
    real = ~K.uniprot_acc.str.startswith("CANARY_").values
    kacc = K.uniprot_acc.values
    kidx = {a: i for i, a in enumerate(kacc)}
    usites = np.unique(pm.site_idx.values)
    col = {s: j for j, s in enumerate(usites)}
    M = np.full((len(kacc), len(usites)), np.nan)
    sel = P[P.s.isin(usites)]
    scol = "score_site" if "score_site" in sel.columns else "score"
    M[sel.k.values, [col[s] for s in sel.s.values]] = sel[scol].values
    res = data.sites.res.values[usites]
    cls = K.residue_class.values
    okc = np.stack([_res_ok(res, c) for c in cls]) & real[:, None]          # kinases x sites
    fam = np.array([f if isinstance(f, str) and f and not f.endswith("Unique") else f"_{a}"
                    for f, a in zip(K.manning_family.values, kacc)])
    n_missing = int((np.isnan(M) & okc).sum())
    if missing == "worst":
        M = np.where(np.isnan(M), -np.inf, M)
    out = []
    for k, s in zip(pm.kinase_acc.values, pm.site_idx.values):
        i, j = kidx[k], col[s]
        cand = okc[:, j].copy()
        others = data.site_kinases.get(s, set()) | data.par.get(k, set())
        for o in others:
            if o in kidx and o != k:
                cand[kidx[o]] = False
        cand[i] = False
        v = M[cand, j]
        t = M[i, j]
        if np.isnan(t):          # only with missing='drop'
            continue
        v = v[~np.isnan(v)]
        ng, ne = float((v > t).sum()), float((v == t).sum())
        # family level: family score = max over remaining candidates (+ k for its own family)
        cf = fam[cand]
        fs = pd.Series(M[cand, j]).groupby(cf).max()
        own = fam[i]
        own_s = max(t, fs.get(own, -np.inf))
        fo = fs.drop(own, errors="ignore").values
        fo = fo[~np.isnan(fo)]
        fg, fe = float((fo > own_s).sum()), float((fo == own_s).sum())
        out.append((k, s, ng, ne, expected_hit_at(ng, ne, 1), expected_hit_at(ng, ne, 5),
                    expected_hit_at(ng, ne, 10), 1.0 / (1 + ng + ne / 2), expected_hit_at(fg, fe, 5),
                    int(data.site_cluster[s])))
    d = pd.DataFrame(out, columns=["kinase_acc", "site_idx", "n_greater", "n_equal", "hit1", "hit5", "hit10", "rr",
                                   "fam_hit5", "cluster"])
    d.attrs["n_missing"] = n_missing
    return d


def _pair_strata_table(ps: pd.DataFrame, data: EvalData) -> pd.DataFrame:
    if not len(ps):
        return pd.DataFrame(columns=["stratum", "kinase_set", "n_kinases", "n_pos", "auroc_hard", "auroc_all",
                                     "top1_recall"])
    ps = ps[~ps.kinase_acc.str.startswith("CANARY_")]
    out = []
    for kname, kset in [("K_primary", set(data.k_primary)), ("K_any", set(data.k_any))]:
        x = ps[ps.kinase_acc.isin(kset)]
        for st, g in x.groupby("stratum"):
            out.append({"stratum": st, "kinase_set": kname, "n_kinases": int(g.auroc_hard.notna().sum()),
                        "n_pos": int(g.n_pos.sum()), "auroc_hard": _nanmean(g.auroc_hard),
                        "auroc_all": _nanmean(g.auroc_all), "top1_recall": _nanmean(g.top1_recall)})
    return pd.DataFrame(out)


def strata(res: EvalResult, cols=("manning_group", "residue_class", "train_bin", "priming", "pikk"),
           kinase_set: str = "K_any") -> pd.DataFrame:
    """Kinase-level strata: macro means of per-kinase metrics within each stratum (DESIGN §6.4)."""
    pk = res.per_kinase[~res.per_kinase.kinase_acc.str.startswith("CANARY_")]
    pk = pk[pk.in_K_any] if kinase_set == "K_any" else pk[pk.in_K_primary]
    mcols = [c for c in ["auroc_hard", "hit5", "top1_recall", "auroc_all", "pauroc01_all"] if c in pk]
    out = []
    for c in cols:
        for v, g in pk.groupby(c):
            r = {"stratum": c, "value": str(v), "n_kinases": len(g)}
            r.update({m: _nanmean(g[m]) for m in mcols})
            out.append(r)
    return pd.DataFrame(out)


# ================================================================================================ canaries
def canary_check(pk_can: pd.DataFrame) -> dict:
    """Acceptance (DESIGN §11.3 / A5, canary v2): mean negative-canary AUROC on E-hard <= 0.57 (E-all reported, not
    gated); positive canary AUROC >= 0.90 on E-hard and E-all."""
    neg = pk_can[pk_can.kinase_acc.str.startswith("CANARY_NEG")]
    pos = pk_can[pk_can.kinase_acc.eq(_canary.POS_ID)]
    out = {"n_neg_canaries": len(neg)}
    for col in ["auroc_hard", "auroc_all"]:
        out[f"neg_mean_{col}"] = _nanmean(neg[col]) if len(neg) else float("nan")
        out[f"pos_{col}"] = float(pos[col].iloc[0]) if len(pos) else float("nan")
    # canary v2 (PI): negatives gated on E-HARD only (E-all reported: generic site-likeness legitimately lifts it)
    neg_ok = (len(neg) > 0) and not (out["neg_mean_auroc_hard"] > NEG_CANARY_MAX)
    pos_ok = all(out[f"pos_{c}"] >= POS_CANARY_MIN for c in ["auroc_hard", "auroc_all"]) if len(pos) else False
    out.update({"neg_pass": bool(neg_ok), "pos_pass": bool(pos_ok), "thresholds": {"neg_max": NEG_CANARY_MAX,
                                                                                 "pos_min": POS_CANARY_MIN}})
    return out


# ================================================================================================ bootstrap
def _cluster_weights(clusters: np.ndarray, B: int, seed: int) -> tuple[np.ndarray, dict]:
    u = np.unique(clusters)
    rng = np.random.default_rng(seed)
    W = np.stack([np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u)) for _ in range(B)]).astype(float)
    return W, {c: j for j, c in enumerate(u)}


def _boot_matrix(res: EvalResult, data_clusters: np.ndarray, B: int, seed: int, kset: list) -> dict:
    """Per-endpoint B x n_kinase matrices of resampled per-kinase metrics."""
    W, cpos = _cluster_weights(data_clusters, B, seed)
    out = {"M1": np.full((B, len(kset)), np.nan), "M3": np.full((B, len(kset)), np.nan),
           "M2": np.full((B, len(kset)), np.nan)}
    sl = res.site_level
    for j, k in enumerate(kset):
        it = res.items.get(k)
        if it is None:
            continue
        cp = np.array([cpos[c] for c in it["pos_c"]], dtype=np.int64)
        ch = np.array([cpos[c] for c in it["hard_c"]], dtype=np.int64)
        for b0 in range(0, B, 1000):
            Wb = W[b0:b0 + 1000]
            if len(cp) and len(ch):
                out["M1"][b0:b0 + 1000, j] = weighted_auroc_matrix(it["pos_s"], it["hard_s"], Wb[:, cp], Wb[:, ch])
            if len(cp):
                wp = Wb[:, cp]
                with np.errstate(invalid="ignore", divide="ignore"):
                    out["M3"][b0:b0 + 1000, j] = (wp * it["top_hits"]).sum(1) / wp.sum(1)
        if len(sl):
            x = sl[sl.kinase_acc.eq(k)]
            if len(x):
                cx = np.array([cpos[c] for c in x.cluster.values], dtype=np.int64)
                wx = W[:, cx]
                with np.errstate(invalid="ignore", divide="ignore"):
                    out["M2"][:, j] = (wx * x.hit5.values).sum(1) / wx.sum(1)
    return out


def _universe_clusters(res: EvalResult) -> np.ndarray:
    cs = [it["pos_c"] for it in res.items.values()] + [it["hard_c"] for it in res.items.values()]
    return np.unique(np.concatenate(cs)) if cs else np.zeros(0, np.int64)


def bootstrap(res: EvalResult, B: int = 1000, seed: int = 0, resample: str = "clusters",
              kinase_set: str = "primary") -> pd.DataFrame:
    """95% CIs for M1-M3 (over K_primary) and M4 (M1 metric over K_any). resample='clusters' (DESIGN §6.8: cluster
    bootstrap over substrate clusters; M3/M2 keep each kinase's top-1% cut-off fixed - documented approximation) or
    'kinases' (sensitivity: resample kinases of the macro average). B: 1,000 for development, 10,000 final (A19)."""
    pk = res.per_kinase.set_index("kinase_acc")
    sets = {"M1": "auroc_hard", "M2": "hit5", "M3": "top1_recall", "M4": "auroc_hard"}
    kp = [k for k in pk.index[pk.in_K_primary]]
    ka = [k for k in pk.index[pk.in_K_any]]
    rows = []
    if resample == "clusters":
        cl = _universe_clusters(res)
        mp = _boot_matrix(res, cl, B, seed, kp)
        ma = _boot_matrix(res, cl, B, seed, ka)
        dist = {"M1": np.nanmean(mp["M1"], 1), "M2": np.nanmean(mp["M2"], 1), "M3": np.nanmean(mp["M3"], 1),
                "M4": np.nanmean(ma["M1"], 1)}
    elif resample == "kinases":
        rng = np.random.default_rng(seed)
        dist = {}
        for e, col in sets.items():
            ks_ = kp if e != "M4" else ka
            v = pk.loc[ks_, col].values.astype(float) if col in pk else np.zeros(0)
            v = v[~np.isnan(v)]
            dist[e] = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(B)]) if len(v) else \
                np.full(B, np.nan)
    else:
        raise ValueError("resample must be 'clusters' or 'kinases'")
    for e, col in sets.items():
        ks_ = kp if e != "M4" else ka
        point = _nanmean(pk.loc[ks_, col]) if col in pk and ks_ else float("nan")
        d = dist[e][~np.isnan(dist[e])]
        rows.append({"endpoint": e, "point": point, "lo95": float(np.percentile(d, 2.5)) if len(d) else np.nan,
                     "hi95": float(np.percentile(d, 97.5)) if len(d) else np.nan,
                     "se": float(d.std(ddof=1)) if len(d) > 1 else np.nan, "B": B, "resample": resample})
    return pd.DataFrame(rows)


def paired_bootstrap(resA: EvalResult, resB: EvalResult, B: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Paired cluster bootstrap of A - B on identical resamples (DESIGN §6.8). Requires a matched universe: identical
    kinase sets and identical positive / hard-negative items per kinase (use matched_universe + missing='drop' or
    complete predictions); raises otherwise. Two-sided bootstrap p-value: 2 * min(P(d<=0), P(d>=0)), (+1)/(B+1)."""
    for k in set(resA.items) | set(resB.items):
        a, b = resA.items.get(k), resB.items.get(k)
        if a is None or b is None or len(a["pos_c"]) != len(b["pos_c"]) or len(a["hard_c"]) != len(b["hard_c"]) \
                or not np.array_equal(a["pos_c"], b["pos_c"]) or not np.array_equal(a["hard_c"], b["hard_c"]) \
                or a["all_idx_hash"] != b["all_idx_hash"]:
            raise ValueError(f"unmatched evaluation items for kinase {k}: compare on a matched universe")
    pa, pb = resA.per_kinase.set_index("kinase_acc"), resB.per_kinase.set_index("kinase_acc")
    kp = list(pa.index[pa.in_K_primary])
    ka = list(pa.index[pa.in_K_any])
    cl = _universe_clusters(resA)
    A1, B1 = _boot_matrix(resA, cl, B, seed, kp), _boot_matrix(resB, cl, B, seed, kp)
    A4, B4 = _boot_matrix(resA, cl, B, seed, ka), _boot_matrix(resB, cl, B, seed, ka)
    spec = {"M1": (A1["M1"], B1["M1"], "auroc_hard", kp), "M2": (A1["M2"], B1["M2"], "hit5", kp),
            "M3": (A1["M3"], B1["M3"], "top1_recall", kp), "M4": (A4["M1"], B4["M1"], "auroc_hard", ka)}
    rows = []
    for e, (ma, mb, col, ks_) in spec.items():
        d = np.nanmean(ma, 1) - np.nanmean(mb, 1)
        d = d[~np.isnan(d)]
        point = _nanmean(pa.loc[ks_, col]) - _nanmean(pb.loc[ks_, col]) if col in pa and col in pb else np.nan
        p = min(1.0, 2 * min(((d <= 0).sum() + 1) / (len(d) + 1), ((d >= 0).sum() + 1) / (len(d) + 1))) \
            if len(d) else np.nan
        rows.append({"endpoint": e, "diff": point, "lo95": float(np.percentile(d, 2.5)) if len(d) else np.nan,
                     "hi95": float(np.percentile(d, 97.5)) if len(d) else np.nan, "p_boot": p, "B": B})
    return pd.DataFrame(rows)


def holm(pvals) -> np.ndarray:
    """Holm step-down adjusted p-values."""
    p = np.asarray(pvals, float)
    o = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (m - r) * p[i]))
        adj[i] = run
    return adj


# ================================================================================================ universe helpers
def matched_universe(predsA: pd.DataFrame, predsB: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Restrict two prediction sets to the (kinase_acc, site_id) pairs both score (PRIOR-1); report coverage."""
    ka = pd.MultiIndex.from_frame(predsA[["kinase_acc", "site_id"]])
    kb = pd.MultiIndex.from_frame(predsB[["kinase_acc", "site_id"]])
    common = ka.intersection(kb)
    a = predsA[ka.isin(common)].reset_index(drop=True)
    b = predsB[kb.isin(common)].reset_index(drop=True)
    rep = {"n_A": len(predsA), "n_B": len(predsB), "n_common": len(common),
           "frac_A_kept": len(common) / max(len(predsA), 1), "frac_B_kept": len(common) / max(len(predsB), 1),
           "kinases_common": int(a.kinase_acc.nunique())}
    return a, b, rep


def coverage(preds: pd.DataFrame, data: EvalData) -> dict:
    """Coverage of the required evaluation pairs (positives + E-all negatives of evaluated kinases) and of sites."""
    P = map_predictions(preds, data)
    have = set(zip(P.k.values, P.s.values))
    kpos = {a: i for i, a in enumerate(data.kinases.uniprot_acc.values)}
    res_arr = data.sites.res.values
    out = {"kinases_scored": int(P.k.nunique()), "sites_scored": int(P.s.nunique()),
           "universe_sites_scored": int(data.universe[np.unique(P.s.values)].sum()),
           "universe_sites": int(data.universe.sum())}
    tot = got = tp = gp = 0
    for k, pos in data.positives.items():
        if k.startswith("CANARY_"):
            continue
        i = kpos[k]
        cand = np.where(data.universe & _res_ok(res_arr, data.kinases.residue_class.values[i]))[0]
        tot += len(cand)
        got += sum((i, s) in have for s in cand)
        tp += len(pos)
        gp += sum((i, s) in have for s in pos)
    out.update({"required_pairs": tot, "required_pairs_scored": got, "frac_required_scored": got / max(tot, 1),
                "positive_pairs": tp, "positive_pairs_scored": gp})
    return out


def coverage_at(prob: pd.DataFrame, taus=(0.5, 0.8), unlabelled_sites: set | None = None) -> dict:
    """A17 coverage: prob = DataFrame[site_id, kinase_acc, p_assign] (site-level assignment probabilities).
    Coverage@tau = #sites whose top kinase has p_assign >= tau (optionally also among `unlabelled_sites`)."""
    top = prob.groupby("site_id").p_assign.max()
    out = {}
    for t in taus:
        out[f"coverage@{t}"] = int((top >= t).sum())
        if unlabelled_sites is not None:
            out[f"coverage@{t}_unlabelled"] = int((top[top.index.isin(unlabelled_sites)] >= t).sum())
    out["n_sites"] = int(len(top))
    return out
