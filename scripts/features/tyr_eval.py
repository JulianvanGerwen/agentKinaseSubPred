"""Sanity-evaluate the F9 tyrosine features on the DEV CV (E-hard, TK kinases only).

Two parts:
  (A) per-feature AUROC alone: for every Y/STY kinase, rank its E-hard candidates by the single
      feature value and take per-kinase AUROC, then the macro mean over the kinases where the
      feature is defined.  Reported for all kinases and for n_pos >= 5 / >= 10.
  (B) a fold-aware LightGBM restricted to TKs: Cantley + F3 + F4 baseline, versus the same plus
      the F9 features.  5-fold CV over the frozen substrate-cluster folds; positives are the
      fold-training DEV evaluation positives, negatives are sampled from fold-training U_ms Y
      sites after the NegativeMasks KP / paralog / Hamming masks (1:50, 50 % E-hard-type).
      Out-of-fold scores are rank-normalised within fold and within kinase (A24) before the
      pooled per-kinase E-hard AUROC is taken.

Outputs -> results/features/tyr/{eval_single_feature.tsv, eval_lgbm_per_kinase.tsv, eval_lgbm.json}
Run: python -I scripts/features/tyr_eval.py [--threads 6]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import index, labels, splits, store  # noqa: E402

RES = ROOT / "results/features/tyr"
RES.mkdir(parents=True, exist_ok=True)
SEED = 0
NEG_RATIO = 50
HARD_FRAC = 0.5

# --------------------------------------------------------------------------- feature declarations
KS_FEATS = {                      # kinase x site
    "f1_cantley_pct": "motif/cantley_unprimed_pct_published",
    "f1_cantley_score": "motif/cantley_unprimed_score",
    "f9_sh2_dock_max": "tyr/sh2_dock_max",
    "f9_sh2_dock_local": "tyr/sh2_dock_local",
    "f9_sh2_dock_near": "tyr/sh2_dock_near",
    "f9_sh2_dock_delta": "tyr/sh2_dock_delta",
    "f9_sh2_dock_spec": "tyr/sh2_dock_spec",
    "f9_rtk_pm": "tyr/rtk_pm",
    "f9_itam_gated": "tyr/itam_gated",
}
KP_FEATS = {                      # kinase x protein
    "f3_ppi_edge_any": "context/ppi/ppi_edge_any",
    "f3_ppi_n_studies": "context/ppi/ppi_n_studies",
    "f3_ppi_rwr_dc_pct": "context/ppi/ppi_rwr_dc_pct",
    "f3_ppi_dist": "context/ppi/ppi_dist",
    "f4_loc_cos_combined": "context/loc/loc_cos_combined",
    "f4_loc_cos_hpa": "context/loc/loc_cos_hpa",
    "f4_loc_cos_hcm": "context/loc/loc_cos_hcm",
}
SITE_VECS = {"tyr/site_sh2": ["f9_n_other_py", "f9_n_other_py_50aa",
                              "f9_adaptor_sh2_max", "f9_adaptor_sh2_local"],
             "tyr/site_itam": ["f9_itam_half", "f9_itam_full", "f9_itim", "f9_itam_order"],
             "context/loc/site_topo_code": ["f4_site_topo_code"]}
KIN_VECS = {"tyr/kinase_desc": ["f9_k_n_sh2", "f9_k_n_sh3", "f9_k_n_ptb",
                                "f9_k_is_tm", "f9_k_is_receptor", "f9_k_has_sh2_matrix"]}
#: Site-density features that double as a site-popularity proxy (B0-pop scores 0.583 on Y kinases).
#: Kept separate so that any F9 gain can be shown not to be popularity.
POP_PROXY = {"f9_n_other_py", "f9_n_other_py_50aa", "f9_adaptor_sh2_max", "f9_adaptor_sh2_local"}
#: Feature-set variants compared in part (B).
VARIANTS = {
    "base":          lambda n: not n.startswith("f9_"),
    "f9_all":        lambda n: True,
    "f9_nopop":      lambda n: n not in POP_PROXY,
    "f9_sh2only":    lambda n: not n.startswith("f9_") or n.startswith("f9_sh2_dock"),
    "f9_poponly":    lambda n: not n.startswith("f9_") or n in POP_PROXY,
    "f9_nosh2":      lambda n: not n.startswith("f9_sh2_dock"),
}


def auc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos, neg])


class Ctx:
    """Shared index / label / mask state for the TK sub-problem."""

    def __init__(self):
        self.ki, self.si = index.load_kinases(), index.load_sites()
        self.prot = index.load_proteins()
        self.kpos = pd.Series(np.arange(len(self.ki)), index=self.ki.uniprot_acc)
        self.spos = pd.Series(np.arange(len(self.si)), index=self.si.site_id)
        self.ppos = pd.Series(np.arange(len(self.prot)), index=self.prot.acc)
        self.site_prot = self.ppos.reindex(self.si.acc).to_numpy()
        self.ums = self.si.in_U_ms.to_numpy().astype(bool)
        self.isY = self.si.res.to_numpy() == "Y"
        self.allp = labels.load_dev("all")
        self.ev = labels.load_dev("eval")
        self.nm = labels.NegativeMasks(self.allp)
        self.par = labels.paralog_sets()
        self.folds = splits.load_cv_folds()            # acc -> fold (DEV proteins)
        self.site_fold = self.folds.reindex(self.si.acc).to_numpy(dtype=float)
        y_cls = self.ki.residue_class.isin(["Y", "STY"])
        self.tk_accs = set(self.ki.uniprot_acc[y_cls])
        self.evY = self.ev[self.ev.res.eq("Y") & self.ev.kinase_acc.isin(self.tk_accs)]
        self.gene = self.ki.set_index("uniprot_acc").gene
        self.family = self.ki.set_index("uniprot_acc").manning_family

    def pool(self, acc):
        """(positives, E-hard negatives) as site-index positions, official universe (U_ms, class Y)."""
        hard = self.nm.hard_candidates(acc)
        hard = hard[self.ums[hard] & self.isY[hard]]
        pi = self.spos.reindex(self.evY.loc[self.evY.kinase_acc.eq(acc), "site_id"]).dropna()
        pi = pi.astype(int).to_numpy()
        neg = np.setdiff1d(np.setdiff1d(hard, self.nm.excluded(acc)), pi)
        return pi, neg


# ------------------------------------------------------------------------------------ (A) singles
def single_feature_sanity(ctx: Ctx) -> pd.DataFrame:
    kacc = sorted(ctx.evY.kinase_acc.unique())
    pools = {a: ctx.pool(a) for a in kacc}
    n_pos = {a: len(pools[a][0]) for a in kacc}
    rows = []

    def record(name, getter, direction=1):
        per = {}
        for a in kacc:
            pi, neg = pools[a]
            if len(pi) == 0 or len(neg) == 0:
                continue
            v = getter(a)
            if v is None:
                continue
            per[a] = auc(direction * v[pi], direction * v[neg])
        per = {k: v for k, v in per.items() if np.isfinite(v)}
        if not per:
            return
        s = pd.Series(per)
        rows.append(dict(feature=name, n_kinases=len(s), macro_auroc=s.mean(),
                         macro_auroc_npos_ge5=s[[a for a in s.index if n_pos[a] >= 5]].mean(),
                         macro_auroc_npos_ge10=s[[a for a in s.index if n_pos[a] >= 10]].mean(),
                         n_kinases_ge5=sum(n_pos[a] >= 5 for a in s.index)))

    for name, path in KS_FEATS.items():
        M, _ = store.load_matrix(path)
        record(name, lambda a, M=M: np.asarray(M[ctx.kpos[a]], float))
    for name, path in KP_FEATS.items():
        M, _ = store.load_matrix(path)
        sgn = -1 if name == "f3_ppi_dist" else 1       # smaller graph distance = closer
        record(name, lambda a, M=M: np.asarray(M[ctx.kpos[a]], float)[ctx.site_prot], sgn)
    for path, cols in SITE_VECS.items():
        V, _ = store.load_matrix(path)
        V = np.asarray(V, float)
        for j, c in enumerate(cols):
            record(c, lambda a, V=V, j=j: V[:, j] if V.ndim > 1 else V)
    return pd.DataFrame(rows).sort_values("macro_auroc_npos_ge5", ascending=False)


# ------------------------------------------------------------------------------------- (B) LightGBM
def gather(ctx, kacc, sidx, mats_ks, mats_kp, site_mat, kin_row, names):
    """Feature block for one kinase over the given site positions."""
    n = len(sidx)
    X = np.empty((n, len(names)), np.float32)
    j = 0
    for _, M in mats_ks:
        X[:, j] = M[ctx.kpos[kacc]][sidx]; j += 1
    for nm_, M in mats_kp:
        X[:, j] = np.asarray(M[ctx.kpos[kacc]])[ctx.site_prot[sidx]]; j += 1
    for col in site_mat:
        X[:, j] = col[sidx]; j += 1
    for v in kin_row[kacc]:
        X[:, j] = v; j += 1
    return X


def lgbm_check(ctx: Ctx, keep, threads: int, seed: int = SEED, min_pos: int = 5):
    import lightgbm as lgb
    rng = np.random.default_rng(seed)

    names_ks = [(n, p) for n, p in KS_FEATS.items() if keep(n)]
    names_kp = [(n, p) for n, p in KP_FEATS.items() if keep(n)]
    mats_ks = [(n, store.load_matrix(p)[0]) for n, p in names_ks]
    mats_kp = [(n, store.load_matrix(p)[0]) for n, p in names_kp]
    site_cols, site_names = [], []
    for path, cols in SITE_VECS.items():
        V = np.asarray(store.load_matrix(path)[0], np.float32)
        for j, c in enumerate(cols):
            if keep(c):
                site_cols.append(V[:, j] if V.ndim > 1 else V)
                site_names.append(c)
    kin_names, kin_row = [], {}
    for path, cols in KIN_VECS.items():
        V = np.asarray(store.load_matrix(path)[0], np.float32)
        use = [j for j, c in enumerate(cols) if keep(c)]
        kin_names += [cols[j] for j in use]
        for a in ctx.ki.uniprot_acc:
            kin_row.setdefault(a, []).extend(V[ctx.kpos[a], use].tolist())
    kin_row = {a: np.array(v, np.float32) for a, v in kin_row.items()}
    names = [n for n, _ in names_ks] + [n for n, _ in names_kp] + site_names + kin_names

    kacc = [a for a in sorted(ctx.evY.kinase_acc.unique()) if len(ctx.pool(a)[0]) >= min_pos]
    pools = {a: ctx.pool(a) for a in kacc}
    sf = ctx.site_fold
    oof = {a: {} for a in kacc}

    for f in sorted(pd.unique(ctx.folds.values)):
        Xs, ys, ks = [], [], []
        tr_mask = (sf != f) & ~pd.isna(sf)
        for a in kacc:
            pi, neg = pools[a]
            pi_tr = pi[tr_mask[pi]]
            if len(pi_tr) == 0:
                continue
            hard_tr = neg[tr_mask[neg]]
            # E-all-type training negatives: fold-training U_ms Y sites after the masks
            pool_all = np.flatnonzero(ctx.ums & ctx.isY & tr_mask)
            pool_all = np.setdiff1d(pool_all, ctx.nm.excluded(a))
            want = NEG_RATIO * len(pi_tr)
            n_hard = min(int(HARD_FRAC * want), len(hard_tr))
            neg_h = rng.choice(hard_tr, n_hard, replace=False) if n_hard else np.zeros(0, np.int64)
            rest = min(want - n_hard, len(pool_all))
            neg_a = rng.choice(pool_all, rest, replace=False) if rest > 0 else np.zeros(0, np.int64)
            sidx = np.concatenate([pi_tr, neg_h, neg_a]).astype(np.int64)
            Xs.append(gather(ctx, a, sidx, mats_ks, mats_kp, site_cols, kin_row, names))
            ys.append(np.r_[np.ones(len(pi_tr)), np.zeros(len(neg_h) + len(neg_a))])
            ks.append(np.full(len(sidx), a))
        X = np.vstack(Xs); y = np.concatenate(ys); kk = np.concatenate(ks)
        w = pd.Series(kk).map(1.0 / pd.Series(kk).value_counts() ** 0.5).to_numpy()
        m = lgb.train(dict(objective="binary", learning_rate=0.05, num_leaves=31,
                           min_data_in_leaf=50, feature_fraction=0.8, bagging_fraction=0.8,
                           bagging_freq=1, num_threads=threads, verbose=-1, seed=seed),
                      lgb.Dataset(X, y, weight=w), num_boost_round=250)
        for a in kacc:
            pi, neg = pools[a]
            ev_sites = np.concatenate([pi[sf[pi] == f], neg[sf[neg] == f]]).astype(np.int64)
            if len(ev_sites) == 0:
                continue
            p = m.predict(gather(ctx, a, ev_sites, mats_ks, mats_kp, site_cols, kin_row, names))
            r = pd.Series(p).rank(pct=True).to_numpy()       # A24: within fold, within kinase
            oof[a].update(dict(zip(ev_sites.tolist(), r.tolist())))


    rows = []
    for a in kacc:
        pi, neg = pools[a]
        sp = np.array([oof[a].get(int(i), np.nan) for i in pi])
        sn = np.array([oof[a].get(int(i), np.nan) for i in neg])
        rows.append(dict(kinase_acc=a, gene=ctx.gene[a], manning_family=ctx.family[a],
                         n_pos=len(pi), n_neg_hard=len(neg), auroc_hard=auc(sp, sn)))
    return pd.DataFrame(rows), names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--skip-lgbm", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    ctx = Ctx()

    single = single_feature_sanity(ctx)
    single.round(4).to_csv(RES / "eval_single_feature.tsv", sep="\t", index=False)
    print("\n=== (A) single-feature E-hard AUROC, TK kinases ===")
    print(single.round(3).to_string(index=False))

    if a.skip_lgbm:
        return
    out, per_k = {}, {}
    for tag in a.variants:
        keep = VARIANTS[tag]
        reps = []
        for s in a.seeds:
            pk, names = lgbm_check(ctx, keep, a.threads, seed=s)
            reps.append(pk.set_index("kinase_acc").auroc_hard.rename(s))
        R = pd.concat(reps, axis=1)
        pk = pk.set_index("kinase_acc")
        pk["auroc_hard"] = R.mean(axis=1)
        pk["auroc_sd_seeds"] = R.std(axis=1)
        per_k[tag] = pk.reset_index()
        out[tag] = dict(n_features=len(names), features=names, seeds=a.seeds,
                        macro_auroc_hard=float(pk.auroc_hard.mean()),
                        macro_npos_ge10=float(pk.loc[pk.n_pos >= 10, "auroc_hard"].mean()),
                        macro_npos_ge20=float(pk.loc[pk.n_pos >= 20, "auroc_hard"].mean()),
                        mean_seed_sd=float(pk.auroc_sd_seeds.mean()), n_kinases=int(len(pk)))
        print(f"  {tag:12s} feats={len(names):3d}  macro={out[tag]['macro_auroc_hard']:.4f}  "
              f"n>=10={out[tag]['macro_npos_ge10']:.4f}  n>=20={out[tag]['macro_npos_ge20']:.4f}  "
              f"seed-sd={out[tag]['mean_seed_sd']:.3f}", flush=True)

    base = per_k["base"].set_index("kinase_acc")
    tabs = {t: per_k[t].set_index("kinase_acc").auroc_hard for t in a.variants}
    m = base[["gene", "manning_family", "n_pos", "n_neg_hard"]].copy()
    for t in a.variants:
        m[f"auroc_{t}"] = tabs[t]
    for t in a.variants:
        if t != "base":
            m[f"delta_{t}"] = m[f"auroc_{t}"] - m["auroc_base"]
            out[t]["delta_macro"] = float(m[f"delta_{t}"].mean())
            out[t]["delta_macro_npos_ge10"] = float(m.loc[m.n_pos >= 10, f"delta_{t}"].mean())
            out[t]["n_kinases_improved"] = int((m[f"delta_{t}"] > 0).sum())
    m = m.sort_values("n_pos", ascending=False)
    m.round(4).to_csv(RES / "eval_lgbm_per_kinase.tsv", sep="\t")
    agg = {"n_kinases": ("gene", "size"), "n_pos": ("n_pos", "sum")}
    agg.update({t: (f"auroc_{t}", "mean") for t in a.variants})
    fam = m.groupby("manning_family").agg(**agg).sort_values("n_pos", ascending=False)
    out["per_family"] = json.loads(fam.round(4).to_json(orient="index"))
    out["runtime_s"] = round(time.time() - t0, 1)
    (RES / "eval_lgbm.json").write_text(json.dumps(out, indent=2))
    print("\n--- per kinase (seed-averaged) ---\n", m.round(3).to_string())
    print("\n--- per Manning family ---\n", fam.round(3).to_string())


if __name__ == "__main__":
    main()
