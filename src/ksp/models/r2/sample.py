"""R2 training-set construction per CV fold (DESIGN §5, §16 A6/A15; ALG-1/2/8).

Positives : outer-training DEV labels (variant L2 by default; U_ms only via load_dev; residue in the kinase's class),
            canary labels included (flag is_canary) - canaries are reported separately, never in metrics.
Negatives : per kinase k, from U_ms sites on outer-TRAINING proteins of k's residue class, minus
            NegativeMasks.excluded(k) (KP(k) incl. identical 15-mers, Hamming<=2 of KP(k), sites labelled to Par(k))
            minus k's positives. Masks are built from TRAINING-fold label rows only (all provenances, + canaries), so
            no eval-fold label information shapes the training set.
            n_neg = ratio * n_pos(k); hard share h: E-hard-type sites (>=1 human-native non-autophos label to a
            kinase outside {k} U Par(k)); if the hard pool is smaller, all hard sites are taken and the rest is
            E-all-type (A15). E-all-type = uniform over the remaining candidates.
Weights   : every row of kinase k gets n_pos(k)^-alpha; normalised to mean 1. Logged: W_k = total weight (ALG-8)
            and rho_k = n_neg_sampled / n_candidates (ALG-1 case-control rate; for prior correction).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ksp import canary, index, labels, splits


def residue_ok(res: np.ndarray, cls: str) -> np.ndarray:
    if cls == "STY":
        return np.isin(res, ["S", "T", "Y"])
    if cls == "ST":
        return np.isin(res, ["S", "T"])
    if cls == "Y":
        return res == "Y"
    raise ValueError(cls)


class FoldData:
    """Everything label-related for one outer CV fold (fold=None: all DEV = final model)."""

    def __init__(self, fold: int | None, variant: str = "L2", exclude_kinases: set | None = None,
                 deploy_loader=None):
        """exclude_kinases: leave-kinase-out - their label rows are removed from positives AND masks.
        deploy_loader: G5 deployment only - callable(variant, with_canaries) -> label table over DEV+TEST proteins
        (ksp.deploy.load_deploy_labels); fold must be None; training sites = U_ms sites on DEV+TEST proteins."""
        self.fold = fold
        excl = set(exclude_kinases or ())
        S = index.load_sites()
        self.S = S
        self.K = canary.kinases_with_canaries()
        self.kpos = {a: i for i, a in enumerate(self.K.uniprot_acc)}
        spos = pd.Series(np.arange(len(S)), index=S.site_id)
        folds = splits.load_cv_folds()
        pf = S.acc.map(folds)
        ums = S.in_U_ms.values
        if deploy_loader is not None:
            if fold is not None:
                raise ValueError("deployment uses fold=None")
            allprot = set(folds.index) | splits.load_test_proteins()
            self.train_sites = np.where(ums & S.acc.isin(allprot).values)[0]
            self.eval_sites = np.zeros(0, np.int64)
            train_prot = allprot
        elif fold is None:
            self.train_sites = np.where(ums & pf.notna().values)[0]
            self.eval_sites = np.zeros(0, np.int64)
            train_prot = set(folds.index)
        else:
            self.train_sites = np.where(ums & pf.notna().values & (pf.values != fold))[0]
            self.eval_sites = np.where(ums & (pf.values == fold))[0]
            train_prot = set(folds.index[folds != fold])
        lab = labels.load_dev(variant, with_canaries=True) if deploy_loader is None else \
            deploy_loader(variant, with_canaries=True)
        self.n_label_rows_loaded = int(len(lab))
        lab = lab[lab.acc.isin(train_prot) & ~lab.kinase_acc.isin(excl)]
        lab = lab[lab.site_id.isin(spos.index)]
        lab = lab.assign(kinase_idx=lab.kinase_acc.map(self.kpos), site_idx=spos[lab.site_id].values)
        if lab.kinase_idx.isna().any():
            raise KeyError(f"labels with kinases outside the index: {lab.kinase_acc[lab.kinase_idx.isna()].unique()[:5]}")
        lab["kinase_idx"] = lab.kinase_idx.astype(np.int64)
        res = S.res.values[lab.site_idx.values]
        cls = self.K.residue_class.values[lab.kinase_idx.values]
        ok = np.array([residue_ok(np.array([r]), c)[0] for r, c in zip(res, cls)])
        self.n_dropped_class = int((~ok).sum())
        lab = lab[ok].drop_duplicates(["kinase_idx", "site_idx"]).reset_index(drop=True)
        self.pos = lab[["kinase_acc", "site_id", "kinase_idx", "site_idx", "is_canary"]].copy()
        # masks from TRAINING-fold rows only (all provenances / mask-only rows) + training canaries
        allp = labels.load_dev("all") if deploy_loader is None else deploy_loader("all", with_canaries=False)
        allp = allp[allp.acc.isin(train_prot) & ~allp.kinase_acc.isin(excl)]
        if deploy_loader is None:
            can = canary.canary_labels()
        else:
            cl = deploy_loader("all", with_canaries=True)
            can = cl[cl.is_canary.astype(bool)]
        can = can[can.acc.isin(train_prot)]
        self.mask_pairs = pd.concat([allp, can], ignore_index=True)
        self.masks = labels.NegativeMasks(self.mask_pairs)

    def sample(self, ratio: int = 50, hard_frac: float = 0.5, alpha: float = 0.5, seed: int = 0) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        S = self.S
        res_all = S.res.values
        is_train = np.zeros(len(S), bool)
        is_train[self.train_sites] = True
        rows = []
        stats = []
        for k, g in self.pos.groupby("kinase_idx"):
            acc = self.K.uniprot_acc.values[k]
            cls = self.K.residue_class.values[k]
            cand = is_train & residue_ok(res_all, cls)
            cand[self.masks.excluded(acc)] = False
            cand[g.site_idx.values] = False
            hard = np.zeros(len(S), bool)
            hard[self.masks.hard_candidates(acc)] = True
            hard &= cand
            n_pos = len(g)
            n_neg = int(ratio * n_pos)
            hp = np.where(hard)[0]
            n_h = min(int(round(hard_frac * n_neg)), len(hp))
            h_sel = rng.choice(hp, n_h, replace=False) if n_h else np.zeros(0, np.int64)
            rest = cand.copy()
            rest[h_sel] = False
            rp = np.where(rest)[0]
            n_e = min(n_neg - n_h, len(rp))
            e_sel = rng.choice(rp, n_e, replace=False) if n_e else np.zeros(0, np.int64)
            neg = np.concatenate([h_sel, e_sel])
            rows.append(pd.DataFrame({"kinase_idx": k, "site_idx": g.site_idx.values, "y": 1, "neg_type": "pos"}))
            rows.append(pd.DataFrame({"kinase_idx": k, "site_idx": neg, "y": 0,
                                      "neg_type": np.r_[np.full(n_h, "hard"), np.full(n_e, "all")]}))
            stats.append({"kinase_idx": k, "kinase_acc": acc, "n_pos": n_pos, "n_hard": n_h, "n_all": n_e,
                          "n_hard_pool": int(hard.sum()), "n_cand": int(cand.sum()),
                          "rho": (n_h + n_e) / max(int(cand.sum()), 1), "is_canary": acc.startswith("CANARY_")})
        d = pd.concat(rows, ignore_index=True)
        st = pd.DataFrame(stats).set_index("kinase_idx")
        w = st.n_pos.astype(float) ** (-alpha)
        d["w"] = d.kinase_idx.map(w).values
        d["w"] *= len(d) / d.w.sum()
        st["W_k"] = d.groupby("kinase_idx").w.sum()
        d.attrs["stats"] = st
        return d
