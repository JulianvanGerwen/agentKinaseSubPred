"""R3 training / scoring (DESIGN §8.4, §5, §16 A6/A13). Owner: R3-neural.

One training run = one outer split:
  CV fold f : training rows = U_ms sites on DEV proteins of folds != f; scored rows = U_ms sites of fold f.
  LKO fold j: training rows = all DEV U_ms sites, labels of kinases in kinase fold j removed (no positives, no
              negatives for them); scored = held-out kinases x all DEV U_ms sites (T2: unseen family, seen proteins).
Early stopping uses ONLY an inner validation split of the training rows (ksp.splits.train_valid_split_inner, grouped by
substrate cluster): macro AUROC(hard negatives) over (non-canary) kinases with >= 3 inner-validation positives.

Each step samples B sites (fraction p_lab from sites carrying >= 1 training positive, rest uniform over training rows)
and scores them against ALL kinases (B x K logits). Loss = balanced masked BCE (positives weighted w_k ~ n_k^-alpha)
+ lambda_nce * filtered listwise InfoNCE over kinases per positive (see model.py).
"""
from __future__ import annotations

import copy
import time

import numpy as np
import pandas as pd
import torch

from ksp import labels, splits
from ksp.models.r3 import data as r3data
from ksp.models.r3.model import R3Net, bce_masked, infonce_filtered

DEFAULT = dict(flank=7, lambda_nce=0.5, alpha=0.5, p_lab=0.5, batch=1024, lr=1e-3, wd=0.01, max_steps=6000,
               eval_every=200, patience=6, desc_dropout=0.2, d_site=192, h=128, d_kin=256, n_layers=2,
               dropout=0.1, esm_dropout=0.1, valid_frac=0.2, seed=0, variant="L2", warmup=200)

_CACHE: dict = {}


def inputs():
    if "si" not in _CACHE:
        _CACHE["si"] = r3data.load_sites()
        _CACHE["ki"] = r3data.load_kinases()
        _CACHE["par"] = labels.paralog_sets()
    return _CACHE["si"], _CACHE["ki"], _CACHE["par"]


class GPUData:
    def __init__(self, si, ki, T, par, flank, dev, fam_keep, grp_keep):
        c0 = r3data.FLANK_MAX - flank
        self.dev = dev
        self.flank = flank
        t0 = time.time()
        self.E = torch.from_numpy(np.ascontiguousarray(si.emb)).to(dev)            # fp16 (R,1280)
        self.load_s = time.time() - t0
        self.rows = torch.from_numpy(si.rows[:, c0:c0 + 2 * flank + 1].copy()).to(dev)
        self.aa = torch.from_numpy(si.res_aa).to(dev)
        self.phos = torch.from_numpy(si.res_phos).to(dev)
        self.prow = torch.from_numpy(si.prot_row).to(dev)
        self.pm = torch.from_numpy(si.prot_mean).to(dev)
        self.T = torch.from_numpy(T).to(dev) if T is not None else None
        K = len(ki.kinases)
        self.k_esm = torch.from_numpy(np.c_[ki.esm, ki.esm_missing]).float().to(dev)
        self.k_cant = torch.from_numpy(np.c_[ki.cantley, ki.cantley_missing, ki.cantley_imputed,
                                             ki.cantley_dist]).float().to(dev)
        g = np.where(grp_keep[ki.group], ki.group, len(ki.group_names))
        f = np.where((ki.family >= 0) & fam_keep[np.maximum(ki.family, 0)], ki.family, len(ki.family_names))
        self.k_grp = torch.from_numpy(g).long().to(dev)
        self.k_fam = torch.from_numpy(f).long().to(dev)
        P = np.zeros((K, K), bool)
        kp = {a: i for i, a in enumerate(ki.kinases.uniprot_acc.values)}
        for a, ss in par.items():
            for b in ss:
                if a in kp and b in kp:
                    P[kp[a], kp[b]] = True
        self.par = torch.from_numpy(P).to(dev)

    def site_batch(self, idx):
        r = self.rows[idx]
        valid = r >= 0
        rr = r.clamp(min=0)
        emb = self.E[rr].float()
        aa = self.aa[rr] * valid
        ph = self.phos[rr] * valid
        ph[:, self.flank] = 0                       # never the centre site's own status
        return emb, aa, ph, valid, self.pm[self.prow[idx]]

    def kin_inputs(self):
        return self.k_esm, self.k_cant, self.k_grp, self.k_fam


def _auroc_t(sp, sn):
    s = torch.cat([sp, sn])
    r = torch.empty_like(s)
    r[torch.argsort(s)] = torch.arange(1, len(s) + 1, device=s.device, dtype=s.dtype)
    n1, n0 = len(sp), len(sn)
    return ((r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)).item()


@torch.no_grad()
def score_rows(model, G, rows_idx, chunk=4096):
    model.eval()
    out = []
    ki = G.kin_inputs()
    for a in range(0, len(rows_idx), chunk):
        idx = torch.as_tensor(rows_idx[a:a + chunk], device=G.dev)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            out.append(model(G.site_batch(idx), ki).float())
    return torch.cat(out)


def inner_val_metric(model, G, val_rows, eval_k):
    L = score_rows(model, G, val_rows)
    T = G.T[torch.as_tensor(val_rows, device=G.dev)]
    hard_row = (T == 1).any(1)
    vals = []
    for k in eval_k:
        p = T[:, k] == 1
        n = (T[:, k] == 0) & hard_row
        if p.sum() >= 3 and n.sum() >= 3:
            vals.append(_auroc_t(L[p, k], L[n, k]))
    return float(np.mean(vals)) if vals else float("nan"), len(vals)


def run(cfg: dict, fold: int | None = None, lko_fold: int | None = None, device: str = "cuda", log=print,
        score: bool = True):
    """Train one R3 model; return (predictions DataFrame[kinase_idx, site_idx, kinase_acc, site_id, score], info)."""
    cfg = {**DEFAULT, **cfg}
    t_start = time.time()
    si, ki, par = inputs()
    K = len(ki.kinases)
    rng = np.random.default_rng(cfg["seed"])
    torch.manual_seed(cfg["seed"])
    if (fold is None) == (lko_fold is None):
        raise ValueError("give exactly one of fold / lko_fold")
    dev_rows = si.cv_fold >= 0
    if fold is not None:
        train_mask = dev_rows & (si.cv_fold != fold)
        score_mask = si.cv_fold == fold
        train_k = None
        score_k = np.ones(K, bool)
    else:
        kf = splits.load_kinase_folds()
        held = set(kf.kinase_acc[kf.lko_fold == lko_fold])
        is_held = ki.kinases.uniprot_acc.isin(held).values
        train_mask = dev_rows
        score_mask = dev_rows
        train_k = ~is_held
        score_k = is_held
    info = {"fold": fold, "lko_fold": lko_fold, "cfg": cfg, "n_train_rows": int(train_mask.sum()),
            "n_score_rows": int(score_mask.sum())}
    tlog: dict = {}
    t0 = time.time()
    T = r3data.build_targets(si, ki, train_mask, train_k, variant=cfg["variant"], par=par, log=tlog)
    info["targets"] = tlog
    info["t_targets_s"] = round(time.time() - t0, 1)
    tr_idx = np.where(train_mask)[0]
    val_flag = splits.train_valid_split_inner(si.site_ids[tr_idx], valid_frac=cfg["valid_frac"], seed=cfg["seed"])
    val_rows = tr_idx[val_flag]
    fit_rows = tr_idx[~val_flag]
    # weights from fit-row positives
    n_k = (T[fit_rows] == 1).sum(0).astype(np.float64)
    trained = n_k > 0
    wk = np.where(trained, np.maximum(n_k, 1) ** (-cfg["alpha"]), 0.0)
    wk = wk / (wk[trained] * n_k[trained]).sum() * n_k[trained].sum()   # mean weight over positives = 1
    lab_rows = fit_rows[(T[fit_rows] == 1).any(1)]
    info.update(n_fit_rows=len(fit_rows), n_val_rows=len(val_rows), n_lab_fit_rows=len(lab_rows),
                n_trained_kinases=int(trained.sum()), n_pos_fit=int(n_k.sum()))
    fam_keep = np.zeros(len(ki.family_names), bool)
    grp_keep = np.zeros(len(ki.group_names), bool)
    fam_keep[ki.family[trained & (ki.family >= 0)]] = True
    grp_keep[ki.group[trained]] = True
    G = GPUData(si, ki, T, par, cfg["flank"], device, fam_keep, grp_keep)
    info["t_gpu_load_s"] = round(G.load_s, 1)
    model = R3Net(si.emb.shape[1], G.k_cant.shape[1], len(ki.group_names), len(ki.family_names), flank=cfg["flank"],
                  d_site=cfg["d_site"], h=cfg["h"], d_kin=cfg["d_kin"], n_layers=cfg["n_layers"],
                  desc_dropout=cfg["desc_dropout"], dropout=cfg["dropout"],
                  esm_dropout=cfg["esm_dropout"]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / cfg["warmup"]))
    wk_t = torch.as_tensor(wk, dtype=torch.float32, device=device)
    eval_k = [k for k in range(K) if trained[k] and not ki.is_canary[k]]
    B, nl = cfg["batch"], int(cfg["batch"] * cfg["p_lab"])
    best, best_state, best_step, bad = -1.0, None, 0, 0
    hist = []
    t_train = time.time()
    for step in range(1, cfg["max_steps"] + 1):
        model.train()
        idx = np.r_[rng.choice(lab_rows, nl), rng.choice(fit_rows, B - nl)]
        idx_t = torch.as_tensor(idx, device=device)
        Tb = G.T[idx_t]
        with torch.autocast("cuda", dtype=torch.bfloat16):
            L = model(G.site_batch(idx_t), G.kin_inputs()).float()
        loss_b = bce_masked(L, Tb, wk_t)
        loss = loss_b
        loss_n = torch.zeros((), device=device)
        if cfg["lambda_nce"] > 0:
            loss_n = infonce_filtered(L, Tb, G.par, wk_t)
            loss = loss + cfg["lambda_nce"] * loss_n
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % cfg["eval_every"] == 0:
            m, nk = inner_val_metric(model, G, val_rows, eval_k)
            hist.append({"step": step, "loss_bce": float(loss_b), "loss_nce": float(loss_n), "val_auroc_hard": m,
                         "val_n_kinases": nk, "t": round(time.time() - t_train, 1)})
            log(f"step {step} bce {float(loss_b):.4f} nce {float(loss_n):.4f} val_auroc_hard {m:.4f} ({nk} k) "
                f"{time.time() - t_train:.0f}s")
            if m > best:
                best, best_step, bad = m, step, 0
                best_state = copy.deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()})
            else:
                bad += 1
                if bad >= cfg["patience"]:
                    break
    info.update(best_val_auroc_hard=best, best_step=best_step, history=hist,
                t_train_s=round(time.time() - t_train, 1))
    model.load_state_dict(best_state)
    if not score:
        del G
        torch.cuda.empty_cache()
        return None, info, model
    t0 = time.time()
    srows = np.where(score_mask)[0]
    S = score_rows(model, G, srows).cpu().numpy()                     # (n, K)
    ok = ki.class_ok[si.res_code[srows]] & score_k[None]
    ri, kk = np.nonzero(ok)
    preds = pd.DataFrame({"kinase_idx": kk.astype(np.int64), "site_idx": si.site_idx[srows[ri]],
                          "score": S[ri, kk].astype(np.float32)})
    kacc = ki.kinases.uniprot_acc.values
    preds["kinase_acc"] = pd.Categorical.from_codes(kk, categories=kacc)
    preds["site_id"] = pd.Categorical.from_codes(ri, categories=si.site_ids[srows])
    info["t_score_s"] = round(time.time() - t0, 1)
    info["n_pred_rows"] = len(preds)
    info["t_total_s"] = round(time.time() - t_start, 1)
    info["peak_gpu_mem_gb"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
    del G
    torch.cuda.empty_cache()
    return preds, info, model
