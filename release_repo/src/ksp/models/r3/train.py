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
               dropout=0.1, esm_dropout=0.1, use_esm=True, use_prot=True, use_phos=True, valid_frac=0.2, seed=0, variant="L2", warmup=200,
               n_seeds=1, esm="650m", pretrain_steps=0, distill_aux=0.0, pretrain_lr=1e-3)

_CACHE: dict = {}


def inputs(esm: str = "650m"):
    if _CACHE.get("esm") != esm:
        _CACHE["si"] = r3data.load_sites(esm)
        _CACHE["esm"] = esm
    if "ki" not in _CACHE:
        _CACHE["ki"] = r3data.load_kinases()
        _CACHE["par"] = labels.paralog_sets()
    return _CACHE["si"], _CACHE["ki"], _CACHE["par"]


def cantley_targets(si, ki, variants=("unprimed", "primed_soft")) -> np.ndarray:
    """Label-free distillation targets: Cantley universe percentile / 100 (own atlas matrix only; NaN where the
    kinase has no own in-class matrix; canaries NaN) for every R3 site x kinase. float16 (N, K, n_variants)."""
    from ksp import store
    nR = int((~ki.is_canary).sum())
    out = np.full((len(si.site_idx), len(ki.kinases), len(variants)), np.nan, np.float16)
    order = np.argsort(si.site_idx)
    cols = si.site_idx[order]
    for j, v in enumerate(variants):
        M, _ = store.load_matrix(f"motif/cantley_{v}_pct_universe")
        X = np.asarray(M[:, cols], np.float32) / 100.0              # (503, N) in sorted-site order
        out[order, :nR, j] = X.T.astype(np.float16)
    return out


class GPUData:
    def __init__(self, si, ki, T, par, flank, dev, fam_keep, grp_keep):
        c0 = si.info["flank_max"] - flank
        if c0 < 0:
            raise ValueError(f"flank {flank} > cached window flank {si.info['flank_max']}")
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
        score: bool = True, test: bool = False, deploy: bool = False, label_loader=None,
        deploy_train_rows: str = "all"):
    """Train one R3 model; return (predictions DataFrame[kinase_idx, site_idx, kinase_acc, site_id, score], info).
    test=True (final model, DESIGN §8.6): train on ALL DEV rows and score the TEST-protein U_ms sites (cv_fold == -1;
    sites only, no TEST labels are read) for real (non-canary) kinases; with lko_fold=j the kinases of LKO fold j are
    removed from training and only they are scored (T3). fold must be None when test=True."""
    cfg = {**DEFAULT, **cfg}
    t_start = time.time()
    si, ki, par = inputs(cfg["esm"])
    K = len(ki.kinases)
    if deploy and (fold is not None or lko_fold is not None or test):
        raise ValueError("deploy=True: no fold / lko_fold / test")
    if test and fold is not None:
        raise ValueError("test=True trains on all DEV: fold must be None")
    if not test and not deploy and (fold is None) == (lko_fold is None):
        raise ValueError("give exactly one of fold / lko_fold")
    in_ums = si.in_ums if si.in_ums is not None else np.ones(len(si.site_ids), bool)
    dev_rows = (si.cv_fold >= 0) & in_ums          # training/eval rows are always U_ms (A6)
    if deploy:
        # DESIGN §8.6 deployment: train on DEV+TEST U_ms rows ('all'; labels via label_loader = ksp.deploy) or on
        # DEV rows only ('dev'; smoke test with DEV labels); score every site of the cache (U_obs) for real kinases
        train_mask = in_ums.copy() if deploy_train_rows == "all" else dev_rows.copy()
        score_mask = np.ones(len(si.site_ids), bool)
        train_k = None
        score_k = ~ki.is_canary
    elif test and lko_fold is None:
        train_mask = dev_rows
        score_mask = si.cv_fold == -1
        train_k = None
        score_k = ~ki.is_canary
    elif fold is not None:
        train_mask = dev_rows & (si.cv_fold != fold)
        score_mask = si.cv_fold == fold
        train_k = None
        score_k = np.ones(K, bool)
    else:
        kf = splits.load_kinase_folds()
        held = set(kf.kinase_acc[kf.lko_fold == lko_fold])
        is_held = ki.kinases.uniprot_acc.isin(held).values
        train_mask = dev_rows
        score_mask = (si.cv_fold == -1) if test else dev_rows
        train_k = ~is_held
        score_k = is_held & ~ki.is_canary
    info = {"fold": fold, "lko_fold": lko_fold, "test": test, "cfg": cfg, "n_train_rows": int(train_mask.sum()),
            "n_score_rows": int(score_mask.sum())}
    tlog: dict = {}
    t0 = time.time()
    T = r3data.build_targets(si, ki, train_mask, train_k, variant=cfg["variant"], par=par, log=tlog,
                             loader=label_loader)
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
    if cfg["pretrain_steps"] > 0 or cfg["distill_aux"] > 0:
        t0 = time.time()
        D = cantley_targets(si, ki)
        G.D = torch.from_numpy(D).to(device)                 # (N, K, 2) fp16, NaN = no target
        info["distill"] = {"n_targets": int(np.isfinite(D).sum()), "t_load_s": round(time.time() - t0, 1),
                           "kinases_with_target": int(np.isfinite(D[:, :, 0]).any(0).sum())}
        del D
    wk_t = torch.as_tensor(wk, dtype=torch.float32, device=device)
    eval_k = [k for k in range(K) if trained[k] and not ki.is_canary[k]]
    seeds = [cfg["seed"] + i for i in range(cfg["n_seeds"])]
    models, per_seed = [], []
    t_train = time.time()
    for sd in seeds:
        model, si_info = _train_one(cfg, sd, si, G, wk_t, lab_rows, fit_rows, val_rows, eval_k, device, log,
                                    pre_rows=tr_idx)
        models.append(model)
        per_seed.append(si_info)
    info.update(seeds=seeds, per_seed=per_seed, history=[h for p in per_seed for h in p["history"]],
                best_val_auroc_hard=float(np.mean([p["best_val_auroc_hard"] for p in per_seed])),
                best_step=[p["best_step"] for p in per_seed], t_train_s=round(time.time() - t_train, 1))
    if len(models) > 1:
        m_ens = _ensemble_val(models, G, val_rows, eval_k)
        info["ensemble_val_auroc_hard"] = m_ens
        log(f"ensemble of {len(models)} seeds: inner val_auroc_hard {m_ens:.4f}")
    if not score:
        del G
        torch.cuda.empty_cache()
        return None, info, models
    t0 = time.time()
    info["t_train_wall_s"] = round(t0 - t_start, 1)      # data/targets/GPU load + all seeds' training
    srows = np.where(score_mask)[0]
    S = sum(score_rows(m, G, srows) for m in models) / len(models)     # mean logit over seeds (n, K)
    S = S.cpu().numpy()
    ok = ki.class_ok[si.res_code[srows]] & score_k[None]
    if deploy:
        info["t_score_s"] = round(time.time() - t0, 1)
        info["n_pred_rows"] = int(ok.sum())
        info["t_total_s"] = round(time.time() - t_start, 1)
        info["peak_gpu_mem_gb"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
        del G
        torch.cuda.empty_cache()
        return {"S": S, "site_idx": si.site_idx[srows], "ok": ok}, info, models
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
    return preds, info, models


def distill_loss(model, G, idx_t):
    """Soft-label BCE between the distillation head and Cantley universe percentiles (both variants) on all
    (site, kinase) entries of the batch that have an own-matrix target. Residue class is implied by the target
    (NaN outside the kinase's class)."""
    with torch.autocast("cuda", dtype=torch.bfloat16):
        Ld = model.distill(G.site_batch(idx_t), G.kin_inputs()).float()       # (B, K, 2)
    Y = G.D[idx_t].float()
    ok = torch.isfinite(Y)
    Y = torch.nan_to_num(Y, 0.0)
    l = torch.nn.functional.binary_cross_entropy_with_logits(Ld, Y, reduction="none")
    return (l * ok).sum() / ok.sum().clamp(min=1)


def _train_one(cfg, seed, si, G, wk_t, lab_rows, fit_rows, val_rows, eval_k, device, log, pre_rows=None):
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    ki = _CACHE["ki"]
    model = R3Net(si.emb.shape[1], G.k_cant.shape[1], len(ki.group_names), len(ki.family_names), flank=cfg["flank"],
                  d_site=cfg["d_site"], h=cfg["h"], d_kin=cfg["d_kin"], n_layers=cfg["n_layers"],
                  desc_dropout=cfg["desc_dropout"], dropout=cfg["dropout"],
                  esm_dropout=cfg["esm_dropout"], use_esm=cfg["use_esm"], use_prot=cfg["use_prot"],
                  use_phos=cfg["use_phos"], d_prot=si.prot_mean.shape[1],
                  d_kin_esm=G.k_esm.shape[1] - 1).to(device)
    B, nl = cfg["batch"], int(cfg["batch"] * cfg["p_lab"])
    hist = []
    t_train = time.time()
    pre_info = {}
    if cfg["pretrain_steps"] > 0:
        # Stage 1: Cantley-atlas distillation on TRAINING-fold U_ms sites only (label-free targets)
        opt = torch.optim.AdamW(model.parameters(), lr=cfg["pretrain_lr"], weight_decay=cfg["wd"])
        n = cfg["pretrain_steps"]
        sched = torch.optim.lr_scheduler.LambdaLR(
            opt, lambda s: min(1.0, (s + 1) / cfg["warmup"]) * max(0.05, 1 - s / n))
        for step in range(1, n + 1):
            model.train()
            idx_t = torch.as_tensor(rng.choice(pre_rows, B), device=device)
            ld = distill_loss(model, G, idx_t)
            opt.zero_grad(set_to_none=True)
            ld.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
            if step % 500 == 0 or step == n:
                log(f"seed {seed} pretrain step {step} distill_bce {float(ld):.4f} {time.time() - t_train:.0f}s")
        pre_info = {"pretrain_final_distill_bce": float(ld), "pretrain_s": round(time.time() - t_train, 1)}
        if eval_k:
            m0, _ = inner_val_metric(model, G, val_rows, eval_k)
            pre_info["val_auroc_hard_after_pretrain"] = m0      # main head untrained here; diagnostic only
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / cfg["warmup"]))
    best, best_state, best_step, bad = -1.0, None, 0, 0
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
        if cfg["distill_aux"] > 0:
            loss = loss + cfg["distill_aux"] * distill_loss(model, G, idx_t)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % cfg["eval_every"] == 0:
            m, nk = inner_val_metric(model, G, val_rows, eval_k)
            hist.append({"seed": seed, "step": step, "loss_bce": float(loss_b), "loss_nce": float(loss_n),
                         "val_auroc_hard": m, "val_n_kinases": nk, "t": round(time.time() - t_train, 1)})
            log(f"seed {seed} step {step} bce {float(loss_b):.4f} nce {float(loss_n):.4f} val_auroc_hard {m:.4f} "
                f"({nk} k) {time.time() - t_train:.0f}s")
            if m > best:
                best, best_step, bad = m, step, 0
                best_state = copy.deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()})
            else:
                bad += 1
                if bad >= cfg["patience"]:
                    break
    model.load_state_dict(best_state)
    model.eval()
    return model, {"seed": seed, "best_val_auroc_hard": best, "best_step": best_step, "history": hist,
                   "t_train_s": round(time.time() - t_train, 1), **pre_info}


@torch.no_grad()
def _ensemble_val(models, G, val_rows, eval_k):
    L = sum(score_rows(m, G, val_rows) for m in models) / len(models)
    T = G.T[torch.as_tensor(val_rows, device=G.dev)]
    hard_row = (T == 1).any(1)
    vals = []
    for k in eval_k:
        p = T[:, k] == 1
        n = (T[:, k] == 0) & hard_row
        if p.sum() >= 3 and n.sum() >= 3:
            vals.append(_auroc_t(L[p, k], L[n, k]))
    return float(np.mean(vals)) if vals else float("nan")
