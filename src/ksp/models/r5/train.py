"""R5 training / scoring. Owner: R5 agent.

Same outer protocol as R3 (src/ksp/models/r3/train.py), so the two are directly comparable:
  CV fold f : training rows = U_ms sites on DEV proteins of folds != f; scored rows = U_ms sites of fold f.
  LKO fold j: training rows = all DEV U_ms sites with the held-out kinases' labels removed; scored = those kinases.
Early stopping uses ONLY an inner-validation split of the training rows (ksp.splits.train_valid_split_inner),
macro AUROC on hard negatives over non-canary kinases with >= 3 inner-validation positives. Because that metric only
ever looks at label-carrying rows (positives, and negatives restricted to rows carrying some positive), R5 scores
only those rows for it -- identical value, much cheaper than R3's full-validation scoring.

Loss = balanced masked BCE + lambda_nce * filtered listwise InfoNCE over kinases per positive
       (+ distill_aux * Cantley-percentile distillation, optionally preceded by a distillation pretraining stage).
"""
from __future__ import annotations

import copy
import time

import numpy as np
import pandas as pd
import torch

from ksp import labels, splits
from ksp.models.r3 import train as r3train
from ksp.models.r5 import data as r5data
from ksp.models.r5.model import R5Net, bce_masked, infonce_filtered

TRUNKS = {"35m": "facebook/esm2_t12_35M_UR50D", "150m": "facebook/esm2_t30_150M_UR50D"}

DEFAULT = dict(
    trunk="35m", lr=3e-4, lr_trunk=2e-5, wd=0.05, alpha=0.5, p_lab=0.25, batch=512, max_steps=3000,
    eval_every=100, patience=10, warmup=100, lambda_nce=0.5, desc_dropout=0.2, dropout=0.1, esm_dropout=0.3,
    d_li=64, h=128, d_kin=256, d_k=192, m_query=8, kin_layers=2,
    use_esm650=True, use_phos=True, use_prot=True, use_li=True, use_pair=True,
    pretrain_steps=0, pretrain_batch=256, pretrain_lr=3e-4, distill_aux=0.0,
    valid_frac=0.2, variant="L2", seed=0, n_seeds=3, esm="650m", score_chunk=1024, emb_on_gpu=True,
)

_CACHE: dict = {}


def inputs(esm: str = "650m"):
    if _CACHE.get("esm") != esm:
        _CACHE["s5"] = r5data.load_sites(esm)
        _CACHE["esm"] = esm
    if "k5" not in _CACHE:
        _CACHE["k5"] = r5data.load_kinases()
        _CACHE["par"] = labels.paralog_sets()
    return _CACHE["s5"], _CACHE["k5"], _CACHE["par"]


class GPUData:
    """Everything the model reads, on the GPU (or pinned CPU for the 9.4 GiB frozen-650M window store)."""

    def __init__(self, s5, k5, T, par, dev, fam_keep, grp_keep, emb_on_gpu=True):
        si, ki = s5.si, k5.ki
        self.dev = dev
        self.flank_cache = si.info["flank_max"]            # cached frozen-650M window half-width (10)
        self.L = s5.tok.shape[1]                           # 33 = CLS + 31 + EOS
        self.centre = 1 + r5data.FLANK                     # token position of the centre residue
        self.e650_lo = 1 + r5data.FLANK - self.flank_cache  # first token position covered by the cached window
        t0 = time.time()
        self.emb_on_gpu = emb_on_gpu
        E = torch.from_numpy(np.ascontiguousarray(si.emb))
        self.E = E.to(dev) if emb_on_gpu else E
        self.load_s = time.time() - t0
        self.rows = torch.from_numpy(si.rows).to(dev)
        self.tok = torch.from_numpy(s5.tok.astype(np.int64)).to(dev)
        self.attn = torch.from_numpy(s5.attn).to(dev)
        self.phos = torch.from_numpy(s5.phos.astype(np.int64)).to(dev)
        self.prow = torch.from_numpy(si.prot_row).to(dev)
        self.pm = torch.from_numpy(si.prot_mean).to(dev)
        self.T = torch.from_numpy(T).to(dev) if T is not None else None
        self.kin_tok = torch.from_numpy(k5.tok).to(dev)
        self.kin_mask = torch.from_numpy(k5.mask).to(dev)
        self.k_esm = torch.from_numpy(np.c_[ki.esm, ki.esm_missing]).float().to(dev)
        self.k_cant = torch.from_numpy(np.c_[ki.cantley, ki.cantley_missing, ki.cantley_imputed,
                                             ki.cantley_dist]).float().to(dev)
        g = np.where(grp_keep[ki.group], ki.group, len(ki.group_names))
        f = np.where((ki.family >= 0) & fam_keep[np.maximum(ki.family, 0)], ki.family, len(ki.family_names))
        self.k_grp = torch.from_numpy(g).long().to(dev)
        self.k_fam = torch.from_numpy(f).long().to(dev)
        K = len(ki.kinases)
        P = np.zeros((K, K), bool)
        kp = {a: i for i, a in enumerate(ki.kinases.uniprot_acc.values)}
        for a, ss in par.items():
            for b in ss:
                if a in kp and b in kp:
                    P[kp[a], kp[b]] = True
        self.par = torch.from_numpy(P).to(dev)

    def site_batch(self, idx):
        B = len(idx)
        r = self.rows[idx]
        valid = r >= 0
        rr = r.clamp(min=0)
        if self.emb_on_gpu:
            win = self.E[rr].float()
        else:
            win = self.E[rr.cpu()].to(self.dev, non_blocking=True).float()
        win = win * valid[..., None]
        e650 = win.new_zeros((B, self.L, win.shape[-1]))
        e650[:, self.e650_lo: self.e650_lo + win.shape[1]] = win
        centre = torch.full((B,), self.centre, device=self.dev, dtype=torch.long)
        return (self.tok[idx], self.attn[idx], self.phos[idx], e650, self.pm[self.prow[idx]], centre)

    def kin_desc(self):
        return self.k_esm, self.k_cant, self.k_grp, self.k_fam


@torch.no_grad()
def kinase_cache(model, G):
    with torch.autocast("cuda", dtype=torch.bfloat16):
        Q, u = model.kinase(G.kin_tok, G.kin_mask, G.kin_desc())
    return Q.float(), u.float()


@torch.no_grad()
def score_rows(model, G, rows_idx, chunk=1024, Qu=None):
    model.eval()
    Q, u = Qu if Qu is not None else kinase_cache(model, G)
    out = []
    for a in range(0, len(rows_idx), chunk):
        idx = torch.as_tensor(rows_idx[a:a + chunk], device=G.dev)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            so = model.site(*G.site_batch(idx))
            out.append(model.head(so, Q.to(torch.bfloat16), u.to(torch.bfloat16)).float())
    return torch.cat(out)


def inner_val_metric(model, G, val_rows, eval_k, chunk=1024):
    L = score_rows(model, G, val_rows, chunk)
    T = G.T[torch.as_tensor(val_rows, device=G.dev)]
    hard_row = (T == 1).any(1)
    vals = []
    for k in eval_k:
        p = T[:, k] == 1
        n = (T[:, k] == 0) & hard_row
        if p.sum() >= 3 and n.sum() >= 3:
            vals.append(r3train._auroc_t(L[p, k], L[n, k]))
    return (float(np.mean(vals)) if vals else float("nan")), len(vals)


def distill_loss(model, G, idx_t, z=None, u=None):
    """Soft-label BCE against the Cantley universe percentiles. During fine-tuning the site/kinase vectors of the
    main forward pass are reused (z, u given), so the auxiliary loss costs one small head, not a second trunk pass."""
    with torch.autocast("cuda", dtype=torch.bfloat16):
        if z is None:
            Ld = model.distill(G.site_batch(idx_t), G.kin_tok, G.kin_mask, G.kin_desc()).float()
        else:
            Ld = model.distill_head(z, u).float()
    Y = G.D[idx_t].float()
    ok = torch.isfinite(Y)
    Y = torch.nan_to_num(Y, 0.0)
    l = torch.nn.functional.binary_cross_entropy_with_logits(Ld, Y, reduction="none")
    return (l * ok).sum() / ok.sum().clamp(min=1)


def _build_model(cfg, s5, k5, G, device):
    from transformers import EsmModel
    trunk = EsmModel.from_pretrained(TRUNKS[cfg["trunk"]], add_pooling_layer=False)
    ki = k5.ki
    return R5Net(trunk, d_cantley=G.k_cant.shape[1], n_groups=len(ki.group_names), n_families=len(ki.family_names),
                 d_kin_esm=G.k_esm.shape[1] - 1, d_prot=s5.si.prot_mean.shape[1], d_esm650=s5.si.emb.shape[1],
                 d_li=cfg["d_li"], h=cfg["h"], d_kin=cfg["d_kin"], m_query=cfg["m_query"],
                 kin_layers=cfg["kin_layers"], d_k=cfg["d_k"], desc_dropout=cfg["desc_dropout"],
                 dropout=cfg["dropout"], esm_dropout=cfg["esm_dropout"], use_esm650=cfg["use_esm650"],
                 use_phos=cfg["use_phos"], use_prot=cfg["use_prot"], use_li=cfg["use_li"],
                 use_pair=cfg["use_pair"], d_kin_tok=k5.tok.shape[2], max_len_kin=k5.tok.shape[1]).to(device)


def _opt(model, cfg, lr):
    trunk = [p for n, p in model.named_parameters() if n.startswith("site.trunk.")]
    rest = [p for n, p in model.named_parameters() if not n.startswith("site.trunk.")]
    return torch.optim.AdamW([{"params": trunk, "lr": cfg["lr_trunk"]}, {"params": rest, "lr": lr}],
                             weight_decay=cfg["wd"])


def run(cfg: dict, fold: int | None = None, lko_fold: int | None = None, device: str = "cuda", log=print,
        score: bool = True):
    cfg = {**DEFAULT, **cfg}
    t_start = time.time()
    s5, k5, par = inputs(cfg["esm"])
    si, ki = s5.si, k5.ki
    K = len(ki.kinases)
    if (fold is None) == (lko_fold is None):
        raise ValueError("give exactly one of fold / lko_fold")
    dev_rows = si.cv_fold >= 0
    if fold is not None:
        train_mask = dev_rows & (si.cv_fold != fold)
        score_mask = si.cv_fold == fold
        train_k, score_k = None, np.ones(K, bool)
    else:
        kf = splits.load_kinase_folds()
        held = set(kf.kinase_acc[kf.lko_fold == lko_fold])
        is_held = ki.kinases.uniprot_acc.isin(held).values
        train_mask = score_mask = dev_rows
        train_k, score_k = ~is_held, is_held
    info = {"fold": fold, "lko_fold": lko_fold, "cfg": cfg, "n_train_rows": int(train_mask.sum()),
            "n_score_rows": int(score_mask.sum()), "data_info": {"sites": s5.info["window"], "kinases": k5.info}}
    tlog: dict = {}
    t0 = time.time()
    T = r5data.build_targets(si, ki, train_mask, train_k, variant=cfg["variant"], par=par, log=tlog)
    info["targets"] = tlog
    info["t_targets_s"] = round(time.time() - t0, 1)
    tr_idx = np.where(train_mask)[0]
    val_flag = splits.train_valid_split_inner(si.site_ids[tr_idx], valid_frac=cfg["valid_frac"], seed=cfg["seed"])
    val_rows = tr_idx[val_flag]
    fit_rows = tr_idx[~val_flag]
    val_rows = val_rows[(T[val_rows] == 1).any(1)]            # only label-carrying rows enter the metric
    n_k = (T[fit_rows] == 1).sum(0).astype(np.float64)
    trained = n_k > 0
    wk = np.where(trained, np.maximum(n_k, 1) ** (-cfg["alpha"]), 0.0)
    wk = wk / (wk[trained] * n_k[trained]).sum() * n_k[trained].sum()
    lab_rows = fit_rows[(T[fit_rows] == 1).any(1)]
    info.update(n_fit_rows=len(fit_rows), n_val_rows=len(val_rows), n_lab_fit_rows=len(lab_rows),
                n_trained_kinases=int(trained.sum()), n_pos_fit=int(n_k.sum()))
    fam_keep = np.zeros(len(ki.family_names), bool)
    grp_keep = np.zeros(len(ki.group_names), bool)
    fam_keep[ki.family[trained & (ki.family >= 0)]] = True
    grp_keep[ki.group[trained]] = True
    G = GPUData(s5, k5, T, par, device, fam_keep, grp_keep, emb_on_gpu=cfg["emb_on_gpu"])
    info["t_gpu_load_s"] = round(G.load_s, 1)
    if cfg["pretrain_steps"] > 0 or cfg["distill_aux"] > 0:
        t0 = time.time()
        D = r3train.cantley_targets(si, ki)
        G.D = torch.from_numpy(D).to(device)
        info["distill"] = {"n_targets": int(np.isfinite(D).sum()), "t_load_s": round(time.time() - t0, 1),
                           "kinases_with_target": int(np.isfinite(D[:, :, 0]).any(0).sum())}
        del D
    wk_t = torch.as_tensor(wk, dtype=torch.float32, device=device)
    eval_k = [k for k in range(K) if trained[k] and not ki.is_canary[k]]
    seeds = [cfg["seed"] + i for i in range(cfg["n_seeds"])]
    models, per_seed = [], []
    t_train = time.time()
    for sd in seeds:
        model, sinfo = _train_one(cfg, sd, s5, k5, G, wk_t, lab_rows, fit_rows, val_rows, eval_k, device, log,
                                  pre_rows=tr_idx)
        models.append(model)
        per_seed.append(sinfo)
    info.update(seeds=seeds, per_seed=per_seed, history=[h for p in per_seed for h in p["history"]],
                best_val_auroc_hard=float(np.mean([p["best_val_auroc_hard"] for p in per_seed])),
                best_step=[p["best_step"] for p in per_seed], t_train_s=round(time.time() - t_train, 1))
    if len(models) > 1:
        L = sum(score_rows(m, G, val_rows, cfg["score_chunk"]) for m in models) / len(models)
        Tv = G.T[torch.as_tensor(val_rows, device=device)]
        hard_row = (Tv == 1).any(1)
        vals = [r3train._auroc_t(L[Tv[:, k] == 1, k], L[(Tv[:, k] == 0) & hard_row, k]) for k in eval_k
                if (Tv[:, k] == 1).sum() >= 3 and ((Tv[:, k] == 0) & hard_row).sum() >= 3]
        info["ensemble_val_auroc_hard"] = float(np.mean(vals)) if vals else float("nan")
        log(f"ensemble of {len(models)} seeds: inner val_auroc_hard {info['ensemble_val_auroc_hard']:.4f}")
    if not score:
        del G
        torch.cuda.empty_cache()
        return None, info, models
    t0 = time.time()
    srows = np.where(score_mask)[0]
    S = sum(score_rows(m, G, srows, cfg["score_chunk"]) for m in models) / len(models)
    S = S.cpu().numpy()
    ok = ki.class_ok[si.res_code[srows]] & score_k[None]
    ri, kk = np.nonzero(ok)
    preds = pd.DataFrame({"kinase_idx": kk.astype(np.int64), "site_idx": si.site_idx[srows[ri]],
                          "score": S[ri, kk].astype(np.float32)})
    preds["kinase_acc"] = pd.Categorical.from_codes(kk, categories=ki.kinases.uniprot_acc.values)
    preds["site_id"] = pd.Categorical.from_codes(ri, categories=si.site_ids[srows])
    info["t_score_s"] = round(time.time() - t0, 1)
    info["n_pred_rows"] = len(preds)
    info["t_total_s"] = round(time.time() - t_start, 1)
    info["peak_gpu_mem_gb"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
    del G
    torch.cuda.empty_cache()
    return preds, info, models


def _train_one(cfg, seed, s5, k5, G, wk_t, lab_rows, fit_rows, val_rows, eval_k, device, log, pre_rows=None):
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    model = _build_model(cfg, s5, k5, G, device)
    B, nl = cfg["batch"], int(cfg["batch"] * cfg["p_lab"])
    hist, pre_info = [], {}
    t_train = time.time()
    if cfg["pretrain_steps"] > 0:
        opt = _opt(model, cfg, cfg["pretrain_lr"])
        n = cfg["pretrain_steps"]
        sched = torch.optim.lr_scheduler.LambdaLR(
            opt, lambda s: min(1.0, (s + 1) / cfg["warmup"]) * max(0.05, 1 - s / n))
        for step in range(1, n + 1):
            model.train()
            idx_t = torch.as_tensor(rng.choice(pre_rows, cfg["pretrain_batch"]), device=device)
            ld = distill_loss(model, G, idx_t)
            opt.zero_grad(set_to_none=True)
            ld.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
            if step % 500 == 0 or step == n:
                log(f"seed {seed} pretrain {step}/{n} distill_bce {float(ld):.4f} {time.time()-t_train:.0f}s")
        pre_info = {"pretrain_final_distill_bce": float(ld), "pretrain_s": round(time.time() - t_train, 1)}
    opt = _opt(model, cfg, cfg["lr"])
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / cfg["warmup"]))
    best, best_state, best_step, bad = -1.0, None, 0, 0
    for step in range(1, cfg["max_steps"] + 1):
        model.train()
        idx = np.r_[rng.choice(lab_rows, nl), rng.choice(fit_rows, B - nl)]
        idx_t = torch.as_tensor(idx, device=device)
        Tb = G.T[idx_t]
        with torch.autocast("cuda", dtype=torch.bfloat16):
            so = model.site(*G.site_batch(idx_t))
            Q, u = model.kinase(G.kin_tok, G.kin_mask, G.kin_desc())
            L = model.head(so, Q, u).float()
        z_site, u_kin = so[1], u
        loss_b = bce_masked(L, Tb, wk_t)
        loss = loss_b
        loss_n = torch.zeros((), device=device)
        if cfg["lambda_nce"] > 0:
            loss_n = infonce_filtered(L, Tb, G.par, wk_t)
            loss = loss + cfg["lambda_nce"] * loss_n
        if cfg["distill_aux"] > 0:
            loss = loss + cfg["distill_aux"] * distill_loss(model, G, idx_t, z=z_site, u=u_kin)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % cfg["eval_every"] == 0:
            m, nk = inner_val_metric(model, G, val_rows, eval_k, cfg["score_chunk"])
            hist.append({"seed": seed, "step": step, "loss_bce": float(loss_b), "loss_nce": float(loss_n),
                         "val_auroc_hard": m, "val_n_kinases": nk, "t": round(time.time() - t_train, 1)})
            log(f"seed {seed} step {step} bce {float(loss_b):.4f} nce {float(loss_n):.4f} "
                f"val_auroc_hard {m:.4f} ({nk} k) {time.time()-t_train:.0f}s")
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
