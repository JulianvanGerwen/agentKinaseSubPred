"""R3 network (DESIGN §8.4; ALG-4/5). Owner: R3-neural.

Site encoder  : per window position t in [-w, w]: Linear(ESM-2 650M per-residue embedding, full-protein context)
                + amino-acid embedding + phospho-channel embedding (1 = that flanking residue is a U_ms phosphosite;
                the centre's own status is never used: forced 0) + position embedding; one extra protein-context
                token (Linear(ESM-2 protein mean)). 2-layer pre-norm transformer; output = [centre ; masked mean].
Kinase encoder: blocks {ESM-2 kinase-domain embedding (+missing flag) | Cantley score-profile PCA (+missing/imputed
                flags) | Manning group | Manning family}, each projected to d_k and zeroed independently per kinase
                with p = desc_dropout during training (descriptor/modality dropout), summed, LayerNorm, MLP.
Head          : bilinear <W z_s, u_k>/sqrt(h) + MLP(z_s * u_k) + kinase bias b(u_k). Scores are logits for every
                (site, kinase) pair of a batch of sites at once (B x K).
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class SiteEncoder(nn.Module):
    def __init__(self, d_esm=1280, d=192, h=128, n_layers=2, n_heads=4, flank=7, n_aa=22, dropout=0.1,
                 esm_dropout=0.1, use_esm=True, use_prot=True, use_phos=True):
        super().__init__()
        self.flank = flank
        self.use_esm, self.use_prot, self.use_phos = use_esm, use_prot, use_phos
        L = 2 * flank + 1
        self.proj = nn.Linear(d_esm, d)
        self.esm_drop = nn.Dropout(esm_dropout)
        self.aa = nn.Embedding(n_aa, d)
        self.phos = nn.Embedding(2, d)
        self.pos = nn.Parameter(torch.zeros(L + 1, d))
        nn.init.normal_(self.pos, std=0.02)
        self.prot = nn.Linear(d_esm, d)
        self.pad = nn.Parameter(torch.zeros(d))
        layer = nn.TransformerEncoderLayer(d, n_heads, 2 * d, dropout, batch_first=True, norm_first=True,
                                           activation="gelu")
        self.enc = nn.TransformerEncoder(layer, n_layers, enable_nested_tensor=False)
        self.norm = nn.LayerNorm(d)
        self.out = nn.Sequential(nn.Linear(2 * d, h), nn.GELU(), nn.Linear(h, h), nn.LayerNorm(h))

    def forward(self, emb, aa, phos, valid, prot):
        """emb (B,L,1280) float, aa (B,L) long, phos (B,L) long (centre already 0), valid (B,L) bool,
        prot (B,1280)."""
        x = self.aa(aa) + self.phos(phos if self.use_phos else torch.zeros_like(phos))
        if self.use_esm:
            x = x + self.proj(self.esm_drop(emb))
        x = torch.where(valid[..., None], x, self.pad.to(x.dtype).expand_as(x))
        p = self.prot(prot)[:, None] if self.use_prot else torch.zeros_like(x[:, :1])
        x = torch.cat([p, x], 1) + self.pos[None]
        kpm = torch.cat([torch.zeros_like(valid[:, :1]), ~valid], 1)
        hdn = self.norm(self.enc(x, src_key_padding_mask=kpm))
        c = hdn[:, 1 + self.flank]
        w = valid.float()[..., None]
        m = (hdn[:, 1:] * w).sum(1) / w.sum(1).clamp(min=1)
        return self.out(torch.cat([c, m], -1))


class KinaseEncoder(nn.Module):
    def __init__(self, d_esm, d_cantley, n_groups, n_families, d=256, h=128, desc_dropout=0.2):
        super().__init__()
        self.p = desc_dropout
        self.esm = nn.Linear(d_esm + 1, d)
        self.cant = nn.Linear(d_cantley, d)
        self.grp = nn.Embedding(n_groups + 1, d, padding_idx=n_groups)
        self.fam = nn.Embedding(n_families + 1, d, padding_idx=n_families)
        self.n_groups, self.n_families = n_groups, n_families
        self.norm = nn.LayerNorm(d)
        self.mlp = nn.Sequential(nn.Linear(d, d), nn.GELU(), nn.Dropout(0.1), nn.Linear(d, h), nn.LayerNorm(h))

    def forward(self, esm, cant, grp, fam):
        blocks = [self.esm(esm), self.cant(cant), self.grp(grp), self.fam(fam)]
        if self.training and self.p > 0:
            K = esm.shape[0]
            keep = (torch.rand(len(blocks), K, 1, device=esm.device) >= self.p).to(blocks[0].dtype)
            blocks = [b * keep[i] for i, b in enumerate(blocks)]
        return self.mlp(self.norm(sum(blocks)))


class R3Net(nn.Module):
    def __init__(self, d_esm, d_cantley, n_groups, n_families, flank=7, d_site=192, h=128, d_kin=256,
                 n_layers=2, desc_dropout=0.2, dropout=0.1, esm_dropout=0.1, use_esm=True, use_prot=True, use_phos=True):
        super().__init__()
        self.site = SiteEncoder(d_esm, d_site, h, n_layers, 4, flank, dropout=dropout, esm_dropout=esm_dropout,
                                use_esm=use_esm, use_prot=use_prot, use_phos=use_phos)
        self.kin = KinaseEncoder(d_esm, d_cantley, n_groups, n_families, d_kin, h, desc_dropout)
        self.W = nn.Linear(h, h, bias=False)
        self.pair = nn.Sequential(nn.Linear(h, 64), nn.GELU(), nn.Linear(64, 1))
        self.kbias = nn.Linear(h, 1)
        self.h = h

    def forward(self, site_inputs, kin_inputs):
        z = self.site(*site_inputs)                 # (B,h)
        u = self.kin(*kin_inputs)                   # (K,h)
        bil = self.W(z) @ u.T / math.sqrt(self.h)   # (B,K)
        prod = z[:, None, :] * u[None, :, :]        # (B,K,h)
        return bil + self.pair(prod).squeeze(-1) + self.kbias(u).squeeze(-1)[None]


def bce_masked(logits, T, wk):
    """Balanced masked BCE: sum_pos w_k*bce / sum_pos w_k + mean_neg bce. T int8 (B,K) in {-1,0,1}."""
    pos = T == 1
    neg = T == 0
    lp = F.softplus(-logits)      # -log sigmoid
    ln = F.softplus(logits)       # -log(1-sigmoid)
    w = wk[None].expand_as(logits)
    lpos = (lp * w * pos).sum() / (w * pos).sum().clamp(min=1e-6)
    lneg = (ln * neg).sum() / neg.sum().clamp(min=1)
    return lpos + lneg


def infonce_filtered(logits, T, par, wk):
    """Listwise softmax over kinases per positive (site i, kinase k): denominator = {k} U {j: T[i,j]==0 (candidate
    negative of matching class, not masked), j not in Par(k)}; other true kinases of i (T==1) and masked entries
    (T==-1) are excluded. par: (K,K) bool."""
    ii, kk = torch.nonzero(T == 1, as_tuple=True)
    if len(ii) == 0:
        return logits.sum() * 0
    rows = logits[ii]                                   # (P,K)
    valid = (T[ii] == 0) & ~par[kk]
    valid[torch.arange(len(kk), device=logits.device), kk] = True
    lse = torch.logsumexp(rows.masked_fill(~valid, float("-inf")), dim=1)
    nll = lse - rows[torch.arange(len(kk), device=logits.device), kk]
    w = wk[kk]
    return (nll * w).sum() / w.sum()
