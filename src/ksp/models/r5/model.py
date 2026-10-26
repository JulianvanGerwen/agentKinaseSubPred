"""R5 network: fine-tuned ESM-2 + ColBERT-style late interaction between substrate-window tokens and
kinase-domain tokens (DESIGN §8.4 family; R5 brief).

Why late interaction and not a full cross-encoder: measured in scripts/models/r5_pilot.py
(results/models/r5/pilot_timing.json). A joint [kinase domain] <sep> [window] cross-encoder on one RTX 4090 runs at
1.4-1.6k pairs/s for ESM-2 35M and 0.41-0.45k pairs/s for 150M, i.e. 3.6-4.3 and 13.2-14.3 GPU-h for the ~21.3M
candidate pairs of one CV fold -- far over the ~1.5 GPU-h/fold budget. Late interaction keeps a *token-level*
kinase-substrate interaction but factorises the encoders, so each side is encoded once per fold and only the
MaxSim head runs per pair (measured well under a GPU-minute per fold).

SITE TOWER (fine-tuned). ESM-2 (35M or 150M) over the +-15 window (31 residues, CLS/EOS added). Three zero-initialised
  additive inputs are injected at the embedding layer, so training starts from the unmodified pretrained trunk:
    * frozen ESM-2 650M per-residue embedding of the same residue taken in FULL-PROTEIN context (cache/r3), which is
      exactly R3's site input - the fine-tuned trunk refines it instead of replacing it;
    * phospho-channel embedding (flanking U_ms phosphosites; the centre's own status is never used);
    * the protein-context vector (ESM-2 650M protein mean), broadcast over the window.
  Outputs: per-token late-interaction vectors S (B, L, d_li) and a pooled site vector z (centre ; masked mean).

KINASE TOWER (trainable on cached inputs). Per-residue ESM-2 650M embeddings of the catalytic domain
  (cache/r5/kin_tok.npy) -> linear -> 2-layer transformer -> m learned query tokens attend over the domain ->
  Q (K, m, d_li). The R3 descriptor blocks (kinase-domain mean embedding, Cantley score-profile PCA, Manning
  group/family, descriptor dropout p) give the global kinase vector u (K, h).

HEAD. score(s, k) = gamma * sum_{j<m} max_t <Q[k,j], S[s,t]> / sqrt(d_li)        <- late interaction
                  + <W z_s, u_k>/sqrt(h) + MLP(z_s * u_k) + b(u_k)               <- R3's bilinear pair head
  The second line makes R5 a strict superset of R3's head, so the late-interaction term is an addition that can be
  ablated (config `li_only` / `no_li`).
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from ksp.models.r3.model import KinaseEncoder, bce_masked, infonce_filtered  # noqa: F401  (re-exported losses)


class SiteTower(nn.Module):
    def __init__(self, trunk, d_esm650=1280, d_prot=1280, d_li=64, h=128, use_esm650=True, use_phos=True,
                 use_prot=True, esm_dropout=0.1):
        super().__init__()
        self.trunk = trunk
        d_m = trunk.config.hidden_size
        self.d_m = d_m
        self.use_esm650, self.use_phos, self.use_prot = use_esm650, use_phos, use_prot
        self.esm650 = nn.Linear(d_esm650, d_m)
        self.prot = nn.Linear(d_prot, d_m)
        self.phos = nn.Embedding(2, d_m)
        for mod in (self.esm650, self.prot):                 # zero init: start from the unmodified trunk
            nn.init.zeros_(mod.weight)
            nn.init.zeros_(mod.bias)
        nn.init.zeros_(self.phos.weight)
        self.esm_drop = nn.Dropout(esm_dropout)
        self.li = nn.Linear(d_m, d_li)
        self.pool = nn.Sequential(nn.Linear(2 * d_m, h), nn.GELU(), nn.Linear(h, h), nn.LayerNorm(h))

    def forward(self, tok, attn, phos, e650, prot, centre):
        # ESM-2 scales word embeddings by (1 - 0.15*0.8) when no <mask> token is present (token_dropout
        # compensation in EsmEmbeddings); the inputs_embeds path skips it, so apply it here to keep the trunk at
        # its pretraining operating point. Position ids are built exactly as EsmEmbeddings would from input_ids.
        pad = self.trunk.config.pad_token_id
        m = (tok != pad).long()
        pos_ids = torch.cumsum(m, dim=1) * m + pad
        emb = self.trunk.embeddings.word_embeddings(tok) * (1.0 - 0.15 * 0.8)
        if self.use_esm650:
            emb = emb + self.esm650(self.esm_drop(e650))
        if self.use_phos:
            emb = emb + self.phos(phos)
        if self.use_prot:
            emb = emb + self.prot(prot)[:, None, :]
        hs = self.trunk(inputs_embeds=emb, attention_mask=attn.long(),
                        position_ids=pos_ids).last_hidden_state          # (B, L, d_m)
        w = attn.to(hs.dtype)[..., None]
        m = (hs * w).sum(1) / w.sum(1).clamp(min=1)
        c = hs[torch.arange(hs.shape[0], device=hs.device), centre]
        return self.li(hs), self.pool(torch.cat([c, m], -1)), attn


class KinaseTower(nn.Module):
    """Cached per-residue domain embeddings -> transformer -> m learned query tokens -> (K, m, d_li)."""

    def __init__(self, d_tok=1280, d_k=192, n_layers=2, n_heads=4, m_query=8, d_li=64, dropout=0.1, max_len=320):
        super().__init__()
        self.inp = nn.Linear(d_tok, d_k)
        self.pos = nn.Parameter(torch.zeros(max_len, d_k))
        nn.init.normal_(self.pos, std=0.02)
        layer = nn.TransformerEncoderLayer(d_k, n_heads, 2 * d_k, dropout, batch_first=True, norm_first=True,
                                           activation="gelu")
        self.enc = nn.TransformerEncoder(layer, n_layers, enable_nested_tensor=False)
        self.norm = nn.LayerNorm(d_k)
        self.q = nn.Parameter(torch.randn(m_query, d_k) * 0.02)
        self.attn = nn.MultiheadAttention(d_k, n_heads, batch_first=True)
        self.out = nn.Linear(d_k, d_li)

    def forward(self, tok, mask):
        x = self.inp(tok) + self.pos[None, : tok.shape[1]]
        hs = self.norm(self.enc(x, src_key_padding_mask=~mask))
        q = self.q[None].expand(tok.shape[0], -1, -1)
        a, _ = self.attn(q, hs, hs, key_padding_mask=~mask, need_weights=False)
        return self.out(a)                                   # (K, m, d_li)


class R5Net(nn.Module):
    def __init__(self, trunk, d_cantley, n_groups, n_families, d_kin_esm, d_prot=1280, d_esm650=1280,
                 d_li=64, h=128, d_kin=256, m_query=8, kin_layers=2, d_k=192, desc_dropout=0.2, dropout=0.1,
                 esm_dropout=0.1, use_esm650=True, use_phos=True, use_prot=True, use_li=True, use_pair=True,
                 d_kin_tok=1280, max_len_kin=320):
        super().__init__()
        self.site = SiteTower(trunk, d_esm650, d_prot, d_li, h, use_esm650, use_phos, use_prot, esm_dropout)
        self.kin_tok = KinaseTower(d_kin_tok, d_k, kin_layers, 4, m_query, d_li, dropout, max_len_kin)
        self.kin = KinaseEncoder(d_kin_esm, d_cantley, n_groups, n_families, d_kin, h, desc_dropout)
        self.W = nn.Linear(h, h, bias=False)
        self.pair = nn.Sequential(nn.Linear(h, 64), nn.GELU(), nn.Linear(64, 1))
        self.kbias = nn.Linear(h, 1)
        self.gamma = nn.Parameter(torch.tensor(1.0))
        self.h, self.d_li = h, d_li
        self.use_li, self.use_pair = use_li, use_pair
        # Cantley-distillation head (2 outputs), as in R3: used for pretraining / the auxiliary loss only.
        self.dW = nn.Linear(h, 2 * h, bias=False)
        self.dpair = nn.Sequential(nn.Linear(h, 64), nn.GELU(), nn.Linear(64, 2))
        self.dbias = nn.Linear(h, 2)

    # --- kinase side ------------------------------------------------------------------------------------------
    def kinase(self, kin_tok, kin_mask, kin_desc):
        return self.kin_tok(kin_tok, kin_mask), self.kin(*kin_desc)        # Q (K,m,d_li), u (K,h)

    # --- scoring ----------------------------------------------------------------------------------------------
    def _maxsim(self, S, attn, Q):
        """S (B,L,d), attn (B,L) bool, Q (K,m,d) -> (B,K) summed MaxSim."""
        sim = torch.einsum("bld,kmd->bkml", S, Q) / math.sqrt(self.d_li)
        sim = sim.masked_fill(~attn[:, None, None, :], float("-inf"))
        return sim.amax(-1).sum(-1)

    def head(self, site_out, Q, u):
        S, z, attn = site_out
        out = 0.0
        if self.use_li:
            out = out + self.gamma * self._maxsim(S, attn, Q)
        if self.use_pair:
            bil = self.W(z) @ u.T / math.sqrt(self.h)
            out = out + bil + self.pair(z[:, None, :] * u[None, :, :]).squeeze(-1)
        return out + self.kbias(u).squeeze(-1)[None]

    def forward(self, site_inputs, kin_tok, kin_mask, kin_desc):
        return self.head(self.site(*site_inputs), *self.kinase(kin_tok, kin_mask, kin_desc))

    def distill_head(self, z, u):
        """(B,h) site vectors x (K,h) kinase vectors -> (B,K,2) Cantley-percentile logits."""
        B = z.shape[0]
        bil = torch.einsum("bjh,kh->bkj", self.dW(z).view(B, 2, self.h), u) / math.sqrt(self.h)
        return bil + self.dpair(z[:, None, :] * u[None, :, :]) + self.dbias(u)[None]

    def distill(self, site_inputs, kin_tok, kin_mask, kin_desc):
        _, z, _ = self.site(*site_inputs)
        return self.distill_head(z, self.kin(*kin_desc))
