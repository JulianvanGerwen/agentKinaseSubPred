"""Phosformer model (paper §2.3): RoBERTa-style encoder (6 layers, 12 heads, 768 d) + head on the central
peptide residue + softmax over 2 classes; focal loss (gamma=2) for fine-tuning.

The encoder is `transformers.RobertaModel` (library). The head is ours: dropout -> Linear(768,768) -> tanh ->
dropout -> Linear(768,2), applied to the final hidden state at token index 6 (paper: "the position corresponds to
the middle residue ... is directed into a feed-forward layer which is followed by the final binary prediction").
Parameter names (`classifier.dense`, `classifier.out_proj`) follow the standard HF RobertaClassificationHead so a
released checkpoint's state dict loads 1:1.
"""
from __future__ import annotations

from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
from transformers import RobertaConfig, RobertaModel

from .tokens import BOS, EOS, PAD, VOCAB_SIZE, CENTER_INDEX


def phosformer_config(**over) -> RobertaConfig:
    """Encoder config from paper §2.3.1 (6 layers x 12 heads x 768). intermediate 3072, GELU, dropout 0.1 are
    RoBERTa defaults (not stated in the paper; match the released config.json). max_position_embeddings=4000 as
    in the released config (paper silent)."""
    kw = dict(vocab_size=VOCAB_SIZE, hidden_size=768, num_hidden_layers=6, num_attention_heads=12,
              intermediate_size=3072, hidden_act="gelu", hidden_dropout_prob=0.1, attention_probs_dropout_prob=0.1,
              max_position_embeddings=4000, type_vocab_size=1, layer_norm_eps=1e-12, initializer_range=0.02,
              pad_token_id=PAD, bos_token_id=BOS, eos_token_id=EOS, position_embedding_type="absolute",
              num_labels=2)
    kw.update(over)
    return RobertaConfig(**kw)


class CenterHead(nn.Module):
    def __init__(self, hidden: int, n_labels: int = 2, dropout: float = 0.1):
        super().__init__()
        self.dense = nn.Linear(hidden, hidden)
        self.dropout = nn.Dropout(dropout)
        self.out_proj = nn.Linear(hidden, n_labels)

    def forward(self, h_center: torch.Tensor) -> torch.Tensor:
        x = self.dropout(h_center)
        x = torch.tanh(self.dense(x))
        x = self.dropout(x)
        return self.out_proj(x)


class PhosformerClassifier(nn.Module):
    def __init__(self, config: RobertaConfig, attn_implementation: str = "sdpa"):
        super().__init__()
        config._attn_implementation = attn_implementation
        self.config = config
        self.roberta = RobertaModel(config, add_pooling_layer=False)
        self.classifier = CenterHead(config.hidden_size, config.num_labels, config.hidden_dropout_prob)

    def forward(self, input_ids, attention_mask):
        h = self.roberta(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        return self.classifier(h[:, CENTER_INDEX])  # logits (B, 2)

    @torch.no_grad()
    def prob(self, input_ids, attention_mask):
        return torch.softmax(self.forward(input_ids, attention_mask).float(), dim=-1)[:, 1]

    # ---- constructors -------------------------------------------------------------------------------------
    @classmethod
    def from_released(cls, path: str | Path, attn_implementation: str = "sdpa") -> "PhosformerClassifier":
        """Load the released fine-tuned checkpoint (external/phosformer/hf_phosformer) — inference only."""
        path = Path(path)
        cfg = RobertaConfig.from_pretrained(path)
        m = cls(cfg, attn_implementation)
        sd = torch.load(path / "pytorch_model.bin", map_location="cpu", weights_only=True)
        sd = {k: v for k, v in sd.items() if not k.endswith("position_ids")}
        missing, unexpected = m.load_state_dict(sd, strict=False)
        if missing or unexpected:
            raise RuntimeError(f"released checkpoint mismatch: missing={missing} unexpected={unexpected}")
        return m

    @classmethod
    def from_mlm(cls, mlm_state: dict, config: RobertaConfig, attn_implementation: str = "sdpa"):
        """Initialise the encoder from our MLM checkpoint (RobertaForMaskedLM state dict); new head."""
        m = cls(config, attn_implementation)
        enc = {k[len("roberta."):]: v for k, v in mlm_state.items() if k.startswith("roberta.")}
        enc = {k: v for k, v in enc.items() if not k.endswith("position_ids") and not k.startswith("pooler.")}
        missing, unexpected = m.roberta.load_state_dict(enc, strict=False)
        missing = [k for k in missing if not k.endswith("position_ids")]
        if missing or unexpected:
            raise RuntimeError(f"MLM checkpoint mismatch: missing={missing} unexpected={unexpected}")
        return m


def focal_loss(logits: torch.Tensor, target: torch.Tensor, gamma: float = 2.0) -> torch.Tensor:
    """Softmax focal loss (Lin et al. 2017), FL(p_t) = -(1-p_t)^gamma log p_t, mean over the batch.
    No alpha weighting (paper gives only gamma)."""
    logp = F.log_softmax(logits.float(), dim=-1)
    logp_t = logp.gather(1, target.view(-1, 1)).squeeze(1)
    p_t = logp_t.exp()
    return (-(1.0 - p_t) ** gamma * logp_t).mean()
