"""Character-level tokenizer for Phosformer (paper §2.3.1-2.3.2).

Vocabulary layout: <s>=0, <pad>=1, </s>=2, <unk>=3, <mask>=4, then printable ASCII chr(33)..chr(122) -> 5..94.
This is the same layout as the released checkpoint's vocab (verified in tests against vocab.json), so released
weights can be scored with this tokenizer and our own checkpoints share the embedding table shape (95).

Pair format (paper §2.3.2: "11-mer flanked by two special tokens, followed by the kinase domain sequence which
is also flanked by two special tokens"); the RoBERTa pair convention realises this as
    <s> p1..p11 </s> </s> k1..kL </s>
so the central peptide residue (p6) is token index 6.
Peptide padding (site within 5 residues of a terminus): '-' or '_' -> <pad> token (paper: "Dash characters ...
were encoded using the <pad> special token"). '.' -> <mask> (used by kinase-masking augmentation).
"""
from __future__ import annotations

import numpy as np

BOS, PAD, EOS, UNK, MASK = 0, 1, 2, 3, 4
VOCAB_SIZE = 95
CENTER_INDEX = 6  # <s> + 5 residues before the site
AA20 = "ACDEFGHIKLMNPQRSTVWY"

_LUT = np.full(256, UNK, dtype=np.int64)
for _c in range(33, 123):
    _LUT[_c] = _c - 28
_LUT[ord("-")] = PAD
_LUT[ord("_")] = PAD
_LUT[ord(".")] = MASK
AA20_IDS = np.array([ord(c) - 28 for c in AA20], dtype=np.int64)


def vocab_dict() -> dict:
    v = {"<s>": BOS, "<pad>": PAD, "</s>": EOS, "<unk>": UNK, "<mask>": MASK}
    v.update({chr(c): c - 28 for c in range(33, 123)})
    return v


def encode_seq(seq: str) -> np.ndarray:
    """Residue string -> token ids (no specials). Lower-case is upper-cased."""
    b = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    return _LUT[b]


def encode_pair(peptide: str, kinase: str) -> np.ndarray:
    if len(peptide) != 11:
        raise ValueError(f"peptide must be an 11-mer, got {len(peptide)}: {peptide!r}")
    return np.concatenate([[BOS], encode_seq(peptide), [EOS, EOS], encode_seq(kinase), [EOS]])


def collate_pairs(id_lists: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Right-pad a list of token-id arrays. Returns (input_ids, attention_mask).

    Note: <pad> tokens *inside* the peptide get attention_mask=1 (as text-encoded tokens) — this is what a
    RoBERTa tokenizer produces for a literal '<pad>' in the text, and RoBERTa's position-id rule
    (positions skip pad ids) applies to them, as in the released model's input pipeline."""
    n = max(len(x) for x in id_lists)
    ids = np.full((len(id_lists), n), PAD, dtype=np.int64)
    att = np.zeros((len(id_lists), n), dtype=np.int64)
    for i, x in enumerate(id_lists):
        ids[i, : len(x)] = x
        att[i, : len(x)] = 1
    return ids, att
