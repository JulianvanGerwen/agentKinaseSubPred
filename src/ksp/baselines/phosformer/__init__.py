"""B3 Phosformer baseline: clean-room reimplementation from Zhou et al. 2023 (Bioinformatics 39:btad046).

Written from the paper (data/paper2.pdf), not from the CC BY-NC-ND reference code. The model is built from
`transformers.RobertaModel` / `RobertaForMaskedLM` (library classes); the classification head, tokenisation,
data pipeline, focal loss and training loops are our own.

Modules:
- tokens:  character tokenizer (vocab layout identical to the released checkpoint so its weights load) and
           RoBERTa-style pair encoding  <s> peptide </s></s> kinase </s>.
- model:   PhosformerClassifier (encoder + head on the central peptide residue, token index 6), focal loss,
           config for the 6-layer/12-head/768-d encoder.
- corpus:  MLM pretraining corpus (their released kinase-domain + substrate FASTA).
- inputs:  our kinase-domain sequences (cache/index/kinases.parquet + UniProt) and 11-mers from window15.
"""
