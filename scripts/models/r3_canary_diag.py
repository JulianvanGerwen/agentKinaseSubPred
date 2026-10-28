"""R3 canary diagnostic: are negative-canary AUROCs > 0.5 explained by LABEL-FREE site properties that the canary
sampler (random proteins, then 1-4 sites per protein) makes informative? Scores each negative canary with label-free,
kinase-agnostic site statistics on the pooled DEV CV universe via ksp.eval:
  inv_prot_ums   = -(number of U_ms sites on the site's protein)
  inv_flank_ums  = -(number of U_ms sites within +-7, centre excluded)   (what the R3 phospho channel sees)
  inv_prot_len   = -(protein length)
A model that encodes protein-level phosphosite density (R3: protein-mean ESM token + phospho channel) can reproduce
these signals legitimately; this is a property of the canary design, not label leakage.
Output: results/models/r3/canary_diagnostic.tsv
"""
import sys
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
import numpy as np, pandas as pd  # noqa: E402
from ksp import eval as kev, index, canary  # noqa: E402
from ksp.labels import _canonical_seqs  # noqa: E402

S = index.load_sites()
ums = S[S.in_U_ms]
n_prot = S.acc.map(ums.groupby("acc").size()).fillna(0).values
key = set(zip(ums.acc, ums.pos))
fl = np.array([sum((a, p + o) in key for o in range(-7, 8) if o) for a, p in zip(S.acc.values, S.pos.values)])
seqs = _canonical_seqs()
plen = S.acc.map({a: len(s) for a, s in seqs.items()}).values
K = canary.kinases_with_canaries()
can_idx = np.where(K.uniprot_acc.str.startswith("CANARY_").values)[0]
data = kev.prepare("cv", fold=None)
rows = []
for name, v in [("inv_prot_ums", -n_prot), ("inv_flank_ums", -fl), ("inv_prot_len", -plen)]:
    parts = []
    for k in can_idx:
        cls = K.residue_class.values[k]
        ok = data.universe & kev._res_ok(S.res.values, cls)
        si = np.where(ok)[0]
        parts.append(pd.DataFrame({"kinase_idx": k, "site_idx": si, "score": v[si].astype(float)}))
    P = pd.concat(parts, ignore_index=True)
    res = kev.evaluate(P, data, kinases=list(K.uniprot_acc.values[can_idx]), missing="drop")
    c = res.canary
    rows.append({"score": name, **{k: c[k] for k in ["neg_mean_auroc_hard", "neg_mean_auroc_all", "pos_auroc_hard",
                                                     "pos_auroc_all"]}})
    print(rows[-1], flush=True)
out = pd.DataFrame(rows)
(ROOT / "results/models/r3").mkdir(parents=True, exist_ok=True)
out.to_csv(ROOT / "results/models/r3/canary_diagnostic.tsv", sep="\t", index=False)
print(out)
