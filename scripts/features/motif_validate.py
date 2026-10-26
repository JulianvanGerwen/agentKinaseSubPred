"""Sanity validation of F1 (Cantley) and F2a (learned PSSM) motif features.

NON-OFFICIAL: uses cache/labels/ks_pairs.parquet L2 rows directly and a TEMPORARY protein-level random
80/20 split (seed 0). These numbers are sanity checks only, not project metrics (no paralog masks, no
substrate clusters, random U_ms negatives).
Writes results/features/motif/{validation_auroc_cantley.tsv, validation_learned.tsv, motif_top_positions.tsv,
validation_runtime.json}.
"""
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp.index import load_kinases, load_sites  # noqa: E402
from ksp.store import load_matrix  # noqa: E402
from ksp.features import motif_cantley as mc  # noqa: E402
from ksp.features import motif_learned as ml  # noqa: E402

RES = ROOT / "results/features/motif"
rng = np.random.default_rng(0)
runtime = {}
S, K = load_sites(), load_kinases()
kpos = {a: i for i, a in enumerate(K.uniprot_acc)}
gene2acc = dict(zip(K.gene, K.uniprot_acc))
res = S.res.to_numpy()
ms = S.in_U_ms.to_numpy(bool)
spos = pd.Index(S.site_id)

# ---- labels (NON-OFFICIAL sanity use) ----
P = pd.read_parquet(ROOT / "cache/labels/ks_pairs.parquet")
P = P[P.L2.fillna(False).astype(bool) & P.kinase_included.fillna(False).astype(bool)]
P = P[P.site_id.isin(spos)][["kinase_acc", "site_id"]].drop_duplicates()
P["si"] = spos.get_indexer(P.site_id)
P["ki"] = P.kinase_acc.map(kpos)
allpos = P.groupby("ki").si.apply(set).to_dict()
any_lab = pd.read_parquet(ROOT / "cache/labels/ks_pairs.parquet", columns=["kinase_acc", "site_id"])
any_lab = any_lab[any_lab.site_id.isin(spos)]
known = {kpos[k]: set(spos.get_indexer(g.site_id)) for k, g in any_lab.groupby("kinase_acc") if k in kpos}
cls = K.residue_class.to_numpy()


def type_mask(k):
    rts = ml.CLASS_RTS[cls[k]]
    return np.isin(res, sum((ml.RT_RES[r] for r in rts), ()))


def sample_negatives(k, pool_mask, n):
    cand = np.where(pool_mask & type_mask(k) & ms)[0]
    excl = known.get(k, set())
    cand = cand[~np.isin(cand, list(excl))] if excl else cand
    return rng.choice(cand, size=min(n, len(cand)), replace=False)


# ---- 1. Cantley AUROC on all L2 labels vs random U_ms negatives ----
t = time.time()
arms = {"score": "unprimed_score", "pct_published": "unprimed_pct_published", "pct_universe": "unprimed_pct_universe",
        "primed_hard_pct_universe": "primed_hard_pct_universe", "primed_soft_pct_universe": "primed_soft_pct_universe",
        "primed_hard_pct_published": "primed_hard_pct_published"}
M = {a: load_matrix(f"motif/cantley_{f}")[0] for a, f in arms.items()}
MI = {a: load_matrix(f"motif/cantley_imputed_{f}")[0] for a, f in arms.items()}


def combined(a, k, idx):
    """own-matrix value where available, else imputed (nearest atlas kinase)."""
    o = np.asarray(M[a][k, idx], float)
    return np.where(np.isfinite(o), o, np.asarray(MI[a][k, idx], float))


def matrix_kind(k, idx):
    o = np.isfinite(np.asarray(M["score"][k, idx], float))
    return "own" if o.all() else ("imputed" if not o.any() else "own+imputed")


rows = []
for k, pos in allpos.items():
    pos = np.array(sorted(pos))
    pos = pos[type_mask(k)[pos]]
    if len(pos) < 5:
        continue
    neg = sample_negatives(k, np.ones(len(S), bool), 2000)
    idx = np.r_[pos, neg]
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    if not np.isfinite(combined("score", k, idx)).all():
        continue
    r = dict(gene=K.gene[k], uniprot_acc=K.uniprot_acc[k], n_pos=len(pos), n_neg=len(neg),
             matrix=matrix_kind(k, pos))
    for a in arms:
        r[a] = roc_auc_score(y, combined(a, k, idx))
    rows.append(r)
ca = pd.DataFrame(rows).sort_values("n_pos", ascending=False)
ca.to_csv(RES / "validation_auroc_cantley.tsv", sep="\t", index=False, float_format="%.4f")
runtime["cantley_auroc_s"] = round(time.time() - t, 1)

# ---- 1b. site-centric (M2-like, sanity): rank of the true kinase among own-matrix kinases of the site's type ----
t = time.time()
own_rows = {rt: np.where(np.isfinite(np.asarray(M["score"][:, np.where(np.isin(res, ml.RT_RES[rt]))[0][0]], float)))[0]
            for rt in ml.RTS}
site_true = P.groupby("si").ki.apply(set)
sc_rows = []
cols = np.array(sorted(site_true.index))
for a in ("score", "pct_published", "pct_universe", "primed_soft_pct_universe"):
    X = np.asarray(M[a][:, cols], float)                      # 503 x n_sites
    for j, sidx in enumerate(cols):
        rt = "Y" if res[sidx] == "Y" else "ST"
        rows_ = own_rows[rt]
        trues = site_true[sidx]
        v = X[rows_, j]
        for k in trues:
            if k not in set(rows_):
                continue
            others = np.array([r_ not in trues or r_ == k for r_ in rows_])   # filtered ranking
            vk = X[k, j]
            rank = 1 + int(((v > vk) & others).sum())
            sc_rows.append(dict(arm=a, ki=k, site=sidx, rank=rank, n_cand=int(others.sum())))
sc = pd.DataFrame(sc_rows)
sc["hit5"], sc["hit10"], sc["rr"] = sc["rank"] <= 5, sc["rank"] <= 10, 1 / sc["rank"]
per_k = sc.groupby(["arm", "ki"])[["hit5", "hit10", "rr"]].mean()
site_centric = per_k.groupby("arm").mean().round(4)
site_centric["n_kinases"] = per_k.groupby("arm").size()
site_centric["n_pairs"] = sc.groupby("arm").size()
site_centric.to_csv(RES / "validation_site_centric_cantley.tsv", sep="\t")
runtime["site_centric_s"] = round(time.time() - t, 1)

# ---- 2. learned PSSM on a temporary protein split ----
acc = P.site_id.str.rsplit("_", n=1).str[0]
prots = np.array(sorted(acc.unique()))
test_prot = set(rng.choice(prots, size=len(prots) // 5, replace=False))
is_te = acc.isin(test_prot).to_numpy()
tr, te = P[~is_te], P[is_te]
models = {}
for name, kw in {"family": dict(prior="family"), "group": dict(prior="group"), "none": dict(prior="none"),
                 "family_phospho": dict(prior="family", phospho=True)}.items():
    t = time.time()
    models[name] = ml.fit_pssm(tr[["kinase_acc", "site_id"]], S, K, **kw)
    runtime[f"fit_pssm_{name}_s"] = round(time.time() - t, 2)
    runtime[f"fit_pssm_{name}_alpha_beta"] = [models[name].alpha, models[name].beta]
t = time.time()
_ = ml.score(models["family"])
runtime["score_all_503x350736_s"] = round(time.time() - t, 1)
te_site_mask = np.isin(S.acc.to_numpy(), list(test_prot))
ntr = tr.groupby("ki").size()
rows = []
for k, g in te.groupby("ki"):
    pos = g.si.to_numpy()
    pos = pos[type_mask(k)[pos]]
    if len(pos) < 5:
        continue
    neg = sample_negatives(k, te_site_mask, 1000)
    if len(neg) < 50:
        continue
    idx = np.r_[pos, neg]
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    r = dict(gene=K.gene[k], n_test_pos=len(pos), n_train=int(ntr.get(k, 0)), n_neg=len(neg))
    for name, m in models.items():
        r[f"pssm_{name}"] = roc_auc_score(y, ml.score_pairs(m, np.full(len(idx), k), idx, S))
    r["cantley_matrix"] = matrix_kind(k, pos)
    for a in ("pct_published", "pct_universe", "primed_soft_pct_universe"):
        v = combined(a, k, idx)
        r[f"cantley_{a}"] = roc_auc_score(y, v) if np.isfinite(v).all() else np.nan
    rows.append(r)
lv = pd.DataFrame(rows).sort_values("n_test_pos", ascending=False)
lv["train_bin"] = pd.cut(lv.n_train, [-1, 0, 9, 49, 10**6], labels=["0", "1-9", "10-49", "50+"])
lv.to_csv(RES / "validation_learned.tsv", sep="\t", index=False, float_format="%.4f")

# ---- 3. top-weighted positions: learned vs Cantley ----
At = {t_: mc.load_atlas(t_) for t_ in mc.KIN_TYPES}
info = mc.kinome_info()
full = {"family": ml.fit_pssm(P[["kinase_acc", "site_id"]], S, K, prior="family"),
        "family_phospho": ml.fit_pssm(P[["kinase_acc", "site_id"]], S, K, prior="family", phospho=True)}
tops = []
for gene, rt in [("PRKACA", "ST"), ("CDK1", "ST"), ("CSNK2A1", "ST"), ("ATM", "ST"), ("GSK3B", "ST"),
                 ("MAPK1", "ST"), ("CSNK1D", "ST"), ("PLK1", "ST"), ("SRC", "Y"), ("EGFR", "Y")]:
    a = gene2acc[gene]
    for name, m in full.items():
        d = ml.top_positions(m, a, rt, 6)
        tops.append(dict(gene=gene, source=f"learned_{name} (all L2, n={int(m.n_train[rt][kpos[a]])})",
                         top=" ".join(f"{r.pos:+d}{r.aa}:{r.log2odds:.2f}" for r in d.itertuples())))
    kt = "ser_thr" if rt == "ST" else "tyrosine"
    mn = info[(info.UNIPROT_ID == a) & (info.TYPE == kt)].MATRIX_NAME
    if len(mn):
        A = At[kt]
        L = A.L[A.names.index(mn.iloc[0])]
        cells = sorted(((L[pi, ai], p, mc.AA23[ai]) for pi, p in enumerate(A.positions) for ai in range(23)), reverse=True)[:6]
        tops.append(dict(gene=gene, source=f"cantley {mn.iloc[0]}", top=" ".join(f"{p:+d}{aa}:{v:.2f}" for v, p, aa in cells)))
pd.DataFrame(tops).to_csv(RES / "motif_top_positions.tsv", sep="\t", index=False)

# ---- 4. crossfit-contract smoke test (each protein its own cluster) ----
try:
    from ksp import crossfit
    t = time.time()
    clusters = pd.Series(np.arange(len(prots)), index=prots)
    q = te.sample(min(2000, len(te)), random_state=0)[["kinase_acc", "site_id"]].reset_index(drop=True)
    fam = crossfit.REGISTRY["f2a_pssm"]().fit(tr[["kinase_acc", "site_id"]], crossfit.FitContext(clusters))
    out = fam.transform(q, crossfit.FitContext(clusters))
    runtime["crossfit_fit_transform_2000_s"] = round(time.time() - t, 1)
    runtime["crossfit_columns"] = list(out.columns)
    runtime["crossfit_nan_frac"] = float(out.isna().mean().mean())
except Exception as e:  # pragma: no cover
    runtime["crossfit_error"] = repr(e)

(RES / "validation_runtime.json").write_text(json.dumps(runtime, indent=1))
print(json.dumps(runtime, indent=1))
print(ca.head(25).round(3).to_string())
print(lv.head(30).round(3).to_string())
num = [c for c in lv.columns if c.startswith(("pssm_", "cantley_p"))]
print(lv.groupby("train_bin", observed=True)[num].mean().round(3).to_string())
print(lv.groupby("cantley_matrix")[num].mean().round(3).to_string(), lv.groupby("cantley_matrix").size())
print(pd.DataFrame(tops).to_string())
print(site_centric.to_string())
print(ca.groupby("matrix")[list(arms)].mean().round(3).to_string(), ca.groupby("matrix").size())
