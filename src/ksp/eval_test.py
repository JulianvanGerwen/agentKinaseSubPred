"""The ONE AND ONLY reader of sealed TEST labels (DESIGN §4.3, AIFM-2, §6, §16-17). Owner: core-infra.

G4 usage (PI): commit FREEZE.json, then
    /workspace/.venv/bin/python -I src/ksp/eval_test.py
Dry run (never touches sealed/, FREEZE.json or results/test/LEDGER.md; toy "TEST" = DEV CV fold 0):
    /workspace/.venv/bin/python -I src/ksp/eval_test.py --dry-run [--fold 0] [--B 1000] [--freeze cfg.json]

Refuses to run (real mode) unless
  * /workspace/FREEZE.json exists, is tracked by git (`git ls-files --error-unmatch`) and unmodified
    (`git diff --quiet HEAD --`) - read-only git calls;
  * FREEZE.json validates against the schema below;
  * every sealed file matches its SHA-256 in results/splits/SEAL.json;
  * this FREEZE.json (by SHA-256) has not been scored before (results/test/LEDGER.md) - single-batch rule.

FREEZE.json schema (version 1)
------------------------------
{
  "version": 1,
  "harness_commit": "<git hash of the frozen harness (A28)>",
  "B": 10000,                                   # paired cluster-bootstrap resamples (A19); default 10000
  "notes": "free text",
  "configs": [                                  # every frozen config AND every baseline scored at G4
    {"name": "final",                           # unique; referenced by comparisons
     "role": "final" | "baseline" | "ablation" | "control",
     "run_id": "<results/runs/<run_id>>",       # training run (manifest with seed, inputs)
     "predictions": "cache/runs/<run_id>/test_predictions.parquet",
         # columns: (kinase_idx, site_idx) or (kinase_acc, site_id), score [, score_site]
         # EXHAUSTIVE over TEST U_ms sites x kinases of matching residue class; a kinase is either scored on all
         # its candidates or not at all (unscored = missing, ranked last, counted).
     "score_site_predictions": null | "<path>", # optional separate file with score_site (A20: M2 score)
     "t3_predictions": null | "<path>",         # optional T3: TEST rows scored by models trained with the
                                                # kinase's LKO family removed (unseen family + unseen protein)
     "runtime": {"train_wall_s": float, "inference_wall_s": float},   # REQUIRED (DESIGN §1)
     "missing": "worst"}                        # optional; 'worst' (default) | 'drop'
  ],
  "comparisons": [
    {"a": "final", "b": "B1_arm_i",  "endpoints": ["M1", "M3"], "primary": true,
     "label": "final vs B1 Cantley (A23 arm i)"},
    {"a": "final", "b": "B1_arm_ii", "endpoints": ["M2"],       "primary": true},
    {"a": "final", "b": "B3_retrained", "endpoints": ["M1", "M2", "M3"], "primary": true},
    {"a": "final", "b": "B4_retrained", "endpoints": ["M1", "M2", "M3"], "primary": true},
    {"a": "final", "b": "B7", "endpoints": ["M1", "M2", "M3", "M4"], "primary": false}
  ]
}
Exactly 9 primary tests (comparison x endpoint) are required (DESIGN §6.8 PRIOR-8); Holm is applied over those 9,
and separately over all secondary tests. Every comparison is computed on the MATCHED universe of (kinase, site)
pairs scored by both configs (PRIOR-1); this also makes M2 matched-kinase (A27: candidate kinases scored by both).
"Beats" = Holm-adjusted p < 0.05 and paired 95% CI of the difference > 0. Resolution: the smallest bootstrap
p-value is 2/(B+1), so Holm over 9 tests can only reach < 0.05 if B >= ~360 (B=10000: floor 2e-4 x 9 = 1.8e-3).

Outputs (results/test/<freeze_sha12>/): <name>.json, <name>_ci.tsv, <name>_strata.tsv, <name>_t3.json,
superiority.tsv (9 primary tests), secondary.tsv, summary.md; one block appended to results/test/LEDGER.md.

Strata (§6.4 + A26), computed by subsetting kinase sets / positive pairs of the frozen ksp.eval outputs (A28):
kinase group, residue class, training-label bins (M4 over K_any), priming-dependent, PIKK, FAM20C (kinase-level);
evidence class, new-since-2021, external curated, hub sites (>= 5 non-paralog labelled kinases) (pair-level);
T3 from t3_predictions. Stratum definitions come from src/ksp/strata.py (canonical; built-in fallback if absent):
kinase-level manning_group, residue_class, train_bin, priming, priming_cantley, pikk, pikk_sq, secretory, fam20c;
pair-level evidence, hub, primed, primed_cantley; plus temporal (ksp.labels.new_since_2021) and external curated.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if not __package__:                       # executed as a script (python src/ksp/eval_test.py): enable `from . import`
    import ksp  # noqa: F401
    __package__ = "ksp"

FREEZE = ROOT / "FREEZE.json"
SEAL = ROOT / "results" / "splits" / "SEAL.json"
SEALED_DIR = ROOT / "sealed"
LEDGER = ROOT / "results" / "test" / "LEDGER.md"
N_PRIMARY = 9
ENDPOINT_COL = {"M1": "auroc_hard", "M2": "hit5", "M3": "top1_recall", "M4": "auroc_hard"}


class FreezeError(RuntimeError):
    pass


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ================================================================================================ guards + schema
def validate_freeze(cfg: dict, root: Path = ROOT, allow_nonstandard_primary: bool = False) -> dict:
    if cfg.get("version", 1) != 1:
        raise FreezeError("unsupported FREEZE.json version")
    cs = cfg.get("configs") or []
    if not cs:
        raise FreezeError("FREEZE.json lists no configs")
    names = [c.get("name") for c in cs]
    if len(set(names)) != len(names) or any(not n for n in names):
        raise FreezeError("config names must be unique and non-empty")
    for c in cs:
        for k in ("predictions", "runtime"):
            if k not in c:
                raise FreezeError(f"config {c['name']}: missing required field {k!r}")
        rt = c["runtime"]
        for k in ("train_wall_s", "inference_wall_s"):
            if not isinstance(rt.get(k), (int, float)):
                raise FreezeError(f"config {c['name']}: runtime.{k} must be a number (DESIGN §1)")
        for k in ("predictions", "score_site_predictions", "t3_predictions"):
            if c.get(k) and not (root / c[k]).exists():
                raise FreezeError(f"config {c['name']}: {k} file {c[k]} does not exist")
    n_primary = 0
    for cp in cfg.get("comparisons", []):
        if cp.get("a") not in names or cp.get("b") not in names:
            raise FreezeError(f"comparison {cp} references unknown configs")
        eps = cp.get("endpoints", ["M1", "M2", "M3"])
        if not set(eps) <= {"M1", "M2", "M3", "M4"}:
            raise FreezeError(f"comparison {cp}: unknown endpoints")
        if cp.get("primary"):
            n_primary += len(eps)
    if n_primary != N_PRIMARY and not allow_nonstandard_primary:
        raise FreezeError(f"{n_primary} primary tests declared; DESIGN §6.8 pre-declares exactly {N_PRIMARY}")
    return cfg


def check_freeze(freeze: Path = FREEZE, require_git: bool = True, allow_nonstandard_primary: bool = False) -> dict:
    if not freeze.exists():
        raise FreezeError(f"{freeze} does not exist: TEST evaluation is only allowed after the freeze (G4)")
    if require_git:
        r = subprocess.run(["git", "ls-files", "--error-unmatch", freeze.name], cwd=freeze.parent,
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise FreezeError("FREEZE.json is not committed (git ls-files does not list it)")
        r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", freeze.name], cwd=freeze.parent)
        if r.returncode != 0:
            raise FreezeError("FREEZE.json has uncommitted changes")
    return validate_freeze(json.loads(freeze.read_text()), allow_nonstandard_primary=allow_nonstandard_primary)


def check_seal(seal: Path = SEAL, sealed_dir: Path = SEALED_DIR) -> dict:
    s = json.loads(Path(seal).read_text())
    for name, meta in s["files"].items():
        f = Path(sealed_dir) / name
        if not f.exists() or _sha(f) != meta["sha256"]:
            raise FreezeError(f"sealed file {name} missing or modified (SHA-256 mismatch with SEAL.json)")
    return s


def check_single_batch(freeze_sha: str, ledger: Path = LEDGER) -> None:
    if Path(ledger).exists() and freeze_sha in Path(ledger).read_text():
        raise FreezeError(f"FREEZE.json {freeze_sha[:12]} was already scored (single-batch rule)")


# ================================================================================================ predictions
def load_predictions(c: dict, key: str = "predictions", root: Path = ROOT) -> pd.DataFrame:
    p = pd.read_parquet(root / c[key])
    cols = [x for x in ["kinase_idx", "site_idx", "kinase_acc", "site_id", "score", "score_site"] if x in p]
    p = p[cols]
    if key == "predictions" and c.get("score_site_predictions"):
        ss = pd.read_parquet(root / c["score_site_predictions"])
        on = ["kinase_idx", "site_idx"] if {"kinase_idx", "site_idx"} <= set(ss) and {"kinase_idx", "site_idx"} \
            <= set(p) else ["kinase_acc", "site_id"]
        p = p.merge(ss[on + ["score_site"]], on=on, how="left")
        if p.score_site.isna().any():
            raise FreezeError(f"config {c['name']}: score_site missing for {p.score_site.isna().sum()} rows")
    return p


def _keys(p: pd.DataFrame, data) -> np.ndarray:
    from . import eval as E
    m = E.map_predictions(p, data)
    return m.k.values * np.int64(len(data.sites)) + m.s.values


def matched_pair(pa: pd.DataFrame, pb: pd.DataFrame, data) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Matched universe (PRIOR-1) on canonical (kinase, site) keys; works for idx- or id-based files."""
    ka, kb = _keys(pa, data), _keys(pb, data)
    common = np.intersect1d(ka, kb)
    a = pa[np.isin(ka, common)].reset_index(drop=True)
    b = pb[np.isin(kb, common)].reset_index(drop=True)
    return a, b, {"n_A": len(pa), "n_B": len(pb), "n_common": int(len(common))}


# ================================================================================================ owner pre-flight
def test_candidate_pairs() -> pd.DataFrame:
    """Public TEST candidate set (no labels): (kinase_idx, site_idx) for every real kinase x every U_ms site on a
    TEST protein (results/splits/test_proteins.tsv) of matching residue class. ~36M rows."""
    from . import eval as E
    from . import index, splits
    S = index.load_sites()
    K = index.load_kinases()
    uni = np.where(S.in_U_ms.values & S.acc.isin(splits.load_test_proteins()).values)[0]
    res = S.res.values[uni]
    parts = [pd.DataFrame({"kinase_idx": np.int32(i), "site_idx": uni[E._res_ok(res, c)].astype(np.int32)})
             for i, c in enumerate(K.residue_class)]
    return pd.concat(parts, ignore_index=True)


def check_predictions(path, kind: str = "test") -> dict:
    """Owner pre-flight for a G4 predictions file (never reads labels): columns, finite scores, no duplicates,
    only TEST U_ms candidates, and the all-or-nothing rule (each kinase scores all its candidates or none).
    kind='t3' relaxes nothing but reports per-kinase coverage the same way."""
    from . import index
    p = pd.read_parquet(path)
    if not {"score"} <= set(p) or not ({"kinase_idx", "site_idx"} <= set(p) or {"kinase_acc", "site_id"} <= set(p)):
        raise ValueError("need (kinase_idx, site_idx) or (kinase_acc, site_id) plus score")
    if "kinase_idx" not in p:
        K, S = index.load_kinases(), index.load_sites()
        p["kinase_idx"] = p.kinase_acc.map({a: i for i, a in enumerate(K.uniprot_acc)})
        p["site_idx"] = p.site_id.map(pd.Series(np.arange(len(S)), index=S.site_id))
        if p[["kinase_idx", "site_idx"]].isna().any().any():
            raise ValueError("unknown kinase_acc/site_id (canaries must not be in TEST files)")
    if not np.isfinite(p.score.values).all():
        raise ValueError("non-finite scores")
    c = test_candidate_pairs()
    n = len(index.load_sites())
    kc = c.kinase_idx.values.astype(np.int64) * n + c.site_idx.values
    kp = p.kinase_idx.values.astype(np.int64) * n + p.site_idx.values
    if len(np.unique(kp)) != len(kp):
        raise ValueError("duplicate (kinase, site) rows")
    extra = ~np.isin(kp, kc)
    need = c.groupby("kinase_idx").size()
    have = p[~extra].groupby("kinase_idx").size().reindex(need.index).fillna(0)
    partial = need.index[(have > 0) & (have < need)].tolist()
    rep = {"rows": int(len(p)), "rows_outside_candidates": int(extra.sum()), "kinases_full": int((have == need).sum()),
           "kinases_unscored": int((have == 0).sum()), "kinases_partial": len(partial),
           "partial_kinase_idx": partial[:20], "ok": bool(extra.sum() == 0 and not partial)}
    return rep


# ================================================================================================ strata
def _strata_module():
    """Canonical stratum definitions from src/ksp/strata.py (strata analyst); None -> built-in fallback below."""
    try:
        from . import strata as S
    except ImportError:
        return None
    return S


def kinase_strata_table(data) -> pd.DataFrame:
    """Kinase-level stratum memberships: kinase_acc, stratum, value (ksp.strata.kinase_sets when available:
    manning_group, residue_class, train_bin, priming, priming_cantley, pikk, pikk_sq, secretory, fam20c)."""
    from . import eval as E
    K = data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY_")]
    S = _strata_module()
    if S is not None:
        T = S.kinase_table(K, train_counts=data.train_counts)
        rows = [(a, st, v) for (st, v), ks in S.kinase_sets(T).items() for a in ks]
        return pd.DataFrame(rows, columns=["kinase_acc", "stratum", "value"])
    tc = K.uniprot_acc.map(data.train_counts).fillna(0)
    rows = []
    for acc, g, cls, fam, gene, t in zip(K.uniprot_acc, K.manning_group.fillna("(none)"), K.residue_class,
                                         K.manning_family, K.gene, tc):
        rows += [(acc, "manning_group", g), (acc, "residue_class", cls),
                 (acc, "train_bin", str(pd.cut([t], E.TRAIN_BINS, labels=E.TRAIN_BIN_LABELS)[0])),
                 (acc, "priming", str(bool(gene in E.PRIMING or fam == "CK1"))),
                 (acc, "pikk", str(bool(fam == "PIKK"))), (acc, "fam20c", str(bool(gene == "FAM20C")))]
    return pd.DataFrame(rows, columns=["kinase_acc", "stratum", "value"])


def _hub_builtin(pm: pd.DataFrame, par: dict) -> dict:
    hub = {}
    for s, g in pm.groupby("site_idx"):
        ks = set(g.kinase_acc)
        seen, n = set(), 0
        for k in ks:
            if k in seen:
                continue
            n += 1
            stack = [k]
            while stack:
                x = stack.pop()
                if x not in seen:
                    seen.add(x)
                    stack += [y for y in par.get(x, ()) if y in ks and y not in seen]
        hub[s] = n >= 5
    return hub


def pair_strata_table(data) -> pd.DataFrame:
    """Pair-level stratum memberships over the split's evaluation positives: kinase_acc, site_idx, stratum, value.
    evidence / hub / primed / primed_cantley from ksp.strata when available; temporal from ksp.labels
    (new_since_2021 flag written at split freeze)."""
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].copy()
    pm["site_id"] = data.sites.site_id.values[pm.site_idx.values]
    out = []
    if "new_since_2021" in pm:
        out.append(pm.assign(stratum="temporal", value=np.where(pm.new_since_2021.astype(bool), "new_since_2021",
                                                                "pre_2021")))
    S = _strata_module()
    if S is not None:
        K = data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY_")]
        if "evidence_native" in pm:
            out.append(pm.assign(stratum="evidence", value=S.evidence_class(pm).values))
        hubs = S.hub_sites(pm[["kinase_acc", "site_id"]], data.par)
        out.append(pm.assign(stratum="hub", value=np.where(pm.site_id.isin(hubs), "hub", "non_hub")))
        for nm, offs in [("primed", S.priming_curated(K)), ("primed_cantley", S.priming_cantley(K))]:
            m = S.primed_mask(pm[["kinase_acc", "site_id"]], offs, sites=data.sites)
            sub = pm[pm.kinase_acc.isin(set(offs))]                     # only kinases with priming offsets
            out.append(sub.assign(stratum=nm, value=np.where(m[pm.kinase_acc.isin(set(offs)).values], "primed",
                                                             "unprimed")))
    else:
        if "evidence_native" in pm:
            out.append(pm.assign(stratum="evidence", value=pm.evidence_native.astype(str)))
        hub = _hub_builtin(pm, data.par)
        out.append(pm.assign(stratum="hub", value=pm.site_idx.map(hub).map({True: "hub", False: "non_hub"})))
    return pd.concat(out, ignore_index=True)[["kinase_acc", "site_idx", "stratum", "value"]]


def _macro(x) -> float:
    x = np.asarray(list(x), float)
    x = x[~np.isnan(x)]
    return float(x.mean()) if len(x) else float("nan")


def strata_results(res, data, name: str) -> pd.DataFrame:
    """All endpoints on every stratum (point estimates). Kinase strata: macro over the stratum's kinases within
    K_eval (M1-M3) / K_any (M4). Pair strata: per kinase, the stratum's positives vs that kinase's own E-hard
    negatives (M1/M4), its site-centric hits (M2) and top-1% hits (M3); macro over kinases with >= 1 such pair.
    External curated: separate positive set (data.ext_positives) vs the same negatives."""
    from . import eval as E
    pk = res.per_kinase.set_index("kinase_acc")
    kp, ka = set(data.k_primary), set(data.k_any)
    rows = []
    for (st, val), g in kinase_strata_table(data).groupby(["stratum", "value"]):
        ks = set(g.kinase_acc) & set(pk.index)
        p_, a_ = sorted(ks & kp), sorted(ks & ka)
        rows.append({"config": name, "level": "kinase", "stratum": st, "value": val, "n_kinases_K_eval": len(p_),
                     "n_kinases_K_any": len(a_),
                     "M1": _macro(pk.loc[p_, "auroc_hard"]) if p_ else np.nan,
                     "M2": _macro(pk.loc[p_, "hit5"]) if p_ and "hit5" in pk else np.nan,
                     "M3": _macro(pk.loc[p_, "top1_recall"]) if p_ else np.nan,
                     "M4": _macro(pk.loc[a_, "auroc_hard"]) if a_ else np.nan})
    ps = pair_strata_table(data)
    sl = res.site_level.set_index(["kinase_acc", "site_idx"]) if len(res.site_level) else None
    for (st, val), g in ps.groupby(["stratum", "value"]):
        per = {}
        for k, gg in g.groupby("kinase_acc"):
            it = res.items.get(k)
            if it is None:
                continue
            pos = data.positives[k]
            if len(it["pos_s"]) != len(pos):
                raise RuntimeError("pair strata need missing='worst' results (positives aligned to data.positives)")
            msk = np.isin(pos, gg.site_idx.values)
            h5 = np.nan
            if sl is not None:
                idx = [(k, s) for s in gg.site_idx.values if (k, s) in sl.index]
                h5 = float(sl.loc[idx, "hit5"].mean()) if idx else np.nan
            per[k] = (E.auroc(it["pos_s"][msk], it["hard_s"]), h5, float(it["top_hits"][msk].mean()), int(msk.sum()))
        p_, a_ = [k for k in per if k in kp], [k for k in per if k in ka]
        rows.append({"config": name, "level": "pair", "stratum": st, "value": val, "n_kinases_K_eval": len(p_),
                     "n_kinases_K_any": len(a_), "n_pairs": int(sum(per[k][3] for k in a_)),
                     "M1": _macro(per[k][0] for k in p_), "M2": _macro(per[k][1] for k in p_),
                     "M3": _macro(per[k][2] for k in p_), "M4": _macro(per[k][0] for k in a_)})
    # external curated (separate positive set)
    ext = [r for r in res.pair_strata.to_dict(orient="records") if r["stratum"] == "external_curated"] \
        if len(res.pair_strata) else []
    for r in ext:
        rows.append({"config": name, "level": "pair", "stratum": "external_curated", "value": r["kinase_set"],
                     "n_kinases_K_eval": r["n_kinases"], "n_kinases_K_any": r["n_kinases"], "n_pairs": r["n_pos"],
                     "M1": r["auroc_hard"], "M2": np.nan, "M3": r["top1_recall"], "M4": np.nan})
    return pd.DataFrame(rows)


# ================================================================================================ superiority
def run_comparisons(cfg: dict, preds: dict, data, B: int, kinases: list | None = None) -> pd.DataFrame:
    """Paired cluster-bootstrap tests on matched universes; Holm over primary tests and separately secondary."""
    from . import eval as E
    rows = []
    for cp in cfg.get("comparisons", []):
        a, b = cp["a"], cp["b"]
        pa, pb, rep = matched_pair(preds[a], preds[b], data)
        ra = E.evaluate(pa, data, missing="drop", kinases=kinases)
        rb = E.evaluate(pb, data, missing="drop", kinases=kinases)
        d = E.paired_bootstrap(ra, rb, B=B, seed=0).set_index("endpoint")
        for ep in cp.get("endpoints", ["M1", "M2", "M3"]):
            r = d.loc[ep]
            rows.append({"comparison": cp.get("label", f"{a} vs {b}"), "a": a, "b": b, "endpoint": ep,
                         "primary": bool(cp.get("primary", False)), "diff": r["diff"], "lo95": r.lo95,
                         "hi95": r.hi95, "p_boot": r.p_boot, "B": B, "n_common_pairs": rep["n_common"],
                         "A_matched": ra.endpoints[_ekey(ep)], "B_matched": rb.endpoints[_ekey(ep)]})
    C = pd.DataFrame(rows)
    if not len(C):
        return C
    C["p_holm"] = np.nan
    for flag in (True, False):
        m = C.primary == flag
        if m.any():
            C.loc[m, "p_holm"] = E.holm(C.loc[m, "p_boot"].values)
    C["beats"] = (C.p_holm < 0.05) & (C.lo95 > 0)
    return C


def _ekey(ep: str) -> str:
    return {"M1": "M1_auroc_hard", "M2": "M2_hit5", "M3": "M3_top1_recall", "M4": "M4_auroc_hard_Kany"}[ep]


# ================================================================================================ batch
def run_batch(cfg: dict, data, outd: Path, B: int, root: Path = ROOT, B_strata: int = 1000,
              kinases: list | None = None) -> dict:
    """Score every config in one batch on `data` (an ksp.eval.EvalData). Returns summary dict; writes outputs."""
    from . import eval as E
    outd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    preds, summ = {}, {}
    for c in cfg["configs"]:
        p = load_predictions(c, root=root)
        preds[c["name"]] = p
        res = E.evaluate(p, data, missing=c.get("missing", "worst"), kinases=kinases)
        ci = E.bootstrap(res, B=B, seed=0)
        res.save_json(outd / f"{c['name']}.json")
        ci.to_csv(outd / f"{c['name']}_ci.tsv", sep="\t", index=False)
        st = strata_results(res, data, c["name"])
        st.to_csv(outd / f"{c['name']}_strata.tsv", sep="\t", index=False)
        e = res.endpoints
        summ[c["name"]] = {k: e.get(_ekey(k)) for k in ("M1", "M2", "M3", "M4")}
        summ[c["name"]].update({"n_missing_pairs": e["n_missing_pairs"], "n_required_pairs": e["n_required_pairs"],
                                "role": c.get("role", ""), "runtime": c["runtime"],
                                "ci": ci.set_index("endpoint")[["lo95", "hi95"]].to_dict(orient="index")})
        if c.get("t3_predictions"):
            p3 = load_predictions(c, key="t3_predictions", root=root)
            r3 = E.evaluate(p3, data, missing=c.get("missing", "worst"), kinases=kinases)
            r3.save_json(outd / f"{c['name']}_t3.json")
            summ[c["name"]]["T3"] = {k: r3.endpoints.get(_ekey(k)) for k in ("M1", "M2", "M3", "M4")}
    C = run_comparisons(cfg, preds, data, B, kinases=kinases)
    if len(C):
        C[C.primary].to_csv(outd / "superiority.tsv", sep="\t", index=False)
        C[~C.primary].to_csv(outd / "secondary.tsv", sep="\t", index=False)
    summ["_comparisons"] = C.to_dict(orient="records")
    summ["_runtime_s"] = round(time.time() - t0, 1)
    _write_summary(summ, outd)
    return summ


def _write_summary(summ: dict, outd: Path) -> None:
    f = lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.4f}"  # noqa: E731
    L = ["| config | role | M1 | M2 Hit@5 | M3 top1% | M4 | T3 M1 | missing/required | train s | infer s |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for n, s in summ.items():
        if n.startswith("_"):
            continue
        L.append(f"| {n} | {s['role']} | {f(s['M1'])} | {f(s['M2'])} | {f(s['M3'])} | {f(s['M4'])} | "
                 f"{f(s.get('T3', {}).get('M1'))} | {s['n_missing_pairs']}/{s['n_required_pairs']} | "
                 f"{s['runtime']['train_wall_s']} | {s['runtime']['inference_wall_s']} |")
    C = pd.DataFrame(summ.get("_comparisons", []))
    if len(C):
        for flag, title in [(True, "Primary superiority tests (Holm over these)"),
                            (False, "Secondary comparisons (Holm separately)")]:
            x = C[C.primary == flag]
            if not len(x):
                continue
            L += ["", f"**{title}**", "", "| comparison | endpoint | diff | 95% CI | p_boot | p_holm | beats |",
                  "|---|---|---|---|---|---|---|"]
            for r in x.itertuples():
                L.append(f"| {r.comparison} | {r.endpoint} | {f(r.diff)} | [{f(r.lo95)}, {f(r.hi95)}] | "
                         f"{f(r.p_boot)} | {f(r.p_holm)} | {bool(r.beats)} |")
    (outd / "summary.md").write_text("\n".join(L) + "\n")
    (outd / "summary.json").write_text(json.dumps(summ, indent=1, default=_jd))


def _jd(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def _append_ledger(ledger: Path, fsha: str, cfg: dict, summ: dict, outd: Path, tag: str = "TEST") -> None:
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger, "a") as fh:
        fh.write(f"\n## {tag} batch {time.strftime('%Y-%m-%d %H:%M:%S %Z')} - FREEZE.json sha256 {fsha}\n\n")
        fh.write(f"Configs: {len(cfg['configs'])}; B={cfg.get('B', 10000)}; runtime {summ['_runtime_s']}s; "
                 f"harness {cfg.get('harness_commit', '?')}; details {outd}\n\n")
        fh.write((outd / "summary.md").read_text())


# ================================================================================================ entry points
def main(freeze: Path = FREEZE, seal: Path = SEAL, sealed_dir: Path = SEALED_DIR, ledger: Path = LEDGER,
         out_root: Path | None = None, require_git: bool = True, data_builder=None, kinases: list | None = None,
         allow_nonstandard_primary: bool = False) -> dict:
    """Real G4 path (all guards). Tests call it with tmp paths and a toy data_builder (dry run)."""
    from . import eval as E
    cfg = check_freeze(freeze, require_git=require_git, allow_nonstandard_primary=allow_nonstandard_primary)
    check_seal(seal, sealed_dir)
    fsha = _sha(freeze)
    check_single_batch(fsha, ledger)
    test_labels = pd.read_parquet(Path(sealed_dir) / "test_labels.parquet")
    data = data_builder(test_labels) if data_builder else E.prepare("test", labels_df=test_labels,
                                                                     with_canaries=False)
    outd = (out_root or Path(ledger).parent) / fsha[:12]
    summ = run_batch(cfg, data, outd, int(cfg.get("B", 10000)), kinases=kinases)
    _append_ledger(Path(ledger), fsha, cfg, summ, outd)
    print(f"scored {len(cfg['configs'])} configs -> {ledger}")
    return summ


def toy_data_builder(fold: int = 0):
    """Dry-run data: DEV CV fold `fold` as a fake TEST (labels from load_dev; never sealed/). K_eval / K_any are
    recomputed on the toy split with the frozen §4.6 rules (>= 3 positives on >= 2 clusters; >= 1)."""
    from . import eval as E

    def build(_labels_ignored=None):
        d = E.prepare("cv", fold=fold, with_canaries=False)
        ke = [k for k, p in d.positives.items() if len(p) >= 3 and len(np.unique(d.site_cluster[p])) >= 2]
        return dataclasses.replace(d, split="test", k_primary=sorted(ke), k_any=sorted(d.positives))
    return build


def write_toy_sealed(sealed_dir: Path, seal: Path, fold: int = 0) -> None:
    """Toy sealed dir + SEAL.json for the dry run: DEV labels of CV fold `fold` (public DEV data)."""
    from . import labels, splits
    sealed_dir.mkdir(parents=True, exist_ok=True)
    d = labels.load_dev("all")
    folds = splits.load_cv_folds()
    d[d.acc.map(folds).eq(fold)].to_parquet(sealed_dir / "test_labels.parquet", index=False)
    seal.write_text(json.dumps({"files": {"test_labels.parquet": {"sha256": _sha(sealed_dir / "test_labels.parquet")}}}))


def dry_run(freeze_cfg: Path | None = None, fold: int = 0, B: int = 1000,
            work: Path = ROOT / "cache" / "eval_test_dryrun") -> dict:
    """End-to-end batch path on DEV fold `fold`; all paths under `work`. With no freeze_cfg, uses
    dry_run_freeze() (existing fold-0 OOF runs as configs)."""
    work.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(freeze_cfg).read_text()) if freeze_cfg else dry_run_freeze(fold, work)
    cfg["B"] = B
    fz = work / "FREEZE.json"
    fz.write_text(json.dumps(cfg, indent=1))
    write_toy_sealed(work / "sealed", work / "SEAL.json", fold)
    led = work / "LEDGER.md"
    if led.exists():
        led.unlink()
    return main(fz, work / "SEAL.json", work / "sealed", led, out_root=work, require_git=False,
                data_builder=toy_data_builder(fold), allow_nonstandard_primary=True)


def dry_run_freeze(fold: int, work: Path) -> dict:
    """FREEZE-like config from existing runs, restricted to CV fold `fold` sites: R2 fold-0 run as 'final',
    B1 arms i/ii, B7, B0-prior (pooled CV runs filtered to the fold), R3 final LKO pooled as T3 of 'final'."""
    from . import index, splits
    runs = ROOT / "cache" / "runs"
    S = index.load_sites()
    sf = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values

    def pick(prefix):
        c = sorted(p.name for p in runs.glob(prefix + "__*") if (p / "predictions.parquet").exists())
        if not c:
            raise FileNotFoundError(f"no run {prefix}")
        return c[-1]

    def sub(run, name):
        p = pd.read_parquet(runs / run / "predictions.parquet")
        if "site_idx" not in p:
            p["site_idx"] = pd.Series(np.arange(len(S)), index=S.site_id)[p.site_id].values
        p = p[sf[p.site_idx.values] == fold]
        if "kinase_idx" in p:                      # TEST files carry no canaries (real kinases only)
            p = p[p.kinase_idx < len(index.load_kinases())]
        if "kinase_acc" in p:
            p = p[~p.kinase_acc.astype(str).str.startswith("CANARY_")]
        out = work / f"{name}.parquet"
        p.to_parquet(out, index=False)
        return str(out.relative_to(ROOT))
    final = pick(f"r2_final_s0_f{fold}")
    cfgs = [{"name": "final_R2", "role": "final", "run_id": final, "predictions": sub(final, "final_R2"),
             "t3_predictions": sub(pick("r3_final_distill_lko_pooled"), "final_T3"),
             "runtime": {"train_wall_s": 0.0, "inference_wall_s": 0.0}}]
    for nm, pre, role in [("B1_arm_i", "simple_b1_published_reeval", "baseline"),
                          ("B1_arm_ii", "simple_b1_universe_reeval", "baseline"),
                          ("B7", "simple_b7_lr_reeval", "baseline"),
                          ("B0_prior", "simple_b0_prior_reeval", "control")]:
        r = pick(pre)
        cfgs.append({"name": nm, "role": role, "run_id": r, "predictions": sub(r, nm),
                     "runtime": {"train_wall_s": 0.0, "inference_wall_s": 0.0}})
    return {"version": 1, "harness_commit": "dry-run", "notes": f"dry run on DEV fold {fold}", "configs": cfgs,
            "comparisons": [{"a": "final_R2", "b": "B1_arm_i", "endpoints": ["M1", "M3"], "primary": True},
                            {"a": "final_R2", "b": "B1_arm_ii", "endpoints": ["M2"], "primary": True},
                            {"a": "final_R2", "b": "B7", "endpoints": ["M1", "M2", "M3"], "primary": True},
                            {"a": "final_R2", "b": "B0_prior", "endpoints": ["M1", "M2", "M3"], "primary": False}]}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--freeze", default=None)
    ap.add_argument("--check-predictions", default=None, help="owner pre-flight for a TEST predictions file")
    a = ap.parse_args()
    if a.check_predictions:
        print(json.dumps(check_predictions(a.check_predictions), indent=1))
    elif a.dry_run:
        s = dry_run(Path(a.freeze) if a.freeze else None, a.fold, a.B)
        print(json.dumps({k: v for k, v in s.items() if not k.startswith("_")}, indent=1, default=_jd))
    else:
        main()
