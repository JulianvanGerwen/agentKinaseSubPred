"""The ONE AND ONLY reader of sealed TEST labels (DESIGN §4.3, AIFM-2, §6, §16-17). Owner: core-infra.

G4 (PI): commit FREEZE.json, then
    /workspace/.venv/bin/python -I src/ksp/eval_test.py --real [--override-completed "<reason>"]
Rehearsal (never touches sealed/, FREEZE.json or results/test/LEDGER.md; toy "TEST" = DEV CV fold 0):
    /workspace/.venv/bin/python -I src/ksp/eval_test.py --dry-run [--fold 0] [--B 1000] [--freeze cfg.json]
Owner pre-flight (no labels):
    /workspace/.venv/bin/python -I src/ksp/eval_test.py --check-predictions <file>
Without --real / --dry-run / --check-predictions the script prints usage and exits. --B/--fold only with --dry-run.

Order of operations in the real batch (main):
  1. FREEZE.json git-tracked + unmodified (read-only `git ls-files` / `git diff`), schema valid (validate_freeze),
     every input file's SHA-256 equal to the value pinned in FREEZE.json (parquets are gitignored: this is the pin).
  2. Real-mode policy: B == 10000; harness_commit resolves to the A28 harness commit 974d520 and
     `git diff --quiet <harness_commit> -- src/ksp/{eval,labels,splits,crossfit,canary}.py`.
  3. Pre-flight of every predictions file (mapped keys in range, no canaries, finite score/score_site, no
     duplicates, only TEST candidates, all-or-nothing per kinase; T3 covers every K_any kinase).
  4. Seal check (sealed files vs results/splits/SEAL.json).
  5. Ledger state machine (results/test/LEDGER.md, machine-readable `<!-- KSP-LEDGER {json} -->` markers):
     - any COMPLETED TEST block -> refuse (single batch), unless --override-completed "<reason>" (logged);
     - an unfinished STARTED block -> rerun only with identical FREEZE sha, input hashes and code hash (logged as
       attempt n); anything else -> refuse.
     A STARTED block (FREEZE sha, input hashes, git HEAD, SHA-256 of every src/ksp/*.py) is appended BEFORE the
     sealed labels are read; a COMPLETED block at the end.
  6. Outputs are written to results/test/<sha12>.partial/ and atomically renamed to results/test/<sha12>/.

FREEZE.json schema (version 1)
------------------------------
{
  "version": 1,
  "harness_commit": "974d520",                  # A28 harness commit (real mode: must resolve to it)
  "B": 10000,                                   # primary paired cluster bootstrap (A19); real mode: exactly 10000
  "B_secondary": 1000,                          # optional; CIs of standalone/strata/sensitivity analyses
  "notes": "free text",
  "configs": [
    {"name": "final", "role": "final",          # exactly ONE config with role 'final'
     "run_id": "<results/runs/<run_id>>",
     "predictions": "cache/runs/<run_id>/test_predictions.parquet", "predictions_sha256": "<64 hex>",
         # (kinase_idx, site_idx) or (kinase_acc, site_id), score [, score_site]; EXHAUSTIVE over TEST candidates
         # (U_ms sites of TEST proteins x kinases of matching residue class), all-or-nothing per kinase.
     "score_site_predictions": null | "<path>", "score_site_predictions_sha256": "<hex>",
         # optional separate score_site file (same key columns as `predictions`). score_site = the cross-kinase
         # score used for M2 (A20; any finite value, e.g. a raw fused score chosen on CV).
     "assignment_prob_predictions": null | "<path>", "assignment_prob_predictions_sha256": "<hex>",
         # optional (inline column `assignment_prob` or separate file): A17 site-level assignment probability
         # P(k | s) in [0, 1]; used ONLY for the A17 reliability analysis and Coverage@tau (not for any endpoint).
     "t3_predictions": null | "<path>", "t3_predictions_sha256": "<hex>",
         # T3 (unseen family + unseen protein): each kinase scored by the model trained with its LKO family
         # (results/splits/kinase_folds.tsv) removed; must cover every K_any kinase fully.
     "runtime": {"train_wall_s": float, "inference_wall_s": float},      # REQUIRED (DESIGN §1)
     "missing": "worst"},                       # optional; only 'worst' is accepted
    {"name": "B1_i", "role": "baseline", "baseline_id": "B1_arm_i", ...},
    {"name": "B1_ii", "role": "baseline", "baseline_id": "B1_arm_ii", ...},
    {"name": "B3", "role": "baseline", "baseline_id": "B3_retrained", ...},
    {"name": "B4", "role": "baseline", "baseline_id": "B4_retrained", ...},
    {"name": "B7", "role": "baseline", "baseline_id": "B7", ...}        # other roles: ablation | control
  ],
  "comparisons": [
    {"a": "final", "b": "B1_i",  "endpoints": ["M1", "M3"], "primary": true},
    {"a": "final", "b": "B1_ii", "endpoints": ["M2"], "primary": true},
    {"a": "final", "b": "B3",    "endpoints": ["M1", "M2", "M3"], "primary": true},
    {"a": "final", "b": "B4",    "endpoints": ["M1", "M2", "M3"], "primary": true},
    {"a": "final", "b": "B7",    "endpoints": ["M1", "M2", "M3", "M4"], "primary": false}
  ],
  "unavailable_primary": [{"baseline_id": "B4_retrained", "reason": "..."}]    # optional contingency (PI)
}
Primary family (DESIGN §6.8, A23, PRIOR-8): exactly the 9 (baseline_id, endpoint) pairs of PRIMARY_SPEC; `a` is the
final config in every primary comparison; no duplicates, no M4. A baseline listed in unavailable_primary keeps its
tests in the family (m = 9) with p = 1, reported as "NOT TESTED (baseline unavailable)". Holm over the 9;
secondary comparisons Holm separately. Every comparison runs on the MATCHED universe of (kinase, site) pairs scored by
both configs (PRIOR-1), which also makes M2 matched-kinase (A27). "Beats" = Holm p < 0.05 AND the (unadjusted)
paired 95% percentile CI lower bound > 0. Resolution: smallest bootstrap p = 2/(B+1).

Outputs (results/test/<sha12>/): per config <name>.json, _ci.tsv (clusters, B_secondary), _ci_kinases.tsv
(kinase-resampling sensitivity, §6.8), _per_kinase.tsv, _site_level.tsv, _strata.tsv, _ext_tiers.tsv,
_reliability.tsv + Coverage@tau (A17, configs with assignment_prob), _t3.json / _t3_per_kinase.tsv; superiority.tsv, secondary.tsv,
summary.md, summary.json. LEDGER.md: STARTED + COMPLETED blocks.
Strata (§6.4 + A26) from src/ksp/strata.py (canonical; built-in fallback), plus temporal (ksp.labels) and external
curated (with ext_tier A/B), computed by subsetting the frozen ksp.eval outputs (A28).
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import re
import shutil
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
A28_COMMIT = "974d520"
FROZEN_FILES = ["src/ksp/eval.py", "src/ksp/labels.py", "src/ksp/splits.py", "src/ksp/crossfit.py",
                "src/ksp/canary.py"]
# code identity for reruns of an unfinished batch: everything the batch imports
IDENTITY_FILES = ["eval.py", "labels.py", "splits.py", "crossfit.py", "canary.py", "strata.py", "eval_test.py",
                  "index.py", "kinases.py", "ids.py", "uniprot.py"]
REAL_B = 10000
N_PRIMARY = 9
PRIMARY_SPEC = {"B1_arm_i": ("M1", "M3"), "B1_arm_ii": ("M2",), "B3_retrained": ("M1", "M2", "M3"),
                "B4_retrained": ("M1", "M2", "M3")}
FILE_KEYS = ("predictions", "score_site_predictions", "assignment_prob_predictions", "t3_predictions")
MARK = re.compile(r"<!-- KSP-LEDGER (\{.*\}) -->")


class FreezeError(RuntimeError):
    pass


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _git(*args, cwd=ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)


# ================================================================================================ schema
def primary_tests(cfg: dict) -> list[dict]:
    """Expanded primary family: one dict per (baseline_id, endpoint) incl. unavailable ones."""
    by = {c["name"]: c for c in cfg["configs"]}
    out = []
    for cp in cfg.get("comparisons", []):
        if cp.get("primary"):
            for ep in cp.get("endpoints", ["M1", "M2", "M3"]):
                out.append({"a": cp["a"], "b": cp["b"], "baseline_id": by[cp["b"]].get("baseline_id"),
                            "endpoint": ep, "available": True})
    for u in cfg.get("unavailable_primary", []):
        for ep in PRIMARY_SPEC.get(u.get("baseline_id"), ()):
            out.append({"a": None, "b": None, "baseline_id": u["baseline_id"], "endpoint": ep, "available": False,
                        "reason": u.get("reason", "")})
    return out


def validate_freeze(cfg: dict, root: Path = ROOT, allow_nonstandard_primary: bool = False,
                    check_files: bool = True) -> dict:
    """Schema + file pins + primary-family structure. allow_nonstandard_primary: tests/dry runs only."""
    if cfg.get("version", 1) != 1:
        raise FreezeError("unsupported FREEZE.json version")
    cs = cfg.get("configs") or []
    if not cs:
        raise FreezeError("FREEZE.json lists no configs")
    names = [c.get("name") for c in cs]
    if len(set(names)) != len(names) or any(not n for n in names):
        raise FreezeError("config names must be unique and non-empty")
    finals = [c["name"] for c in cs if c.get("role") == "final"]
    if len(finals) != 1:
        raise FreezeError(f"exactly one config must have role 'final' (found {finals})")
    for c in cs:
        if "predictions" not in c or "runtime" not in c:
            raise FreezeError(f"config {c['name']}: missing required field 'predictions'/'runtime'")
        for k in ("train_wall_s", "inference_wall_s"):
            if not isinstance(c["runtime"].get(k), (int, float)):
                raise FreezeError(f"config {c['name']}: runtime.{k} must be a number (DESIGN §1)")
        if c.get("missing", "worst") != "worst":
            raise FreezeError(f"config {c['name']}: missing must be 'worst' (got {c.get('missing')!r})")
        for k in FILE_KEYS:
            if not c.get(k):
                continue
            f = root / c[k]
            if not f.exists():
                raise FreezeError(f"config {c['name']}: {k} file {c[k]} does not exist")
            want = c.get(f"{k}_sha256")
            if not (isinstance(want, str) and re.fullmatch(r"[0-9a-f]{64}", want)):
                raise FreezeError(f"config {c['name']}: {k}_sha256 (64 hex) is required")
            if check_files and _sha(f) != want:
                raise FreezeError(f"config {c['name']}: {k} SHA-256 differs from the value pinned in FREEZE.json")
    for cp in cfg.get("comparisons", []):
        if cp.get("a") not in names or cp.get("b") not in names:
            raise FreezeError(f"comparison {cp} references unknown configs")
        if not set(cp.get("endpoints", ["M1", "M2", "M3"])) <= {"M1", "M2", "M3", "M4"}:
            raise FreezeError(f"comparison {cp}: unknown endpoints")
    if allow_nonstandard_primary:
        return cfg
    by = {c["name"]: c for c in cs}
    pt = primary_tests(cfg)
    seen = set()
    for t in pt:
        if t["available"]:
            if t["a"] != finals[0]:
                raise FreezeError(f"primary comparison must have a = the final config ({finals[0]}), got {t['a']}")
            if t["baseline_id"] not in PRIMARY_SPEC:
                raise FreezeError(f"primary comparison b={t['b']} needs baseline_id in {sorted(PRIMARY_SPEC)}")
            if t["endpoint"] not in PRIMARY_SPEC[t["baseline_id"]]:
                raise FreezeError(f"{t['baseline_id']} is not pre-declared for {t['endpoint']} (A23: "
                                  f"{PRIMARY_SPEC[t['baseline_id']]})")
        elif t["baseline_id"] not in PRIMARY_SPEC:
            raise FreezeError(f"unavailable_primary baseline_id {t['baseline_id']!r} unknown")
        if t["endpoint"] == "M4":
            raise FreezeError("M4 is not a primary endpoint")
        key = (t["baseline_id"], t["endpoint"])
        if key in seen:
            raise FreezeError(f"duplicate primary test {key}")
        seen.add(key)
    need = {(b, e) for b, es in PRIMARY_SPEC.items() for e in es}
    if seen != need or len(pt) != N_PRIMARY:
        raise FreezeError(f"primary family must be exactly the {N_PRIMARY} pre-declared tests {sorted(need)}; "
                          f"got {sorted(seen)}")
    for c in cs:
        if c.get("role") == "baseline" and c.get("baseline_id") in PRIMARY_SPEC:
            if any(u.get("baseline_id") == c["baseline_id"] for u in cfg.get("unavailable_primary", [])):
                raise FreezeError(f"{c['baseline_id']} is both a config and unavailable_primary")
    del by
    return cfg


def check_freeze(freeze: Path = FREEZE, require_git: bool = True, allow_nonstandard_primary: bool = False) -> dict:
    if not freeze.exists():
        raise FreezeError(f"{freeze} does not exist: TEST evaluation is only allowed after the freeze (G4)")
    if require_git:
        r = _git("ls-files", "--error-unmatch", freeze.name, cwd=freeze.parent)
        if r.returncode != 0:
            raise FreezeError("FREEZE.json is not committed (git ls-files does not list it)")
        r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", freeze.name], cwd=freeze.parent)
        if r.returncode != 0:
            raise FreezeError("FREEZE.json has uncommitted changes")
    return validate_freeze(json.loads(freeze.read_text()), allow_nonstandard_primary=allow_nonstandard_primary)


def check_real_policy(cfg: dict) -> None:
    """Real mode: B == 10000; harness_commit == A28 commit and frozen harness files unchanged since it."""
    if int(cfg.get("B", 0)) != REAL_B:
        raise FreezeError(f"real mode requires B == {REAL_B} (A19); FREEZE has {cfg.get('B')}")
    hc = str(cfg.get("harness_commit", ""))
    a = _git("rev-parse", "--verify", f"{hc}^{{commit}}")
    b = _git("rev-parse", "--verify", f"{A28_COMMIT}^{{commit}}")
    if a.returncode != 0 or b.returncode != 0 or a.stdout.strip() != b.stdout.strip():
        raise FreezeError(f"harness_commit {hc!r} must be the A28 harness commit {A28_COMMIT}")
    r = subprocess.run(["git", "diff", "--quiet", hc, "--", *FROZEN_FILES], cwd=ROOT)
    if r.returncode != 0:
        raise FreezeError(f"frozen harness files differ from {hc} (A28): {FROZEN_FILES}")


def check_seal(seal: Path = SEAL, sealed_dir: Path = SEALED_DIR) -> dict:
    s = json.loads(Path(seal).read_text())
    for name, meta in s["files"].items():
        f = Path(sealed_dir) / name
        if not f.exists() or _sha(f) != meta["sha256"]:
            raise FreezeError(f"sealed file {name} missing or modified (SHA-256 mismatch with SEAL.json)")
    return s


# ================================================================================================ ledger
def code_hashes() -> dict:
    return {p.name: _sha(p) for p in sorted((ROOT / "src" / "ksp").glob("*.py"))}


def identity_hash(ch: dict) -> str:
    return hashlib.sha256(json.dumps({k: ch.get(k) for k in IDENTITY_FILES}, sort_keys=True).encode()).hexdigest()


def input_hashes(cfg: dict) -> dict:
    return {c[k]: c[f"{k}_sha256"] for c in cfg["configs"] for k in FILE_KEYS if c.get(k)}


def ledger_events(ledger: Path) -> list[dict]:
    if not Path(ledger).exists():
        return []
    return [json.loads(m.group(1)) for m in MARK.finditer(Path(ledger).read_text())]


def _ledger_append(ledger: Path, title: str, ev: dict, body: str = "") -> None:
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger, "a") as fh:
        fh.write(f"\n## {title}\n\n<!-- KSP-LEDGER {json.dumps(ev, sort_keys=True)} -->\n\n{body}")
        fh.flush()
        os.fsync(fh.fileno())


def check_ledger(ledger: Path, fsha: str, inputs: dict, code_id: str, override: str | None = None,
                 tag: str = "TEST") -> int:
    """Return the attempt number for this batch or raise (single-batch + crash-safety rules)."""
    ev = [e for e in ledger_events(ledger) if e.get("tag") == tag]
    if any(e["event"] == "COMPLETED" for e in ev) and not override:
        done = [e["freeze_sha"][:12] for e in ev if e["event"] == "COMPLETED"]
        raise FreezeError(f"a {tag} batch was already COMPLETED ({done}); single-batch rule. Only escape: "
                          f"--override-completed \"<reason>\" (logged)")
    started = [e for e in ev if e["event"] == "STARTED"]
    completed_after = {e.get("batch_id") for e in ev if e["event"] == "COMPLETED"}
    open_ = [e for e in started if e.get("batch_id") not in completed_after]
    if open_:
        last = open_[-1]
        if last["freeze_sha"] != fsha or last["inputs"] != inputs or last["code_identity"] != code_id:
            raise FreezeError("an unfinished STARTED batch exists with a different FREEZE sha / inputs / code; "
                              "a rerun must be identical")
        return int(last["attempt"]) + 1
    return 1


def check_single_batch(freeze_sha: str, ledger: Path = LEDGER) -> None:
    """Refuse if ANY TEST batch was completed (whatever the FREEZE sha)."""
    if any(e["event"] == "COMPLETED" and e.get("tag") == "TEST" for e in ledger_events(ledger)):
        raise FreezeError("a TEST batch was already COMPLETED (single-batch rule)")


# ================================================================================================ predictions
def candidate_pairs(universe: np.ndarray) -> pd.DataFrame:
    """(kinase_idx, site_idx) for every real kinase x universe site of matching residue class."""
    from . import eval as E
    from . import index
    S = index.load_sites()
    K = index.load_kinases()
    uni = np.where(universe)[0]
    res = S.res.values[uni]
    return pd.concat([pd.DataFrame({"kinase_idx": np.int32(i), "site_idx": uni[E._res_ok(res, c)].astype(np.int32)})
                      for i, c in enumerate(K.residue_class)], ignore_index=True)


def test_universe() -> np.ndarray:
    from . import index, splits
    S = index.load_sites()
    return S.in_U_ms.values & S.acc.isin(splits.load_test_proteins()).values


def test_candidate_pairs() -> pd.DataFrame:
    """Public TEST candidate set (no labels): ~25.9M (kinase_idx, site_idx)."""
    return candidate_pairs(test_universe())


def _to_idx(p: pd.DataFrame) -> pd.DataFrame:
    from . import index
    if {"kinase_idx", "site_idx"} <= set(p):
        return p
    if not {"kinase_acc", "site_id"} <= set(p):
        raise FreezeError("need (kinase_idx, site_idx) or (kinase_acc, site_id) columns")
    K, S = index.load_kinases(), index.load_sites()
    p = p.copy()
    p["kinase_idx"] = p.kinase_acc.map({a: i for i, a in enumerate(K.uniprot_acc)})
    p["site_idx"] = p.site_id.map(pd.Series(np.arange(len(S)), index=S.site_id))
    if p[["kinase_idx", "site_idx"]].isna().any().any():
        raise FreezeError("unknown kinase_acc/site_id (canaries are not allowed in TEST files)")
    p["kinase_idx"] = p.kinase_idx.astype(np.int64)
    p["site_idx"] = p.site_idx.astype(np.int64)
    return p


def load_predictions(c: dict, key: str = "predictions", root: Path = ROOT) -> pd.DataFrame:
    """Load a predictions file (+ separate score_site file): canonical (kinase_idx, site_idx, score[, score_site])."""
    raw = pd.read_parquet(root / c[key])
    if list(raw.columns).count("score_site") > 1:
        raise FreezeError(f"config {c['name']}: duplicate score_site columns")
    p = _to_idx(raw)
    if list(raw.columns).count("assignment_prob") > 1:
        raise FreezeError(f"config {c['name']}: duplicate assignment_prob columns")
    p = p[[x for x in ["kinase_idx", "site_idx", "score", "score_site", "assignment_prob"] if x in p]]
    if key != "predictions":
        return p.drop(columns=[x for x in ("assignment_prob",) if x in p])
    for fkey, col in (("score_site_predictions", "score_site"), ("assignment_prob_predictions", "assignment_prob")):
        if not c.get(fkey):
            continue
        if col in p:
            raise FreezeError(f"config {c['name']}: {col} given both inline and as a separate file")
        ss_raw = pd.read_parquet(root / c[fkey])
        if ({"kinase_idx", "site_idx"} <= set(raw)) != ({"kinase_idx", "site_idx"} <= set(ss_raw)):
            raise FreezeError(f"config {c['name']}: predictions and {fkey} use different key columns")
        if list(ss_raw.columns).count(col) != 1:
            raise FreezeError(f"config {c['name']}: {fkey} must have exactly one {col} column")
        ss = _to_idx(ss_raw)[["kinase_idx", "site_idx", col]]
        if ss.duplicated(["kinase_idx", "site_idx"]).any():
            raise FreezeError(f"config {c['name']}: duplicate keys in {fkey}")
        p = p.merge(ss, on=["kinase_idx", "site_idx"], how="left")
        if p[col].isna().any():
            raise FreezeError(f"config {c['name']}: {col} missing for {int(p[col].isna().sum())} rows")
    return p


def check_frame(p: pd.DataFrame, cand: pd.DataFrame, n_sites: int, n_kinases: int, what: str) -> dict:
    """Pre-flight for one canonical predictions frame against a candidate set (raises FreezeError)."""
    if "score" not in p:
        raise FreezeError(f"{what}: no score column")
    for col in ("score", "score_site", "assignment_prob"):
        if col in p and not np.isfinite(p[col].values.astype(float)).all():
            raise FreezeError(f"{what}: non-finite {col}")
    if "assignment_prob" in p and ((p.assignment_prob < 0) | (p.assignment_prob > 1)).any():
        raise FreezeError(f"{what}: assignment_prob must be a probability in [0, 1]")
    if (p.kinase_idx < 0).any() or (p.kinase_idx >= n_kinases).any():
        raise FreezeError(f"{what}: kinase_idx out of range (canaries are not allowed)")
    if (p.site_idx < 0).any() or (p.site_idx >= n_sites).any():
        raise FreezeError(f"{what}: site_idx out of range")
    kp = p.kinase_idx.values.astype(np.int64) * n_sites + p.site_idx.values
    if len(np.unique(kp)) != len(kp):
        raise FreezeError(f"{what}: duplicate (kinase, site) rows")
    kc = cand.kinase_idx.values.astype(np.int64) * n_sites + cand.site_idx.values
    extra = ~np.isin(kp, kc)
    need = cand.groupby("kinase_idx").size()
    have = p[~extra].groupby("kinase_idx").size().reindex(need.index).fillna(0)
    partial = need.index[(have > 0) & (have < need)].tolist()
    return {"rows": int(len(p)), "rows_outside_candidates": int(extra.sum()),
            "kinases_full": int((have == need).sum()), "kinases_unscored": int((have == 0).sum()),
            "kinases_partial": len(partial), "partial_kinase_idx": partial[:20],
            "scored_kinase_idx": need.index[have == need].tolist(),
            "ok": bool(extra.sum() == 0 and not partial)}


def check_predictions(path, universe: np.ndarray | None = None) -> dict:
    """Owner pre-flight for one file (never reads labels)."""
    from . import index
    uni = universe if universe is not None else test_universe()
    p = _to_idx(pd.read_parquet(path))
    r = check_frame(p, candidate_pairs(uni), len(index.load_sites()), len(index.load_kinases()), str(path))
    r.pop("scored_kinase_idx")
    return r


def preflight(cfg: dict, universe: np.ndarray, k_any: list, root: Path = ROOT) -> dict:
    """Pre-flight every predictions / T3 file of the FREEZE (fail fast with a clear message)."""
    from . import index
    K = index.load_kinases()
    kpos = {a: i for i, a in enumerate(K.uniprot_acc)}
    cand = candidate_pairs(universe)
    nS, nK = len(index.load_sites()), len(K)
    rep = {}
    for c in cfg["configs"]:
        p = load_predictions(c, root=root)
        r = check_frame(p, cand, nS, nK, f"config {c['name']} predictions")
        if not r["ok"]:
            raise FreezeError(f"config {c['name']}: predictions not exhaustive per kinase or outside the TEST "
                              f"candidates: {r['kinases_partial']} partial kinases, {r['rows_outside_candidates']} "
                              f"rows outside")
        rep[c["name"]] = {k: v for k, v in r.items() if k != "scored_kinase_idx"}
        if c.get("t3_predictions"):
            p3 = load_predictions(c, key="t3_predictions", root=root)
            r3 = check_frame(p3, cand, nS, nK, f"config {c['name']} t3_predictions")
            if not r3["ok"]:
                raise FreezeError(f"config {c['name']}: T3 predictions not exhaustive per kinase")
            miss = sorted(set(kpos[k] for k in k_any if k in kpos) - set(r3["scored_kinase_idx"]))
            if miss:
                raise FreezeError(f"config {c['name']}: T3 predictions miss {len(miss)} K_any kinases "
                                  f"(e.g. {[K.gene[i] for i in miss[:5]]}); every K_any kinase has an LKO fold")
            rep[c["name"] + ":T3"] = {k: v for k, v in r3.items() if k != "scored_kinase_idx"}
    return rep


def matched_pair(pa: pd.DataFrame, pb: pd.DataFrame, n_sites: int) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Matched universe (PRIOR-1) on canonical (kinase_idx, site_idx) keys."""
    ka = pa.kinase_idx.values.astype(np.int64) * n_sites + pa.site_idx.values
    kb = pb.kinase_idx.values.astype(np.int64) * n_sites + pb.site_idx.values
    common = np.intersect1d(ka, kb)
    a = pa[np.isin(ka, common)].reset_index(drop=True)
    b = pb[np.isin(kb, common)].reset_index(drop=True)
    return a, b, {"n_A": len(pa), "n_B": len(pb), "n_common": int(len(common)),
                  "frac_A_kept": len(common) / max(len(pa), 1), "frac_B_kept": len(common) / max(len(pb), 1),
                  "matched_kinases": int(a.kinase_idx.nunique())}


# ================================================================================================ strata
def _strata_module():
    """Canonical stratum definitions from src/ksp/strata.py; None -> built-in fallback."""
    try:
        from . import strata as S
    except ImportError:
        return None
    return S


def kinase_strata_table(data) -> pd.DataFrame:
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


def pair_strata_table(data, require_columns: bool = True) -> pd.DataFrame:
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")].copy()
    for col in ("evidence_native", "new_since_2021"):
        if require_columns and col not in pm:
            raise FreezeError(f"evaluation positives lack column {col!r} (needed for the §6.4 strata)")
    pm["site_id"] = data.sites.site_id.values[pm.site_idx.values]
    out = [pm.assign(stratum="temporal", value=np.where(pm.new_since_2021.astype(bool), "new_since_2021",
                                                        "pre_2021"))] if "new_since_2021" in pm else []
    S = _strata_module()
    if S is not None:
        K = data.kinases[~data.kinases.uniprot_acc.str.startswith("CANARY_")]
        if "evidence_native" in pm:
            out.append(pm.assign(stratum="evidence", value=S.evidence_class(pm).values))
        hubs = S.hub_sites(pm[["kinase_acc", "site_id"]], data.par)
        out.append(pm.assign(stratum="hub", value=np.where(pm.site_id.isin(hubs), "hub", "non_hub")))
        for nm, offs in [("primed", S.priming_curated(K)), ("primed_cantley", S.priming_cantley(K))]:
            m = S.primed_mask(pm[["kinase_acc", "site_id"]], offs, sites=data.sites)
            sel = pm.kinase_acc.isin(set(offs)).values
            out.append(pm[sel].assign(stratum=nm, value=np.where(m[sel], "primed", "unprimed")))
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


def strata_results(res, data, name: str, require_columns: bool = True) -> pd.DataFrame:
    """Endpoints on every stratum (point estimates, missing='worst' results required)."""
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
    ps = pair_strata_table(data, require_columns)
    sl = res.site_level.set_index(["kinase_acc", "site_idx"]) if len(res.site_level) else None
    for (st, val), g in ps.groupby(["stratum", "value"]):
        per = {}
        for k, gg in g.groupby("kinase_acc"):
            it = res.items.get(k)
            if it is None:
                continue
            pos = data.positives[k]
            if len(it["pos_s"]) != len(pos):
                raise FreezeError("pair strata need missing='worst' results (positives aligned to data.positives)")
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
    # external curated (separate positive set vs the same E-hard negatives): K_eval row -> M1/M3, K_any row -> M4
    ext = [r for r in res.pair_strata.to_dict(orient="records") if r["stratum"] == "external_curated"] \
        if len(res.pair_strata) else []
    for r in ext:
        if r["kinase_set"] == "K_primary":
            rows.append({"config": name, "level": "pair", "stratum": "external_curated", "value": "all",
                         "n_kinases_K_eval": r["n_kinases"], "n_pairs": r["n_pos"], "M1": r["auroc_hard"],
                         "M3": r["top1_recall"]})
        else:
            rows.append({"config": name, "level": "pair", "stratum": "external_curated", "value": "all",
                         "n_kinases_K_any": r["n_kinases"], "n_pairs": r["n_pos"], "M4": r["auroc_hard"]})
    return pd.DataFrame(rows)


def ext_tier_results(res, data, P: pd.DataFrame, labels_df: pd.DataFrame | None, name: str) -> pd.DataFrame:
    """External-curated stratum split by ext_tier (A9: A = SIGNOR/HPRD/UniProt-ECO, B = other): per kinase AUROC of
    the tier's external pairs (not evaluation positives) vs the kinase's E-hard negatives; macro over K_eval (M1)
    and K_any (M4)."""
    from . import eval as E
    if labels_df is None or "ext_tier" not in labels_df:
        return pd.DataFrame()
    S = data.sites
    spos = pd.Series(np.arange(len(S)), index=S.site_id)
    kpos = {a: i for i, a in enumerate(data.kinases.uniprot_acc)}
    x = labels_df[labels_df.ext_tier.isin(["A", "B"]) & labels_df.site_id.isin(spos.index)]
    x = x.assign(site_idx=spos[x.site_id].values)
    x = x[data.universe[x.site_idx.values]]
    sc = pd.Series(P.score.values, index=P.kinase_idx.values.astype(np.int64) * len(S) + P.site_idx.values)
    rows = []
    for tier, g in x.groupby("ext_tier"):
        per = {}
        for k, gg in g.groupby("kinase_acc"):
            it = res.items.get(k)
            if it is None or k not in kpos:
                continue
            s_idx = np.setdiff1d(gg.site_idx.unique(), data.positives.get(k, []))
            v = sc.reindex(kpos[k] * len(S) + s_idx).values
            v = np.where(np.isnan(v), -np.inf, v)
            fl = np.min(np.r_[it["hard_s"], it["pos_s"], 0.0]) - 1
            v = np.where(np.isfinite(v), v, fl)
            if len(v):
                per[k] = (E.auroc(v, it["hard_s"]), len(v))
        kp_ = [k for k in per if k in set(data.k_primary)]
        ka_ = [k for k in per if k in set(data.k_any)]
        rows.append({"config": name, "stratum": "external_curated", "value": f"tier_{tier}",
                     "n_kinases_K_eval": len(kp_), "n_kinases_K_any": len(ka_),
                     "n_pairs": int(sum(per[k][1] for k in ka_)), "M1": _macro(per[k][0] for k in kp_),
                     "M4": _macro(per[k][0] for k in ka_)})
    return pd.DataFrame(rows)


# ================================================================================================ A17 reliability
def reliability(P: pd.DataFrame, data, site_truth: dict, subset: str, name: str, n_bins: int = 10) -> pd.DataFrame:
    """A17: for each labelled site with scored candidates, top kinase = argmax assignment_prob; correct if it is one
    of the site's true kinases. Reliability bins of the top assignment probability; ECE in every row."""
    if "assignment_prob" not in P or not site_truth:
        return pd.DataFrame()
    x = P[P.site_idx.isin(list(site_truth))]
    if not len(x):
        return pd.DataFrame()
    top = x.loc[x.groupby("site_idx").assignment_prob.idxmax()]
    kacc = data.kinases.uniprot_acc.values
    corr = np.array([kacc[k] in site_truth[s] for k, s in zip(top.kinase_idx.values, top.site_idx.values)])
    conf = top.assignment_prob.values
    b = np.minimum((conf * n_bins).astype(int), n_bins - 1)
    rows, ece = [], 0.0
    for i in range(n_bins):
        m = b == i
        if m.any():
            ece += m.sum() / len(conf) * abs(corr[m].mean() - conf[m].mean())
            rows.append({"config": name, "subset": subset, "bin": f"[{i / n_bins:.1f},{(i + 1) / n_bins:.1f})",
                         "n_sites": int(m.sum()), "mean_p": float(conf[m].mean()), "accuracy": float(corr[m].mean())})
    d = pd.DataFrame(rows)
    d["ECE"] = ece
    d["n_total"] = len(conf)
    d["accuracy_overall"] = float(corr.mean())
    return d


def _site_truths(data) -> dict:
    pm = data.pos_meta[~data.pos_meta.kinase_acc.str.startswith("CANARY_")]
    out = {"labelled_TEST_sites": pm.groupby("site_idx").kinase_acc.agg(set).to_dict()}
    if "new_since_2021" in pm:
        nw = pm[pm.new_since_2021.astype(bool)]
        out["new_since_2021_sites"] = nw.groupby("site_idx").kinase_acc.agg(set).to_dict()
    ext = {}
    for k, idx in data.ext_positives.items():
        for s in idx:
            ext.setdefault(int(s), set()).add(k)
    out["external_curated_sites"] = ext
    return out


# ================================================================================================ comparisons
def _ekey(ep: str) -> str:
    return {"M1": "M1_auroc_hard", "M2": "M2_hit5", "M3": "M3_top1_recall", "M4": "M4_auroc_hard_Kany"}[ep]


def run_comparisons(cfg: dict, preds: dict, data, B: int, kinases: list | None = None) -> pd.DataFrame:
    """Paired cluster-bootstrap tests on matched universes; Holm over the primary family (incl. unavailable tests
    at p = 1) and separately over secondary tests."""
    from . import eval as E
    by = {c["name"]: c for c in cfg["configs"]}
    rows = []
    for cp in cfg.get("comparisons", []):
        a, b = cp["a"], cp["b"]
        pa, pb, rep = matched_pair(preds[a], preds[b], len(data.sites))
        ra = E.evaluate(pa, data, missing="drop", kinases=kinases)
        rb = E.evaluate(pb, data, missing="drop", kinases=kinases)
        d = E.paired_bootstrap(ra, rb, B=B, seed=0).set_index("endpoint")
        for ep in cp.get("endpoints", ["M1", "M2", "M3"]):
            r = d.loc[ep]
            rows.append({"comparison": cp.get("label", f"{a} vs {b}"), "a": a, "b": b,
                         "baseline_id": by[b].get("baseline_id", ""), "endpoint": ep,
                         "primary": bool(cp.get("primary", False)), "status": "tested", "diff": r["diff"],
                         "lo95": r.lo95, "hi95": r.hi95, "p_boot": r.p_boot, "B": B,
                         "A_matched": ra.endpoints[_ekey(ep)], "B_matched": rb.endpoints[_ekey(ep)],
                         "A_M1_n_kinases_defined": ra.endpoints["M1_n_kinases_defined"],
                         "B_M1_n_kinases_defined": rb.endpoints["M1_n_kinases_defined"],
                         "matched_kinases": rep["matched_kinases"], "n_common_pairs": rep["n_common"],
                         "frac_A_pairs_kept": rep["frac_A_kept"], "frac_B_pairs_kept": rep["frac_B_kept"],
                         "M2_score_column_A": "score_site" if "score_site" in pa else "score",
                         "M2_score_column_B": "score_site" if "score_site" in pb else "score"})
    for u in cfg.get("unavailable_primary", []):
        for ep in PRIMARY_SPEC.get(u["baseline_id"], ()):
            rows.append({"comparison": f"final vs {u['baseline_id']}", "a": None, "b": None,
                         "baseline_id": u["baseline_id"], "endpoint": ep, "primary": True,
                         "status": f"NOT TESTED (baseline unavailable: {u.get('reason', '')})", "p_boot": 1.0,
                         "B": B})
    C = pd.DataFrame(rows)
    if not len(C):
        return C
    C["p_holm"] = np.nan
    for flag in (True, False):
        m = C.primary == flag
        if m.any():
            C.loc[m, "p_holm"] = E.holm(C.loc[m, "p_boot"].values)
    C["beats"] = (C.status == "tested") & (C.p_holm < 0.05) & (C.lo95 > 0)
    return C


# ================================================================================================ batch
def run_batch(cfg: dict, data, outd: Path, B: int, root: Path = ROOT, B_secondary: int | None = None,
              kinases: list | None = None, labels_df: pd.DataFrame | None = None,
              require_strata_columns: bool = True) -> dict:
    """Score every config in one batch on `data` (ksp.eval.EvalData); write outputs into outd."""
    from . import eval as E
    Bs = int(B_secondary or cfg.get("B_secondary", 1000))
    outd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    preds, summ = {}, {}
    truths = _site_truths(data)
    for c in cfg["configs"]:
        p = load_predictions(c, root=root)
        preds[c["name"]] = p
        res = E.evaluate(p, data, missing="worst", kinases=kinases)
        res.save_json(outd / f"{c['name']}.json")
        ci = E.bootstrap(res, B=Bs, seed=0)
        ci.to_csv(outd / f"{c['name']}_ci.tsv", sep="\t", index=False)
        cik = E.bootstrap(res, B=Bs, seed=0, resample="kinases")
        cik.to_csv(outd / f"{c['name']}_ci_kinases.tsv", sep="\t", index=False)
        res.per_kinase.to_csv(outd / f"{c['name']}_per_kinase.tsv", sep="\t", index=False)
        sl = res.site_level.copy()
        if len(sl):
            sl.insert(1, "site_id", data.sites.site_id.values[sl.site_idx.values])
        sl.to_csv(outd / f"{c['name']}_site_level.tsv", sep="\t", index=False)
        strata_results(res, data, c["name"], require_strata_columns).to_csv(outd / f"{c['name']}_strata.tsv",
                                                                             sep="\t", index=False)
        ext_tier_results(res, data, p, labels_df, c["name"]).to_csv(outd / f"{c['name']}_ext_tiers.tsv", sep="\t",
                                                                    index=False)
        rel = pd.concat([reliability(p, data, t, sub, c["name"]) for sub, t in truths.items()], ignore_index=True)
        rel.to_csv(outd / f"{c['name']}_reliability.tsv", sep="\t", index=False)
        cov = None
        if "assignment_prob" in p:              # A17 Coverage@tau over the TEST U_ms universe (model-based)
            q = p[data.universe[p.site_idx.values]]
            lab = set(data.sites.site_id.values[list(truths["labelled_TEST_sites"])])
            cov = E.coverage_at(pd.DataFrame({"site_id": data.sites.site_id.values[q.site_idx.values],
                                              "kinase_acc": data.kinases.uniprot_acc.values[q.kinase_idx.values],
                                              "p_assign": q.assignment_prob.values}),
                                unlabelled_sites=set(data.sites.site_id.values[data.universe]) - lab)
        e = res.endpoints
        s = {k: e.get(_ekey(k)) for k in ("M1", "M2", "M3", "M4")}
        s.update({"n_missing_pairs": e["n_missing_pairs"], "n_required_pairs": e["n_required_pairs"],
                  "M1_n_kinases_defined": e["M1_n_kinases_defined"], "role": c.get("role", ""),
                  "baseline_id": c.get("baseline_id", ""), "runtime": c["runtime"],
                  "M2_score_column": "score_site" if "score_site" in p else "score",
                  "ci_clusters": ci.set_index("endpoint")[["lo95", "hi95"]].to_dict(orient="index"),
                  "ci_kinases": cik.set_index("endpoint")[["lo95", "hi95"]].to_dict(orient="index"),
                  "A17_ECE": {r.subset: r.ECE for r in rel.drop_duplicates("subset").itertuples()} if len(rel)
                  else None, "A17_coverage": cov})
        if c.get("t3_predictions"):
            p3 = load_predictions(c, key="t3_predictions", root=root)
            r3 = E.evaluate(p3, data, missing="worst", kinases=kinases)
            r3.save_json(outd / f"{c['name']}_t3.json")
            r3.per_kinase.to_csv(outd / f"{c['name']}_t3_per_kinase.tsv", sep="\t", index=False)
            pk3 = r3.per_kinase
            req = pk3.n_pos_required + pk3.n_neg_all
            s["T3"] = {k: r3.endpoints.get(_ekey(k)) for k in ("M1", "M2", "M3", "M4")}
            s["T3"].update({"n_missing_pairs": r3.endpoints["n_missing_pairs"],
                            "n_missing_kinases_K_any": int((pk3.in_K_any & (pk3.n_missing >= req)
                                                            & (req > 0)).sum()),
                            "note": "T3 M2 mixes LKO models (each kinase scored by its own family-held-out model)"})
        summ[c["name"]] = s
    C = run_comparisons(cfg, preds, data, B, kinases=kinases)
    if len(C):
        C[C.primary].to_csv(outd / "superiority.tsv", sep="\t", index=False)
        C[~C.primary].to_csv(outd / "secondary.tsv", sep="\t", index=False)
    summ["_comparisons"] = C.to_dict(orient="records")
    summ["_runtime_s"] = round(time.time() - t0, 1)
    summ["_B"] = {"primary": B, "secondary": Bs}
    _write_summary(summ, outd)
    return summ


def _write_summary(summ: dict, outd: Path) -> None:
    f = lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.4f}"  # noqa: E731
    L = ["**Standalone metrics under missing='worst' (unscored pairs ranked last); NOT head-to-head - see the "
         "matched comparisons below.**", "",
         "| config | role | M1 | M2 Hit@5 (score col) | M3 top1% | M4 | T3 M1 | T3 missing kinases | "
         "missing/required | train s | infer s |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, s in summ.items():
        if n.startswith("_"):
            continue
        t3 = s.get("T3", {})
        L.append(f"| {n} | {s['role']} | {f(s['M1'])} | {f(s['M2'])} ({s['M2_score_column']}) | {f(s['M3'])} | "
                 f"{f(s['M4'])} | {f(t3.get('M1'))} | {t3.get('n_missing_kinases_K_any', '')} | "
                 f"{s['n_missing_pairs']}/{s['n_required_pairs']} | {s['runtime']['train_wall_s']} | "
                 f"{s['runtime']['inference_wall_s']} |")
    if any("T3" in s for n, s in summ.items() if not n.startswith("_")):
        L += ["", "T3: each kinase scored by the model trained with its LKO family removed; T3 M2 mixes LKO models."]
    C = pd.DataFrame(summ.get("_comparisons", []))
    if len(C):
        for flag, title in [(True, f"Primary superiority tests (m = {N_PRIMARY}; Holm over these; "
                                   f"B = {summ['_B']['primary']})"),
                            (False, "Secondary comparisons (Holm separately)")]:
            x = C[C.primary == flag]
            if not len(x):
                continue
            L += ["", f"**{title}**", "",
                  "| comparison | endpoint | status | diff | 95% CI | p_boot | p_holm | beats | matched kinases | "
                  "frac A kept |", "|---|---|---|---|---|---|---|---|---|---|"]
            for r in x.itertuples():
                L.append(f"| {r.comparison} | {r.endpoint} | {r.status} | {f(getattr(r, 'diff', np.nan))} | "
                         f"[{f(getattr(r, 'lo95', np.nan))}, {f(getattr(r, 'hi95', np.nan))}] | {f(r.p_boot)} | "
                         f"{f(r.p_holm)} | {bool(r.beats)} | "
                         f"{'' if pd.isna(getattr(r, 'matched_kinases', np.nan)) else int(r.matched_kinases)} | "
                         f"{f(getattr(r, 'frac_A_pairs_kept', np.nan))} |")
        L += ["", "Footnote: 95% CIs are unadjusted paired percentile intervals; 'beats' requires Holm-adjusted "
              "p < 0.05 AND CI lower bound > 0. Comparisons use the matched universe (PRIOR-1; matched-kinase M2, "
              "A27)."]
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


# ================================================================================================ entry points
def main(freeze: Path = FREEZE, seal: Path = SEAL, sealed_dir: Path = SEALED_DIR, ledger: Path = LEDGER,
         out_root: Path | None = None, mode: str = "real", data_builder=None, kinases: list | None = None,
         allow_nonstandard_primary: bool = False, override_completed: str | None = None,
         universe_builder=None) -> dict:
    """The batch. mode='real': all guards (git, B == 10000, harness commit). mode='dry': same path with the git /
    B / harness policy relaxed (tests and rehearsal; tag DRYRUN in the ledger)."""
    from . import eval as E
    from . import splits
    real = mode == "real"
    if real and (data_builder is not None or allow_nonstandard_primary or kinases is not None):
        raise FreezeError("real mode takes no test hooks")
    tag = "TEST" if real else "DRYRUN"
    cfg = check_freeze(freeze, require_git=real, allow_nonstandard_primary=allow_nonstandard_primary)
    if real:
        check_real_policy(cfg)
    fsha = _sha(freeze)
    inputs = input_hashes(cfg)
    if data_builder is not None:                 # dry run: toy data (no sealed labels needed for the universe)
        data_pre = data_builder()
        uni, k_any = data_pre.universe, data_pre.k_any
    else:
        uni = universe_builder() if universe_builder else test_universe()
        es = splits.load_eval_kinase_sets()
        k_any = es.kinase_acc[es.in_K_any].tolist()
    pre = preflight(cfg, uni, k_any)
    check_seal(seal, sealed_dir)
    ch = code_hashes()
    cid = identity_hash(ch)
    ledger = Path(ledger)
    if override_completed:
        _ledger_append(ledger, f"OVERRIDE of single-batch rule {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
                       {"event": "OVERRIDE", "tag": tag, "freeze_sha": fsha, "reason": override_completed},
                       f"Reason: {override_completed}\n")
    attempt = check_ledger(ledger, fsha, inputs, cid, override=override_completed, tag=tag)
    batch_id = f"{fsha[:12]}-{int(time.time())}" if attempt == 1 else \
        [e for e in ledger_events(ledger) if e["event"] == "STARTED" and e.get("tag") == tag][-1]["batch_id"]
    head = _git("rev-parse", "HEAD").stdout.strip()
    _ledger_append(ledger, f"STARTED {tag} batch {time.strftime('%Y-%m-%d %H:%M:%S %Z')} (attempt {attempt})",
                   {"event": "STARTED", "tag": tag, "batch_id": batch_id, "attempt": attempt, "freeze_sha": fsha,
                    "inputs": inputs, "git_head": head, "code_sha256": ch, "code_identity": cid},
                   f"FREEZE.json sha256 {fsha}; git HEAD {head}; harness {cfg.get('harness_commit')}; "
                   f"{len(cfg['configs'])} configs; pre-flight {json.dumps({k: v['rows'] for k, v in pre.items()})}\n")
    test_labels = pd.read_parquet(Path(sealed_dir) / "test_labels.parquet")
    data = data_builder(test_labels) if data_builder else E.prepare("test", labels_df=test_labels,
                                                                     with_canaries=False)
    base = out_root or ledger.parent
    final_dir, partial = base / fsha[:12], base / f"{fsha[:12]}.partial"
    if final_dir.exists():
        raise FreezeError(f"{final_dir} already exists")
    if partial.exists():
        shutil.rmtree(partial)
    summ = run_batch(cfg, data, partial, int(cfg.get("B", REAL_B)), kinases=kinases, labels_df=test_labels)
    (partial / "preflight.json").write_text(json.dumps(pre, indent=1, default=_jd))
    os.replace(partial, final_dir)
    _ledger_append(ledger, f"COMPLETED {tag} batch {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
                   {"event": "COMPLETED", "tag": tag, "batch_id": batch_id, "freeze_sha": fsha,
                    "attempt": attempt, "outputs": str(final_dir)},
                   f"Configs {len(cfg['configs'])}; B={cfg.get('B')}; runtime {summ['_runtime_s']}s; "
                   f"outputs {final_dir}\n\n" + (final_dir / "summary.md").read_text())
    print(f"scored {len(cfg['configs'])} configs -> {final_dir}; ledger {ledger}")
    return summ


def toy_data_builder(fold: int = 0):
    """Dry-run data: DEV CV fold `fold` as a fake TEST (labels from load_dev; never sealed/). K_eval / K_any are
    recomputed on the toy split with the frozen §4.6 rules."""
    from . import eval as E
    cache = {}

    def build(_labels_ignored=None):
        if "d" not in cache:
            d = E.prepare("cv", fold=fold, with_canaries=False)
            ke = [k for k, p in d.positives.items() if len(p) >= 3 and len(np.unique(d.site_cluster[p])) >= 2]
            cache["d"] = dataclasses.replace(d, split="test", k_primary=sorted(ke), k_any=sorted(d.positives))
        return cache["d"]
    return build


def write_toy_sealed(sealed_dir: Path, seal: Path, fold: int = 0) -> None:
    from . import labels, splits
    sealed_dir.mkdir(parents=True, exist_ok=True)
    d = labels.load_dev("all")
    folds = splits.load_cv_folds()
    d[d.acc.map(folds).eq(fold)].to_parquet(sealed_dir / "test_labels.parquet", index=False)
    seal.write_text(json.dumps({"files": {"test_labels.parquet": {"sha256": _sha(sealed_dir / "test_labels.parquet")}}}))


def pin_files(cfg: dict, root: Path = ROOT) -> dict:
    """Fill *_sha256 for every referenced file (helper for assembling FREEZE.json)."""
    for c in cfg["configs"]:
        for k in FILE_KEYS:
            if c.get(k):
                c[f"{k}_sha256"] = _sha(root / c[k])
    return cfg


def dry_run(freeze_cfg: Path | None = None, fold: int = 0, B: int = 1000,
            work: Path = ROOT / "cache" / "eval_test_dryrun") -> dict:
    """Rehearsal of the full batch path on DEV fold `fold`; every path under `work` (fresh ledger each time)."""
    work.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(freeze_cfg).read_text()) if freeze_cfg else dry_run_freeze(fold, work)
    cfg["B"] = B
    pin_files(cfg)
    fz = work / "FREEZE.json"
    fz.write_text(json.dumps(cfg, indent=1))
    write_toy_sealed(work / "sealed", work / "SEAL.json", fold)
    led = work / "LEDGER.md"
    if led.exists():
        led.unlink()
    sha12 = _sha(fz)[:12]
    for d in (work / sha12, work / f"{sha12}.partial"):
        if d.exists():
            shutil.rmtree(d)
    return main(fz, work / "SEAL.json", work / "sealed", led, out_root=work, mode="dry",
                data_builder=toy_data_builder(fold))


def dry_run_freeze(fold: int, work: Path) -> dict:
    """FREEZE-like config from existing runs restricted to CV fold `fold` (real kinases only): final = R2 fold run,
    B1 arms i/ii, B7, B0-prior; R3 final LKO pooled as T3 of final. B3/B4 have no fold runs -> exercised through
    the unavailable_primary contingency."""
    from . import index, splits
    runs = ROOT / "cache" / "runs"
    S = index.load_sites()
    sf = S.acc.map(splits.load_cv_folds()).fillna(-1).astype(int).values
    nK = len(index.load_kinases())

    def pick(prefix):
        c = sorted(p.name for p in runs.glob(prefix + "__*") if (p / "predictions.parquet").exists())
        if not c:
            raise FileNotFoundError(f"no run {prefix}")
        return c[-1]

    uni = S.in_U_ms.values & (sf == fold)
    need = candidate_pairs(uni).groupby("kinase_idx").size()
    dropped = {}

    def sub(run, name):
        p = _to_idx(pd.read_parquet(runs / run / "predictions.parquet"))
        p = p[(sf[p.site_idx.values] == fold) & uni[p.site_idx.values]]
        p = p[p.kinase_idx < nK][["kinase_idx", "site_idx", "score"]]
        have = p.groupby("kinase_idx").size()
        part = have.index[have.values < need.reindex(have.index).values]
        dropped[name] = int(len(part))          # owner rule: score a kinase fully or not at all
        p = p[~p.kinase_idx.isin(part)]
        out = work / f"{name}.parquet"
        p.to_parquet(out, index=False)
        return str(out.relative_to(ROOT))
    RT = {"train_wall_s": 0.0, "inference_wall_s": 0.0}
    final = pick(f"r2_final_s0_f{fold}")
    cfgs = [{"name": "final_R2", "role": "final", "run_id": final, "predictions": sub(final, "final_R2"),
             "t3_predictions": sub(pick("r3_final_distill_lko_pooled"), "final_T3"), "runtime": RT}]
    for nm, pre, role, bid in [("B1_arm_i", "simple_b1_published_reeval", "baseline", "B1_arm_i"),
                               ("B1_arm_ii", "simple_b1_universe_reeval", "baseline", "B1_arm_ii"),
                               ("B7", "simple_b7_lr_reeval", "baseline", "B7"),
                               ("B0_prior", "simple_b0_prior_reeval", "control", "B0_prior")]:
        r = pick(pre)
        cfgs.append({"name": nm, "role": role, "baseline_id": bid, "run_id": r, "predictions": sub(r, nm),
                     "runtime": RT})
    return {"version": 1, "harness_commit": "dry-run",
            "notes": f"dry run on DEV fold {fold}; partially scored kinases dropped (left unscored): {dropped}",
            "configs": cfgs,
            "comparisons": [{"a": "final_R2", "b": "B1_arm_i", "endpoints": ["M1", "M3"], "primary": True},
                            {"a": "final_R2", "b": "B1_arm_ii", "endpoints": ["M2"], "primary": True},
                            {"a": "final_R2", "b": "B7", "endpoints": ["M1", "M2", "M3"], "primary": False},
                            {"a": "final_R2", "b": "B0_prior", "endpoints": ["M1", "M2", "M3"], "primary": False}],
            "unavailable_primary": [{"baseline_id": "B3_retrained", "reason": "dry run: no fold-0 B3 run"},
                                    {"baseline_id": "B4_retrained", "reason": "dry run: no fold-0 B4 run"}]}


def _cli(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="G4 TEST batch (see module docstring)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--real", action="store_true", help="run the single real TEST batch (PI, G4)")
    g.add_argument("--dry-run", action="store_true", help="rehearsal on DEV fold (never touches sealed/)")
    g.add_argument("--check-predictions", default=None, help="owner pre-flight for a TEST predictions file")
    ap.add_argument("--fold", type=int, default=None)
    ap.add_argument("--B", type=int, default=None)
    ap.add_argument("--freeze", default=None, help="dry-run only: FREEZE-like json")
    ap.add_argument("--override-completed", default=None, metavar="REASON")
    a = ap.parse_args(argv)
    if not (a.real or a.dry_run or a.check_predictions):
        ap.print_usage()
        return 2
    if (a.B is not None or a.fold is not None or a.freeze) and not a.dry_run:
        ap.error("--B/--fold/--freeze are only allowed with --dry-run")
    if a.override_completed and not a.real:
        ap.error("--override-completed only applies to --real")
    if a.check_predictions:
        print(json.dumps(check_predictions(a.check_predictions), indent=1))
    elif a.dry_run:
        s = dry_run(Path(a.freeze) if a.freeze else None, a.fold or 0, a.B or 1000)
        print(json.dumps({k: v for k, v in s.items() if not k.startswith("_")}, indent=1, default=_jd))
    else:
        main(override_completed=a.override_completed)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
