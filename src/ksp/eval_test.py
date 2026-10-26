"""The ONE AND ONLY reader of sealed TEST labels (DESIGN §4.3, AIFM-2). Owner: core-infra.

Usage (PI only, at G4): /workspace/.venv/bin/python -I -m ksp.eval_test   (from /workspace/src) or
                        /workspace/.venv/bin/python -I src/ksp/eval_test.py
Refuses to run unless
  * /workspace/FREEZE.json exists, is tracked by git (`git ls-files --error-unmatch FREEZE.json`) and has no
    uncommitted changes (`git diff --quiet HEAD -- FREEZE.json`) - read-only git calls;
  * every sealed file matches its SHA-256 in results/splits/SEAL.json;
  * this FREEZE.json (by SHA-256) has not been scored before (results/test/LEDGER.md).
FREEZE.json: {"configs": [{"name": str, "predictions": "cache/runs/<run_id>/predictions.parquet",
              "run_id": str, (optional) "score_site_predictions": path}], "B": 10000, "notes": str}
All configurations are scored in ONE batch; per-config metrics -> results/test/<freeze_sha12>/<name>.json; one
ledger block is appended to results/test/LEDGER.md. Paired comparisons (FREEZE.json "comparisons": [[A, B], ...])
use ksp.eval.paired_bootstrap on matched universes with Holm correction across all comparisons x endpoints.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "FREEZE.json"
SEAL = ROOT / "results" / "splits" / "SEAL.json"
SEALED_DIR = ROOT / "sealed"
LEDGER = ROOT / "results" / "test" / "LEDGER.md"


class FreezeError(RuntimeError):
    pass


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check_freeze(freeze: Path = FREEZE) -> dict:
    if not freeze.exists():
        raise FreezeError(f"{freeze} does not exist: TEST evaluation is only allowed after the freeze (G4)")
    r = subprocess.run(["git", "ls-files", "--error-unmatch", freeze.name], cwd=freeze.parent,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise FreezeError("FREEZE.json is not committed (git ls-files does not list it)")
    r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", freeze.name], cwd=freeze.parent)
    if r.returncode != 0:
        raise FreezeError("FREEZE.json has uncommitted changes")
    cfg = json.loads(freeze.read_text())
    if not cfg.get("configs"):
        raise FreezeError("FREEZE.json lists no configs")
    return cfg


def check_seal() -> dict:
    seal = json.loads(SEAL.read_text())
    for name, meta in seal["files"].items():
        f = SEALED_DIR / name
        if not f.exists() or _sha(f) != meta["sha256"]:
            raise FreezeError(f"sealed file {name} missing or modified (SHA-256 mismatch with SEAL.json)")
    return seal


def main(freeze: Path = FREEZE) -> None:
    cfg = check_freeze(freeze)
    check_seal()
    fsha = _sha(freeze)
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    if LEDGER.exists() and fsha in LEDGER.read_text():
        raise FreezeError(f"FREEZE.json {fsha[:12]} was already scored (single-batch rule)")
    sys.path.insert(0, str(ROOT / "src"))
    from ksp import eval as E
    test_labels = pd.read_parquet(SEALED_DIR / "test_labels.parquet")
    data = E.prepare("test", labels_df=test_labels, with_canaries=False)
    B = int(cfg.get("B", 10000))
    outd = LEDGER.parent / fsha[:12]
    outd.mkdir(parents=True, exist_ok=True)
    results, lines = {}, []
    t0 = time.time()
    for c in cfg["configs"]:
        preds = pd.read_parquet(ROOT / c["predictions"])
        if c.get("score_site_predictions"):
            ss = pd.read_parquet(ROOT / c["score_site_predictions"])
            preds = preds.merge(ss[["kinase_acc", "site_id", "score_site"]], on=["kinase_acc", "site_id"],
                                how="left")
        res = E.evaluate(preds, data, missing=c.get("missing", "worst"))
        ci = E.bootstrap(res, B=B, seed=0)
        res.save_json(outd / f"{c['name']}.json")
        ci.to_csv(outd / f"{c['name']}_ci.tsv", sep="\t", index=False)
        results[c["name"]] = (preds, res)
        e = res.endpoints
        lines.append(f"| {c['name']} | {c.get('run_id', '')} | {e['M1_auroc_hard']:.4f} | {e['M2_hit5']:.4f} | "
                     f"{e['M3_top1_recall']:.4f} | {e['M4_auroc_hard_Kany']:.4f} | {e['n_missing_pairs']} |")
    comp_rows = []
    for a, b in cfg.get("comparisons", []):
        pa, pb, _ = E.matched_universe(results[a][0], results[b][0])
        ra, rb = E.evaluate(pa, data, missing="drop"), E.evaluate(pb, data, missing="drop")
        d = E.paired_bootstrap(ra, rb, B=B, seed=0)
        d.insert(0, "comparison", f"{a} - {b}")
        comp_rows.append(d)
    with open(LEDGER, "a") as fh:
        fh.write(f"\n## TEST batch {time.strftime('%Y-%m-%d %H:%M:%S %Z')} - FREEZE.json sha256 {fsha}\n\n")
        fh.write(f"Configs: {len(cfg['configs'])}; B={B}; runtime {time.time() - t0:.0f}s; details {outd.name}/\n\n")
        fh.write("| config | run_id | M1 | M2 Hit@5 | M3 top1% | M4 | missing pairs |\n|---|---|---|---|---|---|---|\n")
        fh.write("\n".join(lines) + "\n")
        if comp_rows:
            C = pd.concat(comp_rows, ignore_index=True)
            C["p_holm"] = E.holm(C.p_boot.values)
            C.to_csv(outd / "comparisons.tsv", sep="\t", index=False)
            fh.write("\nPaired comparisons (Holm over all rows):\n\n" + C.to_string(index=False) + "\n")
    print(f"scored {len(cfg['configs'])} configs -> {LEDGER}")


if __name__ == "__main__":
    main()
