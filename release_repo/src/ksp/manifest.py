"""Run manifests (DESIGN §11.2). Owner: core-infra.

    with RunManifest("r2_lgbm_v1", config=cfg, seed=0, inputs=[path1, path2]) as run:
        ...                                # train / predict
        run.count_fallback("no_cantley_matrix", n)        # every fallback is counted, never silent
        run.save_predictions(preds)        # -> cache/runs/<run_id>/predictions.parquet (+ n_rows, sha256)
        run.save_metrics(res)              # EvalResult or dict -> results/runs/<run_id>/metrics.json
    # on exit: results/runs/<run_id>/manifest.json (git hash, config, seed, input checksums, runtime, n_rows,
    #          fallbacks, status ok/failed + exception text)
run_id = f"{name}__{YYYYmmdd-HHMMSS}__s{seed}" unless given explicitly. Never overwrites an existing run directory.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

import pandas as pd

from ksp.paths import REPO_ROOT, ROOT  # noqa: E402  ROOT = project data/working dir ($KSP_DATA_DIR)
RESULTS = ROOT / "results" / "runs"
CACHE = ROOT / "cache" / "runs"


def _rel(p: Path) -> str:
    try:
        return str(Path(p).relative_to(ROOT))
    except ValueError:
        return str(p)


def git_hash() -> str:
    """Read-only `git rev-parse HEAD` (+ '-dirty' flag is NOT computed: agents never run git status)."""
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True,
                              timeout=20).stdout.strip() or "unknown"
    except Exception as e:  # recorded, not hidden
        return f"unavailable:{type(e).__name__}"


def sha256_file(p, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _checksum(p) -> dict:
    p = Path(p)
    if p.is_dir():
        files = sorted(x for x in p.rglob("*") if x.is_file())
        h = hashlib.sha256()
        for f in files:
            h.update(str(f.relative_to(p)).encode())
            h.update(sha256_file(f).encode())
        return {"path": str(p), "type": "dir", "n_files": len(files), "sha256": h.hexdigest()}
    if not p.exists():
        raise FileNotFoundError(f"manifest input {p} does not exist")
    return {"path": str(p), "type": "file", "bytes": p.stat().st_size, "sha256": sha256_file(p)}


class RunManifest:
    def __init__(self, name: str, config: dict, seed: int, inputs=(), run_id: str | None = None,
                 split: str | None = None, notes: str = ""):
        if "/" in name or " " in name:
            raise ValueError("run name must not contain '/' or spaces")
        self.run_id = run_id or f"{name}__{time.strftime('%Y%m%d-%H%M%S')}__s{seed}"
        self.dir = RESULTS / self.run_id
        self.cache = CACHE / self.run_id
        if self.dir.exists():
            raise FileExistsError(f"run directory {self.dir} exists; runs are never overwritten")
        self.m = {"run_id": self.run_id, "name": name, "git_hash": git_hash(), "config": config, "seed": seed,
                  "split": split, "notes": notes, "inputs": [_checksum(p) for p in inputs],
                  "host": socket.gethostname(), "python": sys.version.split()[0], "platform": platform.platform(),
                  "argv": sys.argv, "pid": os.getpid(), "started": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                  "n_rows": {}, "outputs": {}, "fallbacks": {}}
        self._fallbacks: Counter = Counter()
        self._t0 = time.time()

    # --------------------------------------------------------------------------------------------- recording
    def count_fallback(self, kind: str, n: int = 1) -> None:
        self._fallbacks[kind] += int(n)

    def set(self, key: str, value) -> None:
        self.m[key] = value

    def save_predictions(self, preds: pd.DataFrame, name: str = "predictions") -> Path:
        need = {"score"} | ({"kinase_acc", "site_id"} if "kinase_acc" in preds else {"kinase_idx", "site_idx"})
        if not need <= set(preds.columns):
            raise ValueError(f"predictions need columns {sorted(need)}")
        self.cache.mkdir(parents=True, exist_ok=True)
        f = self.cache / f"{name}.parquet"
        preds.to_parquet(f, index=False)
        self.m["n_rows"][name] = int(len(preds))
        self.m["outputs"][name] = {"path": _rel(f), "sha256": sha256_file(f)}
        return f

    def save_metrics(self, metrics, name: str = "metrics") -> Path:
        self.dir.mkdir(parents=True, exist_ok=True)
        f = self.dir / f"{name}.json"
        if hasattr(metrics, "save_json"):
            metrics.save_json(f)
        else:
            from .eval import _json_default
            f.write_text(json.dumps(metrics, indent=1, default=_json_default))
        self.m["outputs"][name] = {"path": _rel(f), "sha256": sha256_file(f)}
        return f

    def write(self, status: str = "ok", error: str | None = None) -> Path:
        self.dir.mkdir(parents=True, exist_ok=True)
        self.m.update({"status": status, "error": error, "runtime_s": round(time.time() - self._t0, 2),
                       "finished": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "fallbacks": dict(self._fallbacks)})
        f = self.dir / "manifest.json"
        from .eval import _json_default
        f.write_text(json.dumps(self.m, indent=1, default=_json_default))
        return f

    def __enter__(self):
        return self

    def __exit__(self, et, ev, tb):
        if et is None:
            self.write("ok")
        else:
            self.write("failed", "".join(traceback.format_exception(et, ev, tb))[-4000:])
        return False
