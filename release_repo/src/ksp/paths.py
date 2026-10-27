"""Central path configuration for the ksp package and scripts.

Two roots:

- ``REPO_ROOT``: the code checkout (``src/``, ``scripts/``, ``config/``, ``tests/``). Derived from this file's location,
  so it is correct for a source checkout or an editable install (``pip install -e .``).
- ``ROOT``: the project working directory. It holds the inputs that you must obtain yourself (``data/``,
  ``external/``) and everything that the pipeline writes (``cache/``, ``results/``, ``logs/``, the frozen split
  files and ``FREEZE.json``). Set it with the environment variable ``KSP_DATA_DIR``. If the variable is unset,
  ``ROOT`` equals ``REPO_ROOT``; this is the layout in which the study was run.

Only the location of files depends on these values; no computation does.

Other optional environment variables:

- ``KSP_MMSEQS``: path to the MMseqs2 binary (default: ``mmseqs`` on ``PATH``).
- ``KSP_PYTHON``: interpreter for scripts that start sub-processes (default: the running interpreter).
- ``KSP_DATA`` / ``KSP_UNIPROT``: override only the provided-data or UniProt directory (see ``io_psp``, ``uniprot``).
"""
import os
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("KSP_DATA_DIR") or REPO_ROOT).expanduser().resolve()
MMSEQS = os.environ.get("KSP_MMSEQS") or shutil.which("mmseqs") or "mmseqs"
PYTHON = os.environ.get("KSP_PYTHON") or sys.executable


def has_project_data(*rel_paths: str) -> bool:
    """True if every path (relative to ROOT) exists. Used by tests to skip when the licensed data are absent."""
    return all((ROOT / p).exists() for p in rel_paths)
