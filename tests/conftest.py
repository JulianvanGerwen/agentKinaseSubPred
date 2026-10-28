import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ksp.paths import ROOT as DATA_ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

# Inputs that the data-dependent tests read. None of them is redistributed (licensed / derived from PhosphoSitePlus);
# see DATA_SOURCES.md. Paths are relative to KSP_DATA_DIR.
INDEX = ("cache/index/kinases.parquet", "cache/index/sites.parquet")
DEV = ("cache/labels/dev_pairs.parquet",)
CANARY = ("cache/labels/canary_dev.parquet",)
SPLITS = ("results/splits/SEAL.json",)


def missing_data_reason(paths) -> str:
    return (f"needs project data that is not part of this repository (missing under KSP_DATA_DIR={DATA_ROOT}: "
            f"{', '.join(paths)}); see DATA_SOURCES.md")


def pytest_configure(config):
    config.addinivalue_line("markers", "requires_data(*paths): test needs project data files (paths relative to "
                                       "KSP_DATA_DIR); it is skipped when any of them is missing")


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    for m in item.iter_markers("requires_data"):
        missing = [p for p in m.args if not (DATA_ROOT / p).exists()]
        if missing:
            pytest.skip(missing_data_reason(missing))
