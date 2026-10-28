"""Test-seal guard (DESIGN §4.3, §16 A1/A3): no code outside an allowlist may reference sealed TEST material or
raw label sources; feature code never loads labels directly (labels enter only through crossfit fit())."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = ["sealed/", "ks_pairs", "Kinase_Substrate_Dataset", "kinase_label_counts", "enzsub", "omnipath",
            "iptmnet"]
ALLOW = {"scripts/features/phospho_study_lists.py",  # PI 2026-10-26: reads only PMIDs of KS evidence to build the A12 exclusion list; no labels enter features
         "src/ksp/io_psp.py", "src/ksp/labels.py", "scripts/build_labels.py", "src/ksp/eval_test.py",
         "scripts/build_splits.py",
         "src/ksp/deploy.py",  # PI 2026-10-27 00:12: G5 deployment labels (DEV+TEST); refuses unless the TEST batch is COMPLETED in LEDGER.md
         "src/ksp/models/regime/variants.py"}  # PI 2026-10-26 20:23: compares provenance values of load_dev('all') (DEV proteins only); no raw-source read          # build_splits.py: writes the seal (core-infra; runs once at G1)
ALLOW_PREFIX = ("scripts/lit/", "scripts/audit_",
                "scripts/publication/")  # v1.0.0: pre-publication scanners match file NAMES of restricted inputs; they read no labels
# Pre-existing files reported to the PI (docs/core_infra.md) for their owners to fix; remove entries once fixed.
KNOWN_VIOLATORS = {"src/ksp/sites.py", "src/ksp/ids.py", "src/ksp/features/ppi.py", "scripts/make_data_rows.py",
                   "scripts/profile_psp.py"}
FEATURE_FORBIDDEN = ["load_dev(", "dev_pairs", "canary_dev", "load_pairs("]


def _files():
    for d in ["src", "scripts"]:
        for p in (ROOT / d).rglob("*"):
            if p.is_file() and p.suffix in (".py", ".sh", ".R", ".ipynb") and "__pycache__" not in p.parts:
                yield p, str(p.relative_to(ROOT))


def violations():
    out = {}
    for p, rel in _files():
        if rel in ALLOW or rel.startswith(ALLOW_PREFIX):
            continue
        txt = p.read_text(errors="ignore")
        hits = [pat for pat in PATTERNS if re.search(re.escape(pat), txt, flags=re.I)]
        if hits:
            out[rel] = hits
    return out


def test_no_sealed_or_raw_label_access():
    v = violations()
    new = {k: h for k, h in v.items() if k not in KNOWN_VIOLATORS}
    assert not new, f"files referencing sealed/raw label sources outside the allowlist: {new}"


def test_sealed_dir_only_in_eval_test():
    for p, rel in _files():
        if rel in ("src/ksp/eval_test.py", "scripts/build_splits.py", "src/ksp/deploy.py"):  # deploy: G5, post-COMPLETED only
            continue
        if rel.startswith("scripts/publication/"):     # v1.0.0: name patterns for .gitignore/publication checks only
            continue
        assert "sealed/" not in p.read_text(errors="ignore"), f"{rel} references sealed/"


def test_feature_code_gets_labels_only_through_fit():
    feat = ROOT / "src" / "ksp" / "features"
    bad = {}
    for p in feat.rglob("*.py"):
        txt = p.read_text(errors="ignore")
        hits = [x for x in FEATURE_FORBIDDEN if x in txt]
        if hits:
            bad[str(p.relative_to(ROOT))] = hits
    assert not bad, f"feature modules must not load labels directly (use crossfit fit()): {bad}"


if __name__ == "__main__":
    import json
    print(json.dumps(violations(), indent=1))
