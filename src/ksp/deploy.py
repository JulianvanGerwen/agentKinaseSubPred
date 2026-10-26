"""Deployment labels (DESIGN §8.6, G5): the frozen configuration retrained on DEV + TEST for the RELEASED predictions.
Owner: core-infra. A NEW module so that ksp.labels stays frozen (A28).

    from ksp.deploy import load_deploy_labels
    lab = load_deploy_labels("L2")                 # DEV + TEST labels, same schema/semantics as labels.load_dev
    allrows = load_deploy_labels("all")            # every row incl. mask-only rows -> KP / NegativeMasks

Refuses (DeployError) unless
  * FREEZE.json is git-tracked and unmodified (read-only git), and
  * results/test/LEDGER.md holds a COMPLETED TEST block (ksp.eval_test ledger markers) for exactly that FREEZE sha -
    i.e. the single TEST batch is over before any deployment training reads TEST labels;
  * the sealed files match results/splits/SEAL.json.
Semantics are those of ksp.labels.load_dev (variant selection, U_ms restriction by default except variant='all',
is_canary column), applied to DEV rows + sealed TEST rows; DEV and TEST must be protein-disjoint (asserted).

Canaries: the DEFAULT is canary-free (with_canaries=False) - released predictions must never contain canary
pseudo-kinases, and there is no held-out data on which a canary gate could be evaluated after DEV+TEST retraining.
with_canaries=True appends the DEV canary labels (ksp.canary, v2) plus the sealed TEST canary rows (positive canary
only), for pipelines whose code path expects canary kinases to be present; their outputs for CANARY_* kinases must
be dropped before release.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import eval_test as ET
from . import labels


class DeployError(RuntimeError):
    pass


def check_test_completed(freeze: Path = ET.FREEZE, ledger: Path = ET.LEDGER, require_git: bool = True) -> str:
    """Return the FREEZE sha if a COMPLETED TEST batch exists for the committed FREEZE.json, else raise."""
    if not Path(freeze).exists():
        raise DeployError(f"{freeze} does not exist (no FREEZE -> no TEST batch -> no deployment labels)")
    if require_git:
        r = ET._git("ls-files", "--error-unmatch", Path(freeze).name, cwd=Path(freeze).parent)
        if r.returncode != 0:
            raise DeployError("FREEZE.json is not committed")
        r = ET._git("diff", "--quiet", "HEAD", "--", Path(freeze).name, cwd=Path(freeze).parent)
        if r.returncode != 0:
            raise DeployError("FREEZE.json has uncommitted changes")
    fsha = ET._sha(Path(freeze))
    done = [e for e in ET.ledger_events(Path(ledger)) if e.get("event") == "COMPLETED" and e.get("tag") == "TEST"
            and e.get("freeze_sha") == fsha]
    if not done:
        raise DeployError(f"no COMPLETED TEST batch for FREEZE sha {fsha[:12]} in {ledger}: the TEST batch must be "
                          f"over before deployment training may read TEST labels")
    return fsha


def load_deploy_labels(variant: str = "L2", with_canaries: bool = False, ums_only: bool | None = None,
                       freeze: Path = ET.FREEZE, ledger: Path = ET.LEDGER, seal: Path = ET.SEAL,
                       sealed_dir: Path = ET.SEALED_DIR, require_git: bool = True) -> pd.DataFrame:
    """DEV + TEST labels for the deployment model; same schema/semantics as ksp.labels.load_dev()."""
    check_test_completed(freeze, ledger, require_git)
    try:
        ET.check_seal(seal, sealed_dir)
    except ET.FreezeError as e:
        raise DeployError(str(e)) from e
    dev = pd.read_parquet(labels.DEV_PAIRS)
    tst = pd.read_parquet(Path(sealed_dir) / "test_labels.parquet")
    tcan = tst.get("is_canary", pd.Series(False, index=tst.index)).fillna(False).astype(bool)
    t_real, t_can = tst[~tcan], tst[tcan]
    shared = set(dev.acc) & set(t_real.acc)
    if shared:
        raise DeployError(f"DEV and TEST label tables share {len(shared)} proteins (split invariant broken)")
    d = pd.concat([dev, t_real], ignore_index=True)
    d["is_canary"] = False
    d = labels.select_variant(d, variant)
    if ums_only is None:
        ums_only = variant != "all"
    if ums_only:
        d = d[d.in_U_ms.astype(bool)].reset_index(drop=True)
    if with_canaries:
        from . import canary
        c = pd.concat([canary.canary_labels(), t_can], ignore_index=True)
        c["is_canary"] = True
        if variant != "all":
            c = labels.select_variant(c, variant) if variant != "eval" else c[labels.evaluable_positive(c)]
        d = pd.concat([d, c], ignore_index=True)
    return d
