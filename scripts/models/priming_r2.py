#!/usr/bin/env python
"""DESIGN §17 A32(a): R2 final config (r2_final_s0, seed 0) + label-free kinase-specific priming features.

The R2 pipeline (src/ksp/models/r2) is used read-only: this wrapper subclasses assemble.Assembler in-process and adds
three columns after the standard groups (nothing on disk is modified). Config is identical to r2_final_s0_cv5
(groups F1 F0 F2a F3 F4 F6a F6b F6c F5 F7, f2a_sep, desc_v 2, seed 0) plus `priming_features='v1'`.

Pre-check: F1 ALREADY contains the priming-aware Cantley percentile (f1_primed_hard_pct_universe and
f1_primed_soft_pct_universe, + the matching scores; assemble.py F1 block), so the experiment reduces to the
offset-observed features and the phospho-preference strength:
  fp_offset_obs[k, s]   1 if the residue at s + best_offset(k) is S/T and an observed U_ms phosphosite, else 0;
                        best_offset(k) = position of k's maximal pS/pT log2 enrichment in its own S/T Cantley matrix
                        (results/strata/cantley_phospho_preference.tsv, kinase-library 1.8.0).
  fp_offset_score[k, s] sum over positions p of max(0, pref_p(k)) * [s+p is an S/T U_ms site],
                        pref_p(k) = max(log2 pS, log2 pT) at p.
  fp_strength[k]        max_p pref_p(k) (kinase-level).
Kinases without an own S/T matrix (Tyr kinases, non-atlas kinases) and canaries: NaN (as for every kinase x site
store without canary rows). All inputs are label-free (site index window15 + in_U_ms; Cantley matrices). Flank U_ms
status is the same input class as F1 primed_soft/hard (DESIGN §7 F1); the centre site's own flag is not used (A6).
Usage: R2_THREADS=8 R2_GATHER_THREADS=4 python scripts/models/priming_r2.py --folds 0 1 2 3 4
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ksp import index  # noqa: E402
from ksp import strata as ST  # noqa: E402
from ksp.models.r2 import assemble, train  # noqa: E402

POS = ST._ST_POS


def priming_tables():
    """Per index-kinase pS/pT preference (n_kin x 9; NaN rows without own S/T matrix) and per-offset observed-flank
    indicators (9 x n_sites)."""
    K = index.load_kinases()
    S = index.load_sites()
    M = pd.read_csv(ST.KL_MATRICES / "matrices/ser_thr/ser_thr_all_log2_matrices.tsv", sep="\t", index_col=0)
    info = pd.read_csv(ST.KL_MATRICES / "kinase_data/kinome_information.tsv", sep="\t", dtype=str).set_index(
        "MATRIX_NAME")
    L = np.stack([np.maximum(M[[f"{p}s" for p in POS]].to_numpy(float)[:, i],
                             M[[f"{p}t" for p in POS]].to_numpy(float)[:, i]) for i in range(len(POS))], 1)
    acc = info.UNIPROT_ID.reindex(M.index).values
    kpos = {a: i for i, a in enumerate(K.uniprot_acc)}
    pref = np.full((len(K), len(POS)), np.nan)
    for a, row in zip(acc, L):
        i = kpos.get(a)
        if i is None:
            continue
        pref[i] = row if np.isnan(pref[i]).all() else np.maximum(pref[i], row)
    observed = set(S.site_id[S.in_U_ms.values.astype(bool)])
    obs = np.zeros((len(POS), len(S)), np.float32)
    for j, off in enumerate(POS):
        oid = ST.offset_site_ids(S.site_id.values, S.window15.values, off)
        obs[j] = np.fromiter((o is not None and o in observed for o in oid), bool, len(S))
    return pref, obs


class _Rows:
    """Lazy kinase x site payload for the assembler's 'ks' kind: self[k] -> float32 vector over all sites."""

    def __init__(self, fn):
        self.fn = fn

    def __getitem__(self, k):
        return self.fn(k)


_TAB = {}


def _tables():
    if not _TAB:
        pref, obs = priming_tables()
        _TAB.update(pref=pref, obs=obs, best=np.where(np.isnan(pref).all(1), -1,
                                                      np.nanargmax(np.where(np.isnan(pref), -np.inf, pref), 1)))
    return _TAB


class PrimingAssembler(assemble.Assembler):
    def _build(self):
        super()._build()
        T = _tables()
        pref, obs, best = T["pref"], T["obs"], T["best"]
        nS = obs.shape[1]
        nan = np.full(nS, np.nan, np.float32)

        def f_obs(k):
            return nan if best[k] < 0 else obs[best[k]]

        def f_score(k):
            if best[k] < 0:
                return nan
            w = np.maximum(pref[k], 0.0).astype(np.float32)
            return w @ obs
        self._add("fp_offset_obs", "ks", _Rows(f_obs))
        self._add("fp_offset_score", "ks", _Rows(f_score))
        strength = np.where(best >= 0, np.nanmax(np.where(np.isnan(pref), -np.inf, pref), 1), np.nan)
        self._add("fp_strength", "k", strength[:assemble.N_REAL])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    ap.add_argument("--name", default="priming_r2_s0")
    a = ap.parse_args()
    base = json.loads((ROOT / "results/runs/r2_final_s0_cv5__20261026-181435__s0/manifest.json").read_text())["config"]
    cfg = {**train.default_config(), **base, "priming_features": "v1"}
    print("config:", json.dumps(cfg), flush=True)
    assemble.Assembler = PrimingAssembler          # in-process only; train.run builds assemble.Assembler(...)
    for f in a.folds:
        out = train.run(cfg, folds=[f], name=f"{a.name}_f{f}", config_name=f"{a.name}_f{f}", stage="A32a",
                        notes="A32a priming features (fp_offset_obs, fp_offset_score, fp_strength)",
                        register_row=False)
        print("RESULT", json.dumps({"run_id": out["run_id"], "fold": f, "endpoints": out.get("endpoints"),
                                    "canary": out.get("canary")}, default=str), flush=True)


if __name__ == "__main__":
    main()
