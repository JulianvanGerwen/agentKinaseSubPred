#!/usr/bin/env python
"""Re-evaluate every pooled-CV run from its saved predictions with the current ksp.eval (DESIGN §17 A24) and build
results/cv/LEADERBOARD.{md,tsv}. Owner: core-infra.

Discovery: results/runs/<run_id>/manifest.json whose split is CV-like ('cv', 'cv_pooled', or containing 'lko'),
with cache/runs/<run_id>/predictions.parquet, and that is POOLED (no metrics*.json carries an integer 'fold').
Single-fold runs and per-LKO-fold runs (name '..._lko<k>') are skipped (A24 does not change them).
Evaluation: missing='worst' (unscored required pairs ranked last and counted -> coverage column); A24 fold
rank-normalisation ON, except leave-kinase-out pooled runs (one model per kinase): fold_normalise=False.
Writes ONLY results/runs/<run_id>/metrics_a24.json in other agents' run dirs (old metrics kept).
Usage: python scripts/reeval_cv_runs.py [--force] [--only REGEX] [--list] [--board-only]
"""
import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_MAX_THREADS"):
    os.environ.setdefault(v, "8")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp import eval as E  # noqa: E402
from ksp import canary, index, labels, splits  # noqa: E402

RUNS = ROOT / "results" / "runs"
PRED = ROOT / "cache" / "runs"
OUT = ROOT / "results" / "cv"
OUTFILE = "metrics_a24.json"
BINS = ["0", "1-4", "5-19", "20-99", ">=100"]


def family(name: str) -> str:
    for pre, fam in [("simple_b0", "B0 floor"), ("simple_b1", "B1 Cantley"), ("simple_b2", "B2 PSSM"),
                     ("simple_b7", "B7 simple-strong"), ("r2_", "R2 LightGBM"), ("r3_", "R3 neural"), ("r5_", "R5 late-interaction"),
                     ("b3", "B3 Phosformer"), ("b4", "B4 SELPHI2"), ("b5", "B5 NetworKIN"), ("b6", "B6 KolossuS")]:
        if name.startswith(pre):
            return fam
    return "other"


def discover() -> list[dict]:
    out = []
    for m in sorted(RUNS.glob("*/manifest.json")):
        rid = m.parent.name
        try:
            man = json.loads(m.read_text())
        except json.JSONDecodeError:
            continue
        split = str(man.get("split") or "")
        if not (split.startswith("cv") or "lko" in split):
            continue
        pred = PRED / rid / "predictions.parquet"
        if not pred.exists():
            continue
        if re.search(r"_lko\d+$", man.get("name", "")) or re.search(r"lko\d+__", rid):
            continue
        folds = []
        for f in m.parent.glob("metrics*.json"):
            if f.name == OUTFILE:
                continue
            try:
                folds.append(json.loads(f.read_text()).get("fold"))
            except (json.JSONDecodeError, AttributeError):
                pass
        if any(isinstance(x, int) for x in folds):
            continue
        is_lko = "lko" in split or "lko" in rid
        out.append({"run_id": rid, "name": man.get("name", rid), "split": split, "lko": is_lko, "pred": pred,
                    "manifest": man})
    return out


def config_summary(man: dict) -> str:
    c = man.get("config", {})
    if not isinstance(c, dict):
        return str(c)[:120]
    keep = []
    for k, v in c.items():
        if k in ("grid", "fold_logs", "source_runs", "best_params_per_fold") or isinstance(v, (dict, list)) and \
                len(str(v)) > 60:
            continue
        keep.append(f"{k}={v}")
    return "; ".join(keep)[:160]


def reeval(run: dict, data_cv, canary_v2_time: float) -> dict:
    t0 = time.time()
    p = pd.read_parquet(run["pred"])
    cols = [c for c in ["kinase_acc", "site_id", "kinase_idx", "site_idx", "score", "score_site"] if c in p]
    p = p[cols]
    res = E.evaluate(p, data_cv, missing="worst", fold_normalise=not run["lko"], keep_items=False)
    j = res.to_json()
    st = E.strata(res, kinase_set="K_any")
    pk = res.per_kinase
    keep = [c for c in ["kinase_acc", "gene", "n_pos", "n_pos_required", "n_neg_hard", "n_neg_all", "n_missing",
                        "auroc_hard", "auroc_all", "top1_recall", "hit5", "in_K_primary", "in_K_any", "train_bin",
                        "residue_class"] if c in pk]
    j["per_kinase"] = pk[keep].to_dict(orient="records")
    if "kinase_acc" in p:
        can_scored = bool(p.kinase_acc.astype(str).str.startswith("CANARY_").any())
    else:
        can_scored = bool((p.kinase_idx >= len(E.index.load_kinases())).any())
    j["a24"] = {"source_run": run["run_id"], "canary_scored": can_scored, "fold_normalise": not run["lko"], "missing_policy": "worst",
                "evaluated": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "eval_seconds": round(time.time() - t0, 1),
                "n_pred_rows": int(len(p)),
                "canary_labels_version": 2 if _started(run["manifest"]) >= canary_v2_time else 1,
                "strata_K_any": st.to_dict(orient="records")}
    (RUNS / run["run_id"] / OUTFILE).write_text(json.dumps(j, indent=1, default=E._json_default))
    return j


def _started(man: dict) -> float:
    s = man.get("started", "")
    try:
        return time.mktime(time.strptime(s[:19], "%Y-%m-%d %H:%M:%S"))
    except ValueError:
        return 0.0


_UNI = {}


def eall_check(run_id: str) -> dict:
    """Within-kinase E-all coverage: for each kinase with predictions, scored fraction of DEV U_ms sites of its
    residue class. E-all metrics (M3) are valid if <= 5% of scored kinases are partially covered (0-95%): unscored
    whole kinases (e.g. Cantley non-atlas) are fine, label-selected subsets (positives + E-hard only) are not."""
    if not _UNI:
        S = index.load_sites()
        dev = ~S.acc.isin(splits.load_test_proteins()).values & S.in_U_ms.values
        K = canary.kinases_with_canaries()
        _UNI["dev"] = dev
        _UNI["res"] = S.res.values
        _UNI["n_cls"] = {c: int((dev & E._res_ok(S.res.values, c)).sum()) for c in ["ST", "Y", "STY"]}
        _UNI["cls"] = K.residue_class.values
        _UNI["kpos"] = {a: i for i, a in enumerate(K.uniprot_acc)}
        _UNI["spos"] = pd.Series(np.arange(len(S)), index=S.site_id)
    pq = PRED / run_id / "predictions.parquet"
    cols = pd.read_parquet(pq, columns=None).columns if False else None  # noqa: F841
    try:
        p = pd.read_parquet(pq, columns=["kinase_idx", "site_idx"])
        ki, si = p.kinase_idx.values, p.site_idx.values
    except Exception:
        p = pd.read_parquet(pq, columns=["kinase_acc", "site_id"])
        ki = p.kinase_acc.map(_UNI["kpos"]).values
        si = _UNI["spos"][p.site_id].values
    m = _UNI["dev"][si]
    cnt = np.bincount(ki[m].astype(np.int64), minlength=len(_UNI["cls"]))
    need = np.array([_UNI["n_cls"][c] for c in _UNI["cls"]])
    sc = cnt > 0
    frac = cnt[sc] / need[sc]
    part = int(((frac > 0) & (frac < 0.95)).sum())
    return {"n_kinases_scored": int(sc.sum()), "n_partial_kinases": part,
            "eall_valid": bool(part <= 0.05 * max(sc.sum(), 1))}


def board(canary_v2_time: float) -> pd.DataFrame:
    rows = []
    for f in sorted(RUNS.glob(f"*/{OUTFILE}")):
        j = json.loads(f.read_text())
        rid = f.parent.name
        man = json.loads((f.parent / "manifest.json").read_text())
        e, a, c = j["endpoints"], j.get("a24", {}), j.get("canary", {})
        st = pd.DataFrame(a.get("strata_K_any", []))
        r = {"run_id": rid, "family": family(man.get("name", rid)), "config": config_summary(man),
             "split": man.get("split"), "M1": e.get("M1_auroc_hard"), "M1_fold_mean": e.get("M1_fold_mean"),
             "M2_hit5": e.get("M2_hit5"), "M3_top1": e.get("M3_top1_recall"), "M4_Kany": e.get("M4_auroc_hard_Kany"),
             "M5_Kcv": e.get("M5_auroc_hard_Kcv")}
        for b in BINS:
            x = st[(st.stratum == "train_bin") & (st.value == b)] if len(st) else st
            r[f"M5bin_{b}"] = float(x.auroc_hard.iloc[0]) if len(x) else np.nan
            r[f"n_bin_{b}"] = int(x.n_kinases.iloc[0]) if len(x) else 0
        for cl in ["ST", "Y"]:
            x = st[(st.stratum == "residue_class") & (st.value == cl)] if len(st) else st
            r[f"M1_{cl}"] = float(x.auroc_hard.iloc[0]) if len(x) else np.nan
        if c and a.get("canary_scored", True):
            v = a.get("canary_labels_version", 1)
            r["canary_neg_Ehard"] = c.get("neg_mean_auroc_hard")
            r["canary_pos"] = c.get("pos_auroc_hard")
            r["canary_gate"] = ("PASS" if c.get("neg_pass") and c.get("pos_pass") else "FAIL") + \
                               ("" if v == 2 else " (v1-trained)")
        else:
            r["canary_neg_Ehard"] = r["canary_pos"] = np.nan
            r["canary_gate"] = "canaries not scored"
        req, mis = e.get("n_required_pairs", 0), e.get("n_missing_pairs", 0)
        r["coverage"] = 1 - mis / req if req else np.nan
        # E-all validity by within-kinase coverage (cached in metrics_a24.json)
        if "eall_check" not in a and (PRED / rid / "predictions.parquet").exists():
            a["eall_check"] = eall_check(rid)
            j["a24"] = a
            f.write_text(json.dumps(j, indent=1, default=E._json_default))
        ec = a.get("eall_check", {"eall_valid": bool(r["coverage"] >= 0.95), "n_partial_kinases": -1})
        r["eall_valid"] = ec["eall_valid"]
        r["n_partial_kinases"] = ec.get("n_partial_kinases")
        # M3 over K_cv kinases whose required pairs are all scored (matched-coverage view; needs per_kinase)
        pkj = pd.DataFrame(j.get("per_kinase", []))
        if len(pkj) and "n_missing" in pkj:
            kp = pkj[pkj.in_K_primary.astype(bool)]
            fc = kp[kp.n_missing == 0]
            r["M3_fullcov"] = float(fc.top1_recall.mean()) if len(fc) else np.nan
            r["n_fullcov"] = int(len(fc))
            r["eall_valid"] = bool(len(fc) == len(kp))       # every K_cv kinase fully scored -> M3 valid as is
        else:
            r["M3_fullcov"], r["n_fullcov"] = np.nan, np.nan
        if not r["eall_valid"]:
            r["M3_top1"] = np.nan
        cfg = man.get("config", {})
        r["contaminated"] = bool(isinstance(cfg, dict) and cfg.get("contaminated"))
        r["runtime_s"] = man.get("runtime_s")
        r["started"] = man.get("started")
        r["fold_normalised"] = e.get("fold_normalised")
        rows.append(r)
    B = pd.DataFrame(rows)
    if not len(B):
        return B
    B["official"] = ""

    def mark(pattern, label, exclude=r"canary"):
        x = B[B.run_id.str.contains(pattern, regex=True) & ~B.run_id.str.contains(exclude, regex=True)]
        if len(x):
            i = x.sort_values("started").index[-1]
            B.loc[i, "official"] = (B.loc[i, "official"] + "; " if B.loc[i, "official"] else "") + label
    B.loc[B.contaminated, "official"] = "CONTAMINATED (trained on PSP; contrast only)"
    mark(r"^simple_b1_published(_reeval)?__", "B1 official M1/M3 (A23 arm i)")
    mark(r"^simple_b1_universe(_reeval)?__", "B1 official M2 (A23 arm ii)")
    mark(r"^simple_b7_lr(_reeval)?__", "B7 official")
    for b in ["prior", "pop", "miss"]:
        mark(rf"^simple_b0_{b}(_reeval)?__", f"B0-{b} floor")
    return B.sort_values("M1", ascending=False, na_position="last").reset_index(drop=True)


def write_board(B: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    B.to_csv(OUT / "LEADERBOARD.tsv", sep="\t", index=False, float_format="%.4f")
    f = lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"  # noqa: E731
    L = ["# CV leaderboard (pooled OOF, re-evaluated with ksp.eval A24)", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M:%S %Z')} by `scripts/reeval_cv_runs.py` from saved predictions "
         "(metrics_a24.json per run). Sorted by M1 (macro E-hard AUROC over K_cv, fold rank-normalised; LKO pooled "
         "runs not normalised). M5bin = macro E-hard AUROC over K_any by DEV training-label bin; M1_ST/M1_Y = "
         "K_any macro by residue class. Coverage = fraction of required (kinase, U_ms site) pairs scored "
         "(missing = ranked last). M3 is blanked unless every K_cv kinase is fully scored (runs with per-kinase data) "
         "or, for older metrics files, when any scored kinase has partial E-all coverage "
         "(e.g. only positives + E-hard scored): E-all metrics are then not comparable (use ksp.eval.matched_universe); "
         "unscored whole kinases do not invalidate M3. 'M3 full-cov (n)' = M3 over the n K_cv kinases with every "
         "required pair scored (matched-coverage view; available for runs evaluated after 17:40). Canary gate (v2): mean negative E-hard <= 0.57 and positive >= 0.90; "
         "'(v1-trained)' = run trained before the v2 canary labels (gate not meaningful). Official baselines per "
         "A23 (latest run of each arm).", "",
         "| # | run | family | official | M1 | M1 fold-mean | M2 Hit@5 | M3 top1% | M3 full-cov (n) | M4 | M5 0 | 1-4 | 5-19 | 20-99 "
         "| >=100 | M1 ST | M1 Y | canary neg / pos | gate | coverage | runtime s | started |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in B.iterrows():
        L.append(f"| {i + 1} | {r.run_id} | {r.family} | {r.official} | {f(r.M1)} | {f(r.M1_fold_mean)} | "
                 f"{f(r.M2_hit5)} | {f(r.M3_top1)} | "
                 f"{f(r.M3_fullcov)}{'' if pd.isna(r.n_fullcov) else f' ({int(r.n_fullcov)})'} | {f(r.M4_Kany)} | " +
                 " | ".join(f(r[f'M5bin_{b}']) for b in BINS) +
                 f" | {f(r.M1_ST)} | {f(r.M1_Y)} | {f(r.canary_neg_Ehard)} / {f(r.canary_pos)} | {r.canary_gate} | "
                 f"{f(r.coverage)} | {r.runtime_s if r.runtime_s is not None else ''} | {r.started} |")
    L += ["", "Config summaries: see LEADERBOARD.tsv (column `config`)."]
    (OUT / "LEADERBOARD.md").write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--only", default=None, help="regex on run_id")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--board-only", action="store_true")
    ap.add_argument("--shard", default=None, help="i/n: process only every n-th run of the todo list (parallel workers)")
    a = ap.parse_args()
    cfile = labels.DEV_PAIRS.parent / "canary_dev.parquet"
    canary_v2_time = cfile.stat().st_mtime if cfile.exists() else 0.0
    runs = discover()
    if a.only:
        runs = [r for r in runs if re.search(a.only, r["run_id"])]
    todo = [r for r in runs if a.force or not (RUNS / r["run_id"] / OUTFILE).exists()]
    if a.shard:
        i, n = map(int, a.shard.split("/"))
        todo = todo[i::n]
    print(f"{len(runs)} pooled CV runs with predictions; {len(todo)} to (re-)evaluate", flush=True)
    if a.list:
        for r in runs:
            print(("TODO " if r in todo else "done ") + r["run_id"] + ("  [LKO: no fold-normalisation]" if r["lko"]
                                                                        else ""))
        return
    if todo and not a.board_only:
        data = E.prepare("cv", fold=None)
        for r in todo:
            t = time.time()
            try:
                j = reeval(r, data, canary_v2_time)
                print(f"{r['run_id']}: M1 {j['endpoints']['M1_auroc_hard']:.4f} ({time.time() - t:.0f}s)", flush=True)
            except Exception as ex:   # reported, the run is skipped; never silently
                print(f"FAILED {r['run_id']}: {type(ex).__name__}: {ex}", flush=True)
    write_board(board(canary_v2_time))
    print("wrote", OUT / "LEADERBOARD.md")


if __name__ == "__main__":
    main()
