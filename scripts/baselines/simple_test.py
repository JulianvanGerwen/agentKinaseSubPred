"""TEST-batch deliverables for the simple baselines (docs/G4_RUNBOOK.md). Never reads the sealed TEST labels.

    OMP_NUM_THREADS=8 .venv/bin/python -I scripts/baselines/simple_test.py --method b1_arm_i [--t3]

Each method: final model trained on ALL DEV (label-derived ones), scores the exhaustive TEST candidate set
(ksp.eval_test.test_candidate_pairs, 25,922,452 pairs), all-or-nothing per kinase; writes
  cache/runs/<run_id>/test_predictions.parquet (+ test_t3_predictions.parquet with --t3, B2/B7 only)
  results/runs/<run_id>/manifest.json (train_wall_s, inference_wall_s, ...)
  results/baselines/simple/freeze_entries/<method>.json  (FREEZE entry, SHA-256 pinned via eval_test.pin_files,
                                                           + the --check-predictions pre-flight result)
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp import eval as kev  # noqa: E402
from ksp import eval_test as T  # noqa: E402
from ksp import manifest  # noqa: E402
from ksp.baselines import simple as B  # noqa: E402

OUT = ROOT / "results" / "baselines" / "simple"
FE = OUT / "freeze_entries"
SEED = 0
CV_B7_RUN = "simple_b7_lr__20261026-153330__s0"     # CV OOF run used to fit the A17 softmax temperature
B7_CANTLEY = "motif/cantley_unprimed_pct_published"  # B7's Cantley feature = CV best arm (A19 rule; = CV B7 run)

SPEC = {  # method -> (FREEZE name, role, baseline_id, description)
    "b1_arm_i": ("B1_i", "baseline", "B1_arm_i", "Cantley arm (i) as-published percentile (own matrices only)"),
    "b1_arm_ii": ("B1_ii", "baseline", "B1_arm_ii", "Cantley arm (ii) universe-recalibrated percentile"),
    "b2_pssm": ("B2", "baseline", "B2", "learned PSSM (F2a, family prior) on all DEV L2"),
    "b7_lr": ("B7", "baseline", "B7", "simple strong LR on Cantley(i)+B2+PPI+loc, all DEV"),
    "b0_prior": ("B0_prior", "control", "B0_prior", "kinase DEV-label count; score_site = normalised prior"),
    "b0_pop": ("B0_pop", "control", "B0_pop", "site popularity (diagnostic floor)"),
    "b0_miss": ("B0_miss", "control", "B0_miss", "LightGBM on missingness indicators only, all DEV"),
}


def kinase_lko_folds() -> pd.DataFrame:
    return pd.read_csv(ROOT / "results" / "splits" / "kinase_folds.tsv", sep="\t", comment="#")


def run_main(method, cand, data, log):
    """-> dict(score, [score_site], train_s, infer_s)"""
    if method in ("b1_arm_i", "b1_arm_ii"):
        st = {"b1_arm_i": "motif/cantley_unprimed_pct_published",
              "b1_arm_ii": "motif/cantley_unprimed_pct_universe"}[method]
        t0 = time.time()
        v = B.cantley_test(cand, st, log)
        log["store"] = st
        return {"score": v, "train_s": 0.0, "infer_s": time.time() - t0}
    if method == "b2_pssm":
        pct, m, ts, ti = B.b2_test(cand, data, seed=SEED)
        log.update(alpha=m.alpha, beta=m.beta, n_pairs=m.info.get("n_pairs_used"))
        return {"score": pct, "train_s": ts, "infer_s": ti}
    if method == "b7_lr":
        r = B.b7_test(cand, data, B7_CANTLEY, seed=SEED, log=log)
        t0 = time.time()
        temp = B.fit_softmax_temperature(CV_B7_RUN, data)
        t1 = time.time()
        p = B.site_softmax(cand, r["logit"], temp["T"])
        log["softmax_temperature"] = temp
        return {"score": r["lin"], "score_site": p, "train_s": r["train_s"] + (t1 - t0),
                "infer_s": r["infer_s"] + (time.time() - t1)}
    if method == "b0_prior":
        t0 = time.time()
        n = B.fold_train_counts(None, data.kinases.iloc[:B.n_real()])
        nk = n[cand.kinase_idx.values]
        tot = pd.Series(nk).groupby(cand.site_idx.values).transform("sum").values
        log["n_dev_L2_total"] = int(n.sum())
        return {"score": nk.astype(float), "score_site": np.where(tot > 0, nk / np.maximum(tot, 1), 0.0),
                "train_s": 0.0, "infer_s": time.time() - t0}
    if method == "b0_pop":
        t0 = time.time()
        return {"score": B.b0_pop(cand, data), "train_s": 0.0, "infer_s": time.time() - t0}
    if method == "b0_miss":
        v, ts, ti = B.b0_miss_test(cand, data, seed=SEED, log=log)
        return {"score": v, "train_s": ts, "infer_s": ti}
    raise ValueError(method)


def run_t3(method, cand, data, log):
    """Each kinase scored by the model trained without its LKO fold's kinases (DEV labels)."""
    kf = kinase_lko_folds()
    acc = data.kinases.uniprot_acc.values
    lko = pd.Series(kf.lko_fold.values, index=kf.kinase_acc.values)
    kfold = lko.reindex(acc[:B.n_real()]).values
    if np.isnan(kfold.astype(float)).any():
        raise KeyError("kinase without lko_fold")
    score = np.full(len(cand), np.nan)
    site = np.full(len(cand), np.nan) if method == "b7_lr" else None
    ts = ti = 0.0
    for f in sorted(lko.unique()):
        held = set(kf.kinase_acc[kf.lko_fold == f])
        rows = np.where(kfold[cand.kinase_idx.values] == f)[0]
        sub = cand.iloc[rows].reset_index(drop=True)
        lg = {}
        if method == "b2_pssm":
            pct, m, a, b = B.b2_test(sub, data, exclude_kinases=held, seed=SEED)
            lg.update(alpha=m.alpha, beta=m.beta, levels_held=pd.Series(
                np.where(data.sites.res.values[sub.site_idx.values] == "Y", m.level["Y"][sub.kinase_idx.values],
                         m.level["ST"][sub.kinase_idx.values])).value_counts().to_dict())
            score[rows] = pct
        elif method == "b7_lr":
            r = B.b7_test(sub, data, B7_CANTLEY, exclude_kinases=held, seed=SEED, log=lg)
            a, b = r["train_s"], r["infer_s"]
            score[rows] = r["lin"]
            site[rows] = r["logit"]       # softmax over ALL kinases of a site needs all folds: applied after the loop
        else:
            raise ValueError(f"no T3 variant for {method}")
        ts += a
        ti += b
        lg.update(n_kinases_held=len(held), n_rows=len(rows), train_s=round(a, 1), infer_s=round(b, 1))
        log[f"lko_fold{f}"] = lg
    out = {"score": score, "train_s": ts, "infer_s": ti}
    if site is not None:
        temp = B.fit_softmax_temperature(CV_B7_RUN, data)
        out["score_site"] = B.site_softmax(cand, site, temp["T"])
    return out


def frame(cand, r):
    ok = np.isfinite(r["score"])
    p = pd.DataFrame({"kinase_idx": cand.kinase_idx.values[ok].astype(np.int32),
                      "site_idx": cand.site_idx.values[ok].astype(np.int32),
                      "score": r["score"][ok].astype(np.float64)})
    if r.get("score_site") is not None:
        ss = r["score_site"][ok]
        if not np.isfinite(ss).all() or (ss < 0).any() or (ss > 1).any():
            raise ValueError("score_site must be a finite probability in [0, 1]")
        p["score_site"] = ss.astype(np.float64)
    return p, int((~ok).sum())


def preflight(path) -> dict:
    r = subprocess.run([sys.executable, "-I", str(ROOT / "src" / "ksp" / "eval_test.py"), "--check-predictions",
                        str(path)], capture_output=True, text=True, cwd=ROOT)
    txt = r.stdout.strip()
    try:
        j = json.loads(txt[txt.index("{"):])
    except ValueError:
        raise RuntimeError(f"pre-flight output not parseable:\n{txt}\n{r.stderr[-2000:]}")
    if not j.get("ok"):
        raise RuntimeError(f"pre-flight failed: {j}")
    return j


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=sorted(SPEC))
    ap.add_argument("--t3", action="store_true", help="also write test_t3_predictions.parquet (B2/B7 only)")
    a = ap.parse_args()
    FE.mkdir(parents=True, exist_ok=True)
    data = kev.prepare("cv", with_canaries=False)       # DEV labels only (load_dev); TEST labels never touched
    cand = B.test_candidates()
    name, role, bid, desc = SPEC[a.method]
    inputs = [ROOT / "cache" / "labels" / "dev_pairs.parquet", ROOT / "results" / "splits" / "test_proteins.tsv",
              ROOT / "results" / "splits" / "kinase_folds.tsv"]
    with manifest.RunManifest(f"simple_test_{a.method}", config={"method": a.method, "description": desc,
                                                                  "trained_on": "all DEV (L2)", "seed": SEED},
                              seed=SEED, inputs=inputs, split="test") as run:
        log: dict = {}
        r = run_main(a.method, cand, data, log)
        p, n_un = frame(cand, r)
        run.count_fallback("test_candidate_pairs_unscored", n_un)
        f = run.save_predictions(p, "test_predictions")
        del p
        run.set("train_wall_s", round(r["train_s"], 2))
        run.set("inference_wall_s", round(r["infer_s"], 2))
        run.set("log", log)
        pf = preflight(f)
        run.set("preflight", pf)
        entry = {"name": name, "role": role, "baseline_id": bid, "run_id": run.run_id,
                 "predictions": str(f.relative_to(ROOT)), "score_site_predictions": None, "t3_predictions": None,
                 "runtime": {"train_wall_s": round(r["train_s"], 2), "inference_wall_s": round(r["infer_s"], 2)},
                 "missing": "worst", "description": desc}
        if a.t3:
            lg3: dict = {}
            r3 = run_t3(a.method, cand, data, lg3)
            p3, n3 = frame(cand, r3)
            f3 = run.save_predictions(p3, "test_t3_predictions")
            del p3
            pf3 = preflight(f3)
            run.set("t3", {"log": lg3, "train_wall_s": round(r3["train_s"], 2),
                           "inference_wall_s": round(r3["infer_s"], 2), "unscored_pairs": n3, "preflight": pf3})
            entry["t3_predictions"] = str(f3.relative_to(ROOT))
            entry["t3_runtime"] = {"train_wall_s": round(r3["train_s"], 2),
                                   "inference_wall_s": round(r3["infer_s"], 2)}
        cfg = T.pin_files({"configs": [entry]})
        entry = cfg["configs"][0]
        run.set("freeze_entry", entry)
        out = {"freeze_entry": entry, "preflight": pf, "log": log,
               **({"t3_preflight": pf3} if a.t3 else {})}
        (FE / f"{a.method}.json").write_text(json.dumps(out, indent=1, default=kev._json_default))
        print(json.dumps({"method": a.method, "run_id": run.run_id, "preflight": pf,
                          "runtime": entry["runtime"]}, default=kev._json_default))


if __name__ == "__main__":
    main()
