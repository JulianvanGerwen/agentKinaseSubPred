"""Regime driver: one training-label variant x feature mode through the FROZEN final R2 config (v3, seed 0).

Reuses ksp.models.r2 (sample.FoldData, assemble.Assembler, train.train_fold / predict_fold) read-only. Differences
from ksp.models.r2.train.run: (1) the label loader is patched (ksp.models.regime.variants) so FoldData uses the
variant's TRAINING labels; (2) Assembler gets xfit_tag = our refitted store dir when feat='refit'; (3) registry rows
are model='REGIME'. Evaluation is untouched: ksp.eval.prepare('cv') (frozen harness) -> same eval positives / masks.
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from ksp import eval as kev
from ksp import manifest
from ksp.index import ROOT
from ksp.models.regime import variants as V

REGISTRY = ROOT / "results" / "cv" / "REGISTRY.tsv"
OUT = ROOT / "results" / "models" / "regime"
BASE_RUN = "r2_final_s0_cv5__20261026-181435__s0"     # frozen final R2 config v3, seed 0 (matched baseline)


def final_config(seed: int = 0) -> dict:
    m = json.loads((ROOT / "results" / "runs" / BASE_RUN / "manifest.json").read_text())
    cfg = dict(m["config"])
    cfg["seed"] = seed
    return cfg


def run_fold(variant: str, feat: str, fold: int, name: str, seed: int = 0, evaluate_fold: bool = False) -> dict:
    """Train + exhaustive OOF predict one fold; save predictions in a RunManifest. feat in {'L2fixed', 'refit'}."""
    V.patch()
    from ksp.models.r2 import assemble, sample, train
    cfg = final_config(seed)
    cfg["variant"] = variant
    tag = None if feat == "L2fixed" else V.xfit_tag(variant)
    if feat == "refit" and variant == "L2":
        tag = None                                   # L2 refit == R2's own stores (verified identical on fold 0)
    with manifest.RunManifest(name, config={**cfg, "feat": feat, "xfit_tag": tag}, seed=seed, split="cv",
                              notes=f"folds=[{fold}]; regime variant={variant} feat={feat}") as run_:
        t = time.time()
        fd = sample.FoldData(fold, variant)
        A = assemble.Assembler(cfg["groups"], fold, xfit_tag=tag, clean=cfg["clean"], comp_know=cfg["comp_know"],
                               f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
        if tag is not None:          # record which label set the xfit stores were fitted on
            from ksp.models.r2 import xfit
            metas = {fam: xfit.load(fam, fold, tag=tag)[1] for fam in (cfg["f2a_family"], "f3nb", "f6a")}
            bad = {k: m["variant"] for k, m in metas.items() if m["variant"] != variant}
            if bad:
                raise RuntimeError(f"xfit stores fitted on other variants: {bad}")
            run_.set("xfit_n_train_labels", {k: m["n_train_labels"] for k, m in metas.items()})
        lg = {"n_label_pos": int(len(fd.pos)), "n_label_pos_real": int((~fd.pos.is_canary).sum())}
        bst, cols, stats = train.train_fold(fd, A, cfg, lg)
        t3 = time.time()
        P = train.predict_fold(fd, A, bst, cfg)
        lg["t_predict_s"] = round(time.time() - t3, 1)
        lg["t_fold_s"] = round(time.time() - t, 1)
        run_.set("fold_logs", {str(fold): lg})
        run_.save_predictions(P)
        run_.dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"feature": cols, "gain": bst.feature_importance("gain")}).to_csv(
            run_.dir / f"importance_fold{fold}.tsv", sep="\t", index=False)
        if stats is not None:
            stats.to_csv(run_.dir / f"train_stats_fold{fold}.tsv", sep="\t")
        out = {"run_id": run_.run_id, "log": lg}
        if evaluate_fold:
            res = kev.evaluate(to_eval(P, A), kev.prepare("cv", fold=fold), missing="error")
            run_.save_metrics(res)
            out["summary"] = summarise(res)
        print(f"[{name}] fold {fold}: {lg}", flush=True)
        return out


def to_eval(P: pd.DataFrame, A=None) -> pd.DataFrame:
    """kinase_idx/site_idx (kinases_with_canaries / site index order) -> kinase_acc/site_id for ksp.eval."""
    from ksp import canary, index
    K = canary.kinases_with_canaries()
    S = index.load_sites()
    return pd.DataFrame({"kinase_acc": K.uniprot_acc.values[P.kinase_idx.values],
                         "site_id": S.site_id.values[P.site_idx.values], "score": P.score.values,
                         **({"score_site": P.score_site.values} if "score_site" in P else {})})


def summarise(res) -> dict:
    E = res.endpoints
    st = kev.strata(res)
    r = {"M1": E["M1_auroc_hard"], "M1_fold_mean": E.get("M1_fold_mean"), "M1_per_fold": E.get("M1_per_fold"),
         "M2": E["M2_hit5"], "M3": E["M3_top1_recall"], "M4": E["M4_auroc_hard_Kany"]}
    for b in kev.TRAIN_BIN_LABELS:
        x = st[(st.stratum == "train_bin") & (st.value == b)]
        r[f"M5bin_{b}"] = float(x.auroc_hard.iloc[0]) if len(x) else float("nan")
    for cl in ("ST", "Y"):
        x = st[(st.stratum == "residue_class") & (st.value == cl)]
        r[f"M1_{cl}"] = float(x.auroc_hard.iloc[0]) if len(x) else float("nan")
    ps = res.pair_strata
    for _, row in ps[ps.kinase_set == "K_primary"].iterrows():
        r[f"pair_{row.stratum}_M1"] = row.auroc_hard
        r[f"pair_{row.stratum}_nk"] = row.n_kinases
        r[f"pair_{row.stratum}_npos"] = row.n_pos
    c = res.canary
    r.update({"canary_neg_hard": c.get("neg_mean_auroc_hard"), "canary_neg_all": c.get("neg_mean_auroc_all"),
              "canary_pos": c.get("pos_auroc_hard"), "canary_gate": bool(c.get("neg_pass") and c.get("pos_pass"))})
    return r


def register(row: dict) -> int:
    cols = ["timestamp", "model", "run_id", "split", "fold", "config_name", "config_json", "selection_metric", "value",
            "n_configs_tried_model", "notes"]
    r = pd.read_csv(REGISTRY, sep="\t")
    n_prev = int((r.model == "REGIME").sum())
    row = {**row, "model": "REGIME", "n_configs_tried_model": n_prev + 1,
           "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(REGISTRY, "a") as fh:
        fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")
    return n_prev + 1


def pool(name: str, fold_runs: list, config: dict, notes: str = "", bootstrap: bool = False,
         register_row: bool = True, preds: pd.DataFrame | None = None) -> dict:
    """Pool 5 per-fold prediction sets (or a given pooled DataFrame) -> ksp.eval pooled CV (A24), save, register."""
    inputs = []
    if preds is None:
        parts = []
        for rid in fold_runs:
            m = json.loads((ROOT / "results" / "runs" / rid / "manifest.json").read_text())
            p = ROOT / m["outputs"]["predictions"]["path"]
            inputs.append(p)
            parts.append(pd.read_parquet(p))
        preds = pd.concat(parts, ignore_index=True)
    with manifest.RunManifest(name, config=config, seed=config.get("seed", 0), split="cv", inputs=inputs,
                              notes=f"pooled 5-fold OOF of {fold_runs}; {notes}") as run_:
        run_.set("fold_runs", fold_runs)
        P = preds if "kinase_acc" in preds else to_eval(preds)
        run_.save_predictions(P)
        res = kev.evaluate(P, kev.prepare("cv", fold=None), missing="error")
        run_.save_metrics(res)
        run_.dir.mkdir(parents=True, exist_ok=True)
        res.per_kinase.to_csv(run_.dir / "per_kinase.tsv", sep="\t", index=False)
        kev.strata(res).to_csv(run_.dir / "strata.tsv", sep="\t", index=False)
        res.pair_strata.to_csv(run_.dir / "pair_strata.tsv", sep="\t", index=False)
        s = summarise(res)
        if bootstrap:
            bs = kev.bootstrap(res, B=1000, seed=0)
            bs.to_csv(run_.dir / "bootstrap.tsv", sep="\t", index=False)
        (run_.dir / "summary.json").write_text(json.dumps(s, indent=1, default=str))
        if register_row:
            register({"run_id": run_.run_id, "split": "cv", "fold": "all(pooled)", "config_name": name,
                      "config_json": json.dumps({k: v for k, v in config.items() if k != "groups"}, default=str),
                      "selection_metric": "M1_auroc_hard(K_cv) [=M5]", "value": round(s["M1"], 5),
                      "notes": f"M2={s['M2']:.4f} M3={s['M3']:.4f} M4={s['M4']:.4f} fold-mean={s['M1_fold_mean']:.4f}; "
                               f"canary_neg={s['canary_neg_hard']:.3f} pos={s['canary_pos']:.3f}; {notes}"})
        print(f"[{name}] {json.dumps(s, default=str)}", flush=True)
        return {"run_id": run_.run_id, "summary": s, "res": res}
