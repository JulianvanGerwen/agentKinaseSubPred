"""R2 leave-kinase-out T2 (DESIGN §4.5): for kinase fold j (results/splits/kinase_folds.tsv), train the final R2
configuration on all DEV labels of kinases NOT in fold j (positives and negative masks exclude fold-j kinases;
label-derived features from xfit tag lko<j>, fitted without fold-j labels), then score fold-j kinases x all DEV U_ms
sites of matching class. Pool the 5 folds and evaluate like R3's LKO pool: ksp.eval.prepare('cv'),
evaluate(kinases = held-out kinases with positives, missing='drop').
Usage: r2_lko.py --name r2_final_lko --set <config overrides> [--lko 0 1 2 3 4]"""
import argparse, json, time
import numpy as np, pandas as pd
from ksp import eval as kev, manifest, splits, canary
from ksp.index import ROOT
from ksp.models.r2 import assemble, sample, train

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--set", nargs="*", default=[])
ap.add_argument("--lko", nargs="+", type=int, default=[0, 1, 2, 3, 4])
a = ap.parse_args()
cfg = train.default_config()
for kv in a.set:
    k, v = kv.split("=", 1)
    try:
        cfg[k] = json.loads(v)
    except json.JSONDecodeError:
        cfg[k] = v
kf = splits.load_kinase_folds()
K = canary.kinases_with_canaries()
kpos = {acc: i for i, acc in enumerate(K.uniprot_acc)}
with manifest.RunManifest(a.name, config=cfg, seed=cfg["seed"], split="lko_T2", notes=f"lko folds {a.lko}") as run:
    preds, logs = [], {}
    for j in a.lko:
        t = time.time()
        held = set(kf.kinase_acc[kf.lko_fold == j])
        fd = sample.FoldData(None, cfg["variant"], exclude_kinases=held)
        A = assemble.Assembler(cfg["groups"], None, xfit_tag=f"lko{j}", clean=cfg["clean"], comp_know=cfg["comp_know"],
                               f2a_family=cfg.get("f2a_family", "f2a"), desc_v=cfg.get("desc_v", 1))
        if set(fd.pos.kinase_acc) & held:
            raise RuntimeError("held-out kinase labels in training positives")
        lg = {}
        bst, cols, stats = train.train_fold(fd, A, cfg, lg)
        fd.eval_sites = fd.train_sites                    # T2: seen proteins, unseen kinases
        hk = sorted(kpos[x] for x in held if x in kpos)
        P = train.predict_fold(fd, A, bst, cfg, kinases=hk)
        lg.update(n_held=len(held), n_pred_rows=len(P), t_s=round(time.time() - t, 1))
        logs[j] = lg
        preds.append(P)
        print(f"lko {j}: {lg}", flush=True)
    P = pd.concat(preds, ignore_index=True)
    run.set("lko_logs", logs)
    run.save_predictions(P)
    data = kev.prepare("cv", fold=None)
    heldk = sorted(set(K.uniprot_acc.values[np.unique(P.kinase_idx.values)]) & set(data.positives))
    res = kev.evaluate(P, data, kinases=heldk, missing="drop")
    run.save_metrics(res)
    res.per_kinase.to_csv(run.dir / "per_kinase.tsv", sep="\t", index=False)
    kev.strata(res).to_csv(run.dir / "strata.tsv", sep="\t", index=False)
    E = res.endpoints
    print(f"[{a.name}] LKO-T2 M1={E['M1_auroc_hard']:.4f} M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f} "
          f"M4={E['M4_auroc_hard_Kany']:.4f} n_K={E['n_K_primary_evaluated']}", flush=True)
    train.register({"run_id": run.run_id, "split": "lko_T2", "fold": "all(pooled)", "config_name": a.name,
                    "config_json": json.dumps({k: v for k, v in cfg.items() if k != "groups"} | {"groups": "+".join(cfg["groups"])}),
                    "selection_metric": "M1_auroc_hard(K_cv, LKO T2)", "value": round(E["M1_auroc_hard"], 5),
                    "notes": f"LKO T2 of final config (not used for selection); M2={E['M2_hit5']:.4f} M3={E['M3_top1_recall']:.4f}"})
