"""Test, time and sanity-check the label-derived F3 neighbourhood feature (ksp.features.ppi_label).

1. Unit tests on toy inputs (self / same-cluster / autophos exclusion; pair vs dense mode agree).
2. Timing: dense mode, all 503 kinases x 18,702 proteins, training = all L2 pairs.
3. Cross-fitted sanity (NOT an official metric): DEV proteins in the official 5 CV folds (results/splits/cv_folds.tsv);
   feature for proteins in fold f uses training pairs (load_dev L2) whose substrate is outside f. Cluster maps: (a) none (self only), (b) dummy
   gene-family map (HGNC symbol with trailing digits/letters stripped, e.g. CDK1/CDK2 -> CDK) as a stand-in for
   results/splits/substrate_clusters.tsv, (c) the real cluster map if that file exists (columns acc, cluster).
   Enrichment vs RAND/SHUF backgrounds as in context_sanity.py.
Run: OMP_NUM_THREADS=10 .venv/bin/python -I scripts/features/context_ppi_label.py
"""
import json
import re
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts" / "features"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from context_sanity import SEEDS, backgrounds, evaluate_matrix, positives  # noqa: E402
from ksp.features.ppi_label import OUT_NAMES, load_graph_data, substrate_neighbourhood  # noqa: E402
from ksp.index import load_kinases, load_proteins  # noqa: E402

RES = ROOT / "results" / "features" / "context"
CL_FILE = ROOT / "results" / "splits" / "substrate_clusters.tsv"


def unit_tests(g):
    nodes = g.nodes
    # pick a kinase-like seed and some neighbours from the graph
    deg = np.asarray(g.A.sum(1)).ravel()
    hub = int(np.argmax(deg))
    nb = g.A[hub].indices[:4]
    k = "KIN_TEST"
    q = nodes[nb[0]]
    s1, s2, s3 = nodes[nb[1]], nodes[nb[2]], nodes[nb[3]]
    tr = pd.DataFrame({"kinase_acc": [k] * 4, "acc": [q, s1, s2, s3]})
    qp = pd.DataFrame({"kinase_acc": [k], "acc": [q]})
    r0 = substrate_neighbourhood(tr, qp, graph=g).iloc[0]
    assert r0.nb_n_used == 3, r0  # q itself excluded
    r1 = substrate_neighbourhood(tr, qp, exclude_same_cluster={q: "C", s1: "C"}, graph=g).iloc[0]
    assert r1.nb_n_used == 2, r1  # same-cluster s1 excluded
    r2 = substrate_neighbourhood(tr, qp, exclude_same_cluster={q: "C", s1: "C", s2: "C", s3: "C"}, graph=g).iloc[0]
    assert r2.nb_n_used == 0 and np.isnan(r2.nb_max), r2
    # autophos: training pair (k, k) ignored
    tr2 = pd.concat([tr, pd.DataFrame({"kinase_acc": [s1], "acc": [s1]})])
    r3 = substrate_neighbourhood(tr2, pd.DataFrame({"kinase_acc": [s1], "acc": [q]}), graph=g).iloc[0]
    assert r3.nb_n_used == 0, r3
    # manual check of max
    s_idx = [g.pos[a] for a in (s1, s2, s3)]
    assert np.isclose(r0.nb_max, g.DC[g.pos[q], s_idx].max(), atol=1e-5)
    assert np.isclose(r0.nb_mean, g.DC[g.pos[q], s_idx].mean(), atol=1e-5)
    # dense vs pair mode
    dense = substrate_neighbourhood(tr, None, graph=g, kinases=[k], proteins=[q, s1, "NOT_A_NODE"])
    assert np.isclose(dense["nb_max"][0, 0], r0.nb_max) and np.isnan(dense["nb_max"][0, 2])
    return "ok"


def family(sym):
    return re.sub(r"\d+[A-Z]*\d*$", "", sym) or sym


def main():
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    g = load_graph_data()
    t_load = time.time() - t0
    report = {"graph_load_s": round(t_load, 1), "unit_tests": unit_tests(g)}
    kin = load_kinases().uniprot_acc.tolist()
    prot = load_proteins().acc.tolist()
    pi, pj, known, K, P, dev = positives()
    train_all = pd.DataFrame({"kinase_acc": np.array(kin)[pi], "acc": np.array(prot)[pj]})
    t1 = time.time()
    dense = substrate_neighbourhood(train_all, None, graph=g)
    report["dense_all_kinases_s"] = round(time.time() - t1, 1)
    print(report, flush=True)

    # cross-fitted sanity
    from ksp import splits
    cvf = splits.load_cv_folds()
    fold = np.array([cvf.get(a, -1) for a in prot])  # -1 = TEST protein (never queried, never trained on)
    sites = pd.read_parquet(ROOT / "cache" / "index" / "sites.parquet", columns=["acc", "gene"]).drop_duplicates("acc")
    gene = dict(zip(sites.acc, sites.gene.fillna("")))
    cmaps = {"self_only": None, "dummy_gene_family": {a: family(str(gene.get(a, a))) for a in prot}}
    if CL_FILE.exists():
        cmaps["substrate_clusters"] = splits.load_clusters().to_dict()
    bgs = [backgrounds(pi, pj, known, dev, s) for s in SEEDS]
    rows = []
    for cname, cmap in cmaps.items():
        cf = {n: np.full((K, P), np.nan, np.float32) for n in OUT_NAMES}
        for f in range(5):
            tr = train_all[fold[pj] != f]
            cols = np.where(fold == f)[0]
            out = substrate_neighbourhood(tr, None, exclude_same_cluster=cmap, graph=g, kinases=kin,
                                          proteins=[prot[c] for c in cols])
            for n in OUT_NAMES:
                cf[n][:, cols] = out[n]
        for n in OUT_NAMES:
            r = evaluate_matrix(cf[n], pi, pj, bgs)
            r.update(module="ppi_label", feature=f"{n}[cf5,{cname}]")
            rows.append(r)
            print(f"{r['feature']:40s} cov={r['coverage_pos']:.3f} AUCrand={r['auc_rand_all']:.3f} AUCshuf={r['auc_shuf_all']:.3f}", flush=True)
    # in-sample (NOT cross-fitted) for contrast: shows the leak if cross-fitting is skipped
    for n in ["nb_max", "nb_frac_adj"]:
        r = evaluate_matrix(dense[n], pi, pj, bgs)
        r.update(module="ppi_label", feature=f"{n}[IN-SAMPLE,self-excluded]")
        rows.append(r)
    df = pd.DataFrame(rows)
    df = df[["module", "feature"] + [c for c in df.columns if c not in ("module", "feature")]]
    df.to_csv(RES / "sanity_ppi_label.tsv", sep="\t", index=False, float_format="%.4f")
    report["total_s"] = round(time.time() - t0, 1)
    (RES / "ppi_label_report.json").write_text(json.dumps(report, indent=1))
    print(report)


if __name__ == "__main__":
    main()
