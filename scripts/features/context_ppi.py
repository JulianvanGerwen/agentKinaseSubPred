"""Build F3 PPI feature matrices (kinase x protein, aligned to cache/index) from systematic screens only.

Run: OMP_NUM_THREADS=10 .venv/bin/python -I scripts/features/context_ppi.py [--restart 0.5]
Outputs: cache/features/context/ppi/*.npy + descriptor.json (copied to results/features/context/ppi_descriptor.json),
         cache/features/context/ppi/edges.parquet, rwr_full.npy (N x N, float32) + rwr_nodes.txt (for ppi_label),
         results/features/context/id_mapping_rates.json, results/features/context/ppi_graph_stats.json
"""
import argparse
import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse.csgraph as csg  # noqa: E402

from ksp.features.ppi import SOURCES, Graph, load_all_edges  # noqa: E402
from ksp.features.ppi_idmap import MapLog  # noqa: E402
from ksp.index import load_kinases, load_proteins  # noqa: E402
from ksp.store import save_matrix  # noqa: E402

OUT = ROOT / "cache" / "features" / "context" / "ppi"
RES = ROOT / "results" / "features" / "context"
STUDIES = {"bioplex": ["bioplex_293t", "bioplex_hct116"], "huri": ["huri"], "opencell": ["opencell"],
           "buljan": ["buljan"], "hcm": ["hcm"]}


def index_meta():
    man = json.loads((ROOT / "results" / "index_manifest.json").read_text())
    return {"rows": "cache/index/kinases.parquet (uniprot_acc order)", "rows_sha256": man["kinases"]["sha256"],
            "cols": "cache/index/proteins.parquet (acc order)", "cols_sha256": man["proteins"]["sha256"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restart", type=float, default=0.5)
    args = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    kin = load_kinases().uniprot_acc.tolist()
    prot = load_proteins().acc.tolist()
    K, P = len(kin), len(prot)
    log = MapLog()
    edges = load_all_edges(log)
    log.dump(RES / "id_mapping_rates.json")
    edges.to_parquet(OUT / "edges.parquet", index=False)
    kpos = {a: i for i, a in enumerate(kin)}
    ppos = {a: i for i, a in enumerate(prot)}
    feats = {}
    desc = {"module": "F3 PPI (systematic screens only)", "label_derived": False, **index_meta(),
            "restart": args.restart, "sources": SOURCES, "features": {}}

    def save(name, arr, dtype, text, missing, sign=1):
        save_matrix(f"context/ppi/{name}", arr.astype(dtype), "kinase", "protein",
                    meta={"description": text, "missing": missing, "label_derived": False, "module": "F3"})
        desc["features"][name] = {"file": f"cache/features/context/ppi/{name}.npy", "store_name": f"context/ppi/{name}",
                                  "dtype": str(np.dtype(dtype)), "sign": sign,
                                  "shape": list(arr.shape), "description": text, "missing": missing}

    # ---- per-source direct edges + presence
    k_in = np.zeros((K, len(SOURCES)), bool)
    p_in = np.zeros((P, len(SOURCES)), bool)
    any_edge = np.zeros((K, P), np.int8)
    study_count = np.zeros((K, P), np.int8)
    per_src = {}
    for si, s in enumerate(SOURCES):
        e = edges[edges.source == s]
        nodes = set(e.a) | set(e.b)
        k_in[:, si] = [a in nodes for a in kin]
        p_in[:, si] = [a in nodes for a in prot]
        M = np.zeros((K, P), np.int8)
        for a, b in zip(e.a, e.b):
            if a in kpos and b in ppos:
                M[kpos[a], ppos[b]] = 1
            if b in kpos and a in ppos:
                M[kpos[b], ppos[a]] = 1
        per_src[s] = M
        save(f"ppi_edge_{s}", M, np.int8, f"1 if a direct {s} edge links kinase and protein (undirected)",
             "0 = no edge (incl. not screened); see ppi_kinase_in_source / ppi_protein_in_source")
        any_edge |= M
    for st, srcs in STUDIES.items():
        study_count += np.maximum.reduce([per_src[s] for s in srcs])
    save("ppi_edge_any", any_edge, np.int8, "1 if any screen has a direct edge", "0 = none")
    save("ppi_n_studies", study_count, np.int8, "number of independent studies (BioPlex(either line), HuRI, OpenCell, "
         "Buljan, HCM) with a direct edge (0-5)", "0 = none")
    save_matrix("context/ppi/ppi_kinase_in_source", k_in, "kinase", "none", col_ids=SOURCES,
                meta={"description": "kinase is a node of source s"})
    save_matrix("context/ppi/ppi_protein_in_source", p_in, "protein", "none", col_ids=SOURCES,
                meta={"description": "protein is a node of source s"})
    desc["presence"] = {"ppi_kinase_in_source": "cache/features/context/ppi/ppi_kinase_in_source.npy (K x S bool, S in 'sources' order)",
                        "ppi_protein_in_source": "cache/features/context/ppi/ppi_protein_in_source.npy (P x S bool)"}

    # ---- combined graph + RWR
    g = Graph(edges)
    N = len(g.nodes)
    t1 = time.time()
    R = g.rwr_dense(args.restart)  # N x N, column = seed
    print(f"RWR dense {N} nodes in {time.time() - t1:.0f}s", flush=True)
    np.save(OUT / "rwr_full.npy", R)
    (OUT / "rwr_nodes.txt").write_text("\n".join(g.nodes) + "\n")
    bg = R.mean(1)  # mean visiting prob over uniformly chosen seeds (= PageRank with uniform teleport)
    np.save(OUT / "rwr_background.npy", bg)
    k_node = np.array([g.pos.get(a, -1) for a in kin])
    p_node = np.array([g.pos.get(a, -1) for a in prot])
    kin_ok, prot_ok = k_node >= 0, p_node >= 0
    rwr = np.full((K, P), np.nan, np.float32)
    sub = R[np.ix_(p_node[prot_ok], k_node[kin_ok])].T  # kinase seeds x proteins
    rwr[np.ix_(kin_ok, prot_ok)] = sub
    eps = 0.01 / N
    dc = np.log2((rwr + eps) / (bg[np.clip(p_node, 0, None)][None, :] + eps)).astype(np.float32)
    save("ppi_rwr", rwr, np.float32, f"RWR (restart {args.restart}) visiting probability of protein from kinase seed on "
         "the combined screens graph (column-stochastic, sums to 1 over graph nodes)", "NaN = kinase or protein not in graph")
    save("ppi_rwr_dc", dc, np.float32, "degree-corrected RWR: log2((rwr+eps)/(bg+eps)), bg = mean RWR of the protein "
         "over all seeds (hub correction), eps=0.01/N", "NaN = kinase or protein not in graph")
    # percentile of protein within the kinase's RWR row (rank-normalised, within graph proteins)
    pr = np.full((K, P), np.nan, np.float32)
    for i in np.where(kin_ok)[0]:
        v = dc[i, prot_ok]
        pr[i, prot_ok] = (pd.Series(v).rank(pct=True)).to_numpy()
    save("ppi_rwr_dc_pct", pr, np.float16, "within-kinase percentile of ppi_rwr_dc over graph proteins", "NaN as ppi_rwr")
    # shortest path (hops), capped at 5 (= >=5 or unreachable)
    D = csg.shortest_path(g.A, unweighted=True, indices=k_node[kin_ok], directed=False)
    D = np.minimum(D, 5)
    dist = np.full((K, P), np.nan, np.float32)
    dist[np.ix_(kin_ok, prot_ok)] = D[:, p_node[prot_ok]]
    save("ppi_dist", dist, np.float16, "shortest path length (hops) kinase->protein on combined graph, capped at 5",
         "NaN = kinase or protein not in graph", sign=-1)

    # ---- stats
    deg = pd.Series(g.deg, index=g.nodes)
    stats = {"n_nodes": int(N), "n_edges": int(g.A.nnz // 2),
             "edges_per_source": edges.groupby("source").size().to_dict(),
             "nodes_per_source": {s: int(len(set(edges[edges.source == s].a) | set(edges[edges.source == s].b))) for s in SOURCES},
             "kinases_in_graph": int(kin_ok.sum()), "kinases_total": K,
             "kinases_in_source": dict(zip(SOURCES, map(int, k_in.sum(0)))),
             "proteins_in_graph": int(prot_ok.sum()), "proteins_total": P,
             "proteins_in_source": dict(zip(SOURCES, map(int, p_in.sum(0)))),
             "median_degree": float(deg.median()), "max_degree": float(deg.max()),
             "n_components": int(csg.connected_components(g.A)[0]), "runtime_s": round(time.time() - t0)}
    (RES / "ppi_graph_stats.json").write_text(json.dumps(stats, indent=1))
    desc["stats"] = stats
    (OUT / "descriptor.json").write_text(json.dumps(desc, indent=1))
    shutil.copy(OUT / "descriptor.json", RES / "ppi_descriptor.json")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
