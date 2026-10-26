#!/usr/bin/env python
"""Build and FREEZE the splits (DESIGN §4, G1) and seal TEST labels.

Inputs : cache/index/*, cache/labels/ks_pairs.parquet (scripts/build_labels.py), cache/splits/allvsall.m8 (mmseqs
         all-vs-all, recomputed if missing), external/kincore (kinase domain identity).
Outputs: results/splits/{substrate_clusters,test_proteins,cv_folds,kinase_folds,eval_kinase_sets,
         canary_ppi_partners}.tsv, split_report.md, SEAL.json; cache/labels/{dev_pairs,canary_dev}.parquet,
         results/labels/kinase_label_counts_dev.tsv; sealed/{test_labels,ks_pairs_full,...} (moved, not copied).
Usage  : /workspace/.venv/bin/python -I scripts/build_splits.py [--threads 12] [--force]
Refuses to overwrite a frozen split (results/splits/SEAL.json exists) unless --force.
"""
import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ksp import canary, index, labels, sites, splits  # noqa: E402

SEALED = ROOT / "sealed"
SEED = 0


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="optimise + report only; write/move/seal nothing")
    ap.add_argument("--report", action="store_true", help="only (re)write results/splits/split_report.md")
    ap.add_argument("--relabel", action="store_true",
                    help="keep the FROZEN split; rebuild DEV/TEST label tables from a rebuilt ks_pairs.parquet and "
                         "re-seal (verifies that evaluation positives / kinase sets are unchanged)")
    ap.add_argument("--iter-test", type=int, default=400_000)
    ap.add_argument("--iter-cv", type=int, default=800_000)
    a = ap.parse_args()
    if a.relabel:
        relabel()
        return write_report()
    if a.report:
        return write_report()
    t0 = time.time()
    out = splits.SPLIT_DIR if not a.dry_run else splits.CACHE / "dry_run"
    if (out / "SEAL.json").exists() and not a.force:
        sys.exit("split already frozen (results/splits/SEAL.json); use --force only with PI approval")
    out.mkdir(parents=True, exist_ok=True)
    splits.CACHE.mkdir(parents=True, exist_ok=True)
    SEALED.mkdir(exist_ok=True)
    if not labels.FULL_PAIRS.exists():
        sys.exit(f"{labels.FULL_PAIRS} missing (already sealed?) - rerun scripts/build_labels.py first")
    R = {}  # report numbers

    # ------------------------------------------------------------------ index + labels
    K = index.load_kinases()
    S = index.load_sites()
    P = index.load_proteins().acc
    kpos = {acc: i for i, acc in enumerate(K.uniprot_acc)}
    ums = set(S.site_id[S.in_U_ms])
    in_idx = set(S.site_id)
    full = annotate_full(pd.read_parquet(labels.FULL_PAIRS), ums, in_idx)
    L2 = labels.select_variant(full, "L2")
    EV = labels.select_variant(full, "eval")
    R["L2_pairs"] = len(L2)
    R["L2_pairs_not_in_index"] = int((~L2.site_id.isin(in_idx)).sum())
    R["L2_pairs_in_index_not_U_ms"] = int((L2.site_id.isin(in_idx) & ~L2.site_id.isin(ums)).sum())
    R["eval_pairs"] = len(EV)
    R["eval_pairs_not_U_ms(dropped)"] = int((~EV.site_id.isin(ums)).sum())
    L2u = L2[L2.site_id.isin(ums)]
    EVu = EV[EV.site_id.isin(ums)]
    R["L2_pairs_U_ms(used)"], R["eval_pairs_U_ms(used)"] = len(L2u), len(EVu)

    # ------------------------------------------------------------------ canaries (before split; PI r2)
    n2 = L2.groupby("kinase_acc").site_id.nunique()
    size_src = tuple(sorted(int(x) for x in n2[(n2 >= 10) & (n2 <= 100)]))
    can, partners = canary.generate_v1(size_src)   # split freeze used v1; --relabel installs v2
    R["canary_pairs"] = len(can)
    R["canary_sizes"] = can.groupby("kinase_acc").size().to_dict()
    KC = canary.kinases_with_canaries()
    kposc = {acc: i for i, acc in enumerate(KC.uniprot_acc)}

    # ------------------------------------------------------------------ substrate clusters
    m8f = splits.CACHE / "allvsall.m8"
    if not m8f.exists():
        fa = splits.CACHE / "proteins.fasta"
        splits.write_protein_fasta(fa, P.tolist(), sites.default_store(verbose=False).seqs)
        splits.run_mmseqs_allvsall(fa, m8f, threads=a.threads)
    R["mmseqs_m8_sha256"] = sha(m8f)
    cl, crep = splits.substrate_clusters(splits.read_m8(m8f), P, S[["acc", "window15"]], L2.site_id)
    R.update({f"clusters:{k}": v for k, v in crep.items()})
    cmap = pd.Series(cl.cluster_id.values, index=cl.acc.values)
    nC = int(cl.cluster_id.max()) + 1

    def cl_of(site_ids):
        return cmap[splits.site_protein(site_ids)].values

    # per-cluster global quantities: L2 labels(U_ms), U_ms sites, S/T/Y L2 labels, cluster count
    glob = np.zeros((nC, 6))
    np.add.at(glob[:, 0], cl_of(L2u.site_id), 1)
    Su = S[S.in_U_ms]
    np.add.at(glob[:, 1], cmap[Su.acc].values, 1)
    for j, r in enumerate("STY"):
        x = L2u[L2u.res.eq(r)]
        np.add.at(glob[:, 2 + j], cl_of(x.site_id), 1)
    glob[:, 5] = 1
    gnames = ["L2_labels", "U_ms_sites", "S_labels", "T_labels", "Y_labels", "clusters"]
    tol = np.array([0.01, 0.01, 0.02, 0.02, 0.02, 0.02])
    gw = np.array([1.0, 1.0, 0.5, 0.5, 0.5, 0.3])
    # kinase positives used for balancing: evaluable positives (U_ms) of real kinases + all canary labels
    kp_sites = pd.concat([EVu[["kinase_acc", "site_id"]], can[["kinase_acc", "site_id"]]], ignore_index=True)
    kin_idx = kp_sites.kinase_acc.map(kposc).values.astype(np.int64)
    kin_cl = cl_of(kp_sites.site_id)

    # ------------------------------------------------------------------ TEST
    t1 = time.time()
    tgt = np.array([0.8, 0.2])
    a_test = splits.optimise_assignment(np.arange(nC), kin_cl, kin_idx, len(KC), glob, tgt,
                                        np.array([False, True]), tol, seed=SEED, n_iter=a.iter_test, glob_w=gw)
    R["anneal_test_seconds"] = round(time.time() - t1, 1)
    is_test_c = a_test == 1
    cl["split"] = np.where(is_test_c[cl.cluster_id.values], "TEST", "DEV")
    test_prot = set(cl.acc[cl.split.eq("TEST")])

    # ------------------------------------------------------------------ DEV CV
    t1 = time.time()
    dev_c = np.where(~is_test_c)[0]
    remap = -np.ones(nC, dtype=np.int64)
    remap[dev_c] = np.arange(len(dev_c))
    dmask = ~is_test_c[kin_cl]
    a_cv = splits.optimise_assignment(np.arange(len(dev_c)), remap[kin_cl[dmask]], kin_idx[dmask], len(KC),
                                      glob[dev_c], np.full(5, 0.2), np.ones(5, bool), tol, seed=SEED,
                                      n_iter=a.iter_cv, glob_w=gw)
    R["anneal_cv_seconds"] = round(time.time() - t1, 1)
    fold_c = -np.ones(nC, dtype=np.int64)
    fold_c[dev_c] = a_cv
    cl["cv_fold"] = fold_c[cl.cluster_id.values]

    # ------------------------------------------------------------------ write split tables
    ns = S.groupby("acc").size()
    nu = Su.groupby("acc").size()
    nl = L2u.groupby("acc").size()
    cl["n_sites"] = cl.acc.map(ns).fillna(0).astype(int)
    cl["n_U_ms"] = cl.acc.map(nu).fillna(0).astype(int)
    cl["n_L2_labels"] = cl.acc.map(nl).fillna(0).astype(int)
    cl[["acc", "cluster_id", "n_sites", "n_U_ms", "n_L2_labels"]].to_csv(out / "substrate_clusters.tsv", sep="\t",
                                                                         index=False)
    cl.loc[cl.split.eq("TEST"), ["acc", "cluster_id"]].to_csv(out / "test_proteins.tsv", sep="\t", index=False)
    cl.loc[cl.split.eq("DEV"), ["acc", "cluster_id", "cv_fold"]].rename(columns={"cv_fold": "fold"}).to_csv(
        out / "cv_folds.tsv", sep="\t", index=False)
    pd.DataFrame({"acc": partners}).to_csv(out / "canary_ppi_partners.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ label tables + seal
    full_acc = full.acc
    in_prot_index = full_acc.isin(set(P))
    R["label_rows_total"] = len(full)
    R["label_rows_protein_not_in_index(dropped from DEV/TEST; kept in sealed full)"] = int((~in_prot_index).sum())
    is_t = full_acc.isin(test_prot)
    dev = full[in_prot_index & ~is_t].reset_index(drop=True)
    tst = full[in_prot_index & is_t].reset_index(drop=True)
    can_t = can.acc.isin(test_prot)
    can_dev = can[~can_t].reset_index(drop=True)
    tst = pd.concat([tst, can[can_t]], ignore_index=True)
    assert not set(dev.acc) & test_prot
    if not a.dry_run:
        dev.to_parquet(labels.DEV_PAIRS, index=False)
        can_dev.to_parquet(labels.DEV_PAIRS.parent / canary.CANARY_FILE_DEV, index=False)
        tst.to_parquet(SEALED / "test_labels.parquet", index=False)
    # DEV per-kinase counts (public)
    dcnt = (dev.groupby("kinase_acc").agg(**{v: (v, "sum") for v in ["L1", "L2", "L3", "L2_noauto", "EXT"]})
            .reset_index())
    dcnt = K[["uniprot_acc", "gene"]].merge(dcnt, left_on="uniprot_acc", right_on="kinase_acc", how="left") \
        .drop(columns="kinase_acc").fillna(0)
    dcnt.to_csv((out if a.dry_run else ROOT / "results/labels") / "kinase_label_counts_dev.tsv", sep="\t",
                index=False)

    # ------------------------------------------------------------------ kinase folds (LKO)
    ide = labels.kinase_domain_identity()
    kcl = splits.kinase_clusters(K, ide, 0.50)
    nlab = K.uniprot_acc.map(dev[dev.L2].groupby("kinase_acc").site_id.nunique()).fillna(0).values
    kf = splits.kinase_folds(kcl, nlab, 5, SEED)
    kfold = pd.DataFrame({"kinase_acc": K.uniprot_acc, "gene": K.gene, "kinase_cluster": kcl,
                          "manning_group": K.manning_group, "manning_family": K.manning_family,
                          "n_dev_L2": nlab.astype(int), "lko_fold": kf})
    kfold.to_csv(out / "kinase_folds.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ evaluable kinase sets (§4.6)
    ev_t = EVu[EVu.acc.isin(test_prot)]
    ev_d = EVu[~EVu.acc.isin(test_prot)]
    tot_l2 = L2.groupby("kinase_acc").site_id.nunique()

    def stats(e):
        return e.assign(c=cl_of(e.site_id)).groupby("kinase_acc").agg(n=("site_id", "nunique"), nc=("c", "nunique"))
    st_t, st_d = stats(ev_t), stats(ev_d)
    es = K[["uniprot_acc", "gene", "residue_class", "manning_group", "manning_family", "is_pseudokinase"]].rename(
        columns={"uniprot_acc": "kinase_acc"})
    es["n_test_pos"] = es.kinase_acc.map(st_t.n).fillna(0).astype(int)
    es["n_test_clusters"] = es.kinase_acc.map(st_t.nc).fillna(0).astype(int)
    es["n_dev_pos"] = es.kinase_acc.map(st_d.n).fillna(0).astype(int)
    es["n_dev_clusters"] = es.kinase_acc.map(st_d.nc).fillna(0).astype(int)
    es["n_L2_total"] = es.kinase_acc.map(tot_l2).fillna(0).astype(int)
    es["n_dev_L2"] = nlab.astype(int)
    pseudo_ok = ~es.is_pseudokinase | (es.n_L2_total >= 3)
    es["in_K_eval"] = (es.n_test_pos >= 3) & (es.n_test_clusters >= 2) & pseudo_ok
    es["in_K_any"] = (es.n_test_pos >= 1) & pseudo_ok
    es["in_K_cv"] = (es.n_dev_pos >= 5) & (es.n_dev_clusters >= 2) & pseudo_ok
    es["train_bin"] = pd.cut(es.n_dev_L2, [-1, 0, 4, 19, 99, 1e9], labels=["0", "1-4", "5-19", "20-99", ">=100"])
    es.to_csv(out / "eval_kinase_sets.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ bio controls
    bc_path = ROOT / "config" / "bio_controls.tsv"
    if bc_path.exists() and not a.dry_run:
        lines = bc_path.read_text().splitlines()
        com = [x for x in lines if x.startswith("#")]
        bc = pd.read_csv(bc_path, sep="\t", comment="#")
        bc["in_test_split"] = bc.site_id.str.rsplit("_", n=1).str[0].isin(test_prot)
        with open(bc_path, "w") as fh:
            fh.write("\n".join(com) + "\n")
            bc.to_csv(fh, sep="\t", index=False)
        R["bio_controls_in_test"] = f"{int(bc.in_test_split.sum())}/{len(bc)}"

    # ------------------------------------------------------------------ report numbers
    def frac(mask_sites_acc):
        return float(np.mean(mask_sites_acc))
    R["clusters_test"] = int(is_test_c.sum())
    R["clusters_frac_test"] = round(float(is_test_c.mean()), 4)
    R["proteins_test"] = len(test_prot)
    for j, gname in enumerate(gnames):
        R[f"frac_test:{gname}"] = round(float(glob[is_test_c, j].sum() / glob[:, j].sum()), 4)
        for f in range(5):
            R[f"frac_cv{f}:{gname}"] = round(float(glob[fold_c == f, j].sum() / glob[~is_test_c, j].sum()), 4)
    R["test_eval_pos(U_ms)"] = len(ev_t)
    R["test_L2_pairs(U_ms)"] = int(L2u.acc.isin(test_prot).sum())
    R["test_U_ms_sites"] = int(Su.acc.isin(test_prot).sum())
    # per-kinase band compliance on TEST
    nk = EVu.groupby("kinase_acc").site_id.nunique()
    nt = ev_t.groupby("kinase_acc").site_id.nunique().reindex(nk.index).fillna(0)
    big = nk[nk >= 5].index
    fr_k = (nt[big] / nk[big])
    lo, hi = np.floor(0.15 * nk[big]), np.ceil(0.25 * nk[big])
    R["kinases_ge5_eval_pos"] = len(big)
    R["kinases_ge5_in_15_25pct"] = int(((fr_k >= 0.15) & (fr_k <= 0.25)).sum())
    R["kinases_ge5_in_count_band(floor15%,ceil25%)"] = int(((nt[big] >= lo) & (nt[big] <= hi)).sum())
    R["kinases_ge20_in_15_25pct"] = f"{int(((fr_k >= 0.15) & (fr_k <= 0.25))[nk[big] >= 20].sum())}/" \
                                    f"{int((nk[big] >= 20).sum())}"
    R["per_kinase_test_frac_quantiles(ge5)"] = fr_k.quantile([0, 0.1, 0.5, 0.9, 1]).round(3).tolist()
    worst = fr_k.sub(0.2).abs().sort_values(ascending=False).head(10).index
    R["worst_kinases"] = {K.gene[kpos[k]]: f"{int(nt[k])}/{int(nk[k])}" for k in worst}
    # canaries
    cn = can.groupby("kinase_acc").site_id.nunique()
    ct = can[can_t].groupby("kinase_acc").site_id.nunique().reindex(cn.index).fillna(0)
    R["canary_test_frac"] = (ct / cn).round(3).to_dict()
    for s_ in ["in_K_eval", "in_K_any", "in_K_cv"]:
        R[f"|{s_[3:]}|"] = int(es[s_].sum())
        R[f"{s_[3:]}_by_group"] = es[es[s_]].manning_group.fillna("(none)").value_counts().to_dict()
        R[f"{s_[3:]}_by_class"] = es[es[s_]].residue_class.value_counts().to_dict()
    R["kinase_clusters"] = int(len(set(kcl)))
    R["lko_fold_kinases"] = kfold.groupby("lko_fold").size().to_dict()
    R["lko_fold_dev_L2"] = kfold.groupby("lko_fold").n_dev_L2.sum().to_dict()
    R["test_eval_pos_new_since_2021"] = int(tst[labels.evaluable_positive(tst) & tst.new_since_2021 &
                                                tst.in_U_ms].shape[0])
    R["test_ext_pairs(DL-5 filtered, U_ms)"] = int((tst.ext_test & tst.in_U_ms).sum())
    R["dev_label_rows"], R["test_label_rows"] = len(dev), len(tst)

    if a.dry_run:
        R["runtime_s"] = round(time.time() - t0, 1)
        (out / "split_numbers.json").write_text(json.dumps(R, indent=1, default=str))
        print(json.dumps(R, indent=1, default=str))
        return
    # ------------------------------------------------------------------ move pre-split tables into sealed/
    moved = {}
    for src in [labels.FULL_PAIRS, labels.FULL_PAIRS.parent / "ks_unmapped.parquet",
                ROOT / "results/labels/kinase_label_counts.tsv"]:
        if src.exists():
            dst = SEALED / (src.name if src.suffix == ".parquet" and src.stem != "ks_pairs" else
                            ("ks_pairs_full.parquet" if src.stem == "ks_pairs" else src.name))
            shutil.move(str(src), str(dst))
            moved[str(src.relative_to(ROOT))] = str(dst.relative_to(ROOT))
    seal = {"created": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "seed": SEED, "moved_into_sealed": moved,
            "files": {f.name: {"sha256": sha(f), "bytes": f.stat().st_size} for f in sorted(SEALED.glob("*"))
                      if f.is_file()},
            "public": {str(p.relative_to(ROOT)): sha(p) for p in
                       [out / "substrate_clusters.tsv", out / "test_proteins.tsv", out / "cv_folds.tsv",
                        out / "kinase_folds.tsv", out / "eval_kinase_sets.tsv", out / "canary_ppi_partners.tsv",
                        labels.DEV_PAIRS, labels.DEV_PAIRS.parent / canary.CANARY_FILE_DEV]},
            "test_label_rows": len(tst), "test_proteins": len(test_prot)}
    (out / "SEAL.json").write_text(json.dumps(seal, indent=1))
    R["runtime_s"] = round(time.time() - t0, 1)
    (out / "split_numbers.json").write_text(json.dumps(R, indent=1, default=str))
    print(json.dumps(R, indent=1, default=str))


def write_report():
    out = splits.SPLIT_DIR
    R = json.loads((out / "split_numbers.json").read_text())
    seal = json.loads((out / "SEAL.json").read_text())
    es = splits.load_eval_kinase_sets()
    L = ["# Split report (frozen at G1; DESIGN §4, §16)", "",
         f"Seed 0. Frozen {seal['created']}; relabelled {R.get('relabel', {}).get('when', 'never')} "
         "(label tables only; split unchanged). Generated by `scripts/build_splits.py`; numbers in split_numbers.json.",
         "", "## Substrate clusters", "",
         f"- {R['clusters:n_proteins']:,} proteins; mmseqs2 all-vs-all edges (fident>=0.40, min(qcov,tcov)>=0.50): "
         f"{R['clusters:mmseqs_edges']:,}; identical-15-mer groups {R['clusters:identical_window_groups']:,} "
         f"({R['clusters:identical_window_edges']:,} edges).",
         f"- Connected components: **{R['clusters:n_clusters']:,} clusters** ({R['clusters:singleton_clusters']:,} "
         f"singletons; largest {R['clusters:largest_cluster_proteins']} proteins).",
         f"- Largest component label share {R['clusters:max_component_label_frac_before']:.2%} <= 3% -> "
         f"no re-thresholding ({R['clusters:components_rethresholded']} components).", "",
         "## TEST (simulated annealing, soft penalties)", "",
         f"- {R['clusters_test']:,} clusters ({R['clusters_frac_test']:.1%}), {R['proteins_test']:,} proteins, "
         f"{R['test_U_ms_sites']:,} U_ms sites, {R['test_L2_pairs(U_ms)']:,} L2 pairs, "
         f"**{R['test_eval_pos(U_ms)']:,} evaluation positives** (human-native, in class, non-auto, U_ms).",
         "", "| quantity | TEST | CV0 | CV1 | CV2 | CV3 | CV4 |", "|---|---|---|---|---|---|---|"]
    for g in ["L2_labels", "U_ms_sites", "S_labels", "T_labels", "Y_labels", "clusters"]:
        L.append(f"| {g} | {R['frac_test:' + g]:.4f} | " + " | ".join(f"{R[f'frac_cv{f}:{g}']:.4f}" for f in range(5))
                 + " |")
    L += ["", "(CV fractions are of DEV.)", "",
          f"- Per-kinase band (kinases with >=5 evaluation positives: {R['kinases_ge5_eval_pos']}): "
          f"{R['kinases_ge5_in_15_25pct']} in 15-25%; {R['kinases_ge5_in_count_band(floor15%,ceil25%)']} within the "
          f"integer band [floor(0.15n), ceil(0.25n)]; kinases with >=20: {R['kinases_ge20_in_15_25pct']} in 15-25%.",
          f"- Per-kinase TEST fraction quantiles (0/10/50/90/100%): {R['per_kinase_test_frac_quantiles(ge5)']}.",
          f"- Worst kinases (TEST/total): {R['worst_kinases']} - their positives sit in 1-3 clusters.", "",
          "## Dropped from training/evaluation (A6)", "",
          f"- L2 pairs without a site-index row: {R['L2_pairs_not_in_index']}; L2 pairs in index but not U_ms: "
          f"{R['L2_pairs_in_index_not_U_ms']}; evaluation pairs not in U_ms: {R['eval_pairs_not_U_ms(dropped)']} "
          f"(of {R['eval_pairs']}). Used: L2 {R['L2_pairs_U_ms(used)']}, eval {R['eval_pairs_U_ms(used)']}.",
          f"- Label rows on proteins outside the index (kept only in sealed full table): "
          f"{R['label_rows_protein_not_in_index(dropped from DEV/TEST; kept in sealed full)']}.", "",
          "## Evaluable kinase sets (§4.6)", "",
          f"- |K_eval| = **{R['|K_eval|']}** (>=3 TEST positives on >=2 clusters), |K_any| = **{R['|K_any|']}**, "
          f"|K_cv| = **{R['|K_cv|']}** (>=5 DEV positives on >=2 clusters). Pseudokinases need >=3 labels.", "",
          "| group | K_eval | K_any | K_cv |", "|---|---|---|---|"]
    for g in sorted(set(R["K_any_by_group"]) | set(R["K_cv_by_group"])):
        L.append(f"| {g} | {R['K_eval_by_group'].get(g, 0)} | {R['K_any_by_group'].get(g, 0)} | "
                 f"{R['K_cv_by_group'].get(g, 0)} |")
    L += ["", f"By residue class: K_eval {R['K_eval_by_class']}, K_any {R['K_any_by_class']}, K_cv {R['K_cv_by_class']}.",
          "", "Training-count bins (DEV L2) of K_any: " +
          str(es[es.in_K_any].train_bin.value_counts().to_dict()) + "; of K_cv: " +
          str(es[es.in_K_cv].train_bin.value_counts().to_dict()) + ".", "",
          "## Leave-kinase-out", "",
          f"- {R['kinase_clusters']} kinase clusters (>=50% KinCore domain identity U Manning family); 5 folds: "
          f"kinases {R['lko_fold_kinases']}, DEV L2 labels {R['lko_fold_dev_L2']}.", "",
          "## Canaries", "",
          f"- v1 (used for the split freeze; provenance only): sizes {R['canary_sizes']}; TEST share "
          f"{R['canary_test_frac']}.",
          f"- **v{R.get('canary_version', 1)} (current)**: negatives drawn from real DEV-labelled sites "
          f"(protein-clustered); sizes {R.get('canary_v2_sizes', '-')}; TEST share {R.get('canary_v2_test_frac', '-')} "
          "(DEV-only pool, so negative canaries are CV-only; the positive canary keeps its v1 TEST part). "
          "Gate: mean negative AUROC on E-hard <= 0.57; positive >= 0.90.", "",
          "## Secondary TEST strata", ""]
    rl = R.get("relabel", {})
    L += [f"- New since 2021 (absent from 2021 PSP in all species): {rl.get('test_eval_pos_new_since_2021', R.get('test_eval_pos_new_since_2021'))} TEST positives.",
          f"- External curated (DL-5 + A9): {rl.get('test_ext_pairs_U_ms', '?')} TEST pairs in U_ms, tier A "
          f"{rl.get('test_ext_tierA_U_ms', '?')}.",
          f"- Bio controls in TEST: {R.get('bio_controls_in_test', '?')}.", "",
          "## Seal", "", "TEST labels: `sealed/test_labels.parquet` (only ksp.eval_test reads it). SHA-256:", ""]
    for n, m in seal["files"].items():
        L.append(f"- {n}: `{m['sha256']}`")
    (out / "split_report.md").write_text("\n".join(L) + "\n")
    print("wrote", out / "split_report.md")


def annotate_full(full: pd.DataFrame, ums: set, in_idx: set) -> pd.DataFrame:
    full = full.copy()
    full["in_index"] = full.site_id.isin(in_idx)
    full["in_U_ms"] = full.site_id.isin(ums)
    full["is_canary"] = False
    full["new_since_2021"] = labels.new_since_2021(full)
    par = labels.paralog_sets()
    keep = labels.external_test_filter(full, par)
    full["ext_tier"] = labels.ext_test_tiers(full, keep)
    full["ext_test"] = full.ext_tier.ne("")
    return full


def relabel():
    """Rebuild label tables for the frozen split (after a label-build change) and re-seal."""
    t0 = time.time()
    out = splits.SPLIT_DIR
    if not labels.FULL_PAIRS.exists():
        sys.exit(f"{labels.FULL_PAIRS} missing - rerun scripts/build_labels.py first")
    S = index.load_sites()
    P = set(index.load_proteins().acc)
    ums, in_idx = set(S.site_id[S.in_U_ms]), set(S.site_id)
    test_prot = splits.load_test_proteins()
    full = annotate_full(pd.read_parquet(labels.FULL_PAIRS), ums, in_idx)
    # canaries: deterministic; must equal the frozen DEV canaries
    L2 = labels.select_variant(full, "L2")
    n2 = L2.groupby("kinase_acc").site_id.nunique()
    EVp = labels.select_variant(full, "eval")
    dev_pool = EVp[EVp.in_U_ms & ~EVp.acc.isin(test_prot)][["site_id", "acc", "res"]]
    size_src = tuple(sorted(int(x) for x in n2[(n2 >= 10) & (n2 <= 100)]))
    can, _ = canary.generate_all(size_src, dev_pool)
    v1, _ = canary.generate_v1(size_src)
    can_t = can.acc.isin(test_prot)
    cdir = labels.DEV_PAIRS.parent
    if not (cdir / canary.CANARY_FILE_DEV_V1).exists():          # keep v1 for provenance
        old = pd.read_parquet(cdir / canary.CANARY_FILE_DEV)
        if "canary_version" not in old or (old.canary_version != canary.CANARY_VERSION).all():
            shutil.copy(cdir / canary.CANARY_FILE_DEV, cdir / canary.CANARY_FILE_DEV_V1)
    pos_new = set(can[can.kinase_acc.eq(canary.POS_ID)].site_id)
    pos_v1 = set(v1[v1.kinase_acc.eq(canary.POS_ID)].site_id)
    if pos_new != pos_v1:
        sys.exit("positive canary changed - not allowed in --relabel")
    if set(can.kinase_acc) != set(canary.CANARY_IDS):
        sys.exit("canary set incomplete")
    inp = full.acc.isin(P)
    dev = full[inp & ~full.acc.isin(test_prot)].reset_index(drop=True)
    tst = pd.concat([full[inp & full.acc.isin(test_prot)], can[can_t]], ignore_index=True)
    # verify evaluation positives unchanged vs frozen kinase sets
    es = splits.load_eval_kinase_sets()
    EV = labels.select_variant(full, "eval")
    EV = EV[EV.in_U_ms]
    nt = EV[EV.acc.isin(test_prot)].groupby("kinase_acc").site_id.nunique()
    nd = EV[~EV.acc.isin(test_prot)].groupby("kinase_acc").site_id.nunique()
    d_t = (es.kinase_acc.map(nt).fillna(0) != es.n_test_pos).sum()
    d_d = (es.kinase_acc.map(nd).fillna(0) != es.n_dev_pos).sum()
    if d_t or d_d:
        sys.exit(f"evaluation positives changed for {d_t} (TEST) / {d_d} (DEV) kinases - kinase sets are frozen")
    dev.to_parquet(labels.DEV_PAIRS, index=False)
    can[~can_t].reset_index(drop=True).to_parquet(labels.DEV_PAIRS.parent / canary.CANARY_FILE_DEV, index=False)
    tst.to_parquet(SEALED / "test_labels.parquet", index=False)
    shutil.move(str(labels.FULL_PAIRS), str(SEALED / "ks_pairs_full.parquet"))
    um = labels.FULL_PAIRS.parent / "ks_unmapped.parquet"
    if um.exists():
        shutil.move(str(um), str(SEALED / "ks_unmapped.parquet"))
    kc = ROOT / "results/labels/kinase_label_counts.tsv"
    if kc.exists():
        shutil.move(str(kc), str(SEALED / "kinase_label_counts.tsv"))
    seal = json.loads((out / "SEAL.json").read_text())
    hist = seal.setdefault("history", [])
    hist.append({"relabel": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "previous_files": seal["files"]})
    seal["files"] = {f.name: {"sha256": sha(f), "bytes": f.stat().st_size} for f in sorted(SEALED.glob("*"))
                     if f.is_file()}
    seal["public"][str(labels.DEV_PAIRS.relative_to(ROOT))] = sha(labels.DEV_PAIRS)
    cdf = labels.DEV_PAIRS.parent / canary.CANARY_FILE_DEV
    seal["public"][str(cdf.relative_to(ROOT))] = sha(cdf)
    seal["test_label_rows"] = len(tst)
    (out / "SEAL.json").write_text(json.dumps(seal, indent=1))
    R = json.loads((out / "split_numbers.json").read_text())
    cn = can.groupby("kinase_acc").site_id.nunique()
    ct = can[can_t].groupby("kinase_acc").site_id.nunique().reindex(cn.index).fillna(0)
    R["canary_version"] = canary.CANARY_VERSION
    R["canary_v2_sizes"] = cn.astype(int).to_dict()
    R["canary_v2_test_frac"] = (ct / cn).round(3).to_dict()
    R["relabel"] = {"when": hist[-1]["relabel"], "dev_label_rows": len(dev), "test_label_rows": len(tst),
                    "mask_only_rows_dev": int(dev.mask_only.sum()), "mask_only_rows_test": int(tst.mask_only.fillna(False).astype(bool).sum()),
                    "test_eval_pos_new_since_2021": int((labels.evaluable_positive(tst) & tst.new_since_2021.astype(bool)
                                                         & tst.in_U_ms.astype(bool)).sum()),
                    "test_ext_pairs_U_ms": int((tst.ext_test.fillna(False).astype(bool) & tst.in_U_ms.astype(bool)).sum()),
                    "test_ext_tierA_U_ms": int((tst.ext_tier.eq("A") & tst.in_U_ms.astype(bool)).sum()),
                    "runtime_s": round(time.time() - t0, 1)}
    (out / "split_numbers.json").write_text(json.dumps(R, indent=1, default=str))
    print(json.dumps(R["relabel"], indent=1))


if __name__ == "__main__":
    main()
