"""G5 release: fuse the deployment R2 and R3 U_obs scores (same recipe as CV/TEST), attach calibrated probabilities
and A17 assignment probabilities, compute coverage (A17 / §6.6, A31b) against random and Cantley-only baselines, and
write the release files.

Fusion   : z-score each model's raw logits over all its U_obs rows; score = 0.6 z(R2) + 0.4 z(R3);
           score_site = 0.5 z(R2) + 0.5 z(R3)   (weights chosen on DEV CV, r2_fusion_r3__20261026-191901__s0).
prob     : per-kinase tiered calibrator (a_k + b_k * score) fitted on DEV OOF of the CV fusion run.
assign   : softmax over real kinases of matching class at each site of the calibrated logits, T = mean CV fold
           temperature; top kinase -> residue-stratified (ST / Y) isotonic map fitted on all labelled DEV OOF sites;
           other kinases share the remaining mass.  LIMITATION: all maps are fitted on DEV OOF; the deployment model is
           trained on DEV+TEST, for which no OOF predictions exist.
pct      : within-kinase percentile of score over the kinase's U_obs candidate sites (100 = best).
Coverage : Coverage@tau = #U_obs sites whose top assignment_prob >= tau (tau 0.5, 0.8), over all sites and over
           sites with no PSP kinase (no PSP KS row in any species, DEV+TEST); #sites with >= 1 kinase in that kinase's
           top 1% (ceil(0.01 * n_candidates)).
           Baselines: random scores (same calibrators would be meaningless -> raw softmax of N(0,1) logits at the same
           T, and random top-1%); Cantley-only = B1 imputed-universe percentile with ITS OWN DEV-OOF calibrators /
           temperature / ST-Y isotonic map (r2_calib.py on simple_b1_imputed_universe), over the kinases B1 was
           calibrated on.
Outputs  : cache/release/<run_id>/release_scores.parquet (kinase_idx, site_idx, score, score_site, prob, pct,
           assignment_prob; kinase-chunked), results/release/{coverage.tsv, kinase_top200.tsv.gz, site_top5.tsv.gz,
           RELEASE_MANIFEST.json, README_release.md}.
Usage: r2_release.py --r2 <r2_deploy run> --r3 <R3 deploy parquet> [--smoke]"""
import argparse, json, math, time
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from sklearn.isotonic import IsotonicRegression
from ksp import index, labels, manifest, store
from ksp.index import ROOT
from ksp.models.r2 import calib

ap = argparse.ArgumentParser()
ap.add_argument("--r2", required=True); ap.add_argument("--r3", required=True)
ap.add_argument("--cv-run", default="r2_fusion_r3__20261026-191901__s0")
ap.add_argument("--b1-run", default="simple_b1_imputed_universe__20261026-150833__s0")
ap.add_argument("--w", type=float, default=0.6); ap.add_argument("--w-site", type=float, default=0.5)
ap.add_argument("--smoke", action="store_true")
a = ap.parse_args()
t_start = time.time()
S = index.load_sites(); K = index.load_kinases()
nS, nK = len(S), len(K)
res_is_y = S.res.values == "Y"
OUTR = ROOT / "results" / ("release_smoke" if a.smoke else "release")
OUTR.mkdir(parents=True, exist_ok=True)


def maps(run):
    d = ROOT / "results" / "models" / "r2" / "calib" / f"{run}__calibrated"
    C = pd.read_csv(d / "calibrators.tsv", sep="\t", index_col=0)
    T = float(np.mean(list(json.load(open(d / "assignment.json"))["temperature_by_fold"].values())))
    A = pd.read_parquet(ROOT / "cache" / "r2" / f"assign_{run}__calibrated.parquet")
    ry = S.res.values[A.site_idx.values] == "Y"
    iso = {}
    for r, m in (("ST", ~ry), ("Y", ry)):
        L = A[A.labelled.values & m]
        iso[r] = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(L.p_top.values, L.top_is_true.values.astype(float))
    return C, T, iso


def load_sorted(path):
    t = pq.read_table(path, columns=["kinase_idx", "site_idx", "score"])
    k = t.column("kinase_idx").to_numpy().astype(np.int32); s = t.column("site_idx").to_numpy().astype(np.int32)
    v = t.column("score").to_numpy().astype(np.float32)
    o = np.lexsort((s, k))
    return k[o], s[o], v[o]


def assignment(k, s, lc, T, iso, kin_ok=None):
    """Returns per-row assignment_prob (aligned to input rows) and per-site arrays (site, top kinase, p_top_iso)."""
    if kin_ok is not None:
        keep = kin_ok[k]
    else:
        keep = np.ones(len(k), bool)
    idx = np.flatnonzero(keep)
    o = idx[np.lexsort((k[idx], s[idx]))]
    ss = s[o]
    p = calib.site_softmax(ss, lc[o], T)
    b = np.flatnonzero(np.r_[True, ss[1:] != ss[:-1]])
    cnt = np.diff(np.r_[b, len(ss)])
    ptop = np.maximum.reduceat(p, b)
    ry = res_is_y[ss[b]]
    q = np.empty(len(b))
    for r, m in (("ST", ~ry), ("Y", ry)):
        if m.any():
            q[m] = iso[r].predict(ptop[m])
    ptr, qr = np.repeat(ptop, cnt), np.repeat(q, cnt)
    top = p >= ptr
    with np.errstate(divide="ignore", invalid="ignore"):
        ap_ = np.where(top, qr, np.where(ptr < 1, p * (1 - qr) / (1 - ptr), 0.0))
    out = np.full(len(k), np.nan, np.float32)
    out[o] = np.clip(ap_, 0, 1)
    top_k = k[o][b + np.array([np.argmax(p[x:x + c]) for x, c in zip(b, cnt)])]
    return out, ss[b], top_k, q


def top1pct_sites(k, s, v):
    """Sites with >= 1 kinase ranking them in that kinase's top 1% (rows sorted by kinase)."""
    hit = np.zeros(nS, bool)
    bnd = np.flatnonzero(np.r_[True, k[1:] != k[:-1], True])
    for lo, hi in zip(bnd[:-1], bnd[1:]):
        n = hi - lo
        m = max(1, math.ceil(0.01 * n))
        sel = np.argpartition(-v[lo:hi], m - 1)[:m]
        hit[s[lo:hi][sel]] = True
    return hit


def within_pct(k, v):
    out = np.empty(len(v), np.float32)
    bnd = np.flatnonzero(np.r_[True, k[1:] != k[:-1], True])
    for lo, hi in zip(bnd[:-1], bnd[1:]):
        r = pd.Series(v[lo:hi]).rank(method="average").values
        out[lo:hi] = (100.0 * r / (hi - lo)).astype(np.float32)
    return out


r2m = json.load(open(ROOT / "results" / "runs" / a.r2 / "manifest.json"))
r2_path = ROOT / r2m["outputs"]["r2_uobs_scores"]["path"]
with manifest.RunManifest("r2r3_release" + ("_smoke" if a.smoke else ""), config=vars(a), seed=0,
                          split="deploy" if not a.smoke else "smoke", inputs=[r2_path, a.r3],
                          notes="G5 release: fusion + calibration + assignment + coverage") as run:
    out_dir = ROOT / "cache" / "release" / run.run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    k, s, v2 = load_sorted(r2_path)
    k3, s3, v3 = load_sorted(a.r3)
    if not (np.array_equal(k, k3) and np.array_equal(s, s3)):
        raise ValueError(f"R2/R3 deployment rows not aligned: {len(k)} vs {len(k3)}")
    if (k >= nK).any():
        raise ValueError("canary/unknown kinase rows present")
    del k3, s3
    z2 = (v2 - v2.mean(dtype=np.float64)) / v2.std(dtype=np.float64)
    z3 = (v3 - v3.mean(dtype=np.float64)) / v3.std(dtype=np.float64)
    score = (a.w * z2 + (1 - a.w) * z3).astype(np.float32)
    score_site = (a.w_site * z2 + (1 - a.w_site) * z3).astype(np.float32)
    del z2, z3, v3
    C, T, iso = maps(a.cv_run)
    lc = C.a.values[k] + C.b.values[k] * score.astype(np.float64)
    prob = (1 / (1 + np.exp(-lc))).astype(np.float32)
    ap_, site_u, top_k, q = assignment(k, s, lc, T, iso)
    del lc
    pct = within_pct(k, score)
    hit = top1pct_sites(k, s, score)
    t_fuse = time.time() - t0
    # ------------------------------------------------------------------ release scores (kinase-chunked parquet)
    f = out_dir / "release_scores.parquet"
    schema = pa.schema([("kinase_idx", pa.int16()), ("site_idx", pa.int32()), ("score", pa.float32()),
                        ("score_site", pa.float32()), ("prob", pa.float32()), ("pct", pa.float32()),
                        ("assignment_prob", pa.float32())])
    bnd = np.flatnonzero(np.r_[True, k[1:] != k[:-1], True])
    with pq.ParquetWriter(f, schema, compression="zstd") as w:
        for lo, hi in zip(bnd[:-1], bnd[1:]):
            w.write_table(pa.table({"kinase_idx": k[lo:hi].astype(np.int16), "site_idx": s[lo:hi], "score": score[lo:hi],
                                    "score_site": score_site[lo:hi], "prob": prob[lo:hi], "pct": pct[lo:hi],
                                    "assignment_prob": ap_[lo:hi]}, schema=schema))
    run.m["outputs"]["release_scores"] = {"path": str(f.relative_to(ROOT)), "sha256": manifest.sha256_file(f)}
    run.m["n_rows"]["release_scores"] = int(len(k))
    # ------------------------------------------------------------------ coverage
    if a.smoke:
        lab_all = labels.load_dev("all")
    else:
        from ksp.deploy import load_deploy_labels
        lab_all = load_deploy_labels("all")
    psp_sites = set(lab_all.site_id[lab_all.provenance.eq("PSP")])
    no_psp = ~S.site_id.isin(psp_sites).values
    rows = []

    def cov(name, site_u, q, hit, note):
        r = {"method": name, "n_sites": int(nS), "n_sites_no_psp": int(no_psp.sum()), "note": note}
        for tau in (0.5, 0.8):
            m = q >= tau
            r[f"coverage@{tau}"] = int(m.sum())
            r[f"coverage@{tau}_no_psp"] = int((m & no_psp[site_u]).sum())
        r["sites_in_some_kinase_top1pct"] = int(hit.sum())
        r["sites_in_some_kinase_top1pct_no_psp"] = int((hit & no_psp).sum())
        rows.append(r)
    cov("R2+R3 final (deployment)", site_u, q, hit, "assignment via DEV-OOF calibrators + ST/Y isotonic; 503 kinases")
    rng = np.random.default_rng(0)
    rv = rng.standard_normal(len(k))
    pr = calib.site_softmax(s[np.lexsort((k, s))], rv, T)
    ss = np.sort(s, kind="stable")
    bb = np.flatnonzero(np.r_[True, ss[1:] != ss[:-1]])
    cov("random scores", ss[bb], np.maximum.reduceat(pr, bb), top1pct_sites(k, s, rv.astype(np.float32)),
        "raw softmax of N(0,1) logits at the same T (no isotonic map); random top-1%")
    del pr, rv
    # Cantley-only (B1 imputed-universe percentile, own DEV-OOF calibration recipe)
    own = store.load_matrix("motif/cantley_unprimed_pct_universe")[0]
    imp = store.load_matrix("motif/cantley_imputed_unprimed_pct_universe")[0]
    cv_ = np.full(len(k), np.nan, np.float32)
    for lo, hi in zip(bnd[:-1], bnd[1:]):
        kk = int(k[lo])
        x = np.asarray(own[kk][s[lo:hi]], np.float32); y = np.asarray(imp[kk][s[lo:hi]], np.float32)
        cv_[lo:hi] = np.where(np.isnan(x), y, x)
    Cb, Tb, isob = maps(a.b1_run)
    kin_ok = np.zeros(nK, bool); kin_ok[Cb.index.values] = True
    okrow = kin_ok[k] & np.isfinite(cv_)
    lcb = np.full(len(k), -np.inf)
    ia = np.zeros(nK); ib = np.zeros(nK); ia[Cb.index.values] = Cb.a.values; ib[Cb.index.values] = Cb.b.values
    lcb[okrow] = ia[k[okrow]] + ib[k[okrow]] * cv_[okrow]
    _, site_b, _, qb = assignment(k[okrow], s[okrow], lcb[okrow], Tb, isob)
    hb = top1pct_sites(k[np.isfinite(cv_)], s[np.isfinite(cv_)], cv_[np.isfinite(cv_)])
    cov("Cantley-only (B1 imputed-universe pct)", site_b, qb, hb,
        f"assignment over the {int(kin_ok.sum())} kinases B1 was calibrated on (own DEV-OOF calibrators/T/isotonic); "
        f"top-1% over all kinases with an own or imputed matrix")
    CT = pd.DataFrame(rows)
    CT.to_csv(OUTR / "coverage.tsv", sep="\t", index=False)
    print(CT.to_string(), flush=True)
    # ------------------------------------------------------------------ small release files
    genes = K.gene.values; acc = K.uniprot_acc.values; sid = S.site_id.values
    top200 = []
    for lo, hi in zip(bnd[:-1], bnd[1:]):
        m = min(200, hi - lo)
        sel = lo + np.argsort(-score[lo:hi], kind="stable")[:m]
        top200.append(pd.DataFrame({"kinase_acc": acc[k[sel]], "gene": genes[k[sel]], "rank": np.arange(1, m + 1),
                                    "site_id": sid[s[sel]], "score": score[sel], "prob": prob[sel], "pct": pct[sel],
                                    "assignment_prob": ap_[sel]}))
    T200 = pd.concat(top200, ignore_index=True)
    f1 = OUTR / "kinase_top200.tsv.gz"
    T200.round(6).to_csv(f1, sep="\t", index=False)
    o = np.lexsort((-ap_, s))
    so = s[o]
    b = np.flatnonzero(np.r_[True, so[1:] != so[:-1]])
    rk = np.arange(len(so)) - np.repeat(b, np.diff(np.r_[b, len(so)]))
    sel = o[rk < 5]
    T5 = pd.DataFrame({"site_id": sid[s[sel]], "rank": (rk[rk < 5] + 1), "kinase_acc": acc[k[sel]], "gene": genes[k[sel]],
                       "assignment_prob": ap_[sel], "prob": prob[sel], "pct": pct[sel], "score_site": score_site[sel]})
    f2 = OUTR / "site_top5.tsv.gz"
    T5.round(6).to_csv(f2, sep="\t", index=False)
    rel = {"run_id": run.run_id, "r2_deploy_run": a.r2, "r3_deploy_predictions": a.r3,
           "release_scores": {"path": str(f.relative_to(ROOT)), "sha256": run.m["outputs"]["release_scores"]["sha256"],
                              "rows": int(len(k)), "bytes": f.stat().st_size},
           "kinase_top200": {"path": str(f1.relative_to(ROOT)), "sha256": manifest.sha256_file(f1), "bytes": f1.stat().st_size},
           "site_top5": {"path": str(f2.relative_to(ROOT)), "sha256": manifest.sha256_file(f2), "bytes": f2.stat().st_size},
           "coverage": str((OUTR / "coverage.tsv").relative_to(ROOT)),
           "licence": "Derived from PhosphoSitePlus (CC BY-NC-SA 3.0) labels: NON-COMMERCIAL use only, share-alike, attribution to PhosphoSitePlus.",
           "limitations": ["calibrators, temperature and ST/Y isotonic maps were fitted on DEV OOF predictions of the "
                           "DEV-trained CV models; the deployment model is trained on DEV+TEST and has no OOF",
                           "probabilities estimate P(curated PSP-type positive | score) under PU labels: lower bounds",
                           "assignment_prob is conditional on the site being a substrate of one of the 503 modelled kinases; "
                           "no reliable >= 0.8 tier and no Y-site assignment reaches 0.5 (audit r2)"],
           "index_files": ["cache/index/kinases.parquet (kinase_idx)", "cache/index/sites.parquet (site_idx)"],
           "runtime_s": {"fusion_calibration_assignment": round(t_fuse, 1), "total": round(time.time() - t_start, 1)}}
    (OUTR / "RELEASE_MANIFEST.json").write_text(json.dumps(rel, indent=1))
    run.set("release", rel)
    print("RUN", run.run_id, json.dumps({k_: rel[k_] for k_ in ["release_scores", "kinase_top200", "site_top5"]}), flush=True)
