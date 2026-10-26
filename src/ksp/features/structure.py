"""AlphaFold-DB derived per-residue structural features (DESIGN §7 F0), computed for every residue of every
Swiss-Prot canonical human protein and saved for S/T/Y residues as cache/structure/site_structure.parquet
(all residues: cache/structure/residue_structure.parquet).

Input (never extracted to disk, streamed): external/alphafold/UP000005640_9606_HUMAN_v6.tar (AFDB v6 human
proteome; members AF-<acc>-F<n>-model_v6.cif.gz) and optional PAE files external/alphafold/pae/*.json.gz
(scripts/download_afdb_pae.py).  mmCIF parsing is our own minimal loop parser (atom_site, struct_conf,
ma_target_ref_db_details); no downloaded code is executed.

Features per residue (re-implementation of structuremap 0.0.10 / Bludau et al. 2022, PLoS Biol 20:e3001636):
  plddt                  AF pLDDT (B-factor of CA)
  plddt_s5, plddt_s10    mean pLDDT over [i-5, i+5] / [i-10, i+10] (structuremap smooth_score, truncated at ends)
  pse_{12,24}_{70,180}_pae   pPSE: number of residues j != i with PAE(i,j) <= r and d(CAi,CAj) + PAE(i,j) <= r and
                         angle(CAi->CBi, CAi->CAj) <= theta (Gly: pseudo-CB as in structuremap get_gly_vector).
                         theta=180 is the full sphere. NA when no PAE file (fragmented proteins >2700 aa).
  pse_*_nopae            the same with PAE = 0 (available for all modelled residues, incl. fragments)
  pse_24_180_pae_s10     smoothed (±10), and IDR flag `idr_sm` = pse_24_180_pae_s10 <= 34.27 (structuremap tutorial
                         threshold); `high_acc_5` = pse_12_70_pae <= 5 (structuremap 'high accessibility').
  idr_nopae              pse_24_180_nopae_s10 <= T_nopae (T_nopae calibrated to maximise agreement with idr_sm on
                         proteins having PAE; used as fallback in `idr` where idr_sm is NA; source in idr_source)
  disorder_plddt         plddt_s10 < 50
  sasa, rsa              Shrake-Rupley SASA (own vectorised implementation: 100 golden-spiral points per atom,
                         probe 1.4 A, Bondi radii C1.70 N1.55 O1.52 S1.80, occlusion test against the 6 nearest
                         atom centres) over heavy atoms, relative to Tien et al. 2013 theoretical max ASA.
  ss8 / ss3              DSSP secondary structure as provided by AFDB in the mmCIF _struct_conf block
                         (H alpha, G 3-10, I pi, P PPII, E strand, T turn, S bend, C none) and 3-state H/E/C.
  af_frag, af_n_frags    AFDB fragment used (proteins > 2700 aa are split into 1400-aa fragments, step 200; each
                         residue is taken from the fragment in which it is most central)
  af_status              ok | not_in_afdb | seq_mismatch  (features NA unless ok)
Usage: python -I src/ksp/features/structure.py --procs 14
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import tarfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from ksp.features.esm_embed import canonical_swissprot  # noqa: E402

TAR = ROOT / "external/alphafold/UP000005640_9606_HUMAN_v6.tar"
PAE_DIR = ROOT / "external/alphafold/pae"
OUT = ROOT / "cache/structure"
NAME_RE = re.compile(r"AF-([A-Z0-9]+)-F(\d+)-model_v(\d+)\.cif\.gz$")
SS8 = {"HELX_RH_AL_P": "H", "HELX_RH_3T_P": "G", "HELX_RH_PI_P": "I", "HELX_LH_PP_P": "P", "STRN": "E",
       "TURN_TY1_P": "T", "BEND": "S"}
SS3 = {"H": "H", "G": "H", "I": "H", "E": "E"}
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90}
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201, K=236, M=224, F=240,
              P=159, S=155, T=172, W=285, Y=263, V=174)
PSE_GRID = [(12, 70), (12, 180), (24, 70), (24, 180)]
IDR_SM_T = 34.27
PROBE = 1.4


# ----------------------------------------------------------------------------------------------- mmCIF parsing
def _loops(lines):
    """Yield (category, fields, rows) for loop_ blocks; rows are whitespace-split (sufficient for AFDB atom_site
    and struct_conf, which contain no quoted values with spaces)."""
    i, n = 0, len(lines)
    while i < n:
        if lines[i].startswith("loop_"):
            i += 1
            fields = []
            while i < n and lines[i].startswith("_"):
                fields.append(lines[i].strip())
                i += 1
            rows = []
            while i < n and not (lines[i].startswith("#") or lines[i].startswith("loop_") or lines[i].startswith("_")):
                if lines[i].strip():
                    rows.append(lines[i].split())
                i += 1
            cat = fields[0].split(".")[0] if fields else ""
            yield cat, [f.split(".", 1)[1] for f in fields], rows
        else:
            i += 1


def parse_cif(text: str) -> dict:
    lines = text.splitlines()
    kv = {}
    for ln in lines:
        if ln.startswith("_ma_target_ref_db_details.") or ln.startswith("_ma_qa_metric_global."):
            p = ln.split(None, 1)
            if len(p) == 2:
                kv[p[0]] = p[1].strip()
    atoms, ss = None, []
    for cat, fields, rows in _loops(lines):
        if cat == "_atom_site":
            atoms = (fields, rows)
        elif cat == "_struct_conf":
            fi = {f: k for k, f in enumerate(fields)}
            ss = [(int(r[fi["beg_label_seq_id"]]), int(r[fi["end_label_seq_id"]]), r[fi["conf_type_id"]]) for r in rows]
    fields, rows = atoms
    fi = {f: k for k, f in enumerate(fields)}
    seq_id = np.array([int(r[fi["label_seq_id"]]) for r in rows])
    atom = np.array([r[fi["label_atom_id"]] for r in rows])
    elem = np.array([r[fi["type_symbol"]] for r in rows])
    xyz = np.array([[float(r[fi["Cartn_x"]]), float(r[fi["Cartn_y"]]), float(r[fi["Cartn_z"]])] for r in rows])
    b = np.array([float(r[fi["B_iso_or_equiv"]]) for r in rows])
    unp_num = np.array([int(r[fi["pdbx_sifts_xref_db_num"]]) if r[fi["pdbx_sifts_xref_db_num"]] not in "?." else -1
                        for r in rows])
    unp_res = np.array([r[fi["pdbx_sifts_xref_db_res"]] for r in rows])
    return dict(kv=kv, ss=ss, seq_id=seq_id, atom=atom, elem=elem, xyz=xyz, b=b, unp_num=unp_num, unp_res=unp_res)


# ----------------------------------------------------------------------------------------------- geometry
def _rotation(axis, theta_deg):
    """structuremap rotate_vector_around_axis (Euler-Rodrigues), vectorised over axes (n,3)."""
    th = np.radians(theta_deg)
    axis = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    a = np.cos(th / 2.0)
    bcd = -axis * np.sin(th / 2.0)
    b, c, d = bcd[:, 0], bcd[:, 1], bcd[:, 2]
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    R = np.stack([np.stack([aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)], -1),
                  np.stack([2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)], -1),
                  np.stack([2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc], -1)], 1)
    return R


def side_chain_vectors(ca, cb, c, n):
    """Unit CA->CB vectors; for Gly (NaN CB) the structuremap pseudo vector (N rotated -120 deg about CA->C)."""
    u = cb - ca
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    g = np.isnan(cb[:, 0])
    if g.any():
        un = (n[g] - ca[g]) / np.linalg.norm(n[g] - ca[g], axis=1, keepdims=True)
        uc = (c[g] - ca[g]) / np.linalg.norm(c[g] - ca[g], axis=1, keepdims=True)
        R = _rotation(uc, -120.0)
        u[g] = np.einsum("nij,nj->ni", R, un)
    return u


def pse_counts(ca, u, pae=None):
    """pPSE counts for PSE_GRID. ca (L,3), u (L,3) unit side-chain vectors, pae (L,L) or None."""
    L = len(ca)
    D = np.sqrt(((ca[:, None, :] - ca[None, :, :]) ** 2).sum(-1)).astype(np.float32)
    G = (u @ ca.T) - (u * ca).sum(1)[:, None]  # u_i . (ca_j - ca_i)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosang = (G / D).astype(np.float32)
        ang = np.degrees(np.arccos(cosang))  # NaN where |cos|>1 numerically or i==j, as in structuremap
    eye = np.eye(L, dtype=bool)
    P = np.zeros_like(D) if pae is None else pae.astype(np.float32)
    out = {}
    for r, th in PSE_GRID:
        m = (P <= r) & (D + P <= r) & (ang <= th) & ~eye
        out[(r, th)] = m.sum(1).astype(np.int16)
    return out


_SPHERE = None


def _sphere(n=100):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], 1)


def sasa_atoms(xyz, radii, n_points=100, k=6):
    global _SPHERE
    if _SPHERE is None or len(_SPHERE) != n_points:
        _SPHERE = _sphere(n_points)
    R = radii + PROBE
    pts = (xyz[:, None, :] + R[:, None, None] * _SPHERE[None]).reshape(-1, 3)
    owner = np.repeat(np.arange(len(xyz)), n_points)
    tree = cKDTree(xyz)
    kk = min(k + 1, len(xyz))
    d, j = tree.query(pts, k=kk)
    if kk == 1:
        d, j = d[:, None], j[:, None]
    buried = ((d < R[j] - 1e-6) & (j != owner[:, None])).any(1)
    exposed = (~buried).reshape(len(xyz), n_points).sum(1)
    return 4 * np.pi * R ** 2 * exposed / n_points


# ----------------------------------------------------------------------------------------------- per fragment
def process_member(args):
    name, raw = args
    m = NAME_RE.search(name)
    acc, frag, ver = m.group(1), int(m.group(2)), int(m.group(3))
    try:
        d = parse_cif(gzip.decompress(raw).decode())
    except Exception as e:  # malformed file -> reported, NA
        return dict(acc=acc, frag=frag, error=f"parse: {e!r}")
    seq_id, atom = d["seq_id"], d["atom"]
    res_ids, first = np.unique(seq_id, return_index=True)
    L = len(res_ids)
    idx = np.searchsorted(res_ids, seq_id)
    unp_num = d["unp_num"][first]
    unp_res = d["unp_res"][first]
    plddt = np.full(L, np.nan)
    coords = {a: np.full((L, 3), np.nan) for a in ("CA", "CB", "C", "N")}
    for a in coords:
        s = atom == a
        coords[a][idx[s]] = d["xyz"][s]
        if a == "CA":
            plddt[idx[s]] = d["b"][s]
    ss8 = np.array(["C"] * L)
    for beg, end, t in d["ss"]:
        ss8[(res_ids >= beg) & (res_ids <= end)] = SS8.get(t, "C")
    u = side_chain_vectors(coords["CA"], coords["CB"], coords["C"], coords["N"])
    pae = None
    pf = PAE_DIR / f"AF-{acc}-F{frag}-predicted_aligned_error_v{ver}.json.gz"
    pae_status = "none"
    if pf.exists():
        try:
            pj = json.loads(gzip.decompress(pf.read_bytes()))
            pj = pj[0] if isinstance(pj, list) else pj
            pae = np.asarray(pj["predicted_aligned_error"], dtype=np.float32)
            pae_status = "ok" if pae.shape == (L, L) else "shape_mismatch"
            if pae_status != "ok":
                pae = None
        except Exception:
            pae, pae_status = None, "unreadable"
    out = dict(acc=acc, frag=frag, ver=ver, L=L, pos=unp_num, aa=unp_res, plddt=plddt, ss8=ss8,
               align_begin=d["kv"].get("_ma_target_ref_db_details.seq_db_align_begin"),
               align_end=d["kv"].get("_ma_target_ref_db_details.seq_db_align_end"),
               seq_date=d["kv"].get("_ma_target_ref_db_details.seq_db_sequence_version_date"),
               pae_status=pae_status)
    for (r, th), v in pse_counts(coords["CA"], u, None).items():
        out[f"pse_{r}_{th}_nopae"] = v
    if pae is not None:
        for (r, th), v in pse_counts(coords["CA"], u, pae).items():
            out[f"pse_{r}_{th}_pae"] = v
    rad = np.array([RADII.get(e, 1.8) for e in d["elem"]])
    sa = sasa_atoms(d["xyz"], rad)
    sasa = np.bincount(idx, weights=sa, minlength=L)
    out["sasa"] = sasa
    out["rsa"] = sasa / np.array([MAXASA.get(a, np.nan) for a in unp_res], dtype=float)
    return out


def smooth(x, hw):
    """structuremap smooth_score: mean over [i-hw, i+hw] truncated at the ends (NaN-aware)."""
    x = np.asarray(x, float)
    ok = ~np.isnan(x)
    c = np.concatenate([[0], np.cumsum(np.where(ok, x, 0))])
    n = np.concatenate([[0], np.cumsum(ok)])
    i = np.arange(len(x))
    lo, hi = np.maximum(i - hw, 0), np.minimum(i + hw + 1, len(x))
    with np.errstate(invalid="ignore"):
        return (c[hi] - c[lo]) / (n[hi] - n[lo])


# ----------------------------------------------------------------------------------------------- assembly
FEATS_NUM = (["plddt", "sasa", "rsa"] + [f"pse_{r}_{t}_{p}" for p in ("pae", "nopae") for r, t in PSE_GRID])


def assemble(acc, seq, frags):
    """Combine fragment results into a per-residue table for one canonical protein; returns (df, status)."""
    L = len(seq)
    frags = sorted(frags, key=lambda f: f["frag"])
    for f in frags:
        if "error" in f:
            return None, "parse_error"
        p = f["pos"]
        if (p < 1).any() or (p > L).any():
            return None, "seq_mismatch"
        can = np.array(list(seq))[p - 1]
        ok = (can == f["aa"]) | ((can == "U") & np.isin(f["aa"], ["C", "U", "X"]))
        if not ok.all():
            return None, "seq_mismatch"
    covered = np.zeros(L, bool)
    best = np.full(L, -1)
    bestc = np.full(L, -1)
    for k, f in enumerate(frags):
        p = f["pos"] - 1
        covered[p] = True
        cen = np.minimum(p - p.min(), p.max() - p)
        upd = cen > bestc[p]
        best[p[upd]] = k
        bestc[p[upd]] = cen[upd]
    if not covered.all():
        return None, "incomplete"
    cols = {c: np.full(L, np.nan) for c in FEATS_NUM}
    ss8 = np.array(["C"] * L, dtype=object)
    for k, f in enumerate(frags):
        p = f["pos"] - 1
        sel = best[p] == k
        for c in FEATS_NUM:
            if c in f:
                cols[c][p[sel]] = f[c][sel]
        ss8[p[sel]] = f["ss8"][sel]
    df = pd.DataFrame(cols)
    df.insert(0, "pos", np.arange(1, L + 1))
    df.insert(1, "res", list(seq))
    df.insert(0, "acc", acc)
    df["plddt_s5"] = smooth(df.plddt, 5)
    df["plddt_s10"] = smooth(df.plddt, 10)
    df["pse_24_180_pae_s10"] = smooth(df.pse_24_180_pae, 10)
    df["pse_24_180_nopae_s10"] = smooth(df.pse_24_180_nopae, 10)
    df["ss8"] = ss8
    df["ss3"] = [SS3.get(s, "C") for s in ss8]
    df["af_frag"] = [frags[k]["frag"] for k in best]
    df["af_n_frags"] = len(frags)
    df["af_version"] = frags[0]["ver"]
    df["af_seq_date"] = frags[0]["seq_date"]
    return df, "ok"


def main():
    global PAE_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=14)
    ap.add_argument("--limit", type=int, default=0, help="stop after N matching members (testing)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--pae-dir", default=str(PAE_DIR))
    a = ap.parse_args()
    PAE_DIR = Path(a.pae_dir)  # read by forked workers
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    prot = canonical_swissprot()
    seqs = dict(zip(prot.acc, prot.seq))
    t0 = time.time()

    def members():
        n = 0
        with tarfile.open(TAR, "r|") as tf:
            for m in tf:
                mm = NAME_RE.search(m.name)
                if not m.isfile() or not mm or mm.group(1) not in seqs:
                    continue
                yield m.name, tf.extractfile(m).read()
                n += 1
                if a.limit and n >= a.limit:
                    return

    by_acc: dict[str, list] = {}
    n_members = 0
    with Pool(a.procs) as pool:
        for r in pool.imap_unordered(process_member, members(), chunksize=4):
            by_acc.setdefault(r["acc"], []).append(r)
            n_members += 1
            if n_members % 2000 == 0:
                print(f"{n_members} models, {time.time() - t0:.0f}s", flush=True)
    t_parse = time.time() - t0
    print(f"parsed {n_members} models for {len(by_acc)} proteins in {t_parse:.0f}s", flush=True)

    dfs, status = [], {}
    pae_status = {}
    for acc, seq in seqs.items():
        if acc not in by_acc:
            status[acc] = "not_in_afdb"
            continue
        pae_status[acc] = ",".join(sorted({f.get("pae_status", "na") for f in by_acc[acc]}))
        df, st = assemble(acc, seq, by_acc[acc])
        status[acc] = st
        if df is not None:
            dfs.append(df)
    res = pd.concat(dfs, ignore_index=True)

    # IDR / exposure flags. structuremap (with PAE): IDR = smoothed (±10) full-sphere 24 A pPSE <= 34.27;
    # high accessibility = pPSE(12 A, 70 deg) <= 5. AFDB rate-limits PAE downloads (HTTP 429), so PAE is only
    # available for a random calibration subset (cache/structure/pae_calibration_accs.txt). Main flags use the
    # PAE-free pPSE with thresholds calibrated on that subset to best reproduce the PAE-based flags.
    calib_accs = set((out / "pae_calibration_accs.txt").read_text().split()) if (out / "pae_calibration_accs.txt").exists() \
        else set((OUT / "pae_calibration_accs.txt").read_text().split())
    cal = res[res.acc.isin(calib_accs) & res.pse_24_180_pae.notna()].copy()
    cal["idr_sm"] = (cal.pse_24_180_pae_s10 <= IDR_SM_T).astype(int)
    cal["high_acc_5"] = (cal.pse_12_70_pae <= 5).astype(int)

    def best_t(x, y, grid):
        acc_ = np.array([((x <= t).astype(int) == y).mean() for t in grid])
        k = int(np.argmax(acc_))
        t = float(grid[k]); p = (x <= t).astype(int)
        po, pe = acc_[k], p.mean() * y.mean() + (1 - p.mean()) * (1 - y.mean())
        return t, float(acc_[k]), float((po - pe) / (1 - pe))

    t_idr, acc_idr, kap_idr = best_t(cal.pse_24_180_nopae_s10.values, cal.idr_sm.values, np.arange(20, 90, 0.25))
    t_hacc, acc_hacc, kap_hacc = best_t(cal.pse_12_70_nopae.values, cal.high_acc_5.values, np.arange(0, 20, 1.0))
    calib = dict(n_proteins=int(cal.acc.nunique()), n_residues=len(cal),
                 idr=dict(t_nopae=t_idr, accuracy=acc_idr, kappa=kap_idr, frac_idr_pae=float(cal.idr_sm.mean())),
                 high_acc=dict(t_nopae=t_hacc, accuracy=acc_hacc, kappa=kap_hacc, frac_pae=float(cal.high_acc_5.mean())),
                 spearman={c: float(cal[c + "_pae"].corr(cal[c + "_nopae"], method="spearman"))
                           for c in [f"pse_{r}_{t}" for r, t in PSE_GRID]})
    cal.to_parquet(out / "pae_calibration_residues.parquet", index=False)
    res.columns = [c.replace("_nopae", "") if "_nopae" in c else
                   (("p" + c.replace("_pae", "")) if c.startswith("pse_") and "_pae" in c else c) for c in res.columns]
    has = res.ppse_24_180.notna()
    res["idr_sm"] = np.where(has, (res.ppse_24_180_s10 <= IDR_SM_T).astype(float), np.nan)
    res["high_acc_5"] = np.where(has, (res.ppse_12_70 <= 5).astype(float), np.nan)
    res["idr_nopae"] = (res.pse_24_180_s10 <= t_idr).astype("int8")
    res["high_acc_nopae"] = (res.pse_12_70 <= t_hacc).astype("int8")
    res["idr"] = res.idr_sm.fillna(res.idr_nopae).astype("int8")
    res["high_acc"] = res.high_acc_5.fillna(res.high_acc_nopae).astype("int8")
    res["idr_source"] = np.where(has, "pae", "nopae_calibrated")
    res["disorder_plddt"] = (res.plddt_s10 < 50).astype("int8")
    for c in [c for c in res.columns if c.startswith("pse_") and not c.endswith("_s10")]:
        res[c] = res[c].astype("int16")
    for c in [c for c in res.columns if c.startswith("ppse_") and not c.endswith("_s10")]:
        res[c] = res[c].astype("Int16")
    res.to_parquet(out / "residue_structure.parquet", index=False)

    # S/T/Y table on the full canonical universe (NA rows for proteins without usable AFDB model)
    sty = [(acc, i + 1, ch) for acc, seq in seqs.items() for i, ch in enumerate(seq) if ch in "STY"]
    sty = pd.DataFrame(sty, columns=["acc", "pos", "res"])
    sty["site_id"] = sty.acc + "_" + sty.res + sty.pos.astype(str)
    sty = sty.merge(res.drop(columns=["res"]), on=["acc", "pos"], how="left")
    sty["af_status"] = sty.acc.map(status)
    sty = sty[["site_id"] + [c for c in sty.columns if c != "site_id"]]
    sty.to_parquet(out / "site_structure.parquet", index=False)

    st = pd.Series(status).value_counts().to_dict()
    meta = dict(tar=str(TAR.relative_to(ROOT)), n_models=n_members, n_proteins_canonical=len(seqs), status=st,
                pae_status=pd.Series(pae_status).value_counts().to_dict(), pae_calibration=calib, n_sty=len(sty), n_sty_ok=int((sty.af_status == "ok").sum()),
                seconds_parse=t_parse, seconds_total=time.time() - t0,
                mismatch_accs=sorted(k for k, v in status.items() if v == "seq_mismatch"),
                not_in_afdb_accs=sorted(k for k, v in status.items() if v == "not_in_afdb"),
                other_fail=sorted(k for k, v in status.items() if v not in ("ok", "seq_mismatch", "not_in_afdb")))
    (out / "site_structure.meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if not k.endswith("accs") and k != "other_fail"}, indent=1))


if __name__ == "__main__":
    main()
