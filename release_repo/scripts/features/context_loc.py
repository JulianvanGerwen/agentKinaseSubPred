"""Build F4 localisation features (kinase x protein cosine similarities) + site topology vector.

Run: .venv/bin/python -I scripts/features/context_loc.py
Outputs (via ksp.store, under cache/features/context/loc/):
  loc_cos_<src>           K x P float16  cosine of coarse compartment vectors (src in hpa, opencell, hcm, hein, comp_exp, comp_know)
  loc_cos_combined        K x P float16  cosine of the mean unit vector over {hpa, opencell, hcm, hein}
  loc_cos_combined_know   K x P float16  ABLATION: same + COMPARTMENTS knowledge channel
  loc_n_sources           K x P int8     number of primary sources where both kinase and protein have a vector
  loc_hcm_nmf_cos         K x P float16  cosine of native 20-dim HCM NMF profiles
  loc_hein_profile_corr   K x P float16  Pearson correlation of Hein 2025 organelle-IP enrichment profiles
  loc_hpa_min_reliability K x P float16  min(HPA reliability of kinase, of protein) (1/0.8/0.5/0.2)
  vec_<src>_protein (P x 15), vec_<src>_kinase (K x 15) float32 compartment vectors (NaN rows = missing)
  site_topo_code (sites,) int8, site_lumenal (sites,) float16 (1 extracellular/lumenal, 0 cytoplasmic side, NaN unknown/TM/signal)
Also results/features/context/loc_coverage.json and loc_descriptor.json.
"""
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ksp.features import loc as L  # noqa: E402
from ksp.features.ppi_idmap import MapLog  # noqa: E402
from ksp.index import load_kinases, load_proteins, load_sites  # noqa: E402
from ksp.store import save_matrix  # noqa: E402

OUT = ROOT / "cache" / "features" / "context" / "loc"
RES = ROOT / "results" / "features" / "context"
PRIMARY = ["hpa", "opencell", "hcm", "hein"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kin = load_kinases().uniprot_acc.tolist()
    prot = load_proteins().acc.tolist()
    log = MapLog()
    hpa, hpa_rel = L.load_hpa(log)
    hcm, hcm_nat = L.load_hcm(log)
    hein, hein_nat = L.load_hein(log)
    vecs = {"hpa": hpa, "opencell": L.load_opencell(log), "hcm": hcm, "hein": hein,
            "comp_exp": L.load_compartments("experiments", log), "comp_know": L.load_compartments("knowledge", log)}
    log.dump(RES / "id_mapping_rates.json")
    desc = {"module": "F4 localisation", "label_derived": False, "vocab": L.VOCAB, "features": {}}

    def save(name, arr, text, missing, row="kinase", col="protein", ablation=False, sign=1):
        save_matrix(f"context/loc/{name}", arr, row, col, meta={"description": text, "missing": missing,
                                                              "label_derived": False, "ablation_only": ablation})
        desc["features"][name] = {"file": f"cache/features/context/loc/{name}.npy", "store_name": f"context/loc/{name}",
                                  "dtype": str(arr.dtype), "shape": list(arr.shape), "description": text,
                                  "missing": missing, "ablation_only": ablation, "sign": sign}

    cov = {}
    kset, pset = set(kin), set(prot)
    for s, V in vecs.items():
        abl = s == "comp_know"
        C = L.cosine_kxp(V, kin, prot).astype(np.float16)
        save(f"loc_cos_{s}", C, f"cosine similarity of {s} coarse compartment vectors (kinase vs protein)",
             "NaN = kinase or protein has no annotation in this source", ablation=abl)
        for ax, ids in (("kinase", kin), ("protein", prot)):
            M = V.reindex(ids).to_numpy(np.float32)
            save_matrix(f"context/loc/vec_{s}_{ax}", M, ax, "none", col_ids=L.VOCAB,
                        meta={"description": f"{s} compartment vector (VOCAB order); NaN row = missing"})
        cov[s] = {"n_proteins_any": int(len(V)), "kinases": round(len(kset & set(V.index)) / len(kin), 4),
                  "index_proteins": round(len(pset & set(V.index)) / len(prot), 4)}

    def combined(srcs):
        acc = sorted(set().union(*[set(vecs[s].index) for s in srcs]))
        S = np.zeros((len(acc), len(L.VOCAB)))
        pos = {a: i for i, a in enumerate(acc)}
        for s in srcs:
            U = L.unit_rows(vecs[s].to_numpy(float))
            idx = [pos[a] for a in vecs[s].index]
            S[idx] += np.nan_to_num(U)
        return pd.DataFrame(S, index=acc, columns=L.VOCAB)

    save("loc_cos_combined", L.cosine_kxp(combined(PRIMARY), kin, prot).astype(np.float16),
         "cosine of mean unit compartment vector over HPA, OpenCell, HCM, Hein", "NaN = no primary source for one side")
    save("loc_cos_combined_know", L.cosine_kxp(combined(PRIMARY + ["comp_know"]), kin, prot).astype(np.float16),
         "ABLATION: combined incl. COMPARTMENTS knowledge channel (GO/UniProt curation)", "NaN as above", ablation=True)
    n_src = np.zeros((len(kin), len(prot)), np.int8)
    for s in PRIMARY:
        kp = np.array([a in vecs[s].index for a in kin])
        pp = np.array([a in vecs[s].index for a in prot])
        n_src += np.outer(kp, pp).astype(np.int8)
    save("loc_n_sources", n_src, "number of primary localisation sources with both kinase and protein annotated (0-4)",
         "0 = none")
    save("loc_hcm_nmf_cos", L.cosine_kxp(hcm_nat, kin, prot).astype(np.float16),
         "cosine of native Human Cell Map NMF rank profiles (20-dim)", "NaN = not an HCM prey")
    save("loc_hein_profile_corr", L.corr_kxp(hein_nat, kin, prot).astype(np.float16),
         "Pearson correlation of z-scored Hein 2025 organelle-IP enrichment profiles", "NaN = not quantified")
    rk = hpa_rel.reindex(kin).to_numpy(float)
    rp = hpa_rel.reindex(prot).to_numpy(float)
    save("loc_hpa_min_reliability", np.minimum(rk[:, None], rp[None, :]).astype(np.float16),
         "min of HPA reliability (Enhanced 1, Supported .8, Approved .5, Uncertain .2) of kinase and protein",
         "NaN = not in HPA")
    cov["hcm_nmf"] = {"kinases": round(len(kset & set(hcm_nat.index)) / len(kin), 4),
                      "index_proteins": round(len(pset & set(hcm_nat.index)) / len(prot), 4)}
    cov["hein_profile"] = {"kinases": round(len(kset & set(hein_nat.index)) / len(kin), 4),
                           "index_proteins": round(len(pset & set(hein_nat.index)) / len(prot), 4)}

    # ---- site topology
    sites = load_sites()
    code = L.site_topology(sites[["acc", "pos"]])
    save_matrix("context/loc/site_topo_code", code, "site", "none", meta={"codes": L.TOPO_CODES})
    lum = np.full(len(code), np.nan, np.float16)
    lum[np.isin(code, [2, 5])] = 1
    lum[np.isin(code, [0, 1, 6, 7])] = 0
    save_matrix("context/loc/site_lumenal", lum, "site", "none",
                meta={"description": "1 = extracellular/lumenal (TOPO_DOM or after signal peptide of non-TM protein), "
                                     "0 = cytoplasmic/nuclear/mito/perox side, NaN = TM/signal/unknown"})
    desc["site_features"] = {"site_topo_code": {"store_name": "context/loc/site_topo_code", "codes": L.TOPO_CODES},
                             "site_lumenal": {"store_name": "context/loc/site_lumenal"}}
    vc = pd.Series(code).value_counts().sort_index()
    cov["site_topo_code_counts"] = {f"{int(k)}:{L.TOPO_CODES[int(k)]}": int(v) for k, v in vc.items()}
    cov["site_lumenal_frac1"] = float(np.nanmean(lum == 1))
    (RES / "loc_coverage.json").write_text(json.dumps(cov, indent=1))
    desc["coverage"] = cov
    (OUT / "descriptor.json").write_text(json.dumps(desc, indent=1))
    shutil.copy(OUT / "descriptor.json", RES / "loc_descriptor.json")
    print(json.dumps(cov, indent=1))


if __name__ == "__main__":
    main()
