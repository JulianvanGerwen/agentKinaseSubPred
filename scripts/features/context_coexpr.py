"""Build F5 co-expression / co-essentiality kinase x protein matrices (low priority).

Run: OMP_NUM_THREADS=10 .venv/bin/python -I scripts/features/context_coexpr.py
Outputs via ksp.store: context/coexpr/{coexpr_hpa_tissue_spearman, coexpr_depmap_pearson, coexpr_depmap_pearson_pc5}
"""
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from ksp.features import coexpr as C  # noqa: E402
from ksp.features.loc import corr_kxp  # noqa: E402
from ksp.features.ppi_idmap import MapLog  # noqa: E402
from ksp.index import load_kinases, load_proteins  # noqa: E402
from ksp.store import save_matrix  # noqa: E402

OUT = ROOT / "cache" / "features" / "context" / "coexpr"
RES = ROOT / "results" / "features" / "context"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kin = load_kinases().uniprot_acc.tolist()
    prot = load_proteins().acc.tolist()
    log = MapLog()
    desc = {"module": "F5 co-expression / co-dependency", "label_derived": False, "features": {}, "coverage": {}}

    def save(name, M, text, P):
        save_matrix(f"context/coexpr/{name}", M, "kinase", "protein",
                    meta={"description": text, "missing": "NaN = gene not in source / zero variance", "label_derived": False})
        desc["features"][name] = {"file": f"cache/features/context/coexpr/{name}.npy", "store_name": f"context/coexpr/{name}",
                                  "dtype": str(M.dtype), "shape": list(M.shape), "description": text, "sign": 1}
        desc["coverage"][name] = {"kinases": round(float(np.mean([a in P.index for a in kin])), 4),
                                  "index_proteins": round(float(np.mean([a in P.index for a in prot])), 4)}

    T = C.hpa_tissue_profiles(log)
    T = T[T.std(1) > 0]
    save("coexpr_hpa_tissue_spearman", corr_kxp(T, kin, prot).astype(np.float16),
         f"Spearman correlation of HPA RNA tissue consensus profiles ({T.shape[1]} tissues)", T)
    D = C.depmap_profiles(log)
    save("coexpr_depmap_pearson", corr_kxp(D, kin, prot).astype(np.float16),
         f"Pearson correlation of DepMap 26Q1 Chronos gene effects ({D.shape[1]} cell lines)", D)
    D5 = C.remove_top_pcs(D, 5)
    save("coexpr_depmap_pearson_pc5", corr_kxp(D5, kin, prot).astype(np.float16),
         "as coexpr_depmap_pearson after removing the top 5 principal components", D5)
    log.dump(RES / "id_mapping_rates.json")
    (OUT / "descriptor.json").write_text(json.dumps(desc, indent=1))
    shutil.copy(OUT / "descriptor.json", RES / "coexpr_descriptor.json")
    print(json.dumps(desc["coverage"], indent=1))


if __name__ == "__main__":
    main()
