"""Write canonical index files cache/index/{kinases,sites,proteins}.parquet (+ checksums)."""
import hashlib, json, sys
from pathlib import Path
import pandas as pd
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
from ksp.index import INDEX_DIR

INDEX_DIR.mkdir(parents=True, exist_ok=True)
k = pd.read_csv(ROOT / "results/labels/kinase_list.tsv", sep="\t")
k = k[k["include"]].sort_values("uniprot_acc").reset_index(drop=True)
k = k[["uniprot_acc", "gene", "residue_class", "manning_group", "manning_family",
       "manning_subfamily", "is_pseudokinase", "kinase_domain_start", "kinase_domain_end"]]
assert k.uniprot_acc.is_unique
u = pd.read_parquet(ROOT / "cache/sites/universe.parquet")
u = u[u["in_U_obs"]].sort_values("site_id").reset_index(drop=True)
ms = (u.psp_ms_lit.fillna(0) + u.psp_ms_cst.fillna(0) > 0) | u.src_compendium | u.src_cptac | u.src_ochoa
s = u[["site_id", "acc", "res", "pos", "gene", "window15", "window31"]].copy()
s["in_U_ms"] = ms.values
assert s.site_id.is_unique
p = pd.DataFrame({"acc": sorted(s.acc.unique())})
out = {}
for name, df in [("kinases", k), ("sites", s), ("proteins", p)]:
    f = INDEX_DIR / f"{name}.parquet"
    df.to_parquet(f, index=False)
    out[name] = {"n": len(df), "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}
(ROOT / "results/index_manifest.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
