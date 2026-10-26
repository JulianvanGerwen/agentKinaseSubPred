"""One-off migration (2026-10-26): store.py switched from whole-file index hashes to ID-column
hashes. Re-stamp descriptors whose recorded hash equals a known historical index FILE hash whose
ID order is identical to the current index (verified: kinases order/inclusion unchanged; only the
residue_class column changed). Any other hash is left untouched and will fail load_matrix."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from ksp import store
KNOWN = {  # axis -> historical whole-file sha256 values with identical ID order
    "kinase": {"c9e3e1731d2dc94b1eeea6b7cb15ace1ff11d779e92d343a11beb27026699aed"},
    "site": {"a280d3bb1a843c965b148b2067715567c9a515e9f3e18ab486bb97ed0d052975"},
    "protein": {"996d3bee3a3440ff0fefa99af67422a6e08c18007b8b2fa03311c304cba790e5"},
}
import hashlib
KNOWN["kinase"].add(hashlib.sha256((store.INDEX_DIR / "kinases.parquet").read_bytes()).hexdigest())
n = 0
for f in store.STORE.rglob("*.json"):
    d = json.loads(f.read_text())
    if "row_axis" not in d: continue
    ch = False
    for ax, hk in (("row_axis", "row_index_sha256"), ("col_axis", "col_index_sha256")):
        a = d.get(ax)
        if a in KNOWN and d.get(hk) in KNOWN[a]:
            d[hk] = store._index_hash(a); ch = True
    if ch: f.write_text(json.dumps(d, indent=1)); n += 1
print("restamped", n)
