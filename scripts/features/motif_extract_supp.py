"""Stream-extract Cantley supplementary site tables (Johnson 2023 Suppl. Table 3; Yaron-Barir 2024
Suppl. Table 3) from the strict-OOXML xlsx (openpyxl cannot read them) into parquet.
Keeps SITE_+/-7_AA, phos_res, Phosphosite, Uniprot cols and all '<kinase>_percentile' columns.
Usage: python -I scripts/features/motif_extract_supp.py <xlsx> <sheetN> <out.parquet> [max_rows]
Untrusted-data rule: only XML is parsed (xml.etree), nothing is executed."""
import re, sys, zipfile
import xml.etree.ElementTree as ET
import pandas as pd

f, sheet, out = sys.argv[1:4]
max_rows = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
z = zipfile.ZipFile(f)
ss = []
for ev, el in ET.iterparse(z.open("xl/sharedStrings.xml")):
    if el.tag.endswith("}si"):
        ss.append("".join(t.text or "" for t in el.iter() if t.tag.endswith("}t")))
        el.clear()
col_re = re.compile(r"^([A-Z]+)")
header, keep, rows = None, None, []
for ev, el in ET.iterparse(z.open(f"xl/worksheets/{sheet}.xml")):
    if not el.tag.endswith("}row"):
        continue
    vals = {}
    for c in el:
        v = [x for x in c if x.tag.endswith("}v")]
        if not v:
            continue
        v = v[0].text
        if c.get("t") == "s":
            v = ss[int(v)]
        vals[col_re.match(c.get("r")).group(1)] = v
    el.clear()
    if header is None:
        header = vals
        keep = {k: h for k, h in vals.items()
                if h.endswith("_percentile") or h in ("SITE_+/-7_AA", "phos_res", "Phosphosite", "Uniprot Primary Accession", "Database")}
        continue
    rows.append({h: vals.get(k) for k, h in keep.items()})
    if len(rows) >= max_rows:
        break
d = pd.DataFrame(rows)
for c in d.columns:
    if c.endswith("_percentile"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
d.to_parquet(out, index=False)
print(out, d.shape)
