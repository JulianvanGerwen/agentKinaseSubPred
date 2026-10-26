"""Paired cluster bootstrap (ksp.eval.paired_bootstrap) between two pooled R3 CV prediction sets (same universe).
  python -I scripts/models/r3_paired.py --a <5 run ids> --b <5 run ids> --name <out name> [--B 1000]
Output: results/models/r3/paired_<name>.tsv (A - B)."""
import argparse
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd  # noqa: E402
from ksp import eval as kev, manifest  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--a", nargs=5, required=True)
ap.add_argument("--b", nargs=5, required=True)
ap.add_argument("--name", required=True)
ap.add_argument("--B", type=int, default=1000)
a = ap.parse_args()
load = lambda rs: pd.concat([pd.read_parquet(manifest.CACHE / r / "predictions.parquet",
                                             columns=["kinase_idx", "site_idx", "score"]) for r in rs], ignore_index=True)
data = kev.prepare("cv", fold=None)
ra = kev.evaluate(load(a.a), data)
rb = kev.evaluate(load(a.b), data)
d = kev.paired_bootstrap(ra, rb, B=a.B)
d.to_csv(ROOT / f"results/models/r3/paired_{a.name}.tsv", sep="\t", index=False)
print(d.to_string())
