#!/bin/bash
# usage: dl.sh <source_dir> <url> [outname]  -- downloads to $KSP_DATA_DIR/external/<source_dir>/ (default: repo root) (data only, never executed)
d="${KSP_DATA_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/external/$1"; mkdir -p "$d"; u=$2; f=${3:-$(basename "${u%%\?*}")}
curl -L -sS --fail --retry 3 -m 3000 -A "Mozilla/5.0 (research download)" -o "$d/$f" "$u" && echo "OK $1/$f $(stat -c%s "$d/$f")" || { echo "FAIL $1 $u"; rm -f "$d/$f"; }
