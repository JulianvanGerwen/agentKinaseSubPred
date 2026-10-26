#!/bin/bash
# usage: dl.sh <source_dir> <url> [outname]  -- downloads to /workspace/external/<source_dir>/ (data only, never executed)
d=/workspace/external/$1; mkdir -p "$d"; u=$2; f=${3:-$(basename "${u%%\?*}")}
curl -L -sS --fail --retry 3 -m 3000 -A "Mozilla/5.0 (research download)" -o "$d/$f" "$u" && echo "OK $1/$f $(stat -c%s "$d/$f")" || { echo "FAIL $1 $u"; rm -f "$d/$f"; }
