#!/bin/bash
Z=https://zenodo.org/api/records
get(){ d="${KSP_DATA_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/external/$1"; mkdir -p "$d"; for i in 1 2 3; do curl -L -sS --fail --retry 2 -m 3000 -A "Mozilla/5.0 (research download)" -o "$d/$3" "$2" && { if head -c 300 "$d/$3" | grep -q "403 Forbidden"; then rm -f "$d/$3"; sleep 120; else echo "OK $1/$3 $(stat -c%s "$d/$3")"; return; fi; } || sleep 60; done; echo "FAIL $1 $2"; }
get dtc_zenodo_18450948 $Z/18450948/files/DTC_DrugTargetInteractions.RDa/content DTC_DrugTargetInteractions.RDa
get curvecurator_zenodo_8399823 $Z/8399823/files/Kinobeads.zip/content Kinobeads.zip
get kiba_zenodo_5105698 $Z/5105698/files/kiba.zip/content kiba.zip
get hein2025_zenodo_14000085 $Z/14000085/files/Organelle_IP_analyses_and_figures-2024-10-28.zip/content Organelle_IP_analyses_and_figures-2024-10-28.zip
get kinoplex_zenodo_20633591 $Z/20633591/files/Complete_KinoPlex_STY_Structural_Atlas.feather/content Complete_KinoPlex_STY_Structural_Atlas.feather
