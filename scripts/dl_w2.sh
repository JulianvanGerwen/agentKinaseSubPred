#!/bin/bash
# wave 2 downloads: zenodo (via api files/content links) + misc. data only.
Z=https://zenodo.org/api/records
get(){ d=external/$1; mkdir -p "$d"; curl -L -sS --fail --retry 3 -m 3000 -A "Mozilla/5.0 (research download)" -o "$d/$3" "$2" && echo "OK $1/$3 $(stat -c%s "$d/$3")" || { echo "FAIL $1 $2"; rm -f "$d/$3"; }; }
for f in prediction_matrix_overlap.csv prediction_matrix_high_conf_curated.csv kinase-kinase-rels-kinome_kin_20240910.tsv kinase-kinase-rels-kinome_sub_20240910.tsv kinase-citations-2_kin_20240910.tsv positive-set-20240910.tsv kinase-citations-2_sub_20240910.tsv; do get selphi2/zenodo_13821321 $Z/13821321/files/$f/content $f; done
get orphan_top5_zenodo_16881580 "$Z/16881580/files/top5_kinase_predictions_for_orphan_human_phosphosites(50threshold).csv/content" top5_orphan.csv
for f in Table6_ActivationLoopSites.txt Table7_iKiPDB_Modified_Combined_PSP.gmt Table5_Cscore.xlsx Site_Stoichiometry.zip; do get ikipdb_zenodo_21469258 $Z/21469258/files/$f/content $f; done
get ikipdb_zenodo_20037245 $Z/20037245/files/Kinase_Scores_Long_Weighted_OriginalDataset_Gene.csv/content Kinase_Scores_Long_Weighted_OriginalDataset_Gene.csv
get kinasecanvas_zenodo_22941318 $Z/22941318/files/KinaseCanvas.zip/content KinaseCanvas.zip
# --- part 2 (run separately as scripts/dl_w2b.sh)
