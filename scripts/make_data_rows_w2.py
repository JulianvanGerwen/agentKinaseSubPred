#!/usr/bin/env python3
"""Append wave-2 download rows to DATA.md from logs/dl/sha256_wave2.txt (+ git clone commits)."""
import os,collections
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ksp.paths import ROOT  # noqa: E402  project data/working dir ($KSP_DATA_DIR, default: repo root)
Z='https://zenodo.org/api/records/'
src={ # key: (URL, version, license)
'hgnc':('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt','HGNC complete set, current as of 2026-10-26','CC0'),
'buljan2020':('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7427327/supplementaryFiles (Buljan 2020 Mol Cell, doi 10.1016/j.molcel.2020.07.001; mmc3=Table S2 final AP-MS interactome, mmc4=S3 co-purified pairs, mmc7=S6 predicted regulatory phosphosites, mmc2=S1 kinase classification)','Mol Cell 2020 suppl. (PMC7427327)','CC BY-NC-ND (article)'),
'sugiyama2019':('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6642169/supplementaryFiles (Sugiyama 2019 Sci Rep, doi 10.1038/s41598-019-46385-4; MOESM3 = Table S2 198,536 in-vitro kinase-substrate-site rows, MOESM6 = per-kinase counts)','Sci Rep 2019 suppl.','CC BY 4.0'),
'davis2011':('https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.1990/MediaObjects/ (MOESM2-6; Davis 2011 Nat Biotech; .xls OLE2 - xlrd not installed in env, not yet parsed)','Nat Biotech 2011 suppl.','Publisher supplementary; cite'),
'reinecke2024':('https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-023-01459-3/MediaObjects/ (MOESM1-8; Reinecke 2024 Nat Chem Biol "target landscape of 1,000 kinase inhibitors")','Nat Chem Biol 2024 suppl.','Publisher supplementary; cite'),
'kinoplex_natbiotech2026':('https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-026-03239-5/MediaObjects/ (MOESM1-8; "AI-enabled structural atlas decodes kinase specificity", Nat Biotech 2026)','Nat Biotech 2026 suppl.','Publisher supplementary; cite'),
'proteomehd_kustatscher2019':('https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-019-0461-4/MediaObjects/ (MOESM1-9; Kustatscher 2019 Nat Methods ProteomeHD; MOESM3/4 .txt may be the co-regulation tables - verify)','Nat Methods 2019 suppl.','Publisher supplementary; cite'),
'kincore':('https://dunbrack.fccc.edu/kincore/static/downloads/{alignment-files/Human-PK-alignment.{aln,fasta},alignment-files/Human-PK.fasta,fasta_with_labels/PK_labels_PDB.fasta}','KinCore web download as of 2026-10-26 (497 human PK domain seqs; header ">GROUP_GENE/start-end ENTRY GENE ACC")','Not stated (academic); cite Modi & Dunbrack'),
'phosformer/hf_phosformer':('https://huggingface.co/waylandy/phosformer (resolve/main)','HF last modified 2022-11-04; RoBERTa-style Phosformer weights','CC BY-NC-ND 4.0 (LICENSE file in repo) -- NON-COMMERCIAL, NO DERIVATIVES'),
'phosformer/hf_dataset_curated':('https://huggingface.co/datasets/waylandy/phosformer_curated (resolve/main)','HF last modified 2023-01-25; corpus/, curated/, splits/ (train.csv 715 MB)','Not stated on card; repo code is CC BY-NC-ND 4.0 -- treat as non-commercial'),
'phosformer/hf_space_phosformerST':('https://huggingface.co/spaces/gravelcompbio/Phosformer-ST (resolve/main) -- Phosformer-ST weights (ESM2 t30 150M multitask)','HF last modified 2025-04-16','CC BY-NC-ND 4.0 (LICENSE file)'),
'kolossus':('https://pypi.org/project/kolossus/1.0.6/ sdist (kolossus-1.0.6.tar.gz, not installed); weights https://huggingface.co/aparekh2/kolossus (large epoch_59 = 100.7 MB, small epoch_55 = 15.5 MB; README: large uses ESM-2 3B?; small ESM2-650M)','kolossus 1.0.6; HF sha 64088d2bf798be505dc7e03a22a2a09fa684c104 (2025-07-01)','MIT (HF tag); PyPI license field empty'),
'selphi2/zenodo_13821321':(Z+'13821321/files/<name>/content (SELPHI2.0 Zenodo 10.5281/zenodo.13821321; small files only; NOT downloaded: prediction_matrix.csv 8.3 GB, feature_table.tsv 10.7 GB, k_p_*_pssm_scores.tsv 3.9-4.1 GB)','Zenodo 13821321 (2024-09-20)','CC BY 4.0'),
'orphan_top5_zenodo_16881580':(Z+'16881580 "Top5 Kinase Predictions for Orphan Human Phosphosites"','Zenodo 2025-08-14','CC BY 4.0'),
'ikipdb_zenodo_21469258':(Z+'21469258 (phosphoproteome atlas of human cell lines / kinase activity; iKiP-DB tables; NOT downloaded: SpectralLibraries.zip 6.5 GB, Phosphosites_and_ProteinGroups_Tables.zip 525 MB)','Zenodo 2026-07-21','CC BY 4.0'),
'ikipdb_zenodo_20037245':(Z+'20037245 (same atlas; only kinase scores file; NOT downloaded: Phos table 321 MB, Prot table 59 MB)','Zenodo 2026-05-05','CC BY 4.0'),
'kinasecanvas_zenodo_22941318':(Z+'22941318 "KinaseCanvas: integrated resource for curated and novel kinase-substrate interactions" (KinaseCanvas.zip, 1.47 GB; contains literature-derived + predicted KS: do NOT use as independent test data)','Zenodo 2026-09-24','MIT'),
'curvecurator_zenodo_8399823':(Z+'8399823 Kinobeads.zip (CurveCurator example outputs: Klaeger-style kinobeads dose-response, all_binder_data_lfq_annotated.txt, curves.txt); NOT downloaded: CTRP.zip, Simulations.zip','Zenodo 2023-08-04','CC BY 4.0'),
'kiba_zenodo_5105698':(Z+'5105698 kiba.zip (KiBA.xlsx: comparison of Davis/Metz/Klaeger-type kinase inhibitor bioactivity sets; Szwajda/Tang-lab)','Zenodo 2021-07-15','CC BY 4.0'),
'hein2025_zenodo_14000085':(Z+'14000085 (analysis repo of Hein et al. Cell 2025: cluster_annotation_Dec6.csv, enrichment tables NOC_prop, volcano tables; the final per-protein localisation table is in the Cell supplement = 403)','Zenodo 2024-10-28','CC BY 4.0'),
'kinoplex_zenodo_20633591':(Z+'20633591 Complete_KinoPlex_STY_Structural_Atlas.feather (all canonical STY residues x structural features); NOT downloaded: .csv duplicate 3.1 GB, All_ST_bootstrapped_quantiles.feather 3.5 GB','Zenodo 2026-06-10','CC BY 4.0'),
'zenodo_meta':('Zenodo search API JSON (q=recid:ID) for the records above, incl. decryptM 2.0 17533475 (restricted, no files)','2026-10-26','n/a'),
}
rows=collections.defaultdict(list)
for l in open(ROOT / 'logs/dl/sha256_wave2.txt'):
    h,p=l.strip().split('  ',1); p=p[2:]; parts=p.split('/')
    k='/'.join(parts[:2]) if '/'.join(parts[:2]) in src else parts[0]
    rows[k].append((p,h,os.path.getsize(ROOT / 'external' / p),len(k)+1))
out=['\n## Downloaded (external/) -- wave 2, access date 2026-10-26\nAll under external/<dir>/. Full hash list: logs/dl/sha256_wave2.txt. Git clones (data only, not executed): phosformer/github = https://github.com/esbgkannan/phosformer (waylandy/phosformer) commit 6e7c019ba5a4c4f6cf7e44c5359c6115afc14c16, CC BY-NC-ND 4.0; kolossus/github = https://github.com/rohitsinghlab/KolossuS commit e8f3013a9736be07a1ff17a79ca7c9529daec63b; selphi2/selphi_2.0 = https://gitlab.ebi.ac.uk/petsalakilab/selphi_2.0 commit f9822f90f174e09e392d789698d8c5ec33baf1d2, Apache-2.0 (commits: logs/dl/wave2_git_commits.txt). Scripts: scripts/dl_w2.sh, scripts/dl_w2b.sh, scripts/make_data_rows_w2.py.\n']
for d in sorted(rows):
    u,v,lic=src[d]; tot=sum(r[2] for r in rows[d])
    out.append(f'### {d}\n- URL: {u}\n- Version/release: {v}\n- Accessed: 2026-10-26; License: {lic}; total size: {tot/1e6:.1f} MB\n')
    out.append('| file | bytes | SHA-256 |\n|---|---|---|')
    for p,h,s,n in sorted(rows[d]): out.append(f'| {p[len(d)+1:] if p.startswith(d+"/") else p} | {s} | {h} |')
    out.append('')
open(ROOT / 'DATA.md','a').write('\n'.join(out)+'\n')
