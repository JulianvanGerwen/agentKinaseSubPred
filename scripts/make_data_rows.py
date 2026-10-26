#!/usr/bin/env python3
"""Append wave-1 download rows to DATA.md from logs/dl/sha256_wave1.txt (one row per directory-source, file hashes listed)."""
import os,collections
src={ # dir: (URL, version, license)
'string':('https://stringdb-downloads.org/download/ (protein.*.v12.0, 9606)','STRING v12.0','CC BY 4.0'),
'hpa':('https://www.proteinatlas.org/download/ (tsv/subcellular_location, rna_tissue_consensus, rna_celline; proteinatlas.tsv.zip)','HPA current release as of 2026-10-26 (version not stated in file)','CC BY-SA 4.0'),
'compartments':('https://download.jensenlab.org/human_compartment_{knowledge,experiments,integrated}_full.tsv','COMPARTMENTS current as of 2026-10-26','CC BY 4.0'),
'go':('https://current.geneontology.org/annotations/goa_human.gaf.gz ; https://purl.obolibrary.org/obo/go/go-basic.obo','GO current as of 2026-10-26','CC BY 4.0'),
'omnipath':('https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references&format=tsv&license=academic','OmniPath web service as of 2026-10-26','Mixed per-resource; academic licence set (includes non-commercial resources e.g. PSP)'),
'bioplex':('https://bioplex.hms.harvard.edu/data/','BioPlex 3.0 (293T 10K Dec 2019, HCT116 5.5K Dec 2019) + interactionList_v4a (combined)','CC BY 4.0 (check site)'),
'huri':('https://interactome-atlas.org/data/{HuRI,Lit-BM,HI-union}.tsv (ENSG pairs; www. host has a TLS cert mismatch, apex host used)','HuRI as of 2026-10-26','CC BY 4.0 (check site)'),
'opencell':('https://opencell.czbiohub.org/data/datasets/','OpenCell as of 2026-10-26','CC BY 4.0'),
'humancellmap':('https://humancellmap.org/resources/downloads/','Human Cell Map "latest" (2026-10-26)','CC BY 4.0'),
'biogrid':('https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/ (ORGANISM tab3 zip = per-species tab3 files incl. Homo_sapiens; PTMS ptm zip)','BioGRID 5.0.262','MIT'),
'intact':('https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip','IntAct release of 2026-01-14 (MITAB 2.7, all species, 11 GB unzipped)','CC BY 4.0'),
'iptmnet':('https://research.bioinformatics.udel.edu/iptmnet_data/files/current/','iPTMnet 6.2','CC BY-NC-SA 4.0'),
'kinhub':('http://www.kinhub.org/kinases.html (HTML table, to be parsed)','KinHub as of 2026-10-26','Not stated; academic use'),
'kinbase_manning2002':('http://kinase.com/human/kinome/tables/ (Manning et al. 2002 Science SOM Table S1/S2)','Manning 2002','Supplementary data of published paper; cite'),
'depmap':('https://api.figshare.com/v2/articles/31660582 (files 67214582, 67214579); depmap.org portal blocked by Cloudflare challenge','DepMap Public 26Q1 Chronos parameters: gene_effect.csv (= CRISPR gene effect), library_effect.csv','CC BY 4.0'),
'kinase_library':('https://pypi.org/project/kinase-library/1.8.0/ (sdist; extracted as data only into src/, never installed/executed)','kinase-library 1.8.0','CC BY-NC-SA 3.0'),
'cantley_johnson2023':('https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/ (MOESM1-6)','Johnson et al. 2023 Nature suppl. (MOESM3 T1 kinase info, MOESM4 T2 PSSMs, MOESM5 T3 site scores/percentiles [370MB], MOESM6 T4)','Publisher supplementary data; academic use, cite'),
'cantley_yaronbarir2024':('https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/ (MOESM1-9)','Yaron-Barir et al. 2024 Nature suppl. (MOESM4 T2 Tyr PSSMs, MOESM5 T3 annotation, MOESM7 SH2 matrices)','Publisher supplementary data; academic use, cite'),
'ochoa2020':('https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0344-3/MediaObjects/ (MOESM3-7; MOESM5 = functional_score, 116,258 sites; MOESM4 = annotated_phosphoproteome features)','Ochoa et al. 2020 Nat Biotechnol suppl.','Publisher supplementary data; cite'),
'cptac_pancancer_pdc':('https://pdc.cancer.gov/pdc/cptac-pancancer via /pdcapi/file/signedURLFromUuid/<file_id> (open files)','CPTAC pan-cancer harmonized v1 (Li et al. 2023 Cancer Cell)','Open access (NCI PDC / CPTAC data-use policy); cite'),
'cptac_pancan_v1.2':('https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/<COHORT>/ (LinkedOmics mirror)','BCM-harmonised freeze v1.2, 10 cohorts','Open (CPTAC data-use policy); cite'),
'benchmarkin':('https://zenodo.org/records/12566560','Zenodo 10.5281/zenodo.12566560 (2024-06-27; Muller-Dott et al.)','CC BY 4.0'),
}
sizes={}; rows=collections.defaultdict(list)
for l in open('/workspace/logs/dl/sha256_wave1.txt'):
    h,p=l.strip().split('  ',1); p=p[2:]; d=os.path.dirname(p).split('/')[0]
    rows[d].append((p,h,os.path.getsize('/workspace/external/'+p)))
out=['\n## Downloaded (external/) -- wave 1, access date 2026-10-26\nAll under /workspace/external/<dir>/. Per-file SHA-256 below each source; full list in logs/dl/sha256_wave1.txt (kinase_library extracted tree: logs/dl/sha256_kinase_library_extracted.txt).\n']
for d in sorted(rows):
    u,v,lic=src[d]; tot=sum(r[2] for r in rows[d])
    out.append(f'### {d}\n- URL: {u}\n- Version/release: {v}\n- Accessed: 2026-10-26; License: {lic}; total size: {tot/1e6:.1f} MB\n')
    out.append('| file | bytes | SHA-256 |\n|---|---|---|')
    for p,h,s in sorted(rows[d]): out.append(f'| {p.split("/",1)[1]} | {s} | {h} |')
    out.append('')
open('/workspace/DATA.md','a').write('\n'.join(out)+'\n')
