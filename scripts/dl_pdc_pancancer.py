#!/usr/bin/env python3
"""Download selected CPTAC pan-cancer (Li et al. 2023) open files from PDC via its public signed-URL endpoint. Data only."""
import json,sys,os,urllib.request,subprocess
out='/workspace/external/cptac_pancancer_pdc'; os.makedirs(out,exist_ok=True)
def post(query,variables):
    r=urllib.request.Request('https://pdc.cancer.gov/graphql',json.dumps({'query':query,'variables':variables}).encode(),{'Content-Type':'application/json'})
    return json.load(urllib.request.urlopen(r,timeout=120))
files=post('query{getAllUIPancancerFiles{file_id file_name description cohorts}}',{})['data']['getAllUIPancancerFiles']
json.dump(files,open(out+'/pancancer_file_manifest.json','w'),indent=1)
want=sys.argv[1:]
for f in files:
    if not any(f['file_name'].startswith(w) for w in want): continue
    u=json.load(urllib.request.urlopen('https://pdc.cancer.gov/pdcapi/file/signedURLFromUuid/'+f['file_id'],timeout=60))['data']
    dest=f"{out}/{f['file_name']}"
    r=subprocess.run(['curl','-sS','--fail','-L','--retry','3','-o',dest,u]); print(f['file_name'],r.returncode,os.path.getsize(dest) if os.path.exists(dest) else None,flush=True)
