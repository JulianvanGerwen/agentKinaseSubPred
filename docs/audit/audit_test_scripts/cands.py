import sys; sys.path.insert(0,'/workspace/src')
import numpy as np, pandas as pd
from ksp import eval as E, splits
tl=pd.read_parquet('/workspace/sealed/test_labels.parquet')
ks=splits.load_eval_kinase_sets(); Kev=set(ks.kinase_acc[ks.in_K_eval]); Kcv=set(ks.kinase_acc[ks.in_K_cv])
def stats(data, kset, tag):
    pm=data.pos_meta[~data.pos_meta.kinase_acc.str.startswith('CANARY_')]
    K=data.kinases; real=~K.uniprot_acc.str.startswith('CANARY_').values
    kidx={a:i for i,a in enumerate(K.uniprot_acc)}
    cls=K.residue_class.values; res=data.sites.res.values
    nST=int(((cls=='ST')|(cls=='STY'))[real].sum()); nY=int(((cls=='Y')|(cls=='STY'))[real].sum())
    out=[]
    for k,s in zip(pm.kinase_acc.values, pm.site_idx.values):
        if k not in kset: continue
        r=res[s]; base = nY if r=='Y' else nST
        others=(data.site_kinases.get(s,set())|data.par.get(k,set()))-{k}
        okc_o=[o for o in others if o in kidx and real[kidx[o]] and E._res_ok(np.array([r]),cls[kidx[o]])[0]]
        n_true=len([o for o in data.site_kinases.get(s,set()) if o in kidx and real[kidx[o]]])
        out.append((k, r=='Y', base-1-len(okc_o), len(okc_o), n_true))
    d=pd.DataFrame(out,columns=['k','Y','n_comp','n_filtered','n_lab_site'])
    pk=d.groupby('k').mean()
    print(tag, 'pairs',len(d),'kin',d.k.nunique(),'| pair-mean comps %.1f filtered %.2f labelled kin/site %.2f | kinase-macro comps %.1f filtered %.2f fracY %.3f' % (
        d.n_comp.mean(), d.n_filtered.mean(), d.n_lab_site.mean(), pk.n_comp.mean(), pk.n_filtered.mean(), pk.Y.mean()), flush=True)
dt=E.prepare('test',labels_df=tl,with_canaries=False)
stats(dt,Kev,'TEST K_eval')
dc=E.prepare('cv',with_canaries=False)
stats(dc,Kcv,'CV K_cv')
stats(dc,Kev&Kcv,'CV K_eval∩K_cv')
