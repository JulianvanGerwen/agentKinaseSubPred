import json, pandas as pd, numpy as np, sys
sys.path.insert(0,'/workspace/src')
from ksp import splits
SP='/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad'
B='/workspace/results/test/721dc1977366'
ks=splits.load_eval_kinase_sets()
print(ks.columns.tolist()); 
Kev=set(ks.kinase_acc[ks.in_K_eval]); Kcv=set(ks.kinase_acc[ks.in_K_cv])
print('K_eval',len(Kev),'K_cv',len(Kcv),'both',len(Kev&Kcv),'K_eval not in K_cv',len(Kev-Kcv))
runs={'final':'r2_fusion_r3__20261026-191901__s0','R2':'r2_final_ens5__20261026-191213__s-1',
'B1_i':'simple_b1_published_reeval__20261026-155352__s0','B1_ii':'simple_b1_universe_reeval__20261026-155703__s0','B0_prior':'simple_b0_prior_reeval__20261026-155334__s0'}
cv={}
for n,r in runs.items():
    cv[n]=pd.DataFrame(json.load(open(f'results/runs/{r}/metrics_a24.json'))['per_kinase']).set_index('kinase_acc')
cv['R3']=pd.read_csv('results/runs/r3_final_distill_cv_pooled__20261026-173414__s0/per_kinase.tsv',sep='\t').set_index('kinase_acc')
bins=['0','1-4','5-19','20-99','>=100']
tc=pd.Series(ks.n_dev_L2.values,index=ks.kinase_acc)
binof=pd.cut(tc,[-1,0,4,19,99,1e12],labels=bins).astype(str)
print('K_cv bins',binof[list(Kcv)].value_counts().reindex(bins).to_dict())
print('K_eval bins',binof[list(Kev)].value_counts().reindex(bins).to_dict())
print('K_eval&K_cv bins',binof[list(Kev&Kcv)].value_counts().reindex(bins).to_dict())
rows=[]
for n in ['final','R2','R3','B1_ii','B1_i','B0_prior']:
    c=cv[n]; t=pd.read_csv(f'{B}/{n}_per_kinase.tsv',sep='\t').set_index('kinase_acc')
    kc=sorted(Kcv); kb=sorted(Kev&Kcv)
    r={'config':n,
       'CV_M2_Kcv':c.loc[kc,'hit5'].mean(),
       'CV_M2_on_Keval∩Kcv':c.loc[kb,'hit5'].mean(),
       'TEST_M2_Keval':t.loc[sorted(Kev),'hit5'].mean(),
       'TEST_M2_on_Keval∩Kcv':t.loc[kb,'hit5'].mean(),
       'CV_M1_Kcv':c.loc[kc,'auroc_hard'].mean(),'CV_M1_on_Keval∩Kcv':c.loc[kb,'auroc_hard'].mean(),
       'TEST_M1_Keval':t.loc[sorted(Kev),'auroc_hard'].mean(),'TEST_M1_on_Keval∩Kcv':t.loc[kb,'auroc_hard'].mean(),
       'CV_M3_Kcv':c.loc[kc,'top1_recall'].mean(),'CV_M3_on_Keval∩Kcv':c.loc[kb,'top1_recall'].mean(),
       'TEST_M3_Keval':t.loc[sorted(Kev),'top1_recall'].mean(),'TEST_M3_on_Keval∩Kcv':t.loc[kb,'top1_recall'].mean()}
    # CV per-bin M2 reweighted to K_eval composition (bins of K_eval kinases that are in K_cv)
    cb=c.loc[kc,'hit5'].groupby(binof[kc].values).mean()
    w=binof[sorted(Kev)].value_counts(); w=w[w.index.isin(cb.index)]
    r['CV_M2_reweighted_to_Keval_bins']=float((cb[w.index]*w).sum()/w.sum())
    for bb in ['5-19','20-99','>=100']:
        r[f'CV_M2_bin{bb}']=cb.get(bb,np.nan)
        tb=[k for k in sorted(Kev) if binof[k]==bb]; r[f'TEST_M2_bin{bb}']=t.loc[tb,'hit5'].mean()
    # paired per-kinase diff on common kinases
    d=(t.loc[kb,'hit5']-c.loc[kb,'hit5']).dropna()
    r['paired_TEST-CV_M2_mean']=d.mean(); r['paired_se']=d.std()/np.sqrt(len(d)); r['n_paired']=len(d)
    d1=(t.loc[kb,'auroc_hard']-c.loc[kb,'auroc_hard']).dropna()
    r['paired_TEST-CV_M1_mean']=d1.mean(); r['paired_M1_se']=d1.std()/np.sqrt(len(d1))
    rows.append(r)
R=pd.DataFrame(rows).set_index('config').T
pd.set_option('display.width',200)
print(R.round(4).to_string())
R.to_csv(f'{SP}/m2_cv_vs_test.tsv',sep='\t')
