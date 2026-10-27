import sys, pandas as pd
d=sys.argv[1]
cur=pd.read_csv(f"{d}/phosphosites_11mer_kinase_specific.tsv",sep="\t")
cur.columns=[c.replace('|','_') for c in cur.columns]
te=pd.read_csv(f"{d}/test_pos_neg.csv"); va=pd.read_csv(f"{d}/val_pos_neg.csv")
print("curated rows",len(cur),"kinases",cur.kinase_uniprot.nunique())
print("test rows",len(te),"pos",te.label.sum(),"kinases",te.uniprot.nunique())
heldout=pd.concat([te[te.label==1],va[va.label==1]])
ho_pairs=set(zip(heldout.uniprot,heldout.substrate))
cur['pair']=list(zip(cur.kinase_uniprot,cur.substrate_sequence))
train=cur[~cur.pair.isin(ho_pairs)]
print("approx train positive rows",len(train))
tp=te[te.label==1]
trpep=set(train.substrate_sequence)
print("test positive peptides seen as a train POSITIVE (any kinase): %.1f%%"%(100*tp.substrate.isin(trpep).mean()))
# same kinase family, different organism / uniprot
fam=dict(zip(cur.kinase_uniprot,cur.kinase_family)); org=dict(zip(cur.kinase_uniprot,cur.kinase_organism))
trfam=set(zip(train.kinase_family,train.substrate_sequence))
tp=tp.assign(fam=tp.uniprot.map(fam))
print("test positives whose (kinase family, peptide) is a train positive: %.1f%%"%(100*pd.Series([ (f,s) in trfam for f,s in zip(tp.fam,tp.substrate)]).mean()))
# identical kinase domain seq different uniprot? ortholog with same peptide
trks=set(zip(train.kinase_sequence,train.substrate_sequence))
print("test positives with identical (kinase domain seq, peptide) in train: %.1f%%"%(100*pd.Series([(k,s) in trks for k,s in zip(tp.kinase,tp.substrate)]).mean()))
print("test kinase organisms:", tp.uniprot.map(org).value_counts().head(8).to_dict())
tn=te[te.label==0]
print("test hard-negative peptides seen as train positive: %.1f%%"%(100*tn.substrate.isin(trpep).mean()))
