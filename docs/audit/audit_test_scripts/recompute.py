"""Audit recompute of TEST batch 721dc1977366 from pinned predictions + sealed labels (READ-ONLY).
Writes only into the scratchpad."""
import sys, json, time, gc
sys.path.insert(0, '/workspace/src')
import numpy as np, pandas as pd
from pathlib import Path
from ksp import eval as E, eval_test as T

ROOT = Path('/workspace')
OUT = Path('/tmp/claude-588982/-workspace/f2c6b342-58ef-4e11-b665-a8c647ff9a45/scratchpad/recomp')
OUT.mkdir(exist_ok=True)
BATCH = ROOT / 'results/test/721dc1977366'
cfg = json.load(open(ROOT / 'FREEZE.json'))
by = {c['name']: c for c in cfg['configs']}
t0 = time.time()
tl = pd.read_parquet(ROOT / 'sealed/test_labels.parquet')
data = E.prepare('test', labels_df=tl, with_canaries=False)
print('prepared', round(time.time() - t0), data.info, len(data.k_primary), len(data.k_any), flush=True)

KEYS = ['M1_auroc_hard', 'M2_hit5', 'M3_top1_recall', 'M4_auroc_hard_Kany']
rep = {}
names = sys.argv[1].split(',') if len(sys.argv) > 1 else ['final', 'B1_i', 'B1_ii', 'B3', 'B4', 'R2', 'R3', 'B0_prior']
for n in names:
    p = T.load_predictions(by[n], root=ROOT)
    res = E.evaluate(p, data, missing='worst')
    b = json.load(open(BATCH / f'{n}.json'))['endpoints']
    rep[n] = {k: (res.endpoints[k], b[k], abs(res.endpoints[k] - b[k])) for k in KEYS}
    print(n, {k: (round(v[0], 6), round(v[1], 6), f'{v[2]:.1e}') for k, v in rep[n].items()}, flush=True)
    res.per_kinase.to_csv(OUT / f'{n}_per_kinase.tsv', sep='\t', index=False)
    res.site_level.to_csv(OUT / f'{n}_site_level.tsv', sep='\t', index=False)
    if n in ('final', 'B1_i', 'B1_ii'):
        st = T.strata_results(res, data, n)
        st.to_csv(OUT / f'{n}_strata.tsv', sep='\t', index=False)
        bs = pd.read_csv(BATCH / f'{n}_strata.tsv', sep='\t')
        d = np.nanmax(np.abs(st[['M1', 'M2', 'M3', 'M4']].values.astype(float) - bs[['M1', 'M2', 'M3', 'M4']].values.astype(float)))
        print(n, 'strata max abs diff vs batch', d, flush=True)
        et = T.ext_tier_results(res, data, p, tl, n)
        et.to_csv(OUT / f'{n}_ext_tiers.tsv', sep='\t', index=False)
    if n == 'final':
        truths = T._site_truths(data)
        rel = pd.concat([T.reliability(p, data, t, sub, n) for sub, t in truths.items()], ignore_index=True)
        rel.to_csv(OUT / 'final_reliability.tsv', sep='\t', index=False)
        br = pd.read_csv(BATCH / 'final_reliability.tsv', sep='\t')
        print('reliability max abs diff', np.nanmax(np.abs(rel[['mean_p', 'accuracy', 'ECE']].values - br[['mean_p', 'accuracy', 'ECE']].values)), flush=True)
        # top-kinase per labelled site, with residue, for the audit
        lab = truths['labelled_TEST_sites']
        x = p[p.site_idx.isin(list(lab))]
        top = x.loc[x.groupby('site_idx').assignment_prob.idxmax()].copy()
        top['res'] = data.sites.res.values[top.site_idx.values]
        kacc = data.kinases.uniprot_acc.values
        top['correct'] = [kacc[k] in lab[s] for k, s in zip(top.kinase_idx.values, top.site_idx.values)]
        top.to_csv(OUT / 'final_top_assign_labelled.tsv', sep='\t', index=False)
        # label-free: distribution of max assignment prob over all TEST U_ms sites by residue
        q = p[data.universe[p.site_idx.values]]
        mx = q.groupby('site_idx').assignment_prob.max()
        r = data.sites.res.values[mx.index.values]
        print('max p_top by residue (all TEST U_ms sites):', {rr: float(mx.values[r == rr].max()) for rr in ('S', 'T', 'Y')},
              'n>=0.5', int((mx >= 0.5).sum()), 'n>=0.8', int((mx >= 0.8).sum()), flush=True)
    del p, res
    gc.collect()
json.dump({n: {k: list(v) for k, v in d.items()} for n, d in rep.items()}, open(OUT / f'endpoints_{"_".join(names)}.json', 'w'), indent=1)
print('done', round(time.time() - t0), flush=True)
