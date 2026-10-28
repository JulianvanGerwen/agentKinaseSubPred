"""Recompute the 9 primary matched comparisons (point estimates for all; full B=10000 paired bootstrap for the
pairs listed in argv[1]) and compare with superiority.tsv."""
import sys, json, time
sys.path.insert(0, '/workspace/src')
import numpy as np, pandas as pd
from pathlib import Path
from ksp import eval as E, eval_test as T

ROOT = Path('/workspace')
OUT = Path('<scratchpad>/recomp')
cfg = json.load(open(ROOT / 'FREEZE.json'))
by = {c['name']: c for c in cfg['configs']}
boot = set(sys.argv[1].split(',')) if len(sys.argv) > 1 else set()
sup = pd.read_csv(ROOT / 'results/test/721dc1977366/superiority.tsv', sep='\t')
t0 = time.time()
tl = pd.read_parquet(ROOT / 'sealed/test_labels.parquet')
data = E.prepare('test', labels_df=tl, with_canaries=False)
pf = T.load_predictions(by['final'], root=ROOT)
rows = []
for b, eps in [('B1_i', ['M1', 'M3']), ('B1_ii', ['M2']), ('B3', ['M1', 'M2', 'M3']), ('B4', ['M1', 'M2', 'M3'])]:
    pb = T.load_predictions(by[b], root=ROOT)
    pa, pbm, rep = T.matched_pair(pf, pb, len(data.sites))
    ra = E.evaluate(pa, data, missing='drop')
    rb = E.evaluate(pbm, data, missing='drop')
    d = E.paired_bootstrap(ra, rb, B=10000, seed=0).set_index('endpoint') if b in boot else None
    for ep in eps:
        s = sup[(sup.b == b) & (sup.endpoint == ep)].iloc[0]
        r = {'b': b, 'ep': ep, 'A': ra.endpoints[T._ekey(ep)], 'B': rb.endpoints[T._ekey(ep)],
             'diff': ra.endpoints[T._ekey(ep)] - rb.endpoints[T._ekey(ep)], 'batch_diff': s['diff'],
             'matched_kinases': rep['matched_kinases'], 'batch_mk': s.matched_kinases,
             'frac_A': rep['frac_A_kept'], 'batch_frac_A': s.frac_A_pairs_kept,
             'nK_def_A': ra.endpoints['M1_n_kinases_defined'], 'batch_nK': s.A_M1_n_kinases_defined}
        if d is not None:
            r.update({'lo95': d.loc[ep, 'lo95'], 'hi95': d.loc[ep, 'hi95'], 'p': d.loc[ep, 'p_boot'],
                      'batch_lo95': s.lo95, 'batch_hi95': s.hi95, 'batch_p': s.p_boot})
        rows.append(r)
        print(r, flush=True)
    del pb, pa, pbm, ra, rb
R = pd.DataFrame(rows)
R.to_csv(OUT / 'primary_recompute.tsv', sep='\t', index=False)
print('holm on batch p:', E.holm(sup.p_boot.values), 'floor', 2 / 10001, 'done', round(time.time() - t0), flush=True)
