"""(1) Strata of final / B1_i / B1_ii restricted to kinases B1 scores (matched-kinase view of the batch strata; M1/M3/M4
are exact matched values because per-kinase metrics only depend on that kinase's pairs; M2 keeps the standalone
competitor set). (2) A17 reliability of the batch's temporal / external subsets with the site's FULL truth set."""
import sys, json, dataclasses
sys.path.insert(0, '/workspace/src')
import numpy as np, pandas as pd
from pathlib import Path
from ksp import eval as E, eval_test as T

ROOT = Path('/workspace')
OUT = Path('<scratchpad>/recomp')
cfg = json.load(open(ROOT / 'FREEZE.json'))
by = {c['name']: c for c in cfg['configs']}
tl = pd.read_parquet(ROOT / 'sealed/test_labels.parquet')
data = E.prepare('test', labels_df=tl, with_canaries=False)
b1 = pd.read_csv(ROOT / 'results/test/721dc1977366/B1_i_per_kinase.tsv', sep='\t')
unsc = set(b1.kinase_acc[b1.n_missing > 0])
dm = dataclasses.replace(data, k_primary=[k for k in data.k_primary if k not in unsc],
                         k_any=[k for k in data.k_any if k not in unsc])
print('matched K_eval', len(dm.k_primary), 'K_any', len(dm.k_any), flush=True)
out = []
for n in ['final', 'B1_i', 'B1_ii']:
    p = T.load_predictions(by[n], root=ROOT)
    res = E.evaluate(p, data, missing='worst')
    st = T.strata_results(res, dm, n)
    st = st[st.stratum != 'external_curated']          # ext rows come from res.pair_strata (unmatched); skip
    st.loc[len(st)] = {'config': n, 'level': 'all', 'stratum': 'all', 'value': 'all',
                       'n_kinases_K_eval': len(dm.k_primary), 'n_kinases_K_any': len(dm.k_any),
                       'M1': res.per_kinase.set_index('kinase_acc').loc[dm.k_primary, 'auroc_hard'].mean(),
                       'M3': res.per_kinase.set_index('kinase_acc').loc[dm.k_primary, 'top1_recall'].mean(),
                       'M4': res.per_kinase.set_index('kinase_acc').loc[dm.k_any, 'auroc_hard'].mean()}
    out.append(st)
    if n == 'final':
        truths = T._site_truths(data)
        full = truths['labelled_TEST_sites']
        rows = []
        for sub in ['new_since_2021_sites', 'external_curated_sites']:
            t = truths[sub]
            ft = {s: (full.get(s, set()) | t[s]) for s in t}
            for tag, tt in [('batch_truth', t), ('full_site_truth', ft)]:
                r = T.reliability(p, data, tt, sub, n)
                rows.append(r.assign(truth=tag))
        R = pd.concat(rows, ignore_index=True)
        R.to_csv(OUT / 'final_reliability_fulltruth.tsv', sep='\t', index=False)
        print(R.drop_duplicates(['subset', 'truth'])[['subset', 'truth', 'n_total', 'accuracy_overall', 'ECE']], flush=True)
        for (sub, tg), g in R.groupby(['subset', 'truth']):
            hi = g[g.bin.isin(['[0.5,0.6)', '[0.6,0.7)', '[0.7,0.8)'])]
            n_ = hi.n_sites.sum()
            print(sub, tg, '>=0.5: n', n_, 'claimed', round((hi.n_sites * hi.mean_p).sum() / n_, 3),
                  'observed', round((hi.n_sites * hi.accuracy).sum() / n_, 3), flush=True)
    del p, res
S = pd.concat(out, ignore_index=True)
S.to_csv(OUT / 'matched_strata.tsv', sep='\t', index=False)
S['key'] = S.stratum + '=' + S.value.astype(str)
pd.set_option('display.width', 250)
for m in ['M1', 'M3', 'M4', 'M2']:
    print(m)
    print(S.pivot_table(index='key', columns='config', values=m, aggfunc='first').round(3).to_string())
print(S[S.config == 'final'][['key', 'n_kinases_K_eval', 'n_kinases_K_any', 'n_pairs']].to_string())
