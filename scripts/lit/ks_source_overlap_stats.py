#!/usr/bin/env python
"""Reproduce the descriptive numbers in docs/lit/02_ks_annotation_sources.md.

Usage (run with python -I):
  python -I scripts/lit/ks_source_overlap_stats.py --psp data/Kinase_Substrate_Dataset \
      [--omnipath enzsub.tsv] [--uniprot modres_by.tsv]

Inputs (none are hardcoded; fetch them yourself and record them in DATA.md if kept):
  --psp       PhosphoSitePlus Kinase_Substrate_Dataset (3 header lines, TSV, latin1)
  --omnipath  https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references,curation_effort&types=phosphorylation&organisms=9606
  --uniprot   https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606) AND (reviewed:true) AND (ft_mod_res:by)&fields=accession,gene_primary,ft_mod_res&format=tsv
Assumptions (stated, not silent): kinase identity across species is matched by upper-cased PSP GENE
symbol; UniProt 'by X' strings are matched to PSP GENE symbols verbatim (family names such as 'PKA'
are therefore NOT matched); isoform suffixes are stripped when comparing to OmniPath/UniProt, so a few
mismatches may be residue-numbering artefacts.
"""
import argparse, collections, re
import pandas as pd


def load_psp(path):
    d = pd.read_csv(path, sep='\t', skiprows=3, encoding='latin1', dtype=str).fillna('')
    d['K'] = d.GENE.str.upper()
    d['vivo'] = d.IN_VIVO_RXN.str.strip() == 'X'
    d['vitro'] = d.IN_VITRO_RXN.str.strip() == 'X'
    d['res'] = d.SUB_MOD_RSD.str[0]
    return d


def psp_stats(d):
    print('== PSP rows', len(d))
    print(pd.crosstab(d.vivo, d.vitro))
    h = d[(d.KIN_ORGANISM == 'human') & (d.SUB_ORGANISM == 'human')]
    print('human-human rows', len(h)); print(pd.crosstab(h.vivo, h.vitro))
    print('isoform substrate accessions (all/human-human)', d.SUB_ACC_ID.str.contains('-').sum(),
          h.SUB_ACC_ID.str.contains('-').sum())
    vc = h.groupby('K').size().sort_values(ascending=False)
    print('kinases', len(vc), 'top10 share %.3f top20 share %.3f' % (vc.head(10).sum() / vc.sum(), vc.head(20).sum() / vc.sum()))
    print('kinases with >=5 rows', (vc >= 5).sum(), '>=20 rows', (vc >= 20).sum())
    print('top substrates', h.groupby('SUB_ACC_ID').size().sort_values(ascending=False).head(5).to_dict())
    # cross-species pooling
    hmr = d[d.KIN_ORGANISM.isin(['human', 'mouse', 'rat']) & d.SUB_ORGANISM.isin(['human', 'mouse', 'rat'])]
    hsub = d[d.SUB_ORGANISM == 'human']
    hp = set(zip(hsub.K, hsub.SITE_GRP_ID)); hsg = set(hsub.SITE_GRP_ID)
    nonh = hmr[hmr.SUB_ORGANISM != 'human']
    key = pd.Series(list(zip(nonh.K, nonh.SITE_GRP_ID)), index=nonh.index)
    new = nonh[~key.isin(hp)].drop_duplicates(['K', 'SITE_GRP_ID'])
    print('mouse/rat-substrate rows', len(nonh), '-> new (kinase,site-group) pairs', len(new),
          'of which site group has no human row in this file', (~new.SITE_GRP_ID.isin(hsg)).sum())
    # biochemically implausible residue annotations vs evidence class
    f = h.groupby('K').apply(lambda x: pd.Series({'fy': (x.res == 'Y').mean(), 'n': len(x)}))
    tk = f[(f.fy >= 0.8) & (f.n >= 10)].index
    x = h[h.K.isin(tk) & h.res.isin(['S', 'T'])]
    print('S/T rows on Y-dominant kinases', len(x), 'in-vivo-only', (x.vivo & ~x.vitro).sum(),
          '| baseline in-vivo-only fraction %.3f' % (h.vivo & ~h.vitro).mean())
    return h


def omnipath_stats(path, h):
    o = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    pk = set(zip(h.KIN_ACC_ID.str.split('-').str[0], h.SUB_ACC_ID.str.split('-').str[0], h.SUB_MOD_RSD))
    o['key'] = list(zip(o.enzyme, o.substrate, o.residue_type + o.residue_offset))
    o['prim'] = o.sources.map(lambda s: set(x.split('_')[0] for x in s.split(';') if x))
    print('== OmniPath human phosphorylation rows', len(o), 'enzymes', o.enzyme.nunique())
    c = collections.Counter(x for s in o.prim for x in s); print('rows per primary source', c.most_common())
    u = collections.Counter(next(iter(s)) for s in o.prim if len(s) == 1); print('single-source rows', u.most_common())
    inpsp = o.key.isin(pk); print('rows matching provided PSP human-human', inpsp.sum())
    cur = {'SIGNOR', 'HPRD', 'phosphoELM', 'dbPTM', 'NCI-PID', 'Reactome'}
    tm = {'ProtMapper', 'Sparser', 'REACH', 'RLIMS-P', 'BEL-Large-Corpus'}
    n = o[~inpsp]
    cat = collections.Counter()
    for s in n.prim:
        if 'PhosphoSite' in s: cat['PSP-labelled but not in provided file'] += 1
        elif s & cur: cat['curated LT db'] += 1
        elif 'PhosphoNetworks' in s and not s - {'PhosphoNetworks', 'MIMP', 'KEA', 'Li2012'}: cat['PhosphoNetworks only'] += 1
        elif s <= tm: cat['text-mining only'] += 1
        elif 'KEA' in s: cat['KEA (mixed)'] += 1
        else: cat['other'] += 1
    print('non-PSP rows by category', cat.most_common())
    nn = n[n.prim.map(lambda s: 'PhosphoSite' not in s and bool(s & cur))]
    pkin = set(h.KIN_ACC_ID.str.split('-').str[0])
    print('curated-LT non-PSP rows', len(nn), 'enzymes', nn.enzyme.nunique(),
          'rows whose enzyme is a PSP human kinase', nn.enzyme.isin(pkin).sum())
    print('top non-kinase/complex enzymes', nn[~nn.enzyme.isin(pkin)].enzyme_genesymbol.value_counts().head(10).to_dict())


def uniprot_stats(path, h):
    psk = set(zip(h.K, h.SUB_ACC_ID.str.split('-').str[0], h.SUB_MOD_RSD.str[1:]))
    kg = set(h.K); rows = []; ev = collections.Counter()
    for i, l in enumerate(open(path)):
        if i == 0: continue
        acc, g, m = l.rstrip('\n').split('\t')
        for feat in m.split('MOD_RES ')[1:]:
            pos = feat.split(';')[0].strip()
            note = re.search(r'/note="([^"]+)"', feat)
            if not note or not note.group(1).startswith('Phospho'): continue
            mm = re.search(r'; by ([^;"]+)', note.group(1))
            if not mm: continue
            for x in set(re.findall(r'ECO:\d+', feat)): ev[x] += 1
            for k in re.split(r',| and ', mm.group(1)):
                k = k.strip().upper()
                if k: rows.append((k, acc, pos, 'ECO:0000269' in feat))
    d = pd.DataFrame(rows, columns=['k', 'acc', 'pos', 'exp'])
    print('== UniProt phospho "by" kinase-site rows', len(d), 'autocatalysis', (d.k == 'AUTOCATALYSIS').sum())
    print('evidence codes (per site)', ev.most_common(6))
    g = d[(d.k != 'AUTOCATALYSIS') & d.k.isin(kg)]
    inp = pd.Series([(a, b, c) in psk for a, b, c in zip(g.k, g.acc, g.pos)], index=g.index)
    print('gene-symbol-matched rows', len(g), 'in PSP', inp.sum(), 'not in PSP', (~inp).sum(),
          'not in PSP with ECO:0000269', (~inp & g.exp).sum())
    print('top not-in-PSP kinases', g[~inp].k.value_counts().head(8).to_dict())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--psp', required=True); ap.add_argument('--omnipath'); ap.add_argument('--uniprot')
    a = ap.parse_args()
    h = psp_stats(load_psp(a.psp))
    if a.omnipath: omnipath_stats(a.omnipath, h)
    if a.uniprot: uniprot_stats(a.uniprot, h)
