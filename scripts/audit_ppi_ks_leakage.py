#!/usr/bin/env python
"""Audit how much PSP kinase-substrate (KS) knowledge is already encoded in PPI resources.

For human-kinase / human-substrate protein pairs in PhosphoSitePlus (PSP), compute the fraction of
pairs that have evidence in each PPI source / STRING channel, versus two backgrounds:
  (a) RAND : kinase x random STRING-mapped protein (same number per kinase, excluding its PSP substrates)
  (b) SHUF : kinase x PSP substrate of *other* kinases (controls for substrate study-bias / degree)
A large KS/SHUF ratio means the source carries KS-specific information (signal or leakage);
sources curated from literature (BioGRID 'Biochemical Activity', STRING database/textmining) are
expected to be leaky. Unbiased screens (BioPlex, HuRI, OpenCell) should show smaller ratios.

Inputs are read-only from data/ and external/ (paths as CLI args with defaults). Run with python -I.
Usage: python -I scripts/audit_ppi_ks_leakage.py [--seed 0] [--out results/leakage/ppi_ks_leakage.tsv]
Assumptions (stated): protein-level pairs (sites collapsed); isoform suffixes stripped; autophosphorylation
pairs (kinase == substrate) excluded; STRING mapping via UniProt_AC aliases (many-to-many allowed:
a pair has evidence if any mapped ENSP pair does).
"""
import argparse, collections, csv, gzip, io, os, random, sys, zipfile

p = argparse.ArgumentParser()
p.add_argument("--psp", default="/workspace/data/Kinase_Substrate_Dataset")
p.add_argument("--ext", default="/workspace/external")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--out", default="/workspace/results/leakage/ppi_ks_leakage.tsv")
a = p.parse_args()
rng = random.Random(a.seed)
E = a.ext
strip = lambda acc: acc.strip().split("-")[0]

# ---------------- PSP human KS protein pairs ----------------
ks = set(); kin_gene = {}
with open(a.psp, encoding="latin1") as f:
    for _ in range(3):
        f.readline()
    for r in csv.DictReader(f, delimiter="\t"):
        if r["KIN_ORGANISM"] != "human" or r["SUB_ORGANISM"] != "human":
            continue
        k, s = strip(r["KIN_ACC_ID"]), strip(r["SUB_ACC_ID"])
        if k and s and k != s:
            ks.add((k, s)); kin_gene[k] = r["GENE"]
kinases = sorted({k for k, _ in ks}); subs_all = sorted({s for _, s in ks})
subs_of = collections.defaultdict(set)
for k, s in ks:
    subs_of[k].add(s)
print(f"PSP human KS protein pairs: {len(ks)}; kinases {len(kinases)}; substrate proteins {len(subs_all)}", file=sys.stderr)

# ---------------- STRING ----------------
u2e = collections.defaultdict(set); e2g = {}
with gzip.open(f"{E}/string/9606.protein.aliases.v12.0.txt.gz", "rt") as f:
    f.readline()
    for line in f:
        e, al, src = line.rstrip("\n").split("\t")
        if src == "UniProt_AC":
            u2e[al].add(e)
        elif src == "Ensembl_gene":
            e2g[e] = al
proteome = sorted(u for u in u2e)  # STRING-mapped UniProt accessions


def sample_pairs():
    rand, shuf = set(), set()
    for k in kinases:
        n = len(subs_of[k])
        pool_r = [u for u in rng.sample(proteome, min(len(proteome), 5 * n + 50)) if u not in subs_of[k] and u != k][:n]
        rand.update((k, u) for u in pool_r)
        others = [s for s in rng.sample(subs_all, min(len(subs_all), 5 * n + 50)) if s not in subs_of[k] and s != k][:n]
        shuf.update((k, s) for s in others)
    return rand, shuf

rand, shuf = sample_pairs()
sets = {"KS": ks, "RAND": rand, "SHUF": shuf}
need_e = set()
for S in sets.values():
    for k, s in S:
        need_e |= u2e.get(k, set()) | u2e.get(s, set())


def load_string(fn, cols):
    d = {}
    with gzip.open(fn, "rt") as f:
        hdr = f.readline().split()
        idx = [hdr.index(c) for c in cols]
        for line in f:
            x = line.split()
            if x[0] in need_e and x[1] in need_e:
                d[(x[0], x[1])] = [int(x[i]) for i in idx]
    return d

full_cols = ["experimental", "database", "textmining", "coexpression", "combined_score"]
phys_cols = ["experimental", "database", "textmining", "combined_score"]
sfull = load_string(f"{E}/string/9606.protein.links.detailed.v12.0.txt.gz", full_cols)
sphys = load_string(f"{E}/string/9606.protein.physical.links.detailed.v12.0.txt.gz", phys_cols)


def string_best(d, k, s, j):
    best = 0
    for ek in u2e.get(k, ()):
        for es in u2e.get(s, ()):
            v = d.get((ek, es))
            if v:
                best = max(best, v[j])
    return best

# ---------------- BioGRID (human) ----------------
bg = collections.defaultdict(set)  # unordered pair -> set of tags
z = zipfile.ZipFile(f"{E}/biogrid/BIOGRID-ORGANISM-LATEST.tab3.zip")
name = [i.filename for i in z.infolist() if "Homo_sapiens" in i.filename][0]
with z.open(name) as fb:
    rd = csv.reader(io.TextIOWrapper(fb, encoding="utf-8", errors="replace"), delimiter="\t")
    h = next(rd); H = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r[H["Experimental System Type"]] != "physical":
            continue
        A = [x for x in r[H["SWISS-PROT Accessions Interactor A"]].split("|") if x and x != "-"]
        B = [x for x in r[H["SWISS-PROT Accessions Interactor B"]].split("|") if x and x != "-"]
        system = r[H["Experimental System"]]; thr = r[H["Throughput"]]
        tags = {"any_physical"}
        if system == "Biochemical Activity":
            tags.add("biochem_activity")
        else:
            tags.add("physical_excl_biochem")
            if "High Throughput" in thr and "Low Throughput" not in thr:
                tags.add("HTP_only_excl_biochem")
        for x in A:
            for y in B:
                bg[frozenset((x, y))] |= tags

# ---------------- Systematic screens ----------------
def add(dct, x, y, tag):
    dct[frozenset((x, y))].add(tag)

screens = collections.defaultdict(set)
for fn, tag in [("BioPlex_293T_Network_10K_Dec_2019.tsv", "BioPlex3_293T"), ("BioPlex_HCT116_Network_5.5K_Dec_2019.tsv", "BioPlex3_HCT116")]:
    with open(f"{E}/bioplex/{fn}") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            add(screens, strip(r["UniprotA"]), strip(r["UniprotB"]), tag)
g2u = collections.defaultdict(set)
for u, es in u2e.items():
    for e in es:
        if e in e2g:
            g2u[e2g[e]].add(u)
with open(f"{E}/huri/HuRI.tsv") as f:
    for line in f:
        x = line.split()
        if len(x) >= 2:
            for ua in g2u.get(x[0].split(".")[0], ()):
                for ub in g2u.get(x[1].split(".")[0], ()):
                    add(screens, ua, ub, "HuRI")
with open(f"{E}/opencell/opencell-protein-interactions.csv") as f:
    for r in csv.DictReader(f):
        for ua in g2u.get(r["target_ensg_id"], ()):
            for ub in r["interactor_uniprot_ids"].split(";"):
                add(screens, ua, strip(ub), "OpenCell")

# ---------------- Tally ----------------
features = ([("STRING_full>0:" + c, lambda k, s, j=j: string_best(sfull, k, s, j) > 0) for j, c in enumerate(full_cols[:-1])]
            + [("STRING_full_combined>=400", lambda k, s: string_best(sfull, k, s, 4) >= 400)]
            + [("STRING_physical>0:" + c, lambda k, s, j=j: string_best(sphys, k, s, j) > 0) for j, c in enumerate(phys_cols[:-1])]
            + [("STRING_physical_combined>=400", lambda k, s: string_best(sphys, k, s, 3) >= 400)]
            + [("BioGRID:" + t, lambda k, s, t=t: t in bg.get(frozenset((k, s)), ())) for t in ["any_physical", "biochem_activity", "physical_excl_biochem", "HTP_only_excl_biochem"]]
            + [("Screen:" + t, lambda k, s, t=t: t in screens.get(frozenset((k, s)), ())) for t in ["BioPlex3_293T", "BioPlex3_HCT116", "HuRI", "OpenCell"]])

# restrict all sets to pairs where both proteins map to STRING (common denominator)
for nm in sets:
    sets[nm] = {(k, s) for k, s in sets[nm] if k in u2e and s in u2e}
os.makedirs(os.path.dirname(a.out), exist_ok=True)
with open(a.out, "w") as fo:
    fo.write(f"# seed={a.seed}; n_KS={len(sets['KS'])} n_RAND={len(sets['RAND'])} n_SHUF={len(sets['SHUF'])}\n")
    fo.write("feature\tfrac_KS\tfrac_RAND\tfrac_SHUF\tratio_KS_SHUF\n")
    for nm, fn in features:
        fr = {S: sum(fn(k, s) for k, s in P) / max(1, len(P)) for S, P in sets.items()}
        ratio = fr["KS"] / fr["SHUF"] if fr["SHUF"] > 0 else float("inf")
        line = f"{nm}\t{fr['KS']:.4f}\t{fr['RAND']:.4f}\t{fr['SHUF']:.4f}\t{ratio:.1f}"
        fo.write(line + "\n"); print(line)
    # Conditional check: does STRING's experiments channel import BioGRID 'Biochemical Activity' (in vitro kinase assay) records?
    fo.write("# conditional (KS pairs only): frac with STRING_physical experimental>0, by BioGRID evidence class\n")
    groups = collections.defaultdict(list)
    for k, s in sets["KS"]:
        t = bg.get(frozenset((k, s)), set())
        g = ("biogrid_biochem_only" if t == {"any_physical", "biochem_activity"} else
             "biogrid_other_physical" if "physical_excl_biochem" in t else "no_biogrid_record")
        groups[g].append(string_best(sphys, k, s, 0) > 0)
    for g, v in sorted(groups.items()):
        line = f"cond:{g}\tn={len(v)}\tfrac_STRINGphys_exp>0={sum(v)/max(1,len(v)):.4f}"
        fo.write("# " + line + "\n"); print(line)
print("wrote", a.out, file=sys.stderr)
