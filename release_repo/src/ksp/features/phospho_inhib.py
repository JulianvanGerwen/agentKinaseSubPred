"""F6b inhibitor response (not label-derived): condition x kinase inhibition matrix + per-site target deconvolution.

Inhibition matrix A (conditions x kinases, canonical kinase order):
  * Hijazi 2020 (PMID 31959955): the provided binary target matrix data/kinhib_inhib_comb__1.csv (A = 1 for listed
    targets), kinase names resolved with ksp.kinases resolver (unknowns counted and reported).
  * every other compendium condition whose Condition/Control text names a kinase drug: fractional occupancy
    occ = dose / (dose + Kd) per kinase, with Kd from (priority) Klaeger 2017 Kinobeads Kd_app (data/Klaeger Table S3;
    'n.i.' = measured, not bound -> occ 0, no fallback), Reinecke 2024 Kinobeads (tool compounds), CurveCurator
    re-fit of the Klaeger curves (pEC50 of 'down' curves; only for drugs absent from Table S3), Davis 2011 KINOMEscan
    Kd (non-mutant constructs; blank = tested, Kd > 10 uM -> occ 0). Combinations: 1 - prod(1 - occ_i).
    A_row = occ(Condition drugs) - occ(Control drugs) (so 'DMSO vs drug' contrasts get a negative row and
    'drug+stimulus vs stimulus' rows keep only the drug).
  * curated class tokens without a named drug (e.g. 'MEKi', 'p38i', 'ATRi') map to designated targets with an assumed
    occupancy (results/features/phospho/drug_aliases.tsv, column occ; documented assumption).
  Doses parsed from the text (pM/nM/uM/mM, decryptM 'Drug_10nM', '1to2' combination ratios); unknown dose ->
  1 uM (flag). Time restriction: conditions with known treatment time > max_minutes (default 120) are dropped; unknown
  time is kept and flagged.

Per-site regression (GPU, batched over sites): y_s (robust-SD-scaled log2FC over the F6b conditions observed for site
s) ~ A_s beta_s, ridge lambda, no intercept. beta < 0 = the site decreases when k is inhibited. z = beta / se with
se from the ridge posterior covariance sigma^2 (G + lambda I)^-1, sigma^2 moderated towards the median (d0 = 4).
NaN where sum_c |A_ck| over observed conditions < min_exposure. Second feature: simple contrast = mean y over
observed conditions with A_ck >= 0.5 minus mean over observed conditions with |A_ck| < 0.1 (>= 2 / >= 3 conditions).
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .phospho import ROOT, RESULTS, CACHE, load_blocks, load_conditions

DATA = ROOT / "data"
EXT = ROOT / "external"
ALIASES = RESULTS / "drug_aliases.tsv"
HIJAZI_PMID = 31959955


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).replace("Ã¶", "o").replace("ö", "o").lower())


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


# ----------------------------------------------------------------------------------------------- affinity sources
def _xlsx_rows(path, sheet):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return list(wb[sheet].iter_rows(values_only=True))


def load_klaeger():
    """Long table drug, gene, kd_nM (inf = 'n.i.' measured not inhibited), + synonyms list [(drug, syn)]."""
    rows = _xlsx_rows(DATA / "Klaeger" / "Table_S3.xlsx", "Kinobeads Drugmatrix detailed")
    hdr = rows[0]
    drugs = hdr[3:]
    recs = []
    for r in rows[1:]:
        gene = r[0]
        for d, v in zip(drugs, r[3:]):
            if v is None or v == "" or v == "n.d.":
                continue
            if v == "n.i.":
                recs.append((d, gene, math.inf))
                continue
            try:
                recs.append((d, gene, float(v)))
            except ValueError:
                pass
    # cell-line-specific re-profiling (EGFR/ERBB2 drugs in BT-474, MET drugs with Caki-1 lysate): numeric Kd only;
    # the per-(drug, gene) minimum over sheets is used (these lysates express the designated targets better)
    for sheet in ("in BT-474 Drugmatrix by Drug", "with Caki-1 Drugmatrix by Drug"):
        rr = _xlsx_rows(DATA / "Klaeger" / "Table_S3.xlsx", sheet)
        for r in rr[1:]:
            for d, v in zip(rr[0][3:], r[3:]):
                try:
                    recs.append((d, r[0], float(v)))
                except (TypeError, ValueError):
                    pass
    s1 = _xlsx_rows(DATA / "Klaeger" / "aan4368_Table_S1.xlsx", "Inhibitor annotation")
    syn = []
    for r in s1[1:]:
        if r[0] and r[3]:
            syn += [(r[0], x.strip()) for x in re.split(r"[,;]", str(r[3])) if x.strip()]
    return pd.DataFrame(recs, columns=["drug", "gene", "kd_nM"]), syn


def load_reinecke():
    rows = _xlsx_rows(EXT / "reinecke2024" / "41589_2023_1459_MOESM3_ESM.xlsx", "Kinobeads Drugmatrix - all")
    drugs = rows[0][1:]
    recs = [(d, r[0], float(v)) for r in rows[1:] for d, v in zip(drugs, r[1:]) if isinstance(v, (int, float))]
    syn = []
    for sheet, a, b in (("Tool compounds", 0, 1), ("Clinical KIs", 0, 1)):
        for r in _xlsx_rows(EXT / "reinecke2024" / "41589_2023_1459_MOESM2_ESM.xlsx", sheet)[1:]:
            if r[a] and r[b]:
                syn += [(r[a], x.strip()) for x in re.split(r"[,;]", str(r[b])) if x.strip()]
    return pd.DataFrame(recs, columns=["drug", "gene", "kd_nM"]), syn


def load_curvecurator():
    import io
    import zipfile
    z = zipfile.ZipFile(EXT / "curvecurator_zenodo_8399823" / "Kinobeads.zip")
    d = pd.read_csv(io.BytesIO(z.read("Kinobeads/curves.txt")), sep="\t", low_memory=False,
                    usecols=["Drug", "Gene names", "pEC50", "Curve Regulation"])
    d = d[d.Drug.str.contains("_") == False]  # noqa: E712  (drop _inBT474/_withCAKI variants)
    d = d[d["Curve Regulation"] == "down"]
    d = d.assign(gene=d["Gene names"].astype(str).str.split(";").str[0], kd_nM=10 ** (9 - d.pEC50))
    return d.rename(columns={"Drug": "drug"})[["drug", "gene", "kd_nM"]], []


def load_davis():
    import xlrd
    wb = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM5_ESM.xls")
    sh = wb.sheet_by_index(0)
    hdr = sh.row_values(0)
    m3 = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM3_ESM.xls").sheet_by_index(0)
    mutant = {m3.row_values(i)[2]: m3.row_values(i)[3] == "YES" for i in range(1, m3.nrows)}
    recs = []
    for i in range(1, sh.nrows):
        r = sh.row_values(i)
        gene, construct = r[1], r[2]
        if mutant.get(construct, False) or "(" in str(construct) and "phosphoryl" in str(construct):
            continue
        for d, v in zip(hdr[3:], r[3:]):
            if v == "" or v is None:
                recs.append((d, gene, math.inf))
            else:
                try:
                    recs.append((d, gene, float(v)))
                except ValueError:
                    pass
    syn = []
    m4 = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM4_ESM.xls").sheet_by_index(0)
    for i in range(1, m4.nrows):
        r = m4.row_values(i)
        if r[0] and r[1]:
            syn += [(r[0], x.strip()) for x in re.split(r"[,;/]", str(r[1])) if x.strip()]
    d = pd.DataFrame(recs, columns=["drug", "gene", "kd_nM"])
    d = d.groupby(["drug", "gene"], as_index=False).kd_nM.min()  # several constructs (e.g. JAK JH1/JH2) -> best
    return d, syn


def load_kinhib_source():
    """Hijazi binary kinhib matrix as the lowest-priority affinity source for the same drugs used in OTHER studies:
    listed targets get Kd = 0 (occupancy 1 at any dose; ASSUMPTION), unlisted kinases are unknown (no row)."""
    kh = pd.read_csv(DATA / "kinhib_inhib_comb__1.csv", index_col=0)
    recs = [(d, g, 0.0) for g, row in kh.iterrows() for d, v in row.items() if v == 1]
    return pd.DataFrame(recs, columns=["drug", "gene", "kd_nM"]), []


def read_aliases(path: Path = ALIASES) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)


def build_affinity(resolver, verbose=True) -> tuple[pd.DataFrame, UnionFind, pd.DataFrame]:
    """Return (aff: drug_id, kinase_acc, kd_nM, source), union-find over normalised drug names, report."""
    srcs = [("klaeger2017", load_klaeger), ("reinecke2024", load_reinecke), ("curvecurator", load_curvecurator),
            ("davis2011", load_davis), ("hijazi_kinhib", load_kinhib_source)]
    uf = UnionFind()
    tabs, rep = [], []
    for name, fn in srcs:
        d, syn = fn()
        for a, b in syn:
            if norm(a) and norm(b) and len(norm(b)) >= 4:
                uf.union(norm(a), norm(b))
        for x in d.drug.unique():
            uf.find(norm(x))
        d = d.assign(source=name)
        tabs.append(d)
    al = read_aliases()
    for r in al.itertuples():
        if r.maps_to:
            uf.union(norm(r.alias), norm(r.maps_to))
    acc_of, unknown = {}, set()
    allg = pd.concat(tabs).gene.dropna().unique()
    for g in allg:
        gg = str(g).split(";")[0].strip()
        sym, how = resolver.resolve(gg)
        if sym is None:
            unknown.add(gg)
            continue
        acc_of[g] = resolver.kl.set_index("gene").uniprot_acc.get(sym)
    out = []
    covered = {}  # drug_id -> best source rank
    for rank, (d, (name, _)) in enumerate(zip(tabs, srcs)):
        d = d.assign(drug_id=d.drug.map(lambda x: uf.find(norm(x))), kinase_acc=d.gene.map(acc_of))
        n_genes_unres = int(d.kinase_acc.isna().sum())
        d = d[d.kinase_acc.notna()]
        drugs_here = set(d.drug_id)
        # a drug is taken from the highest-priority source that profiled it; lower sources only fill kinases that
        # the higher source did not measure at all (Klaeger n.d. / not on beads)
        d = d.groupby(["drug_id", "kinase_acc"], as_index=False).agg(kd_nM=("kd_nM", "min"), drug=("drug", "first"))
        d["source"] = name
        out.append(d)
        rep.append(dict(source=name, n_drugs=len(drugs_here), n_pairs=len(d), n_finite=int(np.isfinite(d.kd_nM).sum()),
                        n_rows_gene_not_in_kinase_list=n_genes_unres))
        for x in drugs_here:
            covered.setdefault(x, rank)
    aff = pd.concat(out, ignore_index=True)
    aff["rank"] = aff.source.map({n: i for i, (n, _) in enumerate(srcs)})
    aff = aff.sort_values(["drug_id", "kinase_acc", "rank"]).drop_duplicates(["drug_id", "kinase_acc"])
    # CurveCurator is a re-fit of Klaeger: only use it for drugs absent from Klaeger
    kl_drugs = set(aff.loc[aff.source == "klaeger2017", "drug_id"])
    aff = aff[~((aff.source == "curvecurator") & aff.drug_id.isin(kl_drugs))]
    rep = pd.DataFrame(rep)
    rep["n_unresolved_gene_names_total"] = len(unknown)
    if verbose:
        print(rep.to_string())
    return aff.reset_index(drop=True), uf, rep


# ----------------------------------------------------------------------------------------------- condition parsing
_DOSE = re.compile(r"(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*(pM|nM|uM|µM|μM|mM)(?![a-zA-Z])", re.I)
_UNIT_NM = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "μm": 1e3, "mm": 1e6}


def parse_doses(text: str) -> list[float]:
    return [float(v) * _UNIT_NM[u.lower()] for v, u in _DOSE.findall(str(text))]


def find_drugs(text: str, lexicon: set, uf: UnionFind, aliases: pd.DataFrame) -> list[tuple[str, str]]:
    """[(drug_id or 'class:<alias>', matched token)] for drug names in a condition string."""
    t = str(text)
    toks = [x for x in re.split(r"[\s_+,;()/:]+", t) if x]
    al = {norm(a): r for a, r in zip(aliases.alias, aliases.itertuples())}
    hits = []
    seen = set()
    cand = toks + ["".join(p) for p in zip(toks, toks[1:])]
    for tok in cand:
        k = norm(tok)
        if len(k) < 3 or k in seen:
            continue
        if k in al and al[k].targets:
            hits.append((f"class:{k}", tok))
            seen.add(k)
        elif k in lexicon:
            hits.append((uf.find(k), tok))
            seen.add(k)
    # de-duplicate the same drug id
    out, ids = [], set()
    for d, tok in hits:
        if d not in ids:
            out.append((d, tok))
            ids.add(d)
    return out


def drug_doses(text: str, drugs: list, pmid) -> tuple[list[float], str]:
    """Dose (nM) per drug found in text. decryptM 'A+B_1to2_10nM' -> A 10, B 20 nM. 32034161 '_0-1' / '_1' -> 0.1 /
    1 uM (study-specific notation, ASSUMPTION). Unknown -> 1000 nM."""
    t = str(text)
    if not drugs:
        return [], ""
    m = re.search(r"_(\d+)to(\d+)_([\d.e+]+)\s*(nM|uM|µM)", t)
    if m and len(drugs) == 2:
        base = float(m.group(3)) * _UNIT_NM[m.group(4).lower()]
        return [base * float(m.group(1)), base * float(m.group(2))], "ratio"
    ds = parse_doses(t)
    if len(ds) == 1:
        return [ds[0]] * len(drugs), "parsed"
    if len(ds) == len(drugs):
        return ds, "parsed"
    if str(pmid) == "32034161":
        m = re.search(r"_(0-\d+|\d+)(?:uM)?$", t)
        if m:
            v = float(m.group(1).replace("-", "."))
            return [v * 1000.0] * len(drugs), "parsed_32034161_uM"
    return [1000.0] * len(drugs), "default_1uM"


def occupancy(drug_id: str, dose_nM: float, aff_by_drug: dict, aliases: pd.DataFrame, kpos: dict,
              resolver, n_k: int) -> np.ndarray:
    v = np.zeros(n_k, np.float32)
    if drug_id.startswith("class:"):
        a = aliases[aliases.alias.map(norm) == drug_id[6:]].iloc[0]
        occ = dose_nM / (dose_nM + float(a.kd_nM)) if a.kd_nM else float(a.occ or 0.9)
        for g in a.targets.split(";"):
            sym, _ = resolver.resolve(g.strip())
            if sym is None:
                raise KeyError(f"alias target {g} not resolvable")
            acc = resolver.kl.set_index("gene").uniprot_acc[sym]
            if acc in kpos:
                v[kpos[acc]] = occ
        return v
    d = aff_by_drug.get(drug_id)
    if d is None:
        return v
    kd = d.kd_nM.to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        o = np.where(np.isfinite(kd) & (dose_nM > 0), dose_nM / (dose_nM + kd), 0.0)
    for acc, oo in zip(d.kinase_acc, o):
        if acc in kpos:
            v[kpos[acc]] = oo
    return v


def build_inhibition_matrix(resolver, kin: pd.DataFrame, max_minutes: float = 120.0, min_occ: float = 0.05,
                            exclude_pmids=(), exclude_drugs=(), verbose=True):
    """Returns (A conditions x kinases float32, cond table with mapping provenance, report dict)."""
    cond = load_conditions().copy()
    kpos = {a: i for i, a in enumerate(kin.uniprot_acc)}
    n_k = len(kin)
    aff, uf, aff_rep = build_affinity(resolver, verbose=verbose)
    aliases = read_aliases()
    aff_by_drug = {d: g for d, g in aff.groupby("drug_id")}
    lexicon = set(uf.p.keys())
    # Hijazi matrix
    kh = pd.read_csv(DATA / "kinhib_inhib_comb__1.csv", index_col=0)
    hij_unres = []
    H = np.zeros((kh.shape[1], n_k), np.float32)
    for name, row in kh.iterrows():
        sym, how = resolver.resolve(str(name))
        acc = resolver.kl.set_index("gene").uniprot_acc.get(sym) if sym else None
        if acc is None or acc not in kpos:
            hij_unres.append((name, how, int(row.sum())))
            continue
        H[:, kpos[acc]] = np.maximum(H[:, kpos[acc]], row.to_numpy(np.float32))
    hcol = {c: i for i, c in enumerate(kh.columns)}
    excl_drug_ids = {uf.find(norm(x)) for x in exclude_drugs} | {norm(x) for x in exclude_drugs}
    rows, recs = [], []
    for r in cond.itertuples():
        rec = dict(conditionID=r.conditionID, PMID=r.PMID, Condition=r.Condition, Control=r.Control,
                   minutes=r.minutes, block_id=r.block_id)
        if r.PMID == HIJAZI_PMID:
            inh = str(r.Condition).split("_", 1)[1]
            if inh in hcol:
                vec = H[hcol[inh]].copy()
                rec.update(source="hijazi_kinhib", drugs=inh, doses_nM="", dose_how="binary",
                           drug_ids=norm(inh))
                if norm(inh) in excl_drug_ids or uf.find(norm(inh)) in excl_drug_ids:
                    rec["dropped"] = "exclude_drugs"
                    vec[:] = 0
            else:
                vec = np.zeros(n_k, np.float32)
                rec.update(source="hijazi_kinhib", drugs=inh, dose_how="", dropped="no_kinhib_column")
        else:
            vec = np.zeros(n_k, np.float32)
            parts = []
            for side, txt, sign in (("cond", r.Condition, 1.0), ("ctrl", r.Control, -1.0)):
                dr = find_drugs(txt, lexicon, uf, aliases)
                if not dr:
                    continue
                doses, how = drug_doses(txt, dr, r.PMID)
                occ_not = np.ones(n_k, np.float32)
                for (did, tok), dose in zip(dr, doses):
                    if did in excl_drug_ids:
                        rec["dropped"] = "exclude_drugs"
                    o = occupancy(did, dose, aff_by_drug, aliases, kpos, resolver, n_k)
                    occ_not *= 1 - o
                    parts.append(f"{side}:{tok}@{dose:g}nM[{how}]")
                    if how != "default_1uM":
                        rec["dose_parsed"] = True
                vec += sign * (1 - occ_not)
            rec.update(source="affinity" if parts else "", drugs=" | ".join(parts))
        vec[np.abs(vec) < min_occ] = 0
        rec["n_kinases_occ50"] = int((np.abs(vec) >= 0.5).sum())
        rec["max_occ"] = float(np.abs(vec).max())
        if "dropped" not in rec:
            # dose-series members below 50 % occupancy are kept as low-inhibition anchors (A row may be ~0);
            # conditions without a parsed dose need >= 50 % occupancy of some kinase
            if not rec.get("drugs"):
                rec["dropped"] = "no_drug"
            elif rec["max_occ"] < 0.5 and not rec.get("dose_parsed", False):
                rec["dropped"] = "no_kinase_target_occ50"
            elif not (np.isnan(r.minutes) or r.minutes <= max_minutes):
                rec["dropped"] = f"time>{max_minutes:g}min"
            elif int(r.PMID) in {int(p) for p in exclude_pmids}:
                rec["dropped"] = "exclude_pmids"
            else:
                rec["dropped"] = ""
        rec["time_known"] = not np.isnan(r.minutes)
        recs.append(rec)
        rows.append(vec)
    A = np.stack(rows)
    tab = pd.DataFrame(recs)
    keep = (tab.dropped == "").to_numpy()
    rep = dict(n_conditions_total=len(tab), n_used=int(keep.sum()),
               n_used_hijazi=int((keep & (tab.source == "hijazi_kinhib")).sum()),
               n_used_affinity=int((keep & (tab.source == "affinity")).sum()),
               n_used_unknown_time=int((keep & ~tab.time_known).sum()),
               dropped=tab.dropped.value_counts().to_dict(),
               hijazi_kinhib_rows_unresolved=hij_unres, affinity_sources=aff_rep.to_dict("records"))
    return A, tab, rep, keep


# ----------------------------------------------------------------------------------------------- regression
def assemble_Y(cond_ids: list, n_sites: int, exclude_pmids=()):
    """Site x condition matrix (scaled log2FC) for the given conditions, from the compendium blocks.
    Returns (rows int array of sites with >= 1 value, Y float32 (rows x len(cond_ids)))."""
    pos = {c: i for i, c in enumerate(cond_ids)}
    blocks = load_blocks(kinds=("compendium",), exclude_pmids=exclude_pmids)
    allrows = np.unique(np.concatenate([b.rows for b in blocks]))
    rpos = np.full(n_sites, -1, np.int64)
    rpos[allrows] = np.arange(len(allrows))
    Y = np.full((len(allrows), len(cond_ids)), np.nan, np.float32)
    for b in blocks:
        cj = [(j, pos[c]) for j, c in enumerate(b.cols) if c in pos]
        if not cj:
            continue
        src = np.array([a for a, _ in cj])
        dst = np.array([b_ for _, b_ in cj])
        Y[np.ix_(rpos[b.rows], dst)] = b.X[:, src]
    ok = np.isfinite(Y).any(1)
    return allrows[ok], Y[ok]


def ridge_deconvolution(Y: np.ndarray, A: np.ndarray, lam: float = 1.0, min_obs: int = 3, min_exposure: float = 1.0,
                        d0: float = 4.0, device="cuda", chunk: int = 1024, verbose=True):
    """Per-row ridge of Y (sites x C, NaN) on A (C x K). Returns beta, z, contrast (K x sites float32, NaN where the
    kinase is not exposed in the site's observed conditions or site has < min_obs observations)."""
    import torch
    dev = torch.device(device)
    torch.backends.cuda.matmul.allow_tf32 = True
    At = torch.as_tensor(A, dtype=torch.float32, device=dev)
    C, K = At.shape
    absA = At.abs()
    Iinh = (At >= 0.5).float()
    Inot = (absA < 0.1).float()
    eye = torch.eye(K, device=dev)
    n = Y.shape[0]
    beta = np.full((K, n), np.nan, np.float32)
    zz = np.full((K, n), np.nan, np.float32)
    con = np.full((K, n), np.nan, np.float32)
    s2_all = []
    raw = []
    for c0 in range(0, n, chunk):
        Yt = torch.as_tensor(Y[c0:c0 + chunk], device=dev)
        M = (~torch.isnan(Yt)).float()
        Y0 = torch.nan_to_num(Yt, nan=0.0)
        nobs = M.sum(1)
        G = torch.einsum("sc,ck,cl->skl", M, At, At)
        b = Y0 @ At  # s x K  (= A^T y)
        Gl = G + lam * eye
        L, info = torch.linalg.cholesky_ex(Gl)
        bet = torch.cholesky_solve(b[:, :, None], L)[:, :, 0]
        Ginv = torch.cholesky_inverse(L)
        fit = bet @ At.T
        rss = (((Y0 - fit) * M) ** 2).sum(1)
        df = torch.einsum("skl,slk->s", G, Ginv)
        resid_df = (nobs - df).clamp(min=0.0)
        expo = M @ absA
        cI, cN = M @ Iinh, M @ Inot
        mI = (Y0 @ Iinh) / cI.clamp(min=1)
        mN = (Y0 @ Inot) / cN.clamp(min=1)
        contrast = torch.where((cI >= 2) & (cN >= 3), mI - mN, torch.full_like(mI, float("nan")))
        raw.append((c0, bet.cpu(), torch.diagonal(Ginv, dim1=1, dim2=2).cpu(), rss.cpu(), resid_df.cpu(),
                    nobs.cpu(), expo.cpu(), contrast.cpu()))
        s2 = rss / resid_df.clamp(min=1e-6)
        s2_all.append(s2[(resid_df >= 1) & (nobs >= min_obs)].cpu())
        del G, Gl, L, Ginv
    s0 = float(torch.cat(s2_all).median()) if s2_all else 1.0
    for c0, bet, dg, rss, rdf, nobs, expo, contrast in raw:
        s2 = (rss + d0 * s0) / (rdf + d0)
        se = torch.sqrt(s2[:, None] * dg)
        ok = (nobs[:, None] >= min_obs) & (expo >= min_exposure)
        nanv = torch.tensor(float("nan"))
        sl = slice(c0, c0 + bet.shape[0])
        beta[:, sl] = torch.where(ok, bet, nanv).T.numpy()
        zz[:, sl] = torch.where(ok, bet / se, nanv).T.numpy()
        con[:, sl] = torch.where(nobs[:, None] >= min_obs, contrast, nanv).T.numpy()
    if verbose:
        print(f"ridge: {n} sites, {C} conditions, {K} kinases; s0^2 = {s0:.3f}", flush=True)
    return beta, zz, con, s0


def nonseparable_pairs(A: np.ndarray, kin: pd.DataFrame, thr: float = 0.95, min_exposure: float = 1.0):
    """Kinase pairs with near-identical inhibition profiles (cosine >= thr over the used conditions)."""
    expo = np.abs(A).sum(0)
    idx = np.nonzero(expo >= min_exposure)[0]
    X = A[:, idx]
    Xn = X / np.linalg.norm(X, axis=0, keepdims=True)
    S = Xn.T @ Xn
    np.fill_diagonal(S, -1)
    maxcos = np.full(len(kin), np.nan, np.float32)
    maxcos[idx] = S.max(1)
    partner = np.full(len(kin), "", dtype=object)
    partner[idx] = kin.gene.to_numpy()[idx[S.argmax(1)]]
    ii, jj = np.nonzero(np.triu(S >= thr, 1))
    pairs = pd.DataFrame(dict(kinase_a=kin.gene.to_numpy()[idx[ii]], kinase_b=kin.gene.to_numpy()[idx[jj]],
                              acc_a=kin.uniprot_acc.to_numpy()[idx[ii]], acc_b=kin.uniprot_acc.to_numpy()[idx[jj]],
                              cosine=S[ii, jj], exposure_a=expo[idx[ii]], exposure_b=expo[idx[jj]]))
    kin_tab = pd.DataFrame(dict(kinase=kin.gene, uniprot_acc=kin.uniprot_acc, exposure=expo,
                                n_conditions_occ50=(np.abs(A) >= 0.5).sum(0), max_cosine_other=maxcos,
                                most_similar=partner))
    kin_tab["nonseparable"] = kin_tab.max_cosine_other >= thr
    return pairs.sort_values("cosine", ascending=False), kin_tab


def condition_ids_for_drugs(drugs, path: Path | None = None) -> set:
    """Compendium conditionIDs whose Condition/Control text names any of `drugs` (leave-drug-out; matching by the
    F6b parser: union-find drug identity incl. synonyms/aliases, or the Hijazi inhibitor name)."""
    path = path or RESULTS / "f6b_conditions.tsv"
    tab = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    from ksp.kinases import load_resolver  # noqa: F401  (only for parity with build_affinity; not needed here)
    uf = UnionFind()
    for name, fn in (("klaeger2017", load_klaeger), ("reinecke2024", load_reinecke), ("davis2011", load_davis)):
        _, syn = fn()
        for a, b in syn:
            if norm(a) and norm(b) and len(norm(b)) >= 4:
                uf.union(norm(a), norm(b))
    for r in read_aliases().itertuples():
        if r.maps_to:
            uf.union(norm(r.alias), norm(r.maps_to))
    want = {uf.find(norm(d)) for d in drugs}
    out = set()
    for cid, txt in zip(tab.conditionID, tab.drugs):
        toks = re.findall(r"(?:cond|ctrl):([^@|]+)@", txt) or ([txt] if txt else [])
        if any(uf.find(norm(t)) in want for t in toks):
            out.add(cid)
    return out
