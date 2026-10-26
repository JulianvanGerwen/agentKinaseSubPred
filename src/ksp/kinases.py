"""Human protein kinase list, symbol-synonym map and a name resolver.

Inputs (all local; see DATA.md): UniProt pkinfam + reviewed annotations (src/ksp/uniprot.py), UniProt KW-0418 detail
table (external/uniprot/kinase_detail), KinHub table (external/kinhub, Manning group/family), Manning 2002 Table S1
(external/kinbase_manning2002), HGNC complete set (external/hgnc), Cantley-lab kinase tables (external/cantley_*),
manual curation (config/kinase_curation.tsv).

Rules (also in docs/labels_universe.md):
* candidates = pkinfam human  ∪  reviewed KW-0418 entries with a 'Protein kinase' domain  ∪  KinHub  ∪  curation 'add'
  ∪  kinase genes seen in PSP (to document their exclusion);
* include = True unless excluded by curation category (lipid kinase, regulatory subunit, fusion, metabolic enzyme,
  glycan kinase, disputed activity, not a protein kinase, pseudogene, guanylate-cyclase pseudokinase) or by policy
  flags (exclude_pseudokinases / exclude_uncertain, both False by default); PSP-only genes absent from all protein
  kinase sources are excluded as 'not_in_protein_kinase_sources' unless curated otherwise;
* residue_class: curation override (dual-specificity list, WEE1/2) > Manning group TK (or pkinfam 'Tyr') -> 'Y' >
  'ST'; the rule used is recorded in residue_class_rule;
* is_pseudokinase: curation flag, OR the entry lacks UniProt keyword KW-0418 (UniProt does not assign 'Kinase' to
  catalytically dead family members), OR the UniProt caution text states pseudokinase / lack of kinase
  activity in the CAUTION comment (function text is not used: it mentions kinase-dead mutants); only evaluated
  for members of the protein-kinase sources; the evidence is recorded in pseudokinase_evidence;
* kinase domain: UniProt DOMAIN features whose note matches 'kinase' (Protein kinase, Alpha-type protein kinase,
  PI3K/PI4K catalytic, Histidine kinase), excluding 'AGC-kinase C-terminal'; primary domain = the one containing an
  ACT_SITE, else the last listed kinase domain (JAK JH1 is C-terminal).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import uniprot

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "external"
CURATION = ROOT / "config" / "kinase_curation.tsv"
EXCLUDE_CATEGORIES = {"lipid_kinase", "regulatory_subunit", "fusion", "metabolic_enzyme", "glycan_kinase",
                      "disputed_activity", "not_protein_kinase", "pseudogene", "guanylate_cyclase_pseudokinase"}
_DOM = re.compile(r'DOMAIN (\d+)\.\.(\d+); /note="([^"]+)"')
_ACT = re.compile(r"ACT_SITE (\d+)")
_PSEUDO_TXT = re.compile(r"pseudokinase|lacks? (?:the )?(?:protein )?kinase activity|catalytically inactive|"
                         r"no (?:detectable )?(?:intrinsic )?(?:protein )?kinase activity|"
                         r"devoid of (?:protein )?kinase activity", re.I)


# ---------------------------------------------------------------------------------------------- loaders
def load_curation(path=CURATION) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)


def load_kinhub(path=EXT / "kinhub" / "kinhub_kinases.html") -> pd.DataFrame:
    """KinHub kinase table (xName, Manning name, HGNC name, protein name, group, family, subfamily, UniProt)."""
    html = Path(path).read_text(encoding="utf-8", errors="replace")
    recs = []
    for r in re.findall(r"<tr>(.*?)</tr>", html, re.S):
        tds = re.findall(r"<td>(.*?)</td>", r, re.S)
        if len(tds) == 8:
            recs.append([re.sub(r"<[^>]+>", "", t).strip() for t in tds])
    return pd.DataFrame(recs, columns=["xname", "manning_name", "hgnc", "protein_name", "group", "family",
                                       "subfamily", "acc"])


def load_manning_s1(path=EXT / "kinbase_manning2002" / "TableS1.txt") -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", encoding="latin1", dtype=str, keep_default_na=False)
    d.columns = [re.sub(r"\s+", " ", c.replace('"', "")).strip() for c in d.columns]
    return d[["Name", "Group", "Family", "Subfamily", "Pseudogene?"]].rename(
        columns={"Name": "manning_name", "Group": "group", "Family": "family", "Subfamily": "subfamily",
                 "Pseudogene?": "pseudogene"})


def load_hgnc(path=EXT / "hgnc" / "hgnc_complete_set.txt") -> pd.DataFrame:
    cols = ["hgnc_id", "symbol", "name", "locus_group", "prev_symbol", "alias_symbol", "uniprot_ids",
            "ensembl_gene_id", "entrez_id"]
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=cols)


def load_kinase_detail(path=EXT / "uniprot" / "kinase_detail" / "human_kw0418_detail.tsv.gz") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_cantley_names() -> pd.DataFrame:
    """(acc, matrix_name, gene) from Johnson 2023 (S/T) and Yaron-Barir 2024 (Y) Table S1."""
    out = []
    for f, sheet, a, m, g in [
        (EXT / "cantley_johnson2023" / "41586_2022_5575_MOESM3_ESM.xlsx", "Table S1 Data", "Uniprot id",
         "Matrix_name", "Gene"),
        (EXT / "cantley_yaronbarir2024" / "41586_2024_7407_MOESM3_ESM.xlsx", "Table_S1_Data", "UNIPROT_ID",
         "MATRIX_NAME", "GENENAME")]:
        if not Path(f).exists():
            continue
        d = pd.read_excel(f, sheet_name=sheet, dtype=str)
        out.append(pd.DataFrame({"acc": d[a].str.strip(), "matrix_name": d[m].str.strip(), "gene": d[g].str.strip(),
                                 "source": f.parent.name}))
    return pd.concat(out, ignore_index=True).dropna(subset=["acc", "matrix_name"]) if out else pd.DataFrame(
        columns=["acc", "matrix_name", "gene", "source"])


# ---------------------------------------------------------------------------------------------- helpers
def parse_domains(ft_domain: str) -> list[tuple[int, int, str]]:
    return [(int(a), int(b), n) for a, b, n in _DOM.findall(ft_domain or "")]


def kinase_domains(ft_domain: str, act_site: str = "") -> tuple[list, tuple | None]:
    doms = [d for d in parse_domains(ft_domain) if re.search(r"kinase|PI3K/PI4K catalytic", d[2], re.I)
            and "C-terminal" not in d[2]]
    acts = [int(x) for x in _ACT.findall(act_site or "")]
    primary = None
    for d in doms:
        if any(d[0] <= a <= d[1] for a in acts):
            primary = d
            break
    if primary is None and doms:
        primary = doms[-1]
    return doms, primary


def _split_syn(s: str, sep=r"[\s|;,]+") -> list[str]:
    return [x for x in re.split(sep, s or "") if x]


# ---------------------------------------------------------------------------------------------- build
def build_kinase_list(psp_kinases: pd.DataFrame | None = None, exclude_pseudokinases: bool = False,
                      exclude_uncertain: bool = False, curation: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return the kinase table (one row per gene). `psp_kinases`: optional DataFrame with columns gene (upper-case
    human symbol as mapped from PSP), acc (human accession or ''), n_rows — PSP kinase entries are added as rows
    (usually excluded) so that every PSP 'kinase' has a documented decision."""
    cur = load_curation() if curation is None else curation
    pk = uniprot.parse_pkinfam()
    ann = uniprot.read_annotations()
    ann_sp = ann[ann.Reviewed.eq("reviewed")]
    det = load_kinase_detail()
    kh = load_kinhub()
    hg = load_hgnc()

    # accession-indexed reference tables
    a_gene = dict(zip(ann_sp.Entry, ann_sp["Gene Names (primary)"]))
    a_entry = dict(zip(ann_sp.Entry, ann_sp["Entry Name"]))
    a_pname = dict(zip(ann_sp.Entry, ann_sp["Protein names"]))
    a_kw = dict(zip(ann_sp.Entry, ann_sp["Keyword ID"]))
    a_dom = dict(zip(ann_sp.Entry, ann_sp["Domain [FT]"]))
    sym2acc = {}
    for a, g in zip(ann_sp.Entry, ann_sp["Gene Names (primary)"]):
        for s in _split_syn(g, r";\s*"):
            sym2acc.setdefault(s.upper(), a)
    det_idx = det.set_index("Entry")
    hg_by_acc = {}
    for r in hg.itertuples(index=False):
        for u in _split_syn(r.uniprot_ids, r"\|"):
            hg_by_acc.setdefault(u, r)
    hg_by_sym = {r.symbol.upper(): r for r in hg.itertuples(index=False)}

    # candidate accessions with their source flags
    cand: dict[str, dict] = {}

    def add(acc, src):
        if acc and acc in a_gene:
            cand.setdefault(acc, set()).add(src)

    for a in pk.acc:
        add(a, "pkinfam")
    kw = det[det["Domain [FT]"].str.contains('note="Protein kinase', regex=False)]
    for a in kw.Entry:
        add(a, "kw0418_pkdomain")
    for a in kh.acc:
        add(a, "kinhub")
    for g in cur.loc[cur.action.eq("add"), "gene"]:
        add(sym2acc.get(g.upper()), "curation_add")
    psp_extra = []
    if psp_kinases is not None:
        for r in psp_kinases.itertuples(index=False):
            acc = r.acc if r.acc in a_gene else sym2acc.get(r.gene.upper())
            if acc:
                add(acc, "psp")
            else:
                psp_extra.append(r)

    pk_group = dict(zip(pk.acc, pk.group_short))
    kh_by_acc = kh.drop_duplicates("acc").set_index("acc")
    kh_by_hgnc = kh.drop_duplicates("hgnc").set_index("hgnc")
    cur_by = {(g.upper(), a): r for g, a, r in zip(cur.gene, cur.action, cur.itertuples(index=False))}

    rows = []
    for acc, srcs in sorted(cand.items()):
        h = hg_by_acc.get(acc)
        gene = h.symbol if h is not None else a_gene[acc].split(";")[0]
        G = gene.upper()
        r = {"gene": gene, "uniprot_acc": acc, "entry_name": a_entry.get(acc, ""), "protein_name": a_pname.get(acc, ""),
             "hgnc_id": h.hgnc_id if h is not None else "", "uniprot_gene": a_gene[acc],
             "in_pkinfam": "pkinfam" in srcs, "in_kw0418_pkdomain": "kw0418_pkdomain" in srcs,
             "in_kinhub": "kinhub" in srcs, "curation_add": "curation_add" in srcs, "in_psp": "psp" in srcs,
             "has_kw0418": "KW-0418" in a_kw.get(acc, ""), "pkinfam_group": pk_group.get(acc, "")}
        # Manning group / family
        k = kh_by_acc.loc[acc] if acc in kh_by_acc.index else (kh_by_hgnc.loc[gene] if gene in kh_by_hgnc.index else None)
        if k is not None:
            r.update(manning_group=k.group, manning_family=k.family, manning_subfamily=k.subfamily,
                     manning_name=k.manning_name, kinhub_xname=k.xname, group_source="kinhub")
        else:
            c = cur_by.get((G, "group"))
            if c is not None:
                grp, fam = (c.value.split("/", 1) + [""])[:2]
                r.update(manning_group=grp, manning_family=fam, group_source="curation")
            elif r["pkinfam_group"]:
                pg = r["pkinfam_group"]
                r.update(manning_group="TK" if pg == "Tyr" else ("Atypical" if pg.startswith("Atypical") else pg),
                         manning_family=pg.replace("Atypical-", "") if pg.startswith("Atypical") else "",
                         group_source="pkinfam")
            else:
                r.update(manning_group="", manning_family="", group_source="")
        # kinase domain
        act = det_idx.loc[acc, "Active site"] if acc in det_idx.index else ""
        doms, prim = kinase_domains(a_dom.get(acc, ""), act)
        r["kinase_domains"] = ";".join(f"{a}-{b}:{n}" for a, b, n in doms)
        r["kinase_domain_start"] = prim[0] if prim else pd.NA
        r["kinase_domain_end"] = prim[1] if prim else pd.NA
        r["kinase_domain_name"] = prim[2] if prim else ""
        r["has_act_site_in_kd"] = bool(prim) and any(prim[0] <= int(x) <= prim[1] for x in _ACT.findall(act))
        # pseudokinase (only meaningful for protein-kinase-family members)
        ev = ["curation"] if (G, "pseudokinase") in cur_by else []
        if srcs & {"pkinfam", "kw0418_pkdomain", "curation_add"} or (
                "kinhub" in srcs and r.get("manning_group") != "Atypical"):
            if not r["has_kw0418"]:
                ev.append("no_KW-0418")
            if acc in det_idx.index:
                m = _PSEUDO_TXT.search(det_idx.loc[acc, "Caution"])
                if m:
                    ev.append(f"uniprot_caution:{m.group(0).lower()}")
        r["is_pseudokinase"] = bool(ev)
        r["pseudokinase_evidence"] = ";".join(ev)
        r["activity_uncertain"] = (G, "uncertain") in cur_by
        r["uncertain_reason"] = cur_by[(G, "uncertain")].rationale if r["activity_uncertain"] else ""
        rows.append(r)

    for r0 in psp_extra:  # PSP 'kinases' that are not Swiss-Prot human entries (fusions, ...)
        rows.append({"gene": r0.gene, "uniprot_acc": r0.acc, "in_psp": True, "in_pkinfam": False,
                     "in_kw0418_pkdomain": False, "in_kinhub": False, "curation_add": False, "has_kw0418": False,
                     "is_pseudokinase": False, "activity_uncertain": False, "pseudokinase_evidence": "",
                     "manning_group": "", "manning_family": "", "group_source": "", "kinase_domains": "",
                     "kinase_domain_start": pd.NA, "kinase_domain_end": pd.NA})
    kl = pd.DataFrame(rows)

    # residue class
    rc, rule = [], []
    for g, grp, pg in zip(kl.gene.str.upper(), kl.manning_group.fillna(""), kl.pkinfam_group.fillna("")):
        c = cur_by.get((g, "residue_class"))
        if c is not None:
            rc.append(c.value), rule.append(f"curation: {c.rationale}")
        elif grp == "TK" or pg == "Tyr":
            rc.append("Y"), rule.append("Manning group TK -> Y")
        else:
            rc.append("ST"), rule.append("non-TK group -> ST (default)")
    kl["residue_class"], kl["residue_class_rule"] = rc, rule

    # include / exclusion
    inc, cat, why = [], [], []
    for r in kl.itertuples(index=False):
        G = r.gene.upper()
        c = cur_by.get((G, "exclude"))
        is_pk_source = r.in_pkinfam or r.in_kw0418_pkdomain or r.in_kinhub or r.curation_add
        if c is not None:
            inc.append(False), cat.append(c.category), why.append(c.rationale)
        elif not is_pk_source:
            inc.append(False), cat.append("not_in_protein_kinase_sources"),
            why.append("PSP kinase entry not in pkinfam / KW-0418(+PK domain) / KinHub / curation")
        elif exclude_pseudokinases and r.is_pseudokinase:
            inc.append(False), cat.append("pseudokinase"), why.append(r.pseudokinase_evidence)
        elif exclude_uncertain and r.activity_uncertain:
            inc.append(False), cat.append("activity_uncertain"), why.append(r.uncertain_reason)
        else:
            inc.append(True), cat.append(""), why.append("")
    kl["include"], kl["exclusion_category"], kl["exclusion_reason"] = inc, cat, why
    kl["kinase_domain_start"] = kl["kinase_domain_start"].astype("Int64")
    kl["kinase_domain_end"] = kl["kinase_domain_end"].astype("Int64")
    if kl.gene.duplicated().any():
        dup = kl[kl.gene.duplicated(keep=False)]
        raise ValueError(f"duplicate gene symbols in kinase list: {dup[['gene', 'uniprot_acc']].values.tolist()}")
    cols = ["gene", "uniprot_acc", "include", "exclusion_category", "exclusion_reason", "residue_class",
            "residue_class_rule", "is_pseudokinase", "pseudokinase_evidence", "activity_uncertain", "uncertain_reason",
            "manning_group", "manning_family", "manning_subfamily", "manning_name", "kinhub_xname", "group_source",
            "pkinfam_group", "kinase_domain_start", "kinase_domain_end", "kinase_domain_name", "kinase_domains",
            "has_act_site_in_kd", "in_pkinfam", "in_kw0418_pkdomain", "in_kinhub", "curation_add", "in_psp",
            "has_kw0418", "hgnc_id", "uniprot_gene", "entry_name", "protein_name"]
    return kl[cols].sort_values("gene").reset_index(drop=True)


def build_synonyms(kl: pd.DataFrame, psp_names: pd.DataFrame | None = None) -> pd.DataFrame:
    """Synonym table: synonym, synonym_uc, gene, uniprot_acc, source, n_genes (distinct kinase-list genes sharing the
    synonym), approved_symbol_of_other_gene (HGNC symbol of a different gene), use_for_join.
    Sources: gene itself, HGNC prev/alias symbols, UniProt gene names, KinHub xName / Manning name, Cantley matrix
    names, PSP kinase names (`psp_names`: columns name, gene).
    use_for_join = the synonym is the gene's own symbol, OR (it points to exactly one kinase-list gene AND it is not
    the approved HGNC symbol of another gene)."""
    hg = load_hgnc()
    ann = uniprot.read_annotations()
    hg_by_id = hg.set_index("hgnc_id")
    gn = dict(zip(ann.Entry, ann["Gene Names"]))
    recs = []
    for r in kl.itertuples(index=False):
        if not isinstance(r.uniprot_acc, str):
            continue
        recs.append((r.gene, r.gene, r.uniprot_acc, "symbol"))
        if r.hgnc_id and r.hgnc_id in hg_by_id.index:
            h = hg_by_id.loc[r.hgnc_id]
            for s in _split_syn(h.prev_symbol, r"\|"):
                recs.append((s, r.gene, r.uniprot_acc, "hgnc_prev"))
            for s in _split_syn(h.alias_symbol, r"\|"):
                recs.append((s, r.gene, r.uniprot_acc, "hgnc_alias"))
        for s in _split_syn(gn.get(r.uniprot_acc, ""), r"[\s;]+"):
            recs.append((s, r.gene, r.uniprot_acc, "uniprot_gene_names"))
        for col, src in [("kinhub_xname", "kinhub_xname"), ("manning_name", "manning_name")]:
            v = getattr(r, col)
            if isinstance(v, str) and v:
                recs.append((v, r.gene, r.uniprot_acc, src))
        if isinstance(r.entry_name, str) and r.entry_name:
            recs.append((r.entry_name.replace("_HUMAN", ""), r.gene, r.uniprot_acc, "uniprot_entry_name"))
    acc2gene = dict(zip(kl.uniprot_acc, kl.gene))
    cn = load_cantley_names()
    for a, m in zip(cn.acc, cn.matrix_name):
        if a in acc2gene and isinstance(m, str):
            recs.append((m, acc2gene[a], a, "cantley_matrix_name"))
    if psp_names is not None:
        g2a = dict(zip(kl.gene.str.upper(), kl.uniprot_acc))
        g2g = dict(zip(kl.gene.str.upper(), kl.gene))
        for n, g in zip(psp_names.name, psp_names.gene):
            if g.upper() in g2g:
                recs.append((n, g2g[g.upper()], g2a[g.upper()], "psp_kinase_name"))
    syn = pd.DataFrame(recs, columns=["synonym", "gene", "uniprot_acc", "source"])
    syn["synonym"] = syn.synonym.str.strip()
    syn = syn[syn.synonym.ne("")]
    syn["synonym_uc"] = syn.synonym.str.upper()
    syn = (syn.groupby(["synonym_uc", "gene", "uniprot_acc"], as_index=False)
           .agg(synonym=("synonym", "first"), source=("source", lambda s: ";".join(sorted(set(s))))))
    syn["n_genes"] = syn.groupby("synonym_uc").gene.transform("nunique")
    approved = dict(zip(hg.symbol.str.upper(), hg.symbol))
    syn["approved_symbol_of_other_gene"] = [approved.get(s, "") if approved.get(s, "").upper() != g.upper() and
                                            s in approved else "" for s, g in zip(syn.synonym_uc, syn.gene)]
    own = syn.synonym_uc.eq(syn.gene.str.upper())
    syn["use_for_join"] = own | (syn.n_genes.eq(1) & syn.approved_symbol_of_other_gene.eq(""))
    # a name that is the own symbol of one kinase but a synonym of another: only the own-symbol row joins
    own_syms = set(syn.loc[own, "synonym_uc"])
    syn.loc[~own & syn.synonym_uc.isin(own_syms), "use_for_join"] = False
    return syn[["synonym", "synonym_uc", "gene", "uniprot_acc", "source", "n_genes", "approved_symbol_of_other_gene",
                "use_for_join"]].sort_values(["gene", "synonym_uc"]).reset_index(drop=True)


class KinaseResolver:
    """Resolve a kinase name / symbol / UniProt accession to the kinase-list gene symbol.
    resolve(name, acc) -> (gene | None, how), how in {'acc','symbol','synonym','family_name','ambiguous','unknown'}.
    Order: accession; own symbol; curated family/group names (config action 'family_name', e.g. PKA, CK2 -> never a
    gene); unambiguous synonyms (use_for_join); otherwise ambiguous/unknown.
    resolve_note(text) handles free-text names as in UniProt 'by PKB/AKT1', 'by GSK3-beta': whole string, then the
    normalised string ('-alpha'/'-beta' suffix -> A/B), then '/'-separated parts — accepted only if all resolvable
    parts give one and the same gene and no part is a family name that disagrees."""

    def __init__(self, kl: pd.DataFrame, syn: pd.DataFrame, family_names=None):
        self.kl = kl
        self.acc2gene = {a: g for a, g in zip(kl.uniprot_acc, kl.gene) if isinstance(a, str) and a}
        self.sym = {g.upper(): g for g in kl.gene}
        j = syn[syn.use_for_join]
        self.syn = dict(zip(j.synonym_uc, j.gene))
        self.ambig = set(syn.loc[~syn.use_for_join, "synonym_uc"]) - set(self.syn)
        if family_names is None:
            cur = load_curation()
            family_names = cur.loc[cur.action.eq("family_name"), "gene"]
        self.family = {str(f).upper() for f in family_names} - set(self.sym)
        self.included = dict(zip(kl.gene, kl.include))

    def resolve(self, name: str | None = None, acc: str | None = None) -> tuple[str | None, str]:
        if acc:
            base = str(acc).split("-", 1)[0]
            if base in self.acc2gene:
                return self.acc2gene[base], "acc"
        if not name:
            return None, "unknown"
        u = str(name).strip().upper()
        if u in self.sym:
            return self.sym[u], "symbol"
        if u in self.family:
            return None, "family_name"
        if u in self.syn:
            return self.syn[u], "synonym"
        if u in self.ambig:
            return None, "ambiguous"
        return None, "unknown"

    @staticmethod
    def _normalise(name: str) -> str:
        n = re.sub(r"[-\s]alpha$", "A", name, flags=re.I)
        n = re.sub(r"[-\s]beta$", "B", n, flags=re.I)
        n = re.sub(r"[-\s]gamma$", "G", n, flags=re.I)
        n = re.sub(r"[-\s]delta$", "D", n, flags=re.I)
        return n

    def resolve_note(self, name: str) -> tuple[str | None, str]:
        g, how = self.resolve(name)
        if g is not None or how == "family_name":
            return g, how
        n = self._normalise(name)
        if n != name:
            g2, how2 = self.resolve(n)
            if g2 is not None or how2 == "family_name":
                return g2, how2 + "(normalised)"
        if "/" in name:
            parts = [self.resolve(self._normalise(p.strip())) for p in name.split("/") if p.strip()]
            genes = {gp for gp, _ in parts if gp is not None}
            if len(genes) == 1:
                return genes.pop(), "slash_parts"
            if len(genes) > 1:
                return None, "ambiguous"
        return None, how


def read_kinase_list(path=ROOT / "results" / "labels" / "kinase_list.tsv") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"kinase_domain_start": "Int64", "kinase_domain_end": "Int64"},
                       keep_default_na=False, na_values={"kinase_domain_start": [""], "kinase_domain_end": [""]})


def read_synonyms(path=ROOT / "results" / "labels" / "kinase_synonyms.tsv") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).assign(
        use_for_join=lambda d: d.use_for_join.eq("True"), n_genes=lambda d: d.n_genes.astype(int))


def load_resolver(results_dir=ROOT / "results" / "labels") -> KinaseResolver:
    kl = read_kinase_list(Path(results_dir) / "kinase_list.tsv")
    kl["include"] = kl.include.astype(str).eq("True")
    return KinaseResolver(kl, read_synonyms(Path(results_dir) / "kinase_synonyms.tsv"))


# ------------------------------------------------------------------------------------------------ strict bulk resolver
HOMONYMS = ROOT / "config" / "kinase_homonyms.tsv"
COMPLEX = ROOT / "config" / "complex_subunits.tsv"


def load_homonyms(path=HOMONYMS) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")


def load_complex_subunits(path=COMPLEX) -> dict:
    d = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    return {s.upper(): c.split() for s, c in zip(d.subunit, d.catalytic)}


def family_members(resolver: "KinaseResolver | None" = None) -> dict:
    """Curated family name (upper) -> list of member kinase-list genes (from kinase_curation.tsv 'family_name'
    values; range notation like 'PLK1-4' and wildcards like 'MAP3K*' / 'EPHA*' are expanded against the list)."""
    import re as _re
    res = resolver or load_resolver()
    genes = sorted(res.kl.gene[res.kl.include.astype(bool)])
    cur = load_curation()
    out = {}
    for name, val in zip(cur.gene[cur.action.eq("family_name")], cur.value[cur.action.eq("family_name")]):
        mem = []
        for tok in str(val).replace("/", " ").split():
            m = _re.fullmatch(r"([A-Z0-9]*?[A-Z])(\d+)-(\d+)", tok)
            if tok.endswith("*"):
                mem += [g for g in genes if g.startswith(tok[:-1])]
            elif m:
                mem += [f"{m.group(1)}{i}" for i in range(int(m.group(2)), int(m.group(3)) + 1)]
            elif tok in genes:
                mem.append(tok)
        out[str(name).upper()] = [g for g in dict.fromkeys(mem) if g in genes]
    return out


def resolve_kinase_names(names, source: str, resolver: "KinaseResolver | None" = None, accs=None,
                         on_unknown: str = "collect") -> tuple[pd.DataFrame, dict]:
    """THE strict bulk name resolver (DESIGN §16 A4). For each name (optionally with a UniProt accession, which wins)
    returns gene, uniprot_acc, how in {acc, homonym_rule, symbol, synonym, family_name, homonym_drop, ambiguous,
    unknown}, and family members for family-level names (mask-only use). Unknown/ambiguous names are collected and
    counted (on_unknown='collect') or raise (on_unknown='raise'). Per-source homonym rules: config/kinase_homonyms.tsv.
    """
    res = resolver or load_resolver()
    hom = load_homonyms()
    rule = {}
    for n, s, t in zip(hom.name.str.upper(), hom.source, hom.target):
        if s == source or (s == "*" and (n, source) not in rule):
            rule[(n, source)] = t
    fam = family_members(res)
    acc_of = dict(zip(res.kl.gene, res.kl.uniprot_acc))
    names = list(names)
    accs = list(accs) if accs is not None else [None] * len(names)
    rows, counts = [], {}
    for nm, ac in zip(names, accs):
        key = str(nm).strip().upper() if nm is not None else ""
        g, how, members = None, "unknown", ""
        if ac:
            g, how = res.resolve(acc=ac)
        if g is None:
            t = rule.get((key, source))
            if t is not None:
                if t == "DROP":
                    how = "homonym_drop"
                elif t.startswith("FAMILY:"):
                    how, members = "family_name", " ".join(fam.get(t[7:].upper(), []))
                else:
                    g, how = t, "homonym_rule"
            else:
                g, how = res.resolve(nm)
                if g is None and how not in ("family_name", "ambiguous"):
                    g, how = res.resolve_note(nm)          # normalised names (GSK3-beta, PKB/AKT1, ...)
                if how == "family_name":
                    members = " ".join(fam.get(key, fam.get(res._normalise(str(nm)).upper(), [])))
        counts[how] = counts.get(how, 0) + 1
        rows.append((nm, ac, g, acc_of.get(g) if g else None, how, members))
    out = pd.DataFrame(rows, columns=["name", "acc_in", "gene", "uniprot_acc", "how", "family_members"])
    bad = out[out.how.isin(["unknown", "ambiguous"])]
    if len(bad) and on_unknown == "raise":
        raise KeyError(f"{len(bad)} unresolved kinase names from {source}: {sorted(set(bad.name))[:20]}")
    return out, {"source": source, "counts": counts, "unresolved": sorted(set(map(str, bad.name)))}
