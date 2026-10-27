"""Build the manuscript figures from committed summary files only.

Run from anywhere:  .venv/bin/python -I manuscript/figures/make_figures.py

Inputs (all committed, metrics only; no label files, nothing under sealed/):
  results/test/721dc1977366/<config>_ci.tsv        standalone TEST endpoints with cluster-bootstrap CIs
  results/test/721dc1977366/superiority.tsv        the 9 pre-registered paired tests
  results/test/721dc1977366/<config>_per_kinase.tsv per-kinase TEST metrics (AUROC, training-label count)
  results/test/721dc1977366/<config>_strata.tsv    TEST strata (descriptive, no CIs)
  results/test/721dc1977366/final_reliability.tsv  assignment_prob reliability on labelled TEST sites
  results/release/coverage.tsv                     deployment coverage counts
Outputs: Fig1..Fig5 as PNG (300 dpi) and PDF next to this script.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEST = ROOT / "results" / "test" / "721dc1977366"
REL = ROOT / "results" / "release"

# Palette: validated categorical slots 1-3 (dataviz reference palette, light mode) + neutral inks.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
LIGHTGRAY = "#c3c2b7"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
})


def save(fig, name, bbox="tight"):
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=300, bbox_inches=bbox)
    plt.close(fig)


def ci(config):
    d = pd.read_csv(TEST / f"{config}_ci.tsv", sep="\t").set_index("endpoint")
    return d


def hgrid(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def vgrid(ax):
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------- Fig 1: method schematic
LINE = 0.0195  # vertical space per body line in figure units (figure 5.8 in tall, 5.8 pt font)


def box(ax, x, y, w, h, title, lines, fc="#f4f3ef", ec=AXIS, tc=INK, bold_color=None, fs=5.8):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.01",
                       fc=fc, ec=bold_color or ec, lw=1.2 if bold_color else 0.8)
    ax.add_patch(p)
    ax.text(x + 0.01, y + h - 0.012, title, ha="left", va="top", fontsize=fs + 1.2,
            fontweight="bold", color=tc)
    ax.text(x + 0.01, y + h - 0.040, "\n".join(lines), ha="left", va="top",
            fontsize=fs, color=INK2, linespacing=1.3)


def bh(n):
    return 0.050 + n * LINE


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, shrinkA=0, shrinkB=0))


def column(ax, x, w, top, specs, gap=0.018):
    """Stack boxes downward from `top`; return list of (y, h)."""
    out = []
    y = top
    for title, lines, kw in specs:
        h = bh(len(lines))
        y -= h
        box(ax, x, y, w, h, title, lines, **kw)
        out.append((y, h))
        y -= gap
    return out


def fig1():
    fig = plt.figure(figsize=(7.5, 5.8))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    W = 0.30
    X = [0.01, 0.35, 0.69]
    ax.text(0.01, 0.985, "Prediction pipeline", fontsize=8.5, fontweight="bold", color=INK, va="top")
    c1 = column(ax, X[0], W, 0.95, [
        ("Labels (positives)", ["PhosphoSitePlus kinase-substrate pairs",
                                "L2: human + mouse/rat mapped by 15-mer",
                                "  (18,638 pairs, 401 kinases)",
                                "No true negatives: PU setting"], {}),
        ("Candidates and negatives", ["503 kinases, residue class ST / Y / STY",
                                      "Sites: U_ms (seen by mass spectrometry)",
                                      "Negatives of kinase k: U_ms sites minus",
                                      "  known positives, paralog mask Par(k),",
                                      "  15-mers within Hamming 2"], {}),
        ("Features (kinase x site)", ["F0 site: AlphaFold, ESM-2 650M",
                                      "F1 Cantley atlas scores, percentiles",
                                      "F2a learned PSSM*",
                                      "F3 PPI screens; F3-label*",
                                      "F4 localisation   F5 co-expression",
                                      "F6a co-regulation*   F6b inhibitor",
                                      "  deconvolution   F6c CPTAC kinase state",
                                      "F7 kinase descriptors",
                                      "* label-derived, cross-fitted by cluster"], {}),
    ])
    c2 = column(ax, X[1], W, 0.95, [
        ("R2: global LightGBM pair model", ["One model for all (kinase, site) pairs",
                                            "Binary loss, 127 leaves, lr 0.03",
                                            "Weights n_k^-0.5; 1:50 negatives,",
                                            "  50% E-hard type",
                                            "Early stopping on inner clusters",
                                            "5 seeds, mean logit"], dict(bold_color=BLUE)),
        ("R3: neural kinase x site model", ["Site: frozen ESM-2 window (+/-7),",
                                            "  phospho channel, 2-layer transformer",
                                            "Kinase: ESM-2 domain, Cantley profile,",
                                            "  Manning group/family (dropout 0.2)",
                                            "Bilinear + MLP head; BCE + 0.5 InfoNCE",
                                            "Pretrained by distilling label-free",
                                            "  Cantley percentiles; aux loss 0.2",
                                            "5 seeds"], dict(bold_color=BLUE)),
    ])
    c3 = column(ax, X[2], W, 0.95, [
        ("Late fusion", ["z-score each model's logits",
                         "score = 0.6 z(R2) + 0.4 z(R3)  (CV M1)",
                         "score_site = 0.5 z(R2) + 0.5 z(R3)  (CV M2)"],
         dict(fc="#eaf2fb", bold_color=BLUE)),
        ("Calibration (fitted on DEV out-of-fold)", ["Tiered Platt: kinase / family / global",
                                                     "Softmax over kinases at a site (T 0.93)",
                                                     "Isotonic map of the top kinase,",
                                                     "  separately for S/T and Y sites",
                                                     "  -> assignment_prob"], {}),
        ("Outputs", ["score: rank sites for a kinase (M1, M3)",
                     "score_site: rank kinases for a site (M2)",
                     "assignment_prob: coverage, reliability",
                     "Deployment: refit on DEV+TEST,",
                     "  score all U_obs sites"], {}),
    ])
    # flow arrows (gutters only)
    ymid = lambda c, i: c[i][0] + c[i][1] / 2
    for i in range(3):
        arrow(ax, X[0] + W, ymid(c1, i), X[1], ymid(c2, 0 if i < 2 else 1))
    arrow(ax, X[1] + W, ymid(c2, 0), X[2], ymid(c3, 0))
    arrow(ax, X[1] + W, ymid(c2, 1), X[2], ymid(c3, 0) - 0.02)
    for i in range(2):
        arrow(ax, X[2] + W / 2, c3[i][0], X[2] + W / 2, c3[i][0] - 0.018)

    # evaluation band
    top = min(c1[-1][0], c2[-1][0], c3[-1][0]) - 0.05
    ax.plot([0.01, 0.99], [top + 0.025, top + 0.025], color=AXIS, lw=0.8)
    ax.text(0.01, top + 0.012, "Evaluation design", fontsize=8.5, fontweight="bold", color=INK, va="top")
    top -= 0.03
    column(ax, X[0], W, top, [("Splits (frozen before modelling)", [
        "Substrate clusters (mmseqs2 >=40% id, >=50%",
        "  cov; + identical 15-mers): 11,609",
        "TEST: 20.3% of clusters, labels sealed",
        "DEV: 5-fold CV by cluster (T1)",
        "Kinase-family-out: T2 (CV), T3 (TEST)"], {})])
    column(ax, X[1], W, top, [("Endpoints, kinase-equal macro", [
        "M1 AUROC vs E-hard negatives",
        "M2 Hit@5 of the true kinase at a site",
        "M3 recall in the top 1% of E-all",
        "Floors B0-prior / pop / miss; canaries",
        "9 pre-registered tests, Holm-corrected"], {})])
    column(ax, X[2], W, top, [("One sealed TEST batch", [
        "FREEZE.json pins each prediction file",
        "Ledger (STARTED / COMPLETED) refuses",
        "  a second batch unless a logged",
        "  override is used",
        "Batch reproduced by an independent",
        "  audit agent"], dict(fc="#fdf0ea", bold_color=ORANGE))])
    from matplotlib.transforms import Bbox
    save(fig, "Fig1_method", bbox=Bbox([[0, 5.8 * (top - bh(6) - 0.02)], [7.5, 5.8]]))


# ---------------------------------------------------------------- Fig 2: TEST endpoints + primary tests
CONFIGS = [
    ("final", "Final (R2+R3 fusion)"),
    ("R2", "R2 alone"),
    ("R3", "R3 alone"),
    ("B7", "B7 simple-strong LR"),
    ("B1", "B1 Cantley atlas"),
    ("B3", "B3 Phosformer (our retrain)"),
    ("B4", "B4 SELPHI2.0 (our retrain)"),
    ("B2", "B2 learned PSSM"),
    ("B0_pop", "B0-pop (site popularity)"),
]


def endpoint_value(config, ep):
    # B1: arm (i) for M1/M3, arm (ii) for M2 (the stronger arm per endpoint, A23)
    if config == "B1":
        config = "B1_ii" if ep == "M2" else "B1_i"
    d = ci(config)
    return d.loc[ep, "point"], d.loc[ep, "lo95"], d.loc[ep, "hi95"]


def fig2():
    fig = plt.figure(figsize=(7.5, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.45, wspace=0.12,
                          left=0.21, right=0.98, top=0.95, bottom=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    eps = [("M1", "M1: macro AUROC (E-hard)", (0.5, 0.85)),
           ("M2", "M2: site-centric Hit@5", (0.0, 0.42)),
           ("M3", "M3: top-1% recall (E-all)", (0.0, 0.42))]
    y = np.arange(len(CONFIGS))[::-1]
    for ax, (ep, title, xl) in zip(axes, eps):
        for yi, (cfg, _) in zip(y, CONFIGS):
            p, lo, hi = endpoint_value(cfg, ep)
            col = BLUE if cfg == "final" else (LIGHTGRAY if cfg == "B0_pop" else MUTED)
            ax.plot([lo, hi], [yi, yi], color=col, lw=2, solid_capstyle="round")
            ax.plot(p, yi, "o", ms=5.5, color=col, mec="white", mew=1.0, zorder=3)
            if cfg == "final":
                ax.text(hi + 0.015 * (xl[1] - xl[0]), yi, f"{p:.3f}", va="center", fontsize=7, color=INK)
        ax.set_xlim(*xl)
        ax.set_title(title, loc="left", fontsize=8)
        vgrid(ax)
        ax.set_yticks(y)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if ax is axes[0]:
            ax.set_yticklabels([lab for _, lab in CONFIGS])
        else:
            ax.set_yticklabels([])
    # panel D: the 9 pre-registered paired differences
    ax = fig.add_subplot(gs[1, :2])
    s = pd.read_csv(TEST / "superiority.tsv", sep="\t")
    s = s[s["primary"]]
    lab = {"B1_i": "Cantley B1 (i)", "B1_ii": "Cantley B1 (ii)", "B3": "B3 Phosformer", "B4": "B4 SELPHI2.0"}
    rows = [(f"{r['endpoint']} vs {lab[r['b']]}", r["diff"], r["lo95"], r["hi95"], r["p_holm"])
            for _, r in s.iterrows()]
    order = [f"M1 vs {lab['B1_i']}", f"M1 vs {lab['B3']}", f"M1 vs {lab['B4']}",
             f"M2 vs {lab['B1_ii']}", f"M2 vs {lab['B3']}", f"M2 vs {lab['B4']}",
             f"M3 vs {lab['B1_i']}", f"M3 vs {lab['B3']}", f"M3 vs {lab['B4']}"]
    rows = sorted(rows, key=lambda t: order.index(t[0]))
    yy = np.arange(len(rows))[::-1]
    for yi, (name, d, lo, hi, ph) in zip(yy, rows):
        col = ORANGE if lo < 0.02 else BLUE
        ax.plot([lo, hi], [yi, yi], color=col, lw=2, solid_capstyle="round")
        ax.plot(d, yi, "o", ms=5.5, color=col, mec="white", mew=1.0, zorder=3)
        ax.text(0.315, yi, f"{d:+.3f} [{lo:.3f}, {hi:.3f}]   Holm p {ph:.4f}", va="center",
                fontsize=6.6, color=INK2)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.02, 0.30)
    ax.set_xlabel("Paired difference, final minus baseline (matched universe, B = 10,000)")
    vgrid(ax)
    ax.set_title("Pre-registered primary tests (orange: the one modest test, CI lower bound below 0.02)", loc="left", fontsize=8)
    for a, letter in zip(axes + [ax], "ABCD"):
        a.text(-0.03, 1.07, letter, transform=a.transAxes, fontsize=10, fontweight="bold", color=INK,
               ha="right")
    save(fig, "Fig2_test_endpoints")


# ---------------------------------------------------------------- Fig 3: per-kinase + strata
def fig3():
    fin = pd.read_csv(TEST / "final_per_kinase.tsv", sep="\t")
    pop = pd.read_csv(TEST / "B0_pop_per_kinase.tsv", sep="\t")[["kinase_acc", "auroc_hard"]]
    fin = fin[fin["in_K_primary"]].merge(pop, on="kinase_acc", suffixes=("", "_pop"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.5, 3.4), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.55))
    for cls, col, lab in [("ST", BLUE, "S/T kinases"), ("STY", BLUE, None), ("Y", ORANGE, "Tyr kinases")]:
        d = fin[fin["residue_class"] == cls]
        a.scatter(d["n_train_dev_L2"], d["auroc_hard"], s=22, color=col, edgecolor="white",
                  linewidth=0.8, label=lab, zorder=3, marker="o" if cls != "Y" else "s")
    a.axhline(0.5, color=AXIS, lw=0.8)
    a.set_xscale("log")
    a.set_xlabel("Training labels of the kinase (DEV, L2)")
    a.set_ylabel("TEST AUROC, E-hard negatives")
    a.set_ylim(0.1, 1.02)
    hgrid(a)
    a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7,
             handletextpad=0.3)
    n_y = int((fin["residue_class"] == "Y").sum())
    beat = int(((fin["residue_class"] == "Y") & (fin["auroc_hard"] > fin["auroc_hard_pop"])).sum())
    a.set_title(f"Per-kinase TEST AUROC (n = {len(fin)})", loc="left", fontsize=8)
    a.text(0.02, 0.03, f"Tyr kinases above their B0-pop AUROC: {beat} of {n_y}", transform=a.transAxes,
           fontsize=6.8, color=INK2)

    # panel b: strata, final vs Cantley arm (i) vs popularity floor (descriptive, no CIs)
    rows = [("kinase", "residue_class", "ST", "S/T kinases"),
            ("kinase", "residue_class", "Y", "Tyr kinases"),
            ("kinase", "train_bin", "5-19", "5-19 labels"),
            ("kinase", "train_bin", "20-99", "20-99 labels"),
            ("kinase", "train_bin", ">=100", ">=100 labels"),
            ("pair", "primed", "primed", "Primed positives*"),
            ("pair", "primed", "unprimed", "Unprimed positives*"),
            ("pair", "hub", "hub", "Hub sites"),
            ("pair", "hub", "non_hub", "Non-hub sites"),
            ("pair", "evidence", "in_vivo_only", "In vivo only"),
            ("pair", "evidence", "in_vitro_only", "In vitro only"),
            ("pair", "temporal", "new_since_2021", "New since 2021")]
    series = [("final", "Final", BLUE), ("B1_i", "Cantley (i)", ORANGE), ("B0_pop", "B0-pop", LIGHTGRAY)]
    tabs = {c: pd.read_csv(TEST / f"{c}_strata.tsv", sep="\t", dtype={"value": str}) for c, _, _ in series}
    yb = np.arange(len(rows))[::-1].astype(float)
    off = {"final": 0.0, "B1_i": -0.12, "B0_pop": 0.12}
    for c, lab, col in series:
        t = tabs[c]
        vals = []
        for lvl, st, v, _ in rows:
            m = t[(t["level"] == lvl) & (t["stratum"] == st) & (t["value"] == v)]
            vals.append(float(m["M1"].iloc[0]))
        b.scatter(vals, yb + off[c], s=26 if c == "final" else 20, color=col, edgecolor="white",
                  linewidth=0.8, zorder=3 if c == "final" else 2, label=lab,
                  marker={"final": "o", "B1_i": "D", "B0_pop": "s"}[c])
    b.set_yticks(yb)
    b.set_yticklabels([r[3] for r in rows])
    b.tick_params(axis="y", length=0)
    b.spines["left"].set_visible(False)
    b.axvline(0.5, color=AXIS, lw=0.8)
    b.set_xlim(0.40, 0.95)
    b.set_xlabel("TEST M1 within stratum")
    vgrid(b)
    b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=3, fontsize=7,
             handletextpad=0.2, columnspacing=1.0)
    b.set_title("TEST strata (descriptive, no CIs)", loc="left", fontsize=8)
    a.text(-0.16, 1.06, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
    fig.text(0.44, 0.925, "B", fontsize=10, fontweight="bold")
    save(fig, "Fig3_kinases_strata")


# ---------------------------------------------------------------- Fig 4: reliability + coverage
def fig4():
    r = pd.read_csv(TEST / "final_reliability.tsv", sep="\t")
    r = r[r["subset"] == "labelled_TEST_sites"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.5, 3.1), gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.6))
    a.plot([0, 1], [0, 1], color=AXIS, lw=1)
    a.plot(r["mean_p"], r["accuracy"], color=BLUE, lw=2, zorder=2)
    a.scatter(r["mean_p"], r["accuracy"], s=8 + r["n_sites"] / 4, color=BLUE, edgecolor="white",
              linewidth=0.8, zorder=3)
    for i, (_, x) in enumerate(r.iterrows()):
        dy = 0.03 if i == len(r) - 1 else -0.045
        a.text(x["mean_p"] + 0.02, x["accuracy"] + dy, f"n={int(x['n_sites'])}", fontsize=5.8, color=MUTED)
    a.axvspan(0.8, 1.0, color="#f4f3ef", zorder=0)
    a.text(0.9, 0.2, "no site\nreaches 0.8", ha="center", fontsize=6.5, color=INK2)
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.set_xlabel("Claimed assignment_prob of top kinase (bin mean)")
    a.set_ylabel("Observed: top kinase is a true kinase")
    ece = float(r["ECE"].iloc[0])
    acc = float(r["accuracy_overall"].iloc[0])
    a.set_title("Reliability on labelled TEST sites", loc="left", fontsize=8)
    a.text(0.03, 0.92, f"n = {int(r['n_total'].iloc[0]):,} sites; ECE {ece:.3f}; top-1 accuracy {acc:.3f}", transform=a.transAxes, fontsize=6.8,
           color=INK2)
    hgrid(a)

    c = pd.read_csv(REL / "coverage.tsv", sep="\t")
    name = {"R2+R3 final (deployment)": "Final (deployment)",
            "Cantley-only (B1 imputed-universe pct)": "Cantley only",
            "random scores": "Random scores"}
    c["label"] = c["method"].map(name)
    c = c.set_index("label").loc[["Final (deployment)", "Cantley only", "Random scores"]]
    yy = np.arange(len(c))[::-1]
    h = 0.34
    b.barh(yy + h / 2 + 0.02, c["coverage@0.5"] / 1000, height=h, color=BLUE, label="assignment_prob >= 0.5")
    b.barh(yy - h / 2 - 0.02, c["coverage@0.8"] / 1000, height=h, color=ORANGE, label="assignment_prob >= 0.8")
    for yi, (_, row) in zip(yy, c.iterrows()):
        b.text(row["coverage@0.5"] / 1000 + 0.6, yi + h / 2 + 0.02, f"{int(row['coverage@0.5']):,}",
               va="center", fontsize=6.8, color=INK)
        b.text(row["coverage@0.8"] / 1000 + 0.6, yi - h / 2 - 0.02, f"{int(row['coverage@0.8']):,}",
               va="center", fontsize=6.8, color=INK)
    b.set_yticks(yy)
    b.set_yticklabels(c.index)
    b.tick_params(axis="y", length=0)
    b.spines["left"].set_visible(False)
    b.set_xlim(0, 38)
    b.set_xlabel("U_obs sites covered (thousands of 350,736)")
    b.set_title("Coverage counts (model-based, not precision)", loc="left", fontsize=8)
    vgrid(b)
    b.legend(frameon=False, loc="lower right", fontsize=6.8)
    a.text(-0.2, 1.06, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
    b.text(-0.5, 1.06, "B", transform=b.transAxes, fontsize=10, fontweight="bold")
    save(fig, "Fig4_assignment")


# ---------------------------------------------------------------- Fig 5: agent setup
def fig5():
    H = 5.8
    fig = plt.figure(figsize=(7.5, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    hum = ["Wrote the brief (PROMPT.md); supplied data, two papers and",
           "PhosphoSitePlus access; answered data requests (DATAINEED.md)",
           "Did not write code, run analyses or write text"]
    box(ax, 0.25, 0.97 - bh(3), 0.50, bh(3), "Human researcher", hum, fc="#fdf0ea", bold_color=ORANGE)
    pi = ["Owns DESIGN.md (amendments A1-A34), RESEARCH_TREE.md, LESSONS.md,",
          "DELEGATION.md; the only agent that runs git",
          "Delegates with written briefs, file ownership, GPU/core budgets",
          "Decides adoption by pre-declared rules; commissions independent recomputation"]
    ypi = 0.97 - bh(3) - 0.04 - bh(4)
    box(ax, 0.25, ypi, 0.50, bh(4), "Lead agent ('PI')", pi, fc="#eaf2fb", bold_color=BLUE)
    arrow(ax, 0.5, 0.97 - bh(3), 0.5, ypi + bh(4))
    cols = [
        ("Literature and\ndesign review",
         ["4 literature agents", "Round 1: 7 reviewers", "  (external PI, prior", "  work, statistics,",
          "  data leakage, compute,", "  AI failure modes,", "  kinase biology)", "Round 2: 4 reviewers"]),
        ("Implementation\n(one owner per path)",
         ["labels/universe", "core-infra: splits,", "  seal, eval, canaries", "features: motif,",
          "  context, phospho,", "  structure+ESM,", "  knn, Tyr", "models: R2, R3, R5", "baselines: simple,",
          "  B3, B4", "strata, regimes,", "  downstream"]),
        ("Independent checks\n(not the producer)",
         ["CV audit r1", "CV audit r2", "TEST-harness code", "  review (A29)", "TEST batch audit",
          "Post-TEST audit", "REPORT consistency", "  review"]),
        ("Process review",
         ["Milestone retros", "  (author not", "  recorded)", "Final retro by a", "  fresh agent", "LESSONS.md pruned",
          "  to 6 rules"]),
        ("Guards in code",
         ["crossfit.py raises", "  LeakageError", "store.py checks", "  index hashes",
          "Test: only eval_test.py", "  reads sealed/", "FREEZE.json SHA pins", "Ledger refuses a",
          "  2nd TEST batch"]),
    ]
    w, gap = 0.184, 0.0125
    top = ypi - 0.07
    nmax = max(len(c[1]) for c in cols)
    h = bh(nmax) + 0.03
    for i, (title, lines) in enumerate(cols):
        x = 0.01 + i * (w + gap)
        p = FancyBboxPatch((x, top - h), w, h, boxstyle="round,pad=0.004,rounding_size=0.01",
                           fc="#f4f3ef", ec=AXIS, lw=0.8)
        ax.add_patch(p)
        ax.text(x + 0.01, top - 0.012, title, ha="left", va="top", fontsize=7.0, fontweight="bold",
                color=INK, linespacing=1.15)
        ax.text(x + 0.01, top - 0.070, "\n".join(lines), ha="left", va="top", fontsize=5.8,
                color=INK2, linespacing=1.3)
        arrow(ax, 0.5, ypi, x + w / 2, top)
    yb = top - h - 0.035
    tl = ["Model finished at T+6.0 h; no CV run after it improved M1 (0.786)",
          "Single sealed TEST spent at T+11.6 h of 24 (FREEZE 5.3 h before its gate); no second shard",
          "12 of 34 amendments (A23-A34) landed in the last 8.6 h around the freeze; A34 after TEST",
          "Priming-aware Cantley arm known to win on primed sites, but left out of FREEZE"]
    box(ax, 0.01, yb - bh(4), 0.98, bh(4), "What went wrong (from the outside-view retro)", tl,
        fc="white", bold_color=ORANGE)
    ax.text(0.5, yb - bh(4) - 0.02, "Opus 5.5 and Sonnet 5.5 were available; the model run by each agent is not recorded. "
            "Agents wrote only in their own paths and never ran git.",
            ha="center", va="top", fontsize=6.2, color=INK2)
    from matplotlib.transforms import Bbox
    save(fig, "Fig5_agents", bbox=Bbox([[0, H * (yb - bh(4) - 0.05)], [7.5, H]]))


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    print("figures written to", HERE.name)
