"""Figures for the research briefing (general scientific audience).

Reads only committed summary files:
  results/test/721dc1977366/*_ci.tsv, superiority.tsv, *_per_kinase.tsv,
  *_strata.tsv, final_reliability.tsv, results/test/LEDGER.md,
  results/cv/LEADERBOARD.tsv.
Event times in the timeline figure that are not in those files are copied
from RESEARCH_TREE.md, docs/retro/final_retro.md, docs/audit/*.md and
DESIGN.md; each is annotated with its source below.

Run from the repository root:  python -I briefing/figures/make_figures.py
"""
import os
import re
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))
TEST = os.path.join(ROOT, "results", "test", "721dc1977366")
CV = os.path.join(ROOT, "results", "cv", "LEADERBOARD.tsv")

# Reference palette (dataviz skill, light mode), fixed categorical order.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
MAGENTA, GREEN, VIOLET = "#e87ba4", "#008300", "#6250d6"
GRAY, LIGHTGRAY = "#8a8985", "#d9d8d4"
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#ffffff"
BOXFILL = "#f4f3f0"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 14,
    "axes.titlesize": 16, "axes.labelsize": 14,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=150, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print("wrote", os.path.relpath(path, ROOT))


def box(ax, x, y, w, h, title, body, edge=GRAY, fill=BOXFILL, tsize=15, bsize=13):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.015",
                                lw=2, ec=edge, fc=fill, transform=ax.transAxes))
    ax.text(x + 0.012, y + h - 0.03, title, transform=ax.transAxes, fontsize=tsize,
            fontweight="bold", va="top", ha="left", color=INK)
    ax.text(x + 0.012, y + h - 0.03 - 0.085, body, transform=ax.transAxes, fontsize=bsize,
            va="top", ha="left", color=INK2, linespacing=1.45)


def arrow(ax, x0, y0, x1, y1, color=GRAY):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=color, lw=2.2, mutation_scale=20))


def ci(cfg):
    return pd.read_csv(os.path.join(TEST, f"{cfg}_ci.tsv"), sep="\t").set_index("endpoint")


def strata(cfg):
    return pd.read_csv(os.path.join(TEST, f"{cfg}_strata.tsv"), sep="\t")


# ---------------------------------------------------------------- Fig 1
def fig1():
    fig = plt.figure(figsize=(15, 6.6))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    w, h, y = 0.138, 0.62, 0.30
    xs = [0.015 + i * (w + 0.0284) for i in range(6)]
    items = [
        ("1  Known pairs", "PhosphoSitePlus:\nkinase-site pairs\nfrom the literature\n\nhuman + mouse/rat\n(same 15-residue\nsequence)\n\nno true negatives", GRAY),
        ("2  A fair test", "group similar\nproteins together\n\nseal 20% as TEST\n\n'hard' negatives:\nsites of unrelated\nkinases", ORANGE),
        ("3  Evidence", "for every kinase-\nsite pair:\nsequence motif\nprotein language\n  model, structure\ninteractions\nlocalisation\nco-expression\nphosphoproteomics", GRAY),
        ("4  Two models", "LightGBM model\nof kinase-site\npairs (all\nevidence)\n\nneural model,\npre-trained on the\nCantley atlas", BLUE),
        ("5  Combine and\n    calibrate", "weighted average\nof the two\nstandardised\nscores\n\nturn the top score\ninto a probability\n('assignment\nprobability')", BLUE),
        ("6  Outputs", "top 200 sites\nper kinase\n\ntop 5 kinases\nper site\n\n503 kinases x\n350,736 sites", GRAY),
    ]
    for x, (t, b, e) in zip(xs, items):
        box(ax, x, y, w, h, t, b, edge=e, tsize=15, bsize=13)
    for i in range(5):
        arrow(ax, xs[i] + w + 0.001, y + h / 2, xs[i + 1] - 0.001, y + h / 2)
    # evaluation bar
    ax.add_patch(FancyBboxPatch((0.015, 0.04), 0.97, 0.19, boxstyle="round,pad=0.01,rounding_size=0.015",
                                lw=2, ec=ORANGE, fc="#fdf0ea", transform=ax.transAxes))
    ax.text(0.02, 0.185, "How it was judged", transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")
    ax.text(0.02, 0.125,
            "All choices made by 5-fold cross-validation on the development set (DEV)  ->  model frozen  ->  "
            "scored once on the sealed TEST set\n"
            "->  checked by an independent agent.  Compared with the Cantley motif atlas and our own retrained "
            "Phosformer and SELPHI2.0.",
            transform=ax.transAxes, fontsize=13, va="top", color=INK2, linespacing=1.5)
    save(fig, "fig1_pipeline.png")


# ---------------------------------------------------------------- Fig 2
def fig2():
    rng = np.random.default_rng(3)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.4))
    fams = ["cluster A", "cluster B", "cluster C", "cluster D", "cluster E"]
    for k, ax in enumerate(axs):
        ax.set_xlim(0, 10); ax.set_ylim(-0.6, 6.2); ax.axis("off")
        for fi, fam in enumerate(fams):
            yy = 5.2 - fi * 1.15
            ax.text(0.0, yy, fam, fontsize=13, va="center", color=INK2)
            for pi in range(2):  # two homologous proteins per family
                x0 = 2.0 + pi * 3.6
                ax.add_patch(FancyBboxPatch((x0, yy - 0.28), 3.0, 0.56, boxstyle="round,pad=0.02",
                                            fc=BOXFILL, ec=LIGHTGRAY, lw=1.5))
                for si in range(5):
                    sx = x0 + 0.35 + si * 0.58
                    if k == 0:
                        test = rng.random() < 0.3
                    else:
                        test = fi == 3
                    ax.scatter(sx, yy, s=170, marker="o" if not test else "s",
                               color=ORANGE if test else BLUE, edgecolor="white", lw=1.5, zorder=3)
        if k == 0:
            ax.set_title("Random split of kinase-site pairs: test sites\nsit on the same or similar proteins as training sites",
                         loc="left", fontsize=15)
        else:
            ax.set_title("Grouped split (this study): whole protein clusters\ngo to one side; TEST is sealed",
                         loc="left", fontsize=15)
            for fi in (3,):
                yy = 5.2 - fi * 1.15
                ax.add_patch(Rectangle((1.75, yy - 0.45), 7.25, 0.9, fill=False, ec=ORANGE, lw=2.2, ls="--"))
                ax.text(9.15, yy, "sealed\nTEST", fontsize=12, va="center", color=INK)
    axs[0].text(0.0, -0.45, "The model can score well by recognising the protein,\nnot by learning what the kinase prefers.",
                fontsize=13, color=INK2, va="center")
    axs[1].text(0.0, -0.45, "DEV clusters are split again into 5 folds for cross-validation;\nthe TEST clusters were opened once, at the end.",
                fontsize=13, color=INK2, va="center")
    h1 = axs[0].scatter([], [], s=170, color=BLUE, label="training site")
    h2 = axs[0].scatter([], [], s=170, marker="s", color=ORANGE, label="test site")
    fig.legend(handles=[h1, h2], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.04))
    save(fig, "fig2_splits.png")


# ---------------------------------------------------------------- Fig 3
def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(16, 7.0), gridspec_kw=dict(wspace=0.35))
    # M1: ranked list for one kinase; true substrates vs hard negatives
    ax = axs[0]
    order = [1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]  # 1 = true substrate, from top
    pos = [i for i, v in enumerate(order) if v]
    neg = [i for i, v in enumerate(order) if not v]
    auc = np.mean([[1.0 if p < n else 0.0 for n in neg] for p in pos])
    for i, v in enumerate(order):
        ax.add_patch(FancyBboxPatch((0.2, 11 - i - 0.38), 2.4, 0.76, boxstyle="round,pad=0.02",
                                    fc=BLUE if v else LIGHTGRAY, ec="white"))
        ax.text(2.8, 11 - i, "true substrate" if v else "hard negative", va="center", fontsize=12,
                color=INK if v else INK2)
        ax.text(0.0, 11 - i, f"{i + 1}", va="center", ha="right", fontsize=11, color=INK2)
    ax.set_xlim(-0.6, 7); ax.set_ylim(-1.2, 12.6); ax.axis("off")
    ax.set_title("M1: does kinase K rank its own\nsubstrates above other kinases'\nsubstrates?", loc="left", fontsize=15)
    ax.text(-0.5, -1.1, f"AUROC here = {auc:.2f}\n(0.5 = coin flip, 1 = perfect)", fontsize=13, va="center")
    # M2: for a site, ranked kinases
    ax = axs[1]
    ax.set_xlim(-0.6, 7.2); ax.set_ylim(-1.2, 12.6); ax.axis("off")
    ax.set_title("M2 (Hit@5): for a known site, is its\ntrue kinase among the model's\ntop 5 of the few hundred\ncandidate kinases?", loc="left", fontsize=15)
    for col, (true_rank, lab) in enumerate([(3, "site 1: hit"), (7, "site 2: miss")]):
        x0 = 0.0 + col * 3.5
        ax.text(x0, 11.2, lab, fontsize=13, fontweight="bold")
        for r in range(1, 10):
            yy = 10.4 - (r - 1) * 1.05
            ax.add_patch(FancyBboxPatch((x0, yy - 0.38), 2.9, 0.76, boxstyle="round,pad=0.02",
                                        fc=ORANGE if r == true_rank else LIGHTGRAY, ec="white"))
            ax.text(x0 + 0.15, yy, f"{r}. {'true kinase' if r == true_rank else 'kinase'}", va="center",
                    fontsize=11.5, color=INK)
        ax.add_patch(Rectangle((x0 - 0.12, 10.4 - 4 * 1.05 - 0.5), 3.14, 5 * 1.05 + 0.05, fill=False,
                               ec=INK, lw=1.8, ls="--"))
    ax.text(-0.5, -0.9, "Hit@5 = share of a kinase's sites that\nare hits, then averaged over kinases",
            fontsize=13, va="center")
    # M3: top 1% of all sites
    ax = axs[2]
    ax.set_xlim(0, 10); ax.set_ylim(-1.2, 12.6); ax.axis("off")
    ax.set_title("M3: how many of K's true substrates\nare in its top 1% of all ~60,000\ncandidate sites?", loc="left", fontsize=15)
    ax.add_patch(Rectangle((0.5, 1.0), 1.4, 9.8, fc=BOXFILL, ec=LIGHTGRAY, lw=1.5))
    ax.add_patch(Rectangle((0.5, 10.8 - 0.098 * 9.8 - 0.6), 1.4, 0.098 * 9.8 + 0.6, fc="#cde2fb", ec=BLUE, lw=1.8))
    for yy in [10.55, 10.25, 9.95]:
        ax.plot([0.5, 1.9], [yy, yy], color=BLUE, lw=3)
    for yy in [6.2, 3.1]:
        ax.plot([0.5, 1.9], [yy, yy], color=BLUE, lw=3)
    ax.text(2.2, 10.2, "top 1% (about 600 sites)\nholds 3 of the 5 true\nsubstrates: recall = 0.6", fontsize=12.5, va="center")
    ax.text(2.2, 4.6, "the other 2 true substrates\nare ranked lower", fontsize=12.5, va="center", color=INK2)
    ax.text(0.5, 0.0, "all candidate sites, best first\n(not to scale)", fontsize=11.5, color=INK2)
    ax.text(0.2, -0.9, "Close to what an experimentalist\nwould test: a short list per kinase", fontsize=13, va="center")
    save(fig, "fig3_metrics.png")


# ---------------------------------------------------------------- Fig 4
def fig4():
    lb = pd.read_csv(CV, sep="\t").set_index("run_id")
    cvrow = lambda rid, col: float(lb.loc[rid, col])
    methods = [  # label, TEST cfg per endpoint, CV run per endpoint (None = not comparable)
        ("Our model (final)", {"M1": "final", "M2": "final", "M3": "final"},
         {"M1": ("r2_fusion_r3__20261026-191901__s0", "M1"), "M2": ("r2_fusion_r3__20261026-191901__s0", "M2_hit5"),
          "M3": ("r2_fusion_r3__20261026-191901__s0", "M3_top1")}, BLUE),
        ("Cantley motif atlas", {"M1": "B1_i", "M2": "B1_ii", "M3": "B1_i"},
         {"M1": ("simple_b1_published_reeval__20261026-155352__s0", "M1"),
          "M2": ("simple_b1_universe_reeval__20261026-155703__s0", "M2_hit5"), "M3": None}, ORANGE),
        ("SELPHI2.0 (our retrain)", {"M1": "B4", "M2": "B4", "M3": "B4"}, {"M1": None, "M2": None, "M3": None}, AQUA),
        ("Phosformer (our retrain)", {"M1": "B3", "M2": "B3", "M3": "B3"},
         {"M1": ("b3a_phosformer_cv_pooled__20261026-222734__s0", "M1"),
          "M2": ("b3a_phosformer_cv_pooled__20261026-222734__s0", "M2_hit5"),
          "M3": ("b3a_phosformer_cv_pooled__20261026-222734__s0", "M3_top1")}, YELLOW),
        ("Popularity floor", {"M1": "B0_pop", "M2": "B0_pop", "M3": "B0_pop"},
         {"M1": ("simple_b0_pop_reeval__20261026-155346__s0", "M1"),
          "M2": ("simple_b0_pop_reeval__20261026-155346__s0", "M2_hit5"),
          "M3": ("simple_b0_pop_reeval__20261026-155346__s0", "M3_top1")}, GRAY),
    ]
    titles = {"M1": "M1: own substrates above\nother kinases' substrates (AUROC)",
              "M2": "M2: true kinase in top 5\nfor a site (Hit@5)",
              "M3": "M3: true substrates in a\nkinase's top 1% (recall)"}
    xlims = {"M1": (0.5, 0.85), "M2": (0.0, 0.42), "M3": (0.0, 0.42)}
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True, gridspec_kw=dict(wspace=0.12))
    ys = np.arange(len(methods))[::-1]
    for ax, ep in zip(axs, ["M1", "M2", "M3"]):
        for y, (lab, tc, cvr, col) in zip(ys, methods):
            c = ci(tc[ep]).loc[ep]
            ax.plot([c.lo95, c.hi95], [y, y], color=col, lw=3, solid_capstyle="round", zorder=2)
            ax.scatter(c.point, y, s=110, color=col, edgecolor="white", lw=1.5, zorder=3)
            if cvr[ep] is not None:
                ax.scatter(cvrow(*cvr[ep]), y + 0.28, s=80, facecolor="white", edgecolor=col, lw=2, zorder=3)
        if ep == "M1":
            ax.axvline(0.5, color=LIGHTGRAY, lw=1.5, ls=":")
        ax.set_xlim(*xlims[ep]); ax.set_title(titles[ep], loc="left", fontsize=14.5)
        ax.grid(axis="x", color="#ecebe8", lw=1); ax.set_axisbelow(True)
    axs[0].set_yticks(ys); axs[0].set_yticklabels([m[0] for m in methods])
    axs[0].set_ylim(-0.6, len(methods) - 0.3)
    h1 = axs[0].scatter([], [], s=110, color=INK2, label="TEST (sealed), with 95% interval")
    h2 = axs[0].scatter([], [], s=80, facecolor="white", edgecolor=INK2, lw=2, label="cross-validation (DEV)")
    fig.legend(handles=[h1, h2], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.55, -0.07))
    save(fig, "fig4_results.png")


# ---------------------------------------------------------------- Fig 5
def fig5():
    fin = pd.read_csv(os.path.join(TEST, "final_per_kinase.tsv"), sep="\t")
    pop = pd.read_csv(os.path.join(TEST, "B0_pop_per_kinase.tsv"), sep="\t")[["kinase_acc", "auroc_hard"]]
    d = fin.merge(pop, on="kinase_acc", suffixes=("", "_pop"))
    d = d[d.in_K_primary.astype(str) == "True"].copy()
    d["beats"] = d.auroc_hard > d.auroc_hard_pop
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.2), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.28))
    ax = axs[0]
    for cls, col, mk, lab in [("ST", INK2, "o", "Ser/Thr kinases"), ("Y", VIOLET, "^", "Tyr kinases")]:
        s = d[d.residue_class == cls]
        for beats, fc in [(True, col), (False, "white")]:
            t = s[s.beats == beats]
            ax.scatter(t.n_train_dev_L2.clip(lower=1), t.auroc_hard, s=75, marker=mk, facecolor=fc,
                       edgecolor=col, lw=1.8, zorder=3, label=None)
        nb = int(s.beats.sum())
        print(f"fig5: {cls}: {nb} of {len(s)} above own popularity floor")
    s = d[d.residue_class == "STY"]
    ax.scatter(s.n_train_dev_L2.clip(lower=1), s.auroc_hard, s=75, marker="D", facecolor=GRAY, edgecolor=GRAY)
    ax.axhline(0.5, color=LIGHTGRAY, ls=":", lw=1.5)
    ax.set_xscale("log"); ax.set_ylim(0.2, 1.02)
    ax.set_xticks([1, 10, 100, 1000]); ax.set_xticklabels(["≤1", "10", "100", "1000"])
    ax.minorticks_off()
    ax.set_xlabel("known substrates of the kinase in training (DEV labels, log scale)")
    ax.set_ylabel("TEST AUROC (M1) for that kinase")
    ax.set_title("A  Each dot is one kinase (119 on TEST)", loc="left")
    hs = [ax.scatter([], [], s=75, marker="o", color=INK2, label="Ser/Thr kinase"),
          ax.scatter([], [], s=75, marker="^", color=VIOLET, label="Tyr kinase"),
          ax.scatter([], [], s=75, marker="o", facecolor="white", edgecolor=INK2, lw=1.8,
                     label="hollow: not above its own\npopularity floor"),
          ax.scatter([], [], s=75, marker="D", color=GRAY, label="dual-specificity kinase")]
    ax.legend(handles=hs, loc="lower right", frameon=False, fontsize=12)
    ax.grid(axis="y", color="#ecebe8"); ax.set_axisbelow(True)
    # B: strata bars
    ax = axs[1]
    cfgs = [("final", "Our model", BLUE), ("B1_i", "Cantley atlas", ORANGE), ("B4", "SELPHI2.0 (retrain)", AQUA),
            ("B3", "Phosformer (retrain)", YELLOW), ("B0_pop", "Popularity floor", GRAY)]
    groups = [("ST", "Ser/Thr kinases (94)"), ("Y", "Tyr kinases (23)")]
    bw = 0.15
    for j, (cfg, lab, col) in enumerate(cfgs):
        st = strata(cfg)
        st = st[st.stratum == "residue_class"].set_index("value")
        vals = [st.loc[g, "M1"] for g, _ in groups]
        xs = np.arange(2) + (j - 2) * (bw + 0.02)
        ax.bar(xs, [v - 0.4 for v in vals], bottom=0.4, width=bw, color=col, label=lab, edgecolor="white", lw=2)
        for x, v in zip(xs, vals):
            ax.text(x, v + 0.006, f"{v:.3f}", ha="center", va="bottom", fontsize=10.5, color=INK2, rotation=90)
    ax.axhline(0.5, color=INK2, ls=":", lw=1.5)
    ax.text(1.55, 0.505, "chance", fontsize=11, color=INK2, va="bottom", ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels([g[1] for g in groups])
    ax.set_ylim(0.4, 0.95); ax.set_ylabel("TEST M1 (macro AUROC)")
    ax.set_title("B  Ser/Thr versus Tyr kinases", loc="left")
    ax.legend(frameon=False, fontsize=11.5, loc="upper right")
    save(fig, "fig5_per_kinase.png")


# ---------------------------------------------------------------- Fig 6
def fig6():
    lb = pd.read_csv(CV, sep="\t").set_index("run_id")
    runs = [
        ("One global model (default)", "regime_L2_L2fixed_s0_cv5__20261026-203141__s0", BLUE),
        ("Partial pooling (selected setting)", "regime_ppool_lam300_cv5__20261026-205817__s0", BLUE),
        ("One model per kinase group", "regime_pergroup_raw_cv5__20261026-220519__s0", ORANGE),
        ("One model per kinase (LightGBM)", "regime_perkinase_lgb_N20_qmap_cv5__20261026-203214__s0", ORANGE),
        ("One model per kinase (logistic)", "regime_perkinase_lr_N20_qmap_cv5__20261026-231144__s0", ORANGE),
        ("Labels: human only", "regime_L1_refit_s0_cv5__20261026-203555__s0", AQUA),
        ("Labels: in vitro evidence only", "regime_L2_vitro_refit_s0_cv5__20261026-210829__s0", AQUA),
        ("Labels: in vivo evidence only", "regime_L2_vivo_refit_s0_cv5__20261026-211636__s0", AQUA),
    ]
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.8), sharey=True, gridspec_kw=dict(wspace=0.08))
    ys = np.arange(len(runs))[::-1]
    base = runs[0][1]
    for ax, col, title, lim in [(axs[0], "M1", "M1 (AUROC against hard negatives)", (0.73, 0.78)),
                                (axs[1], "M2_hit5", "M2 (true kinase in top 5)", (0.09, 0.21))]:
        ref = float(lb.loc[base, col])
        ax.axvline(ref, color=BLUE, lw=1.5, ls="--", alpha=0.7)
        for y, (lab, rid, c) in zip(ys, runs):
            v = float(lb.loc[rid, col])
            ax.plot([ref, v], [y, y], color=LIGHTGRAY, lw=2, zorder=1)
            ax.scatter(v, y, s=120, color=c, edgecolor="white", lw=1.5, zorder=3)
            ax.text(v, y + 0.3, f"{v:.3f}", ha="center", fontsize=11, color=INK2)
        ax.set_xlim(*lim); ax.set_title(title, loc="left")
        ax.grid(axis="x", color="#ecebe8"); ax.set_axisbelow(True)
    axs[0].set_yticks(ys); axs[0].set_yticklabels([r[0] for r in runs])
    axs[0].set_ylim(-0.6, len(runs) - 0.2)
    hs = [axs[0].scatter([], [], s=120, color=BLUE, label="shared model"),
          axs[0].scatter([], [], s=120, color=ORANGE, label="split model"),
          axs[0].scatter([], [], s=120, color=AQUA, label="label variant (global model)")]
    fig.legend(handles=hs, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.55, -0.07))
    save(fig, "fig6_regimes.png")


# ---------------------------------------------------------------- Fig 7
T0 = datetime(2026, 10, 26, 13, 20)


def th(day, hh, mm):
    """Hours since T0 for a CET clock time on 26 or 27 October."""
    return (datetime(2026, 10, day, hh, mm) - T0).total_seconds() / 3600


def fig7():
    lb = pd.read_csv(CV, sep="\t")
    st = pd.to_datetime(lb.started.str.replace(" CET", ""), errors="coerce").dropna()
    cv_hours = ((st - T0).dt.total_seconds() / 3600).values
    led = open(os.path.join(ROOT, "results", "test", "LEDGER.md")).read()
    t_start = re.search(r"STARTED TEST batch (\S+ \S+)", led).group(1)
    t_done = re.search(r"COMPLETED TEST batch (\S+ \S+)", led).group(1)
    ts = (datetime.fromisoformat(t_start) - T0).total_seconds() / 3600
    td = (datetime.fromisoformat(t_done) - T0).total_seconds() / 3600

    fig, ax = plt.subplots(figsize=(16, 7.6))
    lanes = {"Plan and review": 4, "Build and model (CV runs)": 3, "Independent checks": 2,
             "Sealed TEST": 1, "Problems": 0}
    for name, y in lanes.items():
        ax.axhline(y, color="#ecebe8", lw=1, zorder=0)
    # Plan lane (RESEARCH_TREE.md M1; final_retro.md §1.1)
    ax.barh(4, th(26, 14, 33) - 0, left=0, height=0.32, color=BLUE, alpha=0.85)
    ax.text(0.1, 4.28, "literature, design, 11 reviews (to 14:33)", fontsize=11.5, va="bottom")
    ax.barh(4, th(27, 3, 17) - td, left=td, height=0.32, color=BLUE, alpha=0.5)
    ax.text(td, 4.28, "report drafts", fontsize=11.5, va="bottom")
    ax.barh(4, th(27, 9, 47) - th(27, 9, 27), left=th(27, 9, 27), height=0.32, color=BLUE, alpha=0.85)
    ax.text(th(27, 9, 47), 4.28, "final retro + report v3 (09:47)", fontsize=11.5, va="bottom", ha="right")
    # Build lane: every timed CV run (LEADERBOARD.tsv 'started')
    ax.scatter(cv_hours, 3 + np.random.default_rng(0).uniform(-0.13, 0.13, len(cv_hours)), s=14, color=BLUE, alpha=0.6)
    fm = th(26, 19, 19)
    ax.scatter([fm], [3], s=320, marker="*", color=INK, zorder=4)
    ax.text(fm, 3.28, "final model run, 19:19 (T+6.0 h);\nno later run beat it", fontsize=11.5, ha="center", va="bottom")
    # Checks lane (docs/audit/*.md headers; RESEARCH_TREE.md 5.12)
    checks = [(th(26, 16, 45), th(26, 17, 30), "CV audit 1:\nmetric defect A24"),
              (th(26, 20, 25), th(26, 20, 45), "test-code review 20:25:\nguard bypassable;\nfixes applied 23:30"),
              (th(26, 21, 24), th(26, 22, 15), "CV audit 2"),
              (th(27, 1, 0), th(27, 2, 0), "TEST audit"),
              (th(27, 2, 34), th(27, 3, 16), "post-TEST\naudit")]
    for i, (a, b, lab) in enumerate(checks):
        ax.barh(2, b - a, left=a, height=0.32, color=AQUA)
        ax.text((a + b) / 2, 2.22 if i % 2 == 0 else 1.62, lab, fontsize=10.5, ha="center",
                va="bottom" if i % 2 == 0 else "top")
    # fixes from the test-code review applied at 23:30 (RESEARCH_TREE.md 5.12; final_retro.md §3.3)
    ax.plot([th(26, 20, 45), th(26, 23, 30)], [2, 2], color=AQUA, lw=1.5, ls=":")
    ax.scatter([th(26, 23, 30)], [2], marker="|", s=400, color=AQUA, lw=2.5)
    # TEST lane (FREEZE 00:04: final_retro.md; batch: LEDGER.md; gate 05:20: final_retro.md §1.2)
    ax.barh(1, td - ts, left=ts, height=0.36, color=ORANGE)
    ax.scatter([th(27, 0, 4)], [1], marker="|", s=500, color=INK, lw=2.5)
    ax.text(td + 0.2, 1.0, "FREEZE 00:04\nTEST 00:10-00:57", fontsize=11.5, ha="left", va="center")
    g4 = th(27, 5, 20)
    ax.axvline(g4, color=ORANGE, ls="--", lw=1.8, ymin=0.12, ymax=0.32)
    ax.text(g4 + 0.12, 1.0, "planned TEST gate\n05:20 (T+16 h)", fontsize=11.5, va="center")
    # Problems lane (final_retro.md; supplement §6)
    probs = [(th(26, 16, 55), "6 protocol changes\nin 40 min"),
             (th(26, 17, 50), "out of\nmemory"),
             (th(26, 19, 55), "required regime\ntest started late"),
             (th(26, 21, 11), "stronger priming\narm not frozen"),
             (th(26, 22, 12), "M2 score\nchanged late"),
             (th(26, 23, 0), "out of\nmemory")]
    for i, (x, lab) in enumerate(probs):
        ax.scatter([x], [0], s=90, marker="D", color=ORANGE, edgecolor="white", zorder=3)
        ax.text(x, -0.22 if i % 2 == 0 else 0.2, lab, fontsize=10.2, ha="center", va="top" if i % 2 == 0 else "bottom")
    ax.barh(0, th(27, 9, 32) - th(27, 3, 17), left=th(27, 3, 17), height=0.22, color=LIGHTGRAY)
    ax.text(th(27, 5, 40), 0.18, "about 6 h with no file written (03:17-09:32)", fontsize=11, ha="left", va="bottom")
    ax.set_yticks(list(lanes.values())); ax.set_yticklabels(list(lanes.keys()))
    ticks = list(range(0, 25, 2))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"T+{t}h\n{(13 + t) % 24:02d}:20" for t in ticks], fontsize=11.5)
    ax.set_xlim(-0.3, 24.2); ax.set_ylim(-0.95, 4.9)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("hours since start (clock time, 26-27 October 2026)")
    save(fig, "fig7_timeline.png")
    print(f"fig7: {len(cv_hours)} CV runs; {int((cv_hours >= th(26, 20, 0)).sum())} started at or after 20:00; "
          f"TEST batch T+{ts:.2f}h to T+{td:.2f}h")


# ---------------------------------------------------------------- Fig 8
def fig8():
    r = pd.read_csv(os.path.join(TEST, "final_reliability.tsv"), sep="\t")
    r = r[r.subset == "labelled_TEST_sites"]
    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.plot([0, 1], [0, 1], color=LIGHTGRAY, lw=2, ls="--", zorder=1)
    ax.text(0.86, 0.9, "perfect\ncalibration", fontsize=12, color=INK2, ha="center")
    ax.add_patch(Rectangle((0.8, 0), 0.2, 1, color="#f4f3f0", zorder=0))
    ax.text(0.9, 0.12, "no site\nreaches\n0.8", fontsize=12.5, ha="center", color=INK2)
    ax.plot(r.mean_p, r.accuracy, color=BLUE, lw=2, zorder=2)
    ax.scatter(r.mean_p, r.accuracy, s=r.n_sites * 1.1, color=BLUE, edgecolor="white", lw=1.5, zorder=3)
    for _, row in r.iterrows():
        ax.text(row.mean_p, row.accuracy - 0.06, f"n={row.n_sites}", fontsize=10.5, ha="center", color=INK2)
    hi = r[r.mean_p >= 0.5]
    print(f"fig8: n={int(r.n_sites.sum())}; >=0.5 tier n={int(hi.n_sites.sum())}, claimed "
          f"{(hi.mean_p * hi.n_sites).sum() / hi.n_sites.sum():.4f}, observed "
          f"{(hi.accuracy * hi.n_sites).sum() / hi.n_sites.sum():.4f}; ECE {r.ECE.iloc[0]:.4f}")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("assignment probability claimed by the model")
    ax.set_ylabel("share of sites where the top kinase\nwas a known kinase")
    ax.set_title("Labelled TEST sites (n = 1,458)", loc="left")
    ax.grid(color="#ecebe8"); ax.set_axisbelow(True)
    save(fig, "fig8_calibration.png")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); fig7(); fig8()
