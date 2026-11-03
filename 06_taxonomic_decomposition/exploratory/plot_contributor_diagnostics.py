#!/usr/bin/env python3
"""Exploratory — diagnostic plot for the contributor analyses.
A  Q1/Q2: PERMANOVA R2 of the family composition of each function's transcripts
   (raw, and adjusted for the family composition of the whole transcriptome)
B  Decomposition of the change in each function's share of the bacterial transcriptome
C  Effective number of contributing families relative to the whole transcriptome
Run: python3 plot_contributor_diagnostics.py <outdir>
"""
import sys, pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

O = pathlib.Path(sys.argv[1])
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
GREEN, PURPLE, BLUE, ORANGE, GREY = "#5AAE61", "#6B3FA0", "#1F5FA9", "#D95F02", "#8C8C8C"
ORDER = ["Glycolysis / Gluconeogenesis", "Pentose phosphate pathway", "Citrate cycle (TCA cycle)",
         "Pyruvate metabolism", "Ribosome", "Valine, leucine and isoleucine biosynthesis",
         "Arginine biosynthesis", "Peptidoglycan biosynthesis", "Nitrogen metabolism",
         "Ntr regulon (N limitation)"]
SHORT = {"Valine, leucine and isoleucine biosynthesis": "Val/Leu/Ile biosynthesis",
         "Val/Leu/Ile biosynthesis": "Val/Leu/Ile biosynthesis",
         "Ntr regulon (N limitation)": "Ntr regulon", "Ntr regulon (module)": "Ntr regulon",
         "Glycolysis / Gluconeogenesis": "Glycolysis", "Citrate cycle (TCA cycle)": "TCA cycle",
         "Pentose phosphate pathway": "Pentose phosphate", "Ribosome": "Ribosomal proteins"}
lab = lambda f: SHORT.get(f, f)

q = pd.read_csv(O / "Q1Q2_permanova_by_function.tsv", sep="\t")
q = q[(q.domain == "Bacteria") & (q.set_class != "screen")].copy()
q["lab"] = q.function.map(lab)
# panel B: factorial-model contrasts on all samples (decomposition_contrasts.py, A2), not the replicate bootstrap
V = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else O
Gf = pd.read_csv(V / "A2_factorial_contrasts_all_variants.tsv", sep="\t")
Gf = Gf[(Gf.variant == "baseline") & Gf.function.isin(ORDER)]
A = Gf.pivot_table(index=["function", "contrast"], columns="component", values=["est", "lo", "hi", "q"])
A.columns = [f"{c}" if m == "est" else f"{c}_{m}" for m, c in A.columns]
A = A.reset_index().rename(columns={"total": "rel_change_total"})
A["lab"] = A.function.map(lab)
D = pd.read_csv(O / "B_Deff_by_sample.tsv", sep="\t")
D = D[(D.domain == "Bacteria") & (D.set_class != "screen")].copy()
D["lab"] = D.function.map(lab)
rows = [lab(f) for f in ORDER]

fig = plt.figure(figsize=(183 / 25.4, 150 / 25.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.9], height_ratios=[1, 0.78], wspace=0.55, hspace=0.62)

# ---------------- A
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(rows))[::-1]
terms = [("time", "Time"), ("treatment", "Glucose"), ("interaction", "Glucose × time")]
off = {"time": 0.22, "treatment": 0, "interaction": -0.22}
col = {"time": "#333333", "treatment": PURPLE, "interaction": "#A6761D"}
for t, tl in terms:
    d = q.set_index("lab").reindex(rows)
    for kind, mk in (("raw", "o"), ("adj", "D")):
        r2, qq = d[f"{kind}_R2_{t}"].to_numpy(), d[f"{kind}_q_{t}"].to_numpy()
        ax.scatter(r2, y + off[t], marker=mk, s=16 if kind == "raw" else 11,
                   facecolors=np.where(qq < 0.05, col[t], "white"), edgecolors=col[t], linewidths=0.6, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(rows)
ax.set_xlim(0, 0.75); ax.set_xlabel("PERMANOVA R² (family composition\nof the function's transcripts)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
h = [Line2D([], [], ls="", marker="o", mfc=col[t], mec=col[t], ms=4, label=tl) for t, tl in terms]
h += [Line2D([], [], ls="", marker="o", mfc="none", mec="k", ms=4, label="raw"),
      Line2D([], [], ls="", marker="D", mfc="none", mec="k", ms=3.5, label="vs. whole transcriptome"),
      Line2D([], [], ls="", marker="o", mfc="k", mec="k", ms=4, label="filled: BH q < 0.05")]
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, frameon=False, fontsize=6, handletextpad=0.2, columnspacing=0.8)
ax.set_title("A  Who transcribes each function?", loc="left", fontweight="bold", fontsize=8)

# ---------------- B
CONTR = [("time 24->48 (both trt)", None, "24→48 h"),
         ("time 24->72 (both trt)", None, "24→72 h"),
         ("glucose vs control 24 h", None, "Glucose − control, 24 h")]
sub = gs[0, 1].subgridspec(1, 3, wspace=0.08)
for i, (a, b, title) in enumerate(CONTR):
    ax = fig.add_subplot(sub[0, i])
    d = A[A.contrast == a].set_index("lab").reindex(rows)
    for comp, c, dy in (("turnover", BLUE, 0.17), ("within", ORANGE, -0.17)):
        v, lo, hi = d[comp].to_numpy(), d[f"{comp}_lo"].to_numpy(), d[f"{comp}_hi"].to_numpy()
        ax.barh(y + dy, v, height=0.32, color=c, alpha=0.9)
        ax.errorbar(v, y + dy, xerr=[v - lo, hi - v], fmt="none", ecolor="k", elinewidth=0.5, capsize=0)
        qq = d[f"{comp}_q"].to_numpy()
        for yy, vv, hh, q_ in zip(y + dy, v, hi, qq):
            if q_ < 0.05:
                ax.text(max(hh, vv) + 0.03 if vv >= 0 else min(lo[list(y + dy).index(yy)], vv) - 0.03, yy, "*",
                        ha="left" if vv >= 0 else "right", va="center", fontsize=7)
    ax.scatter(d.rel_change_total, y, marker="|", s=40, color="k", zorder=4, linewidths=1.2)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlim(-0.75, 0.75); ax.set_xticks([-0.5, 0, 0.5])
    ax.set_yticks(y); ax.set_yticklabels(rows if i == 0 else [])
    ax.tick_params(axis="y", length=0 if i else 2)
    ax.set_title(title, fontsize=7)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if i == 1:
        ax.set_xlabel("Relative change in the function's share of bacterial transcripts \n(time contrasts pool both treatments; 95% CI from 18 flasks; * BH q < 0.05)")
    if i == 0:
        ax.text(-1.25, 1.10, "B  Turnover of contributors vs. change within families",
                transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom")
h = [plt.Rectangle((0, 0), 1, 1, color=BLUE, label="contributor turnover"),
     plt.Rectangle((0, 0), 1, 1, color=ORANGE, label="within-family change"),
     Line2D([], [], ls="", marker="|", color="k", ms=6, mew=1.2, label="total")]
fig.axes[-1].legend(handles=h, loc="upper center", bbox_to_anchor=(-0.55, -0.2), ncol=3, frameon=False, fontsize=6)

# ---------------- C
sub = gs[1, :].subgridspec(1, 4, wspace=0.35)
SHOW = ["Peptidoglycan biosynthesis", "Ribosomal proteins", "Glycolysis", "Nitrogen metabolism"]
for i, f in enumerate(SHOW):
    ax = fig.add_subplot(sub[0, i])
    d = D[D.lab == f]
    for tr, c in (("C", GREEN), ("G", PURPLE)):
        m = d[d.treatment == tr].groupby("time").D_eff_ratio.agg(["mean", "sem"]).reindex(["24h", "48h", "72h"])
        ax.errorbar([24, 48, 72], m["mean"], yerr=m["sem"], color=c, marker="o", ms=3, lw=1, capsize=1.5,
                    label="control" if tr == "C" else "glucose")
    ax.axhline(1, color=GREY, lw=0.5, ls="--")
    ax.set_xticks([24, 48, 72]); ax.set_xlabel("Time (h)")
    ax.set_ylim(0.3, 1.15)
    if i == 0:
        ax.set_ylabel("D$_{eff}$ of contributors /\nD$_{eff}$ of whole transcriptome")
        ax.legend(frameon=False, fontsize=6, loc="lower right")
        ax.set_title("C  How many families share the function?", loc="left", fontweight="bold",
                     fontsize=8, pad=14, x=-0.25)
    ax.text(0.5, 1.02, f, transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

for ext in ("png", "pdf", "svg"):
    fig.savefig(O / f"contributor_diagnostics.{ext}", dpi=300, bbox_inches="tight")
print("saved")
