#!/usr/bin/env python3
"""Figure for the growth-machinery signature (sensitivity analysis). Run: python3 plot_signature_sensitivity.py <sigdir>"""
import sys, pathlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

O = pathlib.Path(sys.argv[1])
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.linewidth": 0.5, "svg.fonttype": "none"})
GREEN, PURPLE, GREY, ORANGE = "#5AAE61", "#6B3FA0", "#8C8C8C", "#D95F02"
UNITS = ["Ribosomal proteins", "Chaperones (GroESL / DnaK-DnaJ-GrpE)", "Sec / SRP translocon",
         "F-type ATP synthase", "Cell wall: peptidoglycan", "Cell wall: DAP branch", "Signature without DAP branch, core 4 units"]
SHORT = {"Ribosomal proteins": "Ribosomal proteins", "Chaperones (GroESL / DnaK-DnaJ-GrpE)": "Chaperones",
         "Sec / SRP translocon": "Sec / SRP", "F-type ATP synthase": "F-ATP synthase",
         "Cell wall: peptidoglycan": "Peptidoglycan", "Cell wall: DAP branch": "DAP branch",
         "Signature without DAP branch, core 4 units": "Core signature\n(first 4 units)"}
S1 = pd.read_csv(O / "S1_units_factorial.tsv", sep="\t")
S2 = pd.read_csv(O / "S2_family_consistency.tsv", sep="\t")
S3 = pd.read_csv(O / "S3_flask_index.tsv", sep="\t")
S4 = pd.read_csv(O / "S4_gene_contributions.tsv", sep="\t")

fig, axs = plt.subplots(1, 4, figsize=(183 / 25.4, 82 / 25.4), gridspec_kw=dict(width_ratios=[1.15, 1, 0.9, 1]))
y = np.arange(len(UNITS))[::-1]

ax = axs[0]
for v, dy, fc in (("baseline", 0.15, "k"), ("drop_Bacillaceae", -0.15, "white")):
    d = S1[(S1.variant == v) & (S1.component == "within") & (S1.contrast == "time 24->72 (both trt)")].set_index("unit").reindex(UNITS)
    ax.errorbar(d.rel_change, y + dy, xerr=[d.rel_change - d.lo, d.hi - d.rel_change], fmt="o", ms=3.5,
                mfc=fc, mec="k", ecolor="k", elinewidth=0.6, mew=0.6,
                label="all families" if v == "baseline" else "Bacillaceae removed")
ax.axvline(0, color="k", lw=0.4)
ax.set_yticks(y); ax.set_yticklabels([SHORT[u] for u in UNITS])
ax.set_xlabel("Within-family change, 24→72 h\n(relative, 95% CI, 18 flasks)")
ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.32), ncol=2)
ax.set_title("A  Community level", loc="left", fontweight="bold", fontsize=8)

ax = axs[1]
for i, u in enumerate(UNITS):
    d = S2[S2.unit == u]
    jit = np.random.default_rng(i).uniform(-0.22, 0.22, len(d))
    col = np.where(d.family == "Bacillaceae", ORANGE, np.where(d.q < .05, "k", GREY))
    ax.scatter(d.log2FC_72v24.clip(-3, 1.5), y[i] + jit, s=6, c=col, linewidths=0)
    ax.text(1.65, y[i], f"{(d.log2FC_72v24 < 0).sum()}/{len(d)}", va="center", fontsize=6)
ax.axvline(0, color="k", lw=0.4)
ax.set_yticks(y); ax.set_yticklabels([])
ax.set_xlim(-2.6, 2.3)
ax.set_xlabel("log2 (share of family transcriptome,\n72 h / 24 h), one point per family")
ax.set_title("B  Family level", loc="left", fontweight="bold", fontsize=8)
ax.text(-0.1, -0.36, "orange: Bacillaceae · black: q < 0.05 · n declining / n families", transform=ax.transAxes, fontsize=5.5)

ax = axs[2]
d = S3[(S3.variant == "baseline") & (S3.unit == "Signature without DAP branch, core 4 units")]
pos = {"C24": 0, "G24": 0.4, "C48": 1.2, "G48": 1.6, "C72": 2.4, "G72": 2.8}
for _, r in d.iterrows():
    ax.scatter(pos[r.group], r.within_index_rel, s=12, color=GREEN if r.group[0] == "C" else PURPLE, linewidths=0)
ax.set_xticks([0.2, 1.4, 2.6]); ax.set_xticklabels(["24", "48", "72"]); ax.set_xlim(-0.4, 3.2); ax.set_xlabel("Time (h)")
ax.set_ylabel("Core signature, within-family index\n(relative to the 6-group mean)")
ax.set_title("C  Flasks", loc="left", fontweight="bold", fontsize=8)
ax.scatter([], [], color=GREEN, s=12, label="control"); ax.scatter([], [], color=PURPLE, s=12, label="glucose")
ax.legend(frameon=False, fontsize=6, loc="center right")

ax = axs[3]
for i, u in enumerate(UNITS[:-1]):
    d = S4[(S4.variant == "baseline") & (S4.unit == u)]
    jit = np.random.default_rng(i + 10).uniform(-0.2, 0.2, len(d))
    ax.scatter(d.ko_rel_change.clip(-1, 1), y[i] + jit, s=4 + 120 * d.share_of_unit, facecolors="none",
               edgecolors=np.where(d.q < .05, "k", GREY), linewidths=0.5)
    n = len(d); ndown = int((d.ko_rel_change < 0).sum())
    ax.text(1.12, y[i], f"{ndown}/{n}", va="center", fontsize=6)
ax.axvline(0, color="k", lw=0.4)
ax.set_yticks(y); ax.set_yticklabels([])
ax.set_xlim(-1.08, 1.5)
ax.set_xlabel("Within-family change of each KO, 24→72 h\n(relative; circle area = share of unit)")
ax.set_title("D  Genes", loc="left", fontweight="bold", fontsize=8)

for a in axs:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout(w_pad=0.6)
for ext in ("png", "pdf", "svg"):
    fig.savefig(O / f"machinery_signature_sensitivity.{ext}", dpi=300, bbox_inches="tight")
print("saved")
