#!/usr/bin/env python3
"""Is masking concentrated in growth-related functions?

Uses the factorial-model contrasts (A2_factorial_contrasts_all_variants.tsv, 18 flasks,
BH within contrast x component). A map is 'masked' when both components are significant
(q < 0.05) with opposite signs. The a priori grouping is KEGG level 1/2 (BRITE ko00001);
the 'growth machinery' reading of the within-family decline is a post hoc annotation and
is reported as such, together with the gene content of each map so that redundant maps
(Protein export = Bacterial secretion system = Sec/SRP; Lysine = Monobactam = DAP
pathway) are not counted as independent evidence.

Run: python3 masking_by_category.py <valdir> <KO1-4.txt> <counts_ko_bacteria.tsv>
"""
import sys, pathlib
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V, KO14, CNT = pathlib.Path(sys.argv[1]), sys.argv[2], sys.argv[3]
h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
cat = h.dropna(subset=["pid"]).drop_duplicates("pid").set_index("pid")[["PathwayL1", "PathwayL2"]]
G = pd.read_csv(V / "A2_factorial_contrasts_all_variants.tsv", sep="\t")
cb = pd.read_csv(CNT, sep="\t", index_col=0).sum(axis=1)

# what each map actually is: its three most transcribed KOs and their share of the map
top = {}
for pid, d in h.dropna(subset=["pid"]).drop_duplicates(["pid", "KO"]).groupby("pid"):
    s = cb.reindex(d.KO).dropna().sort_values(ascending=False)
    if s.sum() > 0:
        desc = d.set_index("KO").KoDescription
        top[pid] = "; ".join(f"{desc[k].split(';')[0].split(',')[-1].strip()} {s[k] / s.sum():.0%}" for k in s.index[:3])


def table(variant, contrast):
    b = G[(G.variant == variant) & (G.set_class != "module") & (G.contrast == contrast) & (G.component != "total")]
    w = b.pivot_table(index=["pathway", "function"], columns="component", values=["est", "q"]).reset_index()
    w.columns = ["pathway", "function", "turnover", "within", "q_turnover", "q_within"]
    w = w.join(cat, on="pathway")
    w["masked"] = (w.q_turnover < .05) & (w.q_within < .05) & (np.sign(w.turnover) != np.sign(w.within))
    w["within_down"] = (w.q_within < .05) & (w.within < 0)
    w["within_up"] = (w.q_within < .05) & (w.within > 0)
    w["top_KOs"] = w.pathway.map(top)
    w.insert(0, "contrast", contrast); w.insert(0, "variant", variant)
    return w


rows, tests = [], []
for v in ("baseline", "drop_Bacillaceae"):
    for c in ("time 24->48 (both trt)", "time 24->72 (both trt)"):
        w = table(v, c)
        rows.append(w)
        gip = w.PathwayL1 == "Genetic Information Processing"
        for flag in ("masked", "within_down", "within_up"):
            a, b_ = int(w[gip][flag].sum()), int(gip.sum())
            c_, d_ = int(w[~gip][flag].sum()), int((~gip).sum())
            p = stats.fisher_exact([[a, b_ - a], [c_, d_ - c_]])[1]
            tests.append(dict(variant=v, contrast=c, outcome=flag, GIP=f"{a}/{b_}", other=f"{c_}/{d_}", fisher_p=p))
        kw = stats.kruskal(*[g.within.values for _, g in w.groupby("PathwayL1") if len(g) >= 4])
        rho, pr = stats.spearmanr(w.turnover, w.within)
        tests.append(dict(variant=v, contrast=c, outcome="within ~ KEGG level 1 (Kruskal-Wallis)", fisher_p=kw.pvalue))
        tests.append(dict(variant=v, contrast=c, outcome=f"Spearman turnover vs within (rho = {rho:.2f})", fisher_p=pr))
W = pd.concat(rows, ignore_index=True)
W.to_csv(V / "M1_masking_by_map.tsv", sep="\t", index=False)
T = pd.DataFrame(tests)
T.to_csv(V / "M2_masking_category_tests.tsv", sep="\t", index=False)
print(T.round(4).to_string())
print(W[(W.masked | W.within_down) & (W.contrast == "time 24->72 (both trt)")]
      [["variant", "function", "PathwayL2", "turnover", "within", "masked", "top_KOs"]].round(2).to_string())

# ---------------- figure
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.linewidth": 0.5, "svg.fonttype": "none"})
COL = {"Genetic Information Processing": "#6B3FA0", "Metabolism": "#5AAE61",
       "Environmental Information Processing": "#1F5FA9", "Cellular Processes": "#D95F02"}
LAB = {"Ribosome": "Ribosomal proteins", "RNA degradation": "RNA degradation\n(GroEL/DnaK 68 %)",
       "Protein export": "Sec translocon", "Bacterial secretion system": None,
       "Photosynthesis": "F-type ATP synthase\n('Photosynthesis')", "Peptidoglycan biosynthesis": "Peptidoglycan",
       "Lysine biosynthesis": "DAP / lysine", "Monobactam biosynthesis": None,
       "Carbon fixation in photosynthetic organisms": "C fixation (Calvin/glycolytic)",
       "Cationic antimicrobial peptide (CAMP) resistance": "CAMP resistance (dlt)"}
fig, axs = plt.subplots(1, 2, figsize=(183 / 25.4, 85 / 25.4), sharey=True)
for ax, v, title in zip(axs, ("baseline", "drop_Bacillaceae"), ("All families", "Bacillaceae removed")):
    w = W[(W.variant == v) & (W.contrast == "time 24->72 (both trt)")]
    for l1, d in w.groupby("PathwayL1"):
        c = COL.get(l1, "#9C9C9C")
        sig = (d.q_within < .05) | (d.q_turnover < .05)
        ax.scatter(d.turnover[~sig], d.within[~sig], s=10, facecolors="white", edgecolors=c, linewidths=0.6)
        ax.scatter(d.turnover[sig], d.within[sig], s=12, color=c, label=l1, linewidths=0)
    OFF = {"Ribosome": (-60, 3), "RNA degradation": (6, -8), "Protein export": (6, 4),
           "Photosynthesis": (6, -10), "Peptidoglycan biosynthesis": (5, -3), "Lysine biosynthesis": (5, -3),
           "Cationic antimicrobial peptide (CAMP) resistance": (-45, -24)}
    for _, r in w[w.within_down].iterrows():
        lab = LAB.get(r.function, r.function)
        if lab:
            ax.annotate(lab, (r.turnover, r.within), xytext=OFF.get(r.function, (4, -2)),
                        textcoords="offset points", fontsize=5.5)
    mk = w[w.masked]
    ax.scatter(mk.turnover, mk.within, s=34, facecolors="none", edgecolors="k", linewidths=0.6)
    lim = 0.65
    ax.plot([-lim, lim], [lim, -lim], color="#9C9C9C", lw=0.5, ls="--")
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlim(-lim, lim); ax.set_ylim(-0.45, 0.6)
    ax.set_xlabel("Contributor turnover, 24→72 h"); ax.set_title(title, fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axs[0].set_ylabel("Within-family change, 24→72 h")
h_, l_ = axs[0].get_legend_handles_labels()
seen = dict(zip(l_, h_))
fig.legend(seen.values(), seen.keys(), loc="lower center", ncol=4, frameon=False, fontsize=6, bbox_to_anchor=(0.5, -0.04))
fig.text(0.5, 0.98, "Filled: either component BH q < 0.05. Dashed line: turnover = −within (full cancellation). "
         "Black rings: masked maps (both components q < 0.05, opposite signs). Labels: significant within-family decline.", ha="center", fontsize=6)
fig.tight_layout(rect=(0, 0.07, 1, 0.95))
for ext in ("png", "pdf", "svg"):
    fig.savefig(V / f"masking_by_category.{ext}", dpi=300, bbox_inches="tight")
print("saved")
