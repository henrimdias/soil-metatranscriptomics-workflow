"""Module trajectories (treatment x time) with the family of origin of each gene.

(a) Module score per sample = mean over the module's KOs of the per-KO z-score
    of log2(DESeq2-normalised count + 1) across the 18 samples (community
    normalisation, R DESeq2 size factors). Equal weight per KO.
(b) Family of origin: the KO x family x sample table (counts_ko_taxonomy_full)
    restricted to bacterial families; module counts per family per sample are
    divided by the community size factor, so the stacked contributions add up
    to the module's normalised expression. Families shown: those with >= 5% of
    the module's reads in any group; the rest pooled.
Also writes the share of each family in total bacterial mRNA per sample.
"""
import numpy as np
import pandas as pd
from common import *

mods = pd.read_csv(pathlib.Path(__file__).resolve().parent / "gene_sets" / "modules_definition.tsv", sep="\t")
raw = pd.read_csv(IN / "counts_ko_bacteria.tsv", sep="\t", index_col=0)[SAMPLES]
nrm = pd.read_csv(DE / "normalized_counts_bacteria.tsv", sep="\t", index_col=0)
nrm = nrm[SAMPLES] if set(SAMPLES) <= set(nrm.columns) else nrm
sf = (raw.loc[nrm.index] / nrm.replace(0, np.nan)).median()
sf = sf[SAMPLES]
pd.DataFrame({"size_factor": sf}).to_csv(OUT / "community_size_factors_bacteria.tsv", sep="\t")

# ---------------- (a) z-score trajectories
lg = np.log2(nrm[SAMPLES] + 1)
z = lg.sub(lg.mean(1), axis=0).div(lg.std(1, ddof=1).replace(0, np.nan), axis=0)
grp = pd.Series({s: group_of(s) for s in SAMPLES})
rows = []
for mod, g in mods.groupby("module"):
    kos = [k for k in g.KO.unique() if k in z.index]
    score = z.loc[kos].mean()
    for G in GROUPS:
        v = score[grp[grp == G].index]
        rows.append(dict(module=mod, group=G, n_KO=len(kos), mean_z=v.mean(),
                         se=v.std(ddof=1) / np.sqrt(len(v))))
    for sub, gg in g.groupby("submodule"):
        if mod != "Sporulation":
            continue
        kk = [k for k in gg.KO if k in z.index]
        sc = z.loc[kk].mean()
        for G in GROUPS:
            v = sc[grp[grp == G].index]
            rows.append(dict(module=f"Sporulation: {sub}", group=G, n_KO=len(kk),
                             mean_z=v.mean(), se=v.std(ddof=1) / np.sqrt(len(v))))
traj = pd.DataFrame(rows)
traj.to_csv(OUT / "module_trajectories_z.tsv", sep="\t", index=False)

# ---------------- (b) family of origin
fam_b = pd.read_csv(UP / "inputs/"
                    "counts_family_bacteria.tsv", sep="\t", index_col=0)
long = pd.read_csv(UP / "inputs/"
                   "counts_ko_taxonomy_full.tsv", sep="\t")
long = long[long.family.isin(fam_b.index)]
long.to_parquet(OUT / "_long_bacteria.parquet")

share = fam_b[SAMPLES] / fam_b[SAMPLES].sum()
share.to_csv(OUT / "family_share_of_bacterial_mRNA.tsv", sep="\t")

frows = []
for mod, g in mods.groupby("module"):
    sub = long[long.KO.isin(set(g.KO))]
    m = sub.pivot_table(index="family", columns="sample_id", values="count",
                        aggfunc="sum", fill_value=0).reindex(columns=SAMPLES, fill_value=0)
    m = m / sf  # normalised to the community library
    gm = m.T.groupby(grp).mean().T[GROUPS]
    pct = gm / gm.sum() * 100
    keep = pct.index[(pct >= 5).any(axis=1)]
    for fam in list(keep) + ["Other families"]:
        vals = gm.loc[fam] if fam != "Other families" else gm.drop(index=keep).sum()
        for G in GROUPS:
            frows.append(dict(module=mod, family=fam, group=G, norm_expr=vals[G],
                              pct_of_module=100 * vals[G] / gm[G].sum()))
fam = pd.DataFrame(frows)
fam.to_csv(OUT / "module_family_of_origin.tsv", sep="\t", index=False)

print(traj[~traj.module.str.startswith("Sporulation:")].pivot(
    index="module", columns="group", values="mean_z")[GROUPS].round(2).to_string())
top = fam[fam.family == "Bacillaceae"].pivot(index="module", columns="group",
                                               values="pct_of_module")[GROUPS]
print("\n% of module expression from Bacillaceae:\n", top.round(1).to_string())
print("\nBacillaceae share of bacterial mRNA by group (%):")
print((share.loc["Bacillaceae"].groupby(grp).mean()[GROUPS] * 100).round(2).to_string())
