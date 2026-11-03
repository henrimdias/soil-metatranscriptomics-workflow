#!/usr/bin/env python3
"""
Data for the differential-expression plots (counts per contrast, clustered
heatmaps of DE KOs, taxonomic attribution of DE KOs).

  de_counts_by_contrast.tsv          DE KO counts per contrast x domain x direction
  de_heatmap_bacteria_z.tsv          EVERY bacterial KO significant in >=1 contrast,
                                     row z-score over the six treatment x time groups
  de_heatmap_bacteria_clusters.tsv   cutree assignment for each row, in plot order
  de_heatmap_cluster_profiles.tsv    mean z profile of each cluster
  de_heatmap_dendrogram.tsv          dendrogram segments for drawing
  de_cluster_pathways_emapper.tsv    pathway enrichment within each cluster (eggNOG-mapper
                                     KEGG_Pathway sets; 04/de_cluster_pathway_enrichment.py
                                     repeats it with the official KEGG map membership)
  de_heatmap_fungi_z_wide.tsv        every fungal KO significant anywhere
  de_heatmap_fungi_annotation.tsv    their functional annotation
  de_family_attribution_plot.tsv     copied from de_family_attribution.tsv

WHY THE CLUSTERING LIVES HERE AND NOT IN R
The row order of a 1073-row heatmap is a RESULT, not a drawing detail: the
clusters get named in the text. Computing it here means the same ordering feeds
the plot, the cluster table and the enrichment test, and it is reproducible
without opening the plotting script.
"""
import os  # paths are read from environment variables (see README)
import pathlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, dendrogram
from scipy.spatial.distance import pdist
from scipy.stats import fisher_exact

R = pathlib.Path(os.environ.get('MTX_ROOT', '.')) / "results/deseq2"
OLD = pathlib.Path(os.path.join(os.environ.get('MTX_OUT', 'out'), 'de_summary'))
OUT = pathlib.Path(os.path.join(os.environ.get('MTX_OUT', 'out'), 'de_visualization')); OUT.mkdir(parents=True, exist_ok=True)

CONTRAST_ORDER = [
    "glucose_G24_vs_C24", "glucose_G48_vs_C48", "glucose_G72_vs_C72",
    "time_C48_vs_C24", "time_C72_vs_C24", "time_C72_vs_C48",
    "time_G48_vs_G24", "time_G72_vs_G24", "time_G72_vs_G48",
]
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLE_GROUP = {f"{pre}_{tp}_{r}": g
                for g, pre, tp in [("C24", "C", "24h"), ("C48", "C", "48h"),
                                   ("C72", "C", "72h"), ("G24", "G", "24h"),
                                   ("G48", "G", "48h"), ("G72", "G", "72h")]
                for r in (1, 2, 3)}
K = 5                                   # 5 profiles stay distinct; 6 splits a tail

de = pd.read_csv(R / "de_all_contrasts_long.tsv", sep="\t", low_memory=False)
de9 = de[de.contrast.isin(CONTRAST_ORDER)]
sig = de9[de9.significant == True]                                   # noqa: E712


def zmatrix(domain, kos):
    nc = pd.read_csv(R / f"normalized_counts_{domain.lower()}.tsv",
                     sep="\t", index_col=0)
    nc = nc.loc[[k for k in kos if k in nc.index]]
    l = np.log2(nc + 1)
    grp = l.T.groupby(l.columns.map(SAMPLE_GROUP)).mean().T[GROUPS]
    z = grp.sub(grp.mean(axis=1), axis=0) \
           .div(grp.std(axis=1).replace(0, np.nan), axis=0)
    return z.dropna(how="all")


# ---------------------------------------------------------------- 5A
for f in ("de_counts_by_contrast.tsv", "de_heatmap_fungi_z_wide.tsv", "de_heatmap_fungi_z.tsv",
          "de_heatmap_fungi_annotation.tsv"):
    (OUT / f).write_bytes((OLD / f).read_bytes())

# ---------------------------------------------------------------- 5B
kos_b = sorted(sig[sig.domain == "Bacteria"].KO.unique())
zb = zmatrix("Bacteria", kos_b)
print(f"[5B] bacterial KOs significant in >=1 contrast: {len(kos_b)}")
print(f"[5B] with a usable z profile (non-zero variance): {len(zb)}")

X = zb.to_numpy()
link = linkage(pdist(X, metric="euclidean"), method="ward")
order = leaves_list(link)                       # top-to-bottom order
zb_ord = zb.iloc[order]

# how many clusters? report the merge heights so the choice is inspectable
heights = link[-12:, 2][::-1]
print("[5B] top merge heights:", np.round(heights, 1))
cl = fcluster(link, t=K, criterion="maxclust")
cl_ord = cl[order]

# rename clusters 1..K in top-to-bottom order of first appearance, so the strip
# reads down the figure instead of in scipy's internal numbering
seen, remap = [], {}
for c in cl_ord:
    if c not in remap:
        remap[c] = len(seen) + 1
        seen.append(c)
cl_ord = np.array([remap[c] for c in cl_ord])

pd.DataFrame({"KO": zb_ord.index, "row": np.arange(1, len(zb_ord) + 1),
              "cluster": cl_ord}).to_csv(OUT / "de_heatmap_bacteria_clusters.tsv",
                                         sep="\t", index=False)
zb_ord.to_csv(OUT / "de_heatmap_bacteria_z.tsv", sep="\t")

# ------ dendrogram segments, exported rather than recomputed in R -------------
# The tree drawn beside the heatmap MUST be the tree the clusters came from.
# Recomputing it in R (hclust ward.D2) would almost certainly agree, but "almost
# certainly" is not a guarantee worth taking for a figure. scipy hands back the
# drawing coordinates directly: leaf j sits at icoord = 10*j + 5 in leaf order,
# dcoord is the merge height. Each merge is a bracket of three segments.
dd = dendrogram(link, no_plot=True)
seg = []
for xs, ys in zip(dd["icoord"], dd["dcoord"]):
    for k in range(3):
        seg.append(dict(leaf0=(xs[k] - 5) / 10, height0=ys[k],
                        leaf1=(xs[k + 1] - 5) / 10, height1=ys[k + 1]))
seg = pd.DataFrame(seg)
seg["row0"] = seg.leaf0 + 1            # 1 = top row of the heatmap
seg["row1"] = seg.leaf1 + 1
seg[["row0", "height0", "row1", "height1"]].to_csv(
    OUT / "de_heatmap_dendrogram.tsv", sep="\t", index=False)
print(f"[5B] dendrogram: {len(seg)} segments, max height {seg.height1.max():.1f}")

prof = zb_ord.copy(); prof["cluster"] = cl_ord
prof = prof.groupby("cluster")[GROUPS].mean()
prof.insert(0, "n", pd.Series(cl_ord).value_counts().sort_index().values)
prof.round(3).to_csv(OUT / "de_heatmap_cluster_profiles.tsv", sep="\t")
print("[5B] cluster sizes:", dict(zip(prof.index, prof.n)))

# ------ KEGG pathway enrichment inside each cluster (Fisher, BH within cluster)
pw = pd.read_csv(os.path.join(os.environ.get('MTX_ROOT', '.'), 'results/ko_matrices/ko_to_pathway.tsv'), sep="\t")
pw = (pw.assign(pathway=pw.KEGG_Pathway.astype(str).str.split(","))
        .explode("pathway"))
pw = pw[pw.pathway.str.startswith("ko")]        # drop the redundant map##### twins
# pathway names from the KEGG BRITE ko00001 hierarchy (reference/kegg/KO1-4.txt)
ko14 = pd.read_csv(os.path.join(os.environ.get('MTX_REF', os.path.join(os.environ.get('MTX_ROOT', '.'), 'reference')), 'kegg', 'KO1-4.txt'), sep="\t")
ko14["pid"] = ko14.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
names = (ko14.dropna(subset=["pid"]).drop_duplicates("pid").set_index("pid").Pathway
         .str.replace(r" \[PATH:.*", "", regex=True).to_dict())
memb = pw[pw.KO.isin(zb_ord.index)].groupby("pathway").KO.apply(set)
universe = set(zb_ord.index)
rows = []
for c in sorted(set(cl_ord)):
    members = set(zb_ord.index[cl_ord == c])
    for p, kos in memb.items():
        a = len(kos & members)
        if a < 3:
            continue
        b, cc = len(members) - a, len(kos & universe) - a
        d = len(universe) - a - b - cc
        _, pv = fisher_exact([[a, b], [cc, d]], alternative="greater")
        rows.append(dict(cluster=c, pathway=p, pathway_name=names.get(p, p),
                         n_in_cluster=a, n_in_universe=len(kos & universe),
                         p=pv))
enr = pd.DataFrame(rows)
if len(enr):
    enr["padj"] = enr.groupby("cluster").p.transform(
        lambda s: np.minimum.accumulate(
            (s.sort_values(ascending=False) * len(s)
             / np.arange(len(s), 0, -1)).clip(upper=1))[s.index])
    enr.sort_values(["cluster", "p"]).to_csv(OUT / "de_cluster_pathways_emapper.tsv",
                                             sep="\t", index=False)
    top = enr.sort_values(["cluster", "p"]).groupby("cluster").head(3)
    print("[5B] top pathway per cluster:")
    for _, r in top.iterrows():
        print(f"      c{r.cluster}  {r.pathway_name[:44]:<44} "
              f"{r.n_in_cluster}/{r.n_in_universe}  p={r.p:.2g}")

# ---------------------------------------------------------------- 6A, 6B
(OUT / "de_family_attribution_plot.tsv").write_bytes((OLD / "de_family_attribution.tsv").read_bytes())
(OUT / "de_family_global_share_plot.tsv").write_bytes(
    (OLD / "de_family_global_share.tsv").read_bytes())
print("[attribution]  family tables carried over unchanged")
