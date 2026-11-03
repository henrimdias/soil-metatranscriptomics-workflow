"""Pathway over-representation in the row clusters of the DE heatmap, with OFFICIAL
KEGG map membership (ko00001) instead of the eggNOG-mapper KEGG_Pathway union sets. Fisher, one-sided,
universe = all clustered KOs, BH within cluster, maps with >= 3 hits."""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from common import *
h = pd.read_csv(KO14, sep="\t"); h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
h = h.dropna(subset=["pid"])
pname = h.drop_duplicates("pid").set_index("pid").Pathway.str.replace(r" \[PATH:.*", "", regex=True)
cl = pd.read_csv(UP / "results/de_visualization/de_heatmap_bacteria_clusters.tsv", sep="\t")
U = set(cl.KO); mem = h[h.KO.isin(U)].groupby("pid").KO.apply(set)
rows = []
for c, g in cl.groupby("cluster"):
    S = set(g.KO)
    for p, M in mem.items():
        a = len(S & M)
        if a < 3: continue
        b = len(S) - a; cc = len(M - S); d = len(U) - len(S) - cc
        rows.append(dict(cluster=c, pathway=p, name=pname[p], hits=a, map_in_1073=len(M),
                         odds=fisher_exact([[a, b], [cc, d]], alternative="greater")[0],
                         p=fisher_exact([[a, b], [cc, d]], alternative="greater")[1]))
r = pd.DataFrame(rows)
r["padj"] = r.groupby("cluster").p.transform(lambda x: multipletests(x, method="fdr_bh")[1])
r = r.sort_values(["cluster", "p"])
r.to_csv(OUT / "de_cluster_pathways_kegg.tsv", sep="\t", index=False)
old = pd.read_csv(UP / "results/de_visualization/de_cluster_pathways_emapper.tsv", sep="\t")
print(old.head(3).to_string())
print(r[r.p < 0.05].groupby("cluster").head(4).round(4).to_string(index=False))
