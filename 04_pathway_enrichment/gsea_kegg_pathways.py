"""GSEA of KEGG pathways with OFFICIAL map membership.

Why: ko_to_pathway.tsv / the KEGG_Pathway column of de_all_contrasts_long.tsv
come from eggNOG-mapper hits. A hit annotated with several KOs carries the
UNION of their pathways, and that union was attached to every KO of the hit.
Result: ~half the members of a typical pathway set in eggNOG-derived sets are not in
that KEGG map (e.g. 'Bacterial chemotaxis' had 109 KOs, 26 of them real
ko02030 members). Here the sets are rebuilt from the ko00001 hierarchy
(KO -> [PATH:koXXXXX]); the GSEA engine is in common.py.
"""
import sys
import pandas as pd
from common import *

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
DOMS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["Bacteria", "Fungi"]

h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
h = h.dropna(subset=["pid"])
pname = h.drop_duplicates("pid").set_index("pid").Pathway.str.replace(r" \[PATH:.*", "", regex=True)
official = h.groupby("pid").KO.apply(set).to_dict()

de = pd.read_csv(DE / "de_all_contrasts_long.tsv", sep="\t", low_memory=False)
allres = []
for dom in DOMS:
    counts = pd.read_csv(IN / f"counts_ko_{dom.lower()}.tsv", sep="\t", index_col=0)
    ko2rep, rep2mem = equivalence_groups(counts)
    uni = set(counts.index)
    sets = {p: m & uni for p, m in official.items() if len(m & uni) > 0}
    d = de[de.domain == dom]
    stats = {ct: g.set_index("KO").stat for ct, g in d.groupby("contrast")}
    print(f"{dom}: {len(uni)} KOs, {len(rep2mem)} features, {len(sets)} official maps", flush=True)
    r = run_collection(stats, sets, ko2rep, NPERM, 20261009, tag=dom)
    r.insert(0, "domain", dom)
    allres.append(r)
    r.to_csv(OUT / f"gsea_official_{dom}.tsv", sep="\t", index=False)

a = pd.concat(allres, ignore_index=True)
a = a.rename(columns={"set": "pathway"})
a["pathway_name"] = a.pathway.map(pname)
a["display_kept"] = a.pathway.str[2:].astype(int) < 4000
a.to_csv(OUT / "gsea_official_all.tsv", sep="\t", index=False)
print("done", len(a), "tests,", int((a.padj < 0.05).sum()), "padj<0.05")
