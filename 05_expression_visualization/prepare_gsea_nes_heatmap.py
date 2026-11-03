#!/usr/bin/env python3
"""Build gsea_nes_selected.tsv (NES heatmap of selected pathways) from the GSEA run
with the official KEGG KO-to-map sets (04_pathway_enrichment/gsea_kegg_pathways.py).

Rows: the pathways to display for Bacteria and Fungi, in display order (edit PATHWAYS).
Columns: glucose effect at 24, 48 and 72 h and the treatment x time interaction
(G72 - C72) - (G24 - C24). Marks: * padj < 0.05, ° 0.05 <= padj < 0.10.


Run: python3 prepare_gsea_nes_heatmap.py <gsea_official_all.tsv> <out.tsv>
"""
import sys
import pandas as pd

GSEA, OUT = sys.argv[1:3]
PATHWAYS = {
    "Bacteria": ["ko00030", "ko00910", "ko02030", "ko02040", "ko00550", "ko03010",
                 "ko00010", "ko00620", "ko00020"],
    "Fungi": ["ko03010", "ko00010", "ko00620", "ko00020", "ko00280", "ko00650", "ko00410"],
}
COLS = [("glucose_G24_vs_C24", "G24/C24"), ("glucose_G48_vs_C48", "G48/C48"),
        ("glucose_G72_vs_C72", "G72/C72"), ("interaction_72_vs_24", "Δ effect\n72−24 h")]

g = pd.read_csv(GSEA, sep="\t")
rows = []
for dom, pws in PATHWAYS.items():
    for pw in pws:
        for ct, lab in COLS:
            r = g[(g.domain == dom) & (g.pathway == pw) & (g.contrast == ct)]
            if len(r) != 1:
                sys.exit(f"{dom} {pw} {ct}: {len(r)} rows")
            r = r.iloc[0]
            mark = "*" if r.padj < 0.05 else ("°" if r.padj < 0.10 else "")
            rows.append(dict(domain=dom, pathway=pw, pathway_name=r.pathway_name, col=lab,
                             NES=round(r.NES, 2), padj=r.padj, mark=mark, size=float(r["size"])))
pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
print(len(rows), "rows written to", OUT)
