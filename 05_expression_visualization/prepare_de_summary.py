#!/usr/bin/env python3
"""
Summarise the DESeq2 output for plotting (DE counts, heatmap inputs, taxonomic attribution).

Writes small tidy TSVs that the R plotting script consumes, so the plotting code
contains no data wrangling and the plotted numbers are provably the ones
in de_all_contrasts_long.tsv.

de_counts_by_contrast.tsv          DE KO counts per contrast x domain x direction
de_top50_bacteria_z.tsv            top-50 bacterial KOs, row z-score of mean log2(norm+1)
de_heatmap_fungi_z.tsv             every fungal KO significant anywhere, same transform
de_heatmap_fungi_annotation.tsv    its description + pathway codes
de_family_attribution.tsv          family composition of the significant KOs per contrast
"""
import os  # paths are read from environment variables (see README)
import pathlib

import numpy as np
import pandas as pd

R = pathlib.Path(os.environ.get('MTX_ROOT', '.')) / "results/deseq2"
OUT = pathlib.Path(os.path.join(os.environ.get('MTX_OUT', 'out'), 'de_summary')); OUT.mkdir(exist_ok=True)

CONTRAST_ORDER = [
    "glucose_G24_vs_C24", "glucose_G48_vs_C48", "glucose_G72_vs_C72",
    "time_C48_vs_C24", "time_C72_vs_C24", "time_C72_vs_C48",
    "time_G48_vs_G24", "time_G72_vs_G24", "time_G72_vs_G48",
]
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLE_GROUP = {}
for g, pre, tp in [("C24", "C", "24h"), ("C48", "C", "48h"), ("C72", "C", "72h"),
                   ("G24", "G", "24h"), ("G48", "G", "48h"), ("G72", "G", "72h")]:
    for r in (1, 2, 3):
        SAMPLE_GROUP[f"{pre}_{tp}_{r}"] = g

de = pd.read_csv(R / "de_all_contrasts_long.tsv", sep="\t", low_memory=False)
de = de[de.contrast.isin(CONTRAST_ORDER)]          # interaction panels not shown
sig = de[de.significant == True]                    # noqa: E712

# ---------------------------------------------------------------- A
rows = []
for dom in ("Bacteria", "Fungi"):
    for ct in CONTRAST_ORDER:
        x = sig[(sig.domain == dom) & (sig.contrast == ct)]
        lab = de[de.contrast == ct].contrast_label.iloc[0]
        fam = de[de.contrast == ct].contrast_family.iloc[0]
        rows.append(dict(domain=dom, contrast=ct, contrast_label=lab,
                         contrast_family=fam, direction="Upregulated",
                         n=int((x.log2FoldChange > 0).sum())))
        rows.append(dict(domain=dom, contrast=ct, contrast_label=lab,
                         contrast_family=fam, direction="Downregulated",
                         n=int((x.log2FoldChange < 0).sum())))
A = pd.DataFrame(rows)
A.to_csv(OUT / "de_counts_by_contrast.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- B and C
def zmatrix(domain, kos):
    nc = pd.read_csv(R / f"normalized_counts_{domain.lower()}.tsv", sep="\t", index_col=0)
    nc = nc.loc[[k for k in kos if k in nc.index]]
    l = np.log2(nc + 1)
    grp = l.T.groupby(l.columns.map(SAMPLE_GROUP)).mean().T[GROUPS]
    z = grp.sub(grp.mean(axis=1), axis=0).div(grp.std(axis=1).replace(0, np.nan), axis=0)
    return z.dropna(how="all")

# B: rank bacterial KOs by their best padj across the nine contrasts
bb = (sig[sig.domain == "Bacteria"].groupby("KO").padj.min()
      .sort_values().head(50).index.tolist())
zb = zmatrix("Bacteria", bb)
zb.reset_index().melt(id_vars="KO", var_name="group", value_name="z") \
  .to_csv(OUT / "de_top50_bacteria_z.tsv", sep="\t", index=False)
zb.to_csv(OUT / "de_top50_bacteria_z_wide.tsv", sep="\t")

# C: every fungal KO significant in any of the nine contrasts
fk = sorted(sig[sig.domain == "Fungi"].KO.unique())
zf = zmatrix("Fungi", fk)
zf.reset_index().melt(id_vars="KO", var_name="group", value_name="z") \
  .to_csv(OUT / "de_heatmap_fungi_z.tsv", sep="\t", index=False)
zf.to_csv(OUT / "de_heatmap_fungi_z_wide.tsv", sep="\t")

ann = (de[(de.domain == "Fungi") & (de.KO.isin(fk))]
       .drop_duplicates("KO")[["KO", "gene_name", "KEGG_Pathway", "top_family"]])
best = sig[sig.domain == "Fungi"].groupby("KO").agg(
    best_padj=("padj", "min"),
    contrasts=("contrast_label", lambda s: "; ".join(sorted(set(s))))).reset_index()
ann = ann.merge(best, on="KO").sort_values("best_padj")
ann.to_csv(OUT / "de_heatmap_fungi_annotation.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- D
# Which families get a name instead of falling into the grey bucket.
#
# A GLOBAL ">2% of the domain's calls" threshold is a no-op here: it returns
# exactly the same five families as the old top-5 rule, because Micromonosporaceae
# (2.37%) and Nitrobacteraceae (2.16%) are the last above it and the sixth family
# is already at 1.90%. The threshold only does useful work applied WITHIN a
# contrast, which is what lets families that matter to the small glucose
# contrasts be named even though the huge temporal contrasts dominate the totals.
#
# Rule: a family is named if, in any contrast with >= MIN_N significant KOs, it
# accounts for > PCT_MIN of that contrast's KOs AND is the top family for at
# least MIN_CALLS of them (so no family is named on the strength of one gene).
# Fallback: in a domain with very few families overall, name them all - Fungi has
# only three and no contrast reaches MIN_N, so the rule alone would grey out the
# entire panel.
PCT_MIN, MIN_N, MIN_CALLS, SMALL_DOMAIN = 2.0, 5, 2, 8

d_rows = []
for dom in ("Bacteria", "Fungi"):
    s = sig[sig.domain == dom]
    named = set()
    for ct in CONTRAST_ORDER:
        x = s[s.contrast == ct]
        if len(x) < MIN_N:
            continue
        g = x.groupby("top_family").KO.size()
        p = 100 * g / len(x)
        named |= set(p[(p > PCT_MIN) & (g >= MIN_CALLS)].index)
    if s.top_family.nunique() <= SMALL_DOMAIN:
        named = set(s.top_family.unique())

    for ct in CONTRAST_ORDER:
        x = s[s.contrast == ct]
        n = len(x)
        if n == 0:
            d_rows.append(dict(domain=dom, contrast=ct, family="No DE KOs",
                               n=0, pct=np.nan, n_total=0))
            continue
        vc = x.top_family.apply(lambda f: f if f in named else "Other families").value_counts()
        for fam, k in vc.items():
            d_rows.append(dict(domain=dom, contrast=ct, family=fam,
                               n=int(k), pct=100 * k / n, n_total=n))
D = pd.DataFrame(d_rows)
D.to_csv(OUT / "de_family_attribution.tsv", sep="\t", index=False)

# global share of every named family, for the ordering and for reporting
glob = (sig.groupby(["domain", "top_family"]).KO.size().rename("n_calls")
        .reset_index())
glob["pct_of_domain"] = glob.groupby("domain").n_calls.transform(lambda x: 100 * x / x.sum())
glob.sort_values(["domain", "n_calls"], ascending=[True, False]).to_csv(
    OUT / "de_family_global_share.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- report
print("PANEL A — counts (what will be printed above each bar)")
tot = A.groupby(["domain", "contrast_label"], sort=False).n.sum().reset_index()
for dom in ("Bacteria", "Fungi"):
    t = tot[tot.domain == dom]
    print(f"  {dom:9s} " + "  ".join(f"{r.contrast_label}={r.n}" for r in t.itertuples()))
print(f"\nPANEL B — {zb.shape[0]} bacterial KOs x {zb.shape[1]} groups")
print(f"PANEL C — {zf.shape[0]} fungal KOs (all significant anywhere)")
print(ann[["KO", "gene_name", "best_padj", "top_family"]].to_string(index=False))
print("\nPANEL D — families used")
for dom in ("Bacteria", "Fungi"):
    fams = [f for f in D[D.domain == dom].family.unique() if f != "No DE KOs"]
    print(f"  {dom}: {fams}")
print(f"\nwrote to {OUT}")
