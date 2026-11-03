#!/usr/bin/env python3
"""
KO annotation rate per domain — the numerator/denominator table.

🔴 The one thing that must not be got wrong here: the `count` column of
counts_ko_taxonomy_full.tsv is **KO ASSIGNMENTS, NOT READS**. In
13b_parse_ko_taxonomy_full.py the read counter is incremented once
(`ko_hits += 1`) but the table is written inside `for ko in kos:`, so a read
whose UniRef90 hit maps to several KOs contributes one row per KO. Globally that
inflates the count column by a factor of 13,432,047 / 12,372,058 = 1.0857 (8.57%).
Reporting the raw column as "reads" would overstate the annotation rate.

Denominator = taxonomically classified mRNA reads, taken from the same family
tables written by 12a_parse_taxonomy_counts.py (counts_family_{bacteria,fungi,archaea}).
Those tables include `unclassified <order> family` placeholders (~9% of reads in
every domain) — those ARE family-level assignments within a classified lineage,
not unclassified reads, so they belong in the denominator.

Outputs ko_annotation_rate_by_domain.tsv
"""
import os  # paths are read from environment variables (see README)
import pandas as pd

U = os.path.join(os.environ.get('MTX_ROOT', '.'), '')
E = U + "inputs/"
T = U + "inputs/"
S = U + "inputs/"

stats = pd.read_csv(E + "ko_taxonomy_full_stats.tsv", sep="\t")
TOTAL_HITS = int(stats.total_hits.sum())     # every DIAMOND/UniRef90 hit
READS_W_KO = int(stats.ko_hits.sum())        # reads with >= 1 KO (counted once)

long = pd.read_csv(E + "counts_ko_taxonomy_full.tsv", sep="\t",
                   usecols=["KO", "family", "count"],
                   dtype={"KO": "category", "family": "category", "count": "int64"})
ASSIGNMENTS = int(long["count"].sum())
INFLATION = ASSIGNMENTS / READS_W_KO

paths = {"Bacteria": T + "counts_family_bacteria.tsv",
         "Fungi": T + "counts_family_fungi.tsv",
         "Archaea": S + "counts_family_archaea.tsv"}
famset, taxreads = {}, {}
for dom, p in paths.items():
    d = pd.read_csv(p, sep="\t", index_col=0)
    famset[dom] = set(d.index)
    taxreads[dom] = int(d.values.sum())
assert not (famset["Bacteria"] & famset["Fungi"] & famset["Archaea"])

f2d = {f: d for d, fs in famset.items() for f in fs}
long["domain"] = long.family.astype(str).map(f2d).fillna("unclassified / other")
ko_assign = long.groupby("domain", observed=True)["count"].sum()

rows = []
for dom in ("Bacteria", "Fungi", "Archaea"):
    a = int(ko_assign.get(dom, 0))
    r = a / INFLATION                       # assignments -> reads, global factor
    rows.append(dict(domain=dom,
                     tax_classified_mRNA_reads=taxreads[dom],
                     ko_assignments=a,
                     est_ko_assigned_reads=round(r),
                     ko_annotation_rate_pct=round(100 * r / taxreads[dom], 2),
                     naive_rate_assignments_pct=round(100 * a / taxreads[dom], 2)))
out = pd.DataFrame(rows)

unmatched_reads = TOTAL_HITS - sum(taxreads.values())
print(f"total DIAMOND hits                 : {TOTAL_HITS:,}")
print(f"reads with >= 1 KO                 : {READS_W_KO:,}  ({100*READS_W_KO/TOTAL_HITS:.2f}%)")
print(f"KO assignments (count column)      : {ASSIGNMENTS:,}")
print(f"multi-KO inflation factor          : {INFLATION:.4f}")
print(f"classified to family (3 domains)   : {sum(taxreads.values()):,}  "
      f"({100*sum(taxreads.values())/TOTAL_HITS:.2f}%)")
print(f"unclassified / other eukaryotes    : {unmatched_reads:,}\n")
print(out.to_string(index=False))
out.to_csv(os.path.join(os.environ.get('MTX_OUT', 'out'), 'ko_annotation_rate_by_domain.tsv'), sep="\t", index=False)
