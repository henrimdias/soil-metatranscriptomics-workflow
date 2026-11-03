#!/usr/bin/env python3
"""
Build the per-domain KO x 18 raw-count matrices for the DESeq2 scope.

INPUT
  Enrichement_analysis_topActive/counts_ko_taxonomy_full.tsv  (~108 MiB, LONG format:
  sample_id, treatment, timepoint, KO, family, count)

The table carries `family` but NO superkingdom column, so the domain split is done
by joining against the per-domain family tables written by 12a_parse_taxonomy_counts.py:
  inputs/counts_family_bacteria.tsv  -> Bacteria
  inputs/counts_family_fungi.tsv     -> Fungi
  inputs/counts_family_archaea.tsv   -> Archaea (DROPPED)
The three lists are verified disjoint. Anything not in them (unclassified,
non-fungal eukaryotes) is dropped and reported.

FILTER (from the scope): total count >= 10 AND detected in >= 3 of the 18 samples.
The >= 3 rule is what keeps a KO expressed only in one treatment x timepoint cell
(n = 3) from being discarded before it is ever tested.

OUTPUTS (small, committed back to the Mac; DESeq2 reads these, never the 108 MiB file)
  counts_ko_bacteria.tsv / counts_ko_fungi.tsv   KO x 18 raw counts, filtered
  ko_top_family_bacteria.tsv / ..._fungi.tsv     per KO: family contributing most reads
  ko_matrix_build_stats.tsv                      row counts at every filtering step
"""
import os  # paths are read from environment variables (see README)
import pathlib

import pandas as pd

U = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
E = U / "inputs"
T = U / "inputs"
S = U / "inputs"
OUT = U / "results/ko_matrices"
OUT.mkdir(parents=True, exist_ok=True)

MIN_TOTAL = 10
MIN_SAMPLES = 3

meta = pd.read_csv(U / "inputs/sample_metadata.tsv", sep="\t")
SAMPLES = list(meta.sample_id)
assert len(SAMPLES) == 18

fam = {}
for dom, path in [("Bacteria", T / "counts_family_bacteria.tsv"),
                  ("Fungi", T / "counts_family_fungi.tsv"),
                  ("Archaea", S / "counts_family_archaea.tsv")]:
    fam[dom] = set(pd.read_csv(path, sep="\t").iloc[:, 0])
assert not (fam["Bacteria"] & fam["Fungi"]), "family lists overlap"
assert not (fam["Bacteria"] & fam["Archaea"])
assert not (fam["Fungi"] & fam["Archaea"])
fam2dom = {f: d for d, fs in fam.items() for f in fs}

print("reading the long table ...")
long = pd.read_csv(E / "counts_ko_taxonomy_full.tsv", sep="\t",
                   dtype={"sample_id": "category", "KO": "category",
                          "family": "category", "count": "int64"},
                   usecols=["sample_id", "KO", "family", "count"])
print(f"  {len(long):,} rows | {long.KO.nunique():,} KOs | {long.family.nunique():,} families")
assert set(long.sample_id.unique()) == set(SAMPLES), "sample ids do not match the metadata"

long["domain"] = long.family.astype(str).map(fam2dom).fillna("other")
by_dom = long.groupby("domain", observed=True)["count"].sum()
print("\nreads by domain assignment:")
for d, v in by_dom.sort_values(ascending=False).items():
    print(f"  {d:10s} {v:>12,}  ({100*v/by_dom.sum():5.2f}%)")

stats = []
for dom in ("Bacteria", "Fungi"):
    sub = long[long.domain == dom]

    wide = (sub.pivot_table(index="KO", columns="sample_id", values="count",
                            aggfunc="sum", fill_value=0, observed=True)
              .reindex(columns=SAMPLES, fill_value=0)
              .astype("int64"))
    wide.index.name = "KO"
    n_all = len(wide)

    keep = (wide.sum(axis=1) >= MIN_TOTAL) & ((wide > 0).sum(axis=1) >= MIN_SAMPLES)
    filt = wide[keep]

    # how much each criterion removes on its own, for the record
    only_total = int((~(wide.sum(axis=1) >= MIN_TOTAL)).sum())
    only_prev = int((~((wide > 0).sum(axis=1) >= MIN_SAMPLES)).sum())

    # the family contributing the most reads to each surviving KO
    top = (sub[sub.KO.isin(filt.index)]
           .groupby(["KO", "family"], observed=True)["count"].sum()
           .reset_index()
           .sort_values("count", ascending=False)
           .drop_duplicates("KO")
           .rename(columns={"family": "top_family", "count": "top_family_reads"}))
    tot = filt.sum(axis=1).rename("total_reads").reset_index()
    top = top.merge(tot, on="KO")
    top["top_family_pct"] = (100 * top.top_family_reads / top.total_reads).round(2)
    top["domain"] = dom
    top = top.sort_values("KO")[["KO", "domain", "top_family", "top_family_reads",
                                 "total_reads", "top_family_pct"]]

    tag = dom.lower()
    filt.to_csv(OUT / f"counts_ko_{tag}.tsv", sep="\t")
    top.to_csv(OUT / f"ko_top_family_{tag}.tsv", sep="\t", index=False)

    stats.append(dict(domain=dom, families=sub.family.nunique(),
                      ko_before_filter=n_all, ko_after_filter=int(keep.sum()),
                      removed_by_total_lt10=only_total,
                      removed_by_prevalence_lt3=only_prev,
                      reads_before=int(wide.values.sum()),
                      reads_after=int(filt.values.sum()),
                      min_library=int(filt.sum(axis=0).min()),
                      max_library=int(filt.sum(axis=0).max())))
    print(f"\n{dom}: {n_all:,} KOs -> {keep.sum():,} after filter "
          f"({100*keep.sum()/n_all:.1f}% kept, {100*filt.values.sum()/wide.values.sum():.2f}% of reads)")
    print(f"  library sizes {filt.sum(axis=0).min():,} .. {filt.sum(axis=0).max():,}")

st = pd.DataFrame(stats)
st.to_csv(OUT / "ko_matrix_build_stats.tsv", sep="\t", index=False)
print("\n" + st.to_string(index=False))
print(f"\nwrote to {OUT}")
