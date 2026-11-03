#!/usr/bin/env python3
# =============================================================================
# Script:
# 13b_parse_ko_taxonomy_full.py
# Variant of 13_parse_ko_counts.py that keeps the family of origin of every
# KO assignment. Two changes relative to that script:
#
#   1. The `if ko in CARBON_KOS:` gate that restricted the family-level
#      breakdown (ko_tax_counts) to the 5 curated carbon-metabolism blocks is
#      REMOVED. All annotated KOs now get a family breakdown, not just the
#      curated subset -- needed for family-resolved analyses (taxonomic
#      attribution of DE KOs, taxonomic decomposition, family-level enrichment).
#   2. KEGG_Pathway (eggNOG annotation column 13, 0-indexed 12) is parsed
#      alongside KO (column 12, 0-indexed 11) in the same pass over
#      EMAPPER_FILE, producing a KO -> KEGG pathway lookup table. This avoids
#      a separate KEGG REST API call, since eggNOG-mapper already carries
#      this annotation per hit.
#
# Everything else (Steps 1-2, DIAMOND parsing logic, output for counts_ko.tsv
# and counts_ko_carbon.tsv) is UNCHANGED from the original script, so this can
# be run as a drop-in replacement without affecting the existing filtered
# outputs (counts_ko_carbon.tsv, counts_ko_taxonomy.tsv from script 13 are
# untouched on disk; this script writes to NEW filenames, see Outputs below).
#
# Author (script):
# Henrique M. Dias
#
# Affiliation:
# South Dakota State University
#
# Date:
# 2026
#
# Inputs (identical to 13_parse_ko_counts.py):
#   EMAPPER_FILE : eggNOG-mapper annotations (uniref90_hits.emapper.annotations)
#   TAXID_FILE   : UniRef90_ID -> TaxID (uniref90_hits_taxid.tsv)
#   LINEAGE_FILE : TaxID -> lineage (taxid_lineage.tsv)
#   DIAMOND_DIR  : per-sample DIAMOND blastx outputs (*.uniref90.tsv.gz)
#   META_FILE    : sample metadata (sample_metadata.tsv)
#
# Outputs (NEW filenames -- does not overwrite script 13's outputs):
#   counts_ko_taxonomy_full.tsv : long-format KO x family x sample counts,
#                                 ALL annotated KOs (not just carbon blocks)
#   ko_to_pathway.tsv           : KO -> KEGG_Pathway lookup (one row per KO,
#                                 pathway field is the raw comma-separated
#                                 eggNOG string, e.g. "ko00010,map00010")
#   ko_taxonomy_full_stats.tsv  : per-sample annotation stats (same shape as
#                                 ko_parsing_stats.tsv from script 13)
#
# Usage:
#   Submit via SLURM (reuse the run_13_parse_ko_counts.sh wrapper, just point
#   it at this script instead) -- same resource profile as script 13
#   (emapper file is 1.6GB; ~32GB RAM; 1-2h runtime), since it re-reads the
#   same inputs. No new heavy dependency is introduced.
#
# Dependencies:
#   Python 3.x, pandas (included in metatx conda env)
# =============================================================================

import os
import gzip
import pandas as pd
from collections import defaultdict

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR     = os.environ.get('PROJ', '.')
EMAPPER_FILE = os.path.join(BASE_DIR, "results/emapper/uniref90_hits.emapper.annotations")
TAXID_FILE   = os.path.join(BASE_DIR, "db/uniref90/uniref90_hits_taxid.tsv")
LINEAGE_FILE = os.path.join(BASE_DIR, "db/uniref90/taxid_lineage.tsv")
DIAMOND_DIR  = os.path.join(BASE_DIR, "results/diamond/uniref90")
META_FILE    = os.path.join(BASE_DIR, "sample_metadata.tsv")
OUT_DIR      = os.path.join(BASE_DIR, "results/counts")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: Load TaxID → lineage (domain + family) ───────────────────────────
# Identical to 13_parse_ko_counts.py Step 1.
print("Step 1: Loading taxonkit lineage...")
taxid_to_info = {}
skipped = 0

with open(LINEAGE_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            skipped += 1
            continue
        taxid, lineage, rank, kingdom, phylum, clas, order, family, genus = parts
        if family in ("__", ""):
            skipped += 1
            continue
        if ";Bacteria;" in lineage or lineage.startswith("Bacteria"):
            domain = "Bacteria"
        elif ";Eukaryota;" in lineage or lineage.startswith("Eukaryota"):
            domain = "Eukaryota"
        elif ";Archaea;" in lineage or lineage.startswith("Archaea"):
            domain = "Archaea"
        else:
            skipped += 1
            continue
        taxid_to_info[taxid] = (domain, phylum, family)

print(f"  TaxIDs with lineage: {len(taxid_to_info):,} | Skipped: {skipped:,}")

# ── Step 2: Load UniRef90_ID → TaxID ─────────────────────────────────────────
# Identical to 13_parse_ko_counts.py Step 2.
print("\nStep 2: Loading UniRef90_ID → TaxID mapping...")
uniref_to_taxinfo = {}

with open(TAXID_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        uniref_id, taxid = parts
        if taxid in taxid_to_info:
            uniref_to_taxinfo[uniref_id] = taxid_to_info[taxid]

print(f"  UniRef90 IDs with taxonomy: {len(uniref_to_taxinfo):,}")

# ── Step 3: Load eggNOG-mapper annotations → UniRef90_ID → KOs + pathways ────
# CHANGED: also parse KEGG_Pathway (column 13, 0-indexed 12) in the same pass.
print("\nStep 3: Loading eggNOG-mapper KO + KEGG_Pathway annotations...")
uniref_to_kos = {}       # uniref_id -> list of KO terms
ko_to_pathway = {}       # KO -> raw KEGG_Pathway string (first one seen wins)

with open(EMAPPER_FILE) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 13:
            continue
        uniref_id   = parts[0]
        ko_field    = parts[11]
        pathway_field = parts[12]
        if ko_field == "-" or ko_field == "":
            continue
        kos = [k.replace("ko:", "") for k in ko_field.split(",")]
        uniref_to_kos[uniref_id] = kos
        if pathway_field not in ("-", ""):
            for ko in kos:
                if ko not in ko_to_pathway:
                    ko_to_pathway[ko] = pathway_field

print(f"  UniRef90 IDs with KO annotation: {len(uniref_to_kos):,}")
print(f"  Unique KOs with a KEGG_Pathway : {len(ko_to_pathway):,}")

# ── Step 4: Parse DIAMOND outputs and build FULL taxonomy breakdown ─────────
# CHANGED: the `if ko in CARBON_KOS` gate is removed -- every annotated KO
# now contributes to ko_tax_counts, not just the 5 curated blocks.
print("\nStep 4: Parsing DIAMOND outputs (full KO x family breakdown)...")
metadata = pd.read_csv(META_FILE, sep="\t")

ko_tax_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
stats = []

for _, row in metadata.iterrows():
    sample_id    = row["sample_id"]
    file_prefix  = row["file_prefix"]
    diamond_file = os.path.join(DIAMOND_DIR, f"{file_prefix}.uniref90.tsv.gz")

    if not os.path.exists(diamond_file):
        print(f"  WARNING: {diamond_file} not found, skipping.")
        continue

    total_hits = 0
    ko_hits    = 0

    with gzip.open(diamond_file, "rt") as fh:
        for line in fh:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sseqid = parts[1]
            total_hits += 1

            kos = uniref_to_kos.get(sseqid)
            if not kos:
                continue

            ko_hits += 1
            taxinfo = uniref_to_taxinfo.get(sseqid)
            family  = taxinfo[2] if taxinfo else "unclassified"

            for ko in kos:
                ko_tax_counts[sample_id][ko][family] += 1   # no CARBON_KOS gate

    stats.append({
        "sample_id"  : sample_id,
        "treatment"  : row["treatment"],
        "timepoint"  : row["timepoint"],
        "total_hits" : total_hits,
        "ko_hits"    : ko_hits,
        "pct_ko"     : round(ko_hits / total_hits * 100, 2) if total_hits > 0 else 0,
    })
    print(f"  {sample_id}: {total_hits:,} hits | KO: {ko_hits:,} ({stats[-1]['pct_ko']}%)")

# ── Step 5: Build and save outputs ────────────────────────────────────────────
print("\nStep 5: Building full KO x family long-format table...")

sample_ids = [row["sample_id"] for _, row in metadata.iterrows()
              if os.path.exists(os.path.join(DIAMOND_DIR, f"{row['file_prefix']}.uniref90.tsv.gz"))]

tax_rows = []
for sample_id in sample_ids:
    sample_row = metadata[metadata["sample_id"] == sample_id].iloc[0]
    for ko, family_counts in ko_tax_counts[sample_id].items():
        for family, count in family_counts.items():
            tax_rows.append({
                "sample_id": sample_id,
                "treatment": sample_row["treatment"],
                "timepoint": sample_row["timepoint"],
                "KO"       : ko,
                "family"   : family,
                "count"    : count
            })
tax_df = pd.DataFrame(tax_rows)

pathway_df = pd.DataFrame(
    [{"KO": ko, "KEGG_Pathway": pw} for ko, pw in ko_to_pathway.items()]
)

stats_df = pd.DataFrame(stats)

# Save -- NEW filenames, does not touch script 13's existing outputs
tax_out      = os.path.join(OUT_DIR, "counts_ko_taxonomy_full.tsv")
pathway_out  = os.path.join(OUT_DIR, "ko_to_pathway.tsv")
stats_out    = os.path.join(OUT_DIR, "ko_taxonomy_full_stats.tsv")

tax_df.to_csv(tax_out,         sep="\t", index=False)
pathway_df.to_csv(pathway_out, sep="\t", index=False)
stats_df.to_csv(stats_out,     sep="\t", index=False)

print(f"\nOutputs saved:")
print(f"  Full KO x family breakdown : {tax_out}  ({len(tax_df):,} rows)")
print(f"  KO -> KEGG_Pathway lookup  : {pathway_out}  ({len(pathway_df):,} KOs)")
print(f"  Stats                     : {stats_out}")
print("\nDone.")
