#!/usr/bin/env python3
# =============================================================================
# Script:
# Build per-sample read count matrices aggregated by taxonomic family,
# by joining DIAMOND blastx UniRef90 hits with taxonkit lineage annotations.
# Generates separate count matrices for Bacteria, Fungi, and Archaea for
# comparison with 16S and ITS amplicon data respectively.
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
# Description:
# Three-step pipeline:
#   Step 1 — Load TaxID → lineage (from script 11, taxonkit output)
#            Domain (Bacteria/Eukaryota/Archaea) is parsed from the lineage
#            string (col 2) because taxonkit {k} format returns
#            "unclassified X superkingdom" for all entries — NCBI taxonomy
#            has no formal "kingdom" rank for bacteria.
#   Step 2 — Load UniRef90_ID → TaxID (from script 10)
#            Build lookup: UniRef90_ID → domain, phylum, family
#   Step 3 — Parse DIAMOND UniRef90 output per sample, map reads to family,
#            count reads per family per sample, build count matrices
#            separately for Bacteria, Fungi, and Archaea
#
# Filtering:
#   - TaxIDs with incomplete lineage (NF != 9) are excluded
#   - Entries with missing family ("__" or empty) are excluded
#   - Entries with no recognizable domain in lineage are excluded
#   - Bacteria: lineage contains ";Bacteria;"
#   - Fungi: phylum in {Ascomycota, Basidiomycota, Chytridiomycota,
#     Mucoromycota, Glomeromycota, Mortierellomycota}
#   - Archaea: lineage contains ";Archaea;"
#
# Inputs:
#   TAXID_FILE   : UniRef90_ID → TaxID (uniref90_hits_taxid.tsv)
#   LINEAGE_FILE : TaxID → lineage columns (taxid_lineage.tsv)
#   DIAMOND_DIR  : per-sample DIAMOND blastx outputs (*.uniref90.tsv.gz)
#   META_FILE    : sample metadata (sample_metadata.tsv)
#
# Outputs:
#   counts_family_bacteria.tsv : family x samples count matrix (bacteria)
#   counts_family_fungi.tsv    : family x samples count matrix (fungi)
#   counts_family_archaea.tsv  : family x samples count matrix (archaea)
#   taxonomy_parsing_stats.tsv : per-sample annotation stats
#
# Usage:
#   conda activate metatx
#   python 12a_parse_taxonomy_counts.py
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
TAXID_FILE   = os.path.join(BASE_DIR, "db/uniref90/uniref90_hits_taxid.tsv")
LINEAGE_FILE = os.path.join(BASE_DIR, "db/uniref90/taxid_lineage.tsv")
DIAMOND_DIR  = os.path.join(BASE_DIR, "results/diamond/uniref90")
META_FILE    = os.path.join(BASE_DIR, "sample_metadata.tsv")
OUT_DIR      = os.path.join(BASE_DIR, "results/counts")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Fungal phyla to include ───────────────────────────────────────────────────
FUNGAL_PHYLA = {
    "Ascomycota", "Basidiomycota", "Chytridiomycota",
    "Mucoromycota", "Glomeromycota", "Mortierellomycota"
}

# ── Step 1: Load TaxID → lineage ─────────────────────────────────────────────
# NOTE: taxonkit reformat {k} returns "unclassified X superkingdom" for all
# entries — NCBI taxonomy has no formal "kingdom" rank for bacteria.
# Domain (Bacteria/Eukaryota/Archaea) is parsed from the lineage string
# (col 2) instead. Columns 5-9 contain phylum through genus.
print("Step 1: Loading taxonkit lineage file...")
taxid_to_info = {}  # taxid -> (domain, phylum, family)
skipped = 0

with open(LINEAGE_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            skipped += 1
            continue
        taxid, lineage, rank, kingdom, phylum, clas, order, family, genus = parts
        # Skip entries with missing family
        if family in ("__", ""):
            skipped += 1
            continue
        # Determine domain from lineage string
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

print(f"  TaxIDs with complete lineage: {len(taxid_to_info):,}")
print(f"  Skipped (incomplete/unclassified): {skipped:,}")

# ── Step 2: Load UniRef90_ID → TaxID, build UniRef90_ID → family ─────────────
print("\nStep 2: Building UniRef90_ID → family lookup...")
uniref_to_domain = {}
uniref_to_phylum = {}
uniref_to_family = {}

with open(TAXID_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        uniref_id, taxid = parts
        if taxid in taxid_to_info:
            domain, phylum, family = taxid_to_info[taxid]
            uniref_to_domain[uniref_id] = domain
            uniref_to_phylum[uniref_id] = phylum
            uniref_to_family[uniref_id] = family

print(f"  UniRef90 IDs with family annotation: {len(uniref_to_family):,}")

# ── Step 3: Parse DIAMOND outputs and count by family ────────────────────────
print("\nStep 3: Parsing DIAMOND outputs...")
metadata = pd.read_csv(META_FILE, sep="\t")

bact_counts = defaultdict(lambda: defaultdict(int))
fungi_counts = defaultdict(lambda: defaultdict(int))
arch_counts  = defaultdict(lambda: defaultdict(int))
stats = []

for _, row in metadata.iterrows():
    sample_id    = row["sample_id"]
    file_prefix  = row["file_prefix"]
    diamond_file = os.path.join(DIAMOND_DIR, f"{file_prefix}.uniref90.tsv.gz")

    if not os.path.exists(diamond_file):
        print(f"  WARNING: {diamond_file} not found, skipping.")
        continue

    total_hits = 0
    bact_hits  = 0
    fungi_hits = 0
    arch_hits  = 0

    with gzip.open(diamond_file, "rt") as fh:
        for line in fh:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sseqid = parts[1]
            total_hits += 1

            family = uniref_to_family.get(sseqid)
            domain = uniref_to_domain.get(sseqid, "")
            phylum = uniref_to_phylum.get(sseqid, "")

            if family is None:
                continue

            # Bacteria
            if domain == "Bacteria":
                bact_counts[sample_id][family] += 1
                bact_hits += 1
            # Fungi
            elif domain == "Eukaryota" and any(fp in phylum for fp in FUNGAL_PHYLA):
                fungi_counts[sample_id][family] += 1
                fungi_hits += 1
            # Archaea
            elif domain == "Archaea":
                arch_counts[sample_id][family] += 1
                arch_hits += 1

    stats.append({
        "sample_id"    : sample_id,
        "treatment"    : row["treatment"],
        "timepoint"    : row["timepoint"],
        "total_hits"   : total_hits,
        "bacteria_hits": bact_hits,
        "fungi_hits"   : fungi_hits,
        "archaea_hits" : arch_hits,
        "pct_bacteria" : round(bact_hits / total_hits * 100, 2) if total_hits > 0 else 0,
        "pct_fungi"    : round(fungi_hits / total_hits * 100, 2) if total_hits > 0 else 0,
        "pct_archaea"  : round(arch_hits / total_hits * 100, 2) if total_hits > 0 else 0,
    })
    print(f"  {sample_id}: {total_hits:,} hits | bacteria: {bact_hits:,} ({stats[-1]['pct_bacteria']}%) | fungi: {fungi_hits:,} ({stats[-1]['pct_fungi']}%) | archaea: {arch_hits:,} ({stats[-1]['pct_archaea']}%)")

# ── Build count matrices ──────────────────────────────────────────────────────
print("\nBuilding count matrices...")

sample_ids = [row["sample_id"] for _, row in metadata.iterrows()
              if os.path.exists(os.path.join(DIAMOND_DIR, f"{row['file_prefix']}.uniref90.tsv.gz"))]

# Bacteria matrix
all_bact_families = sorted(set(f for s in sample_ids for f in bact_counts[s].keys()))
bact_matrix = pd.DataFrame(
    {s: [bact_counts[s].get(f, 0) for f in all_bact_families] for s in sample_ids},
    index=all_bact_families
)
bact_matrix.index.name = "family"

# Fungi matrix
all_fungi_families = sorted(set(f for s in sample_ids for f in fungi_counts[s].keys()))
fungi_matrix = pd.DataFrame(
    {s: [fungi_counts[s].get(f, 0) for f in all_fungi_families] for s in sample_ids},
    index=all_fungi_families
)
fungi_matrix.index.name = "family"

# Archaea matrix
all_arch_families = sorted(set(f for s in sample_ids for f in arch_counts[s].keys()))
arch_matrix = pd.DataFrame(
    {s: [arch_counts[s].get(f, 0) for f in all_arch_families] for s in sample_ids},
    index=all_arch_families
)
arch_matrix.index.name = "family"

# Stats table
stats_df = pd.DataFrame(stats)

# ── Save outputs ──────────────────────────────────────────────────────────────
bact_out  = os.path.join(OUT_DIR, "counts_family_bacteria.tsv")
fungi_out = os.path.join(OUT_DIR, "counts_family_fungi.tsv")
arch_out  = os.path.join(OUT_DIR, "counts_family_archaea.tsv")
stats_out = os.path.join(OUT_DIR, "taxonomy_parsing_stats.tsv")

bact_matrix.to_csv(bact_out,  sep="\t")
fungi_matrix.to_csv(fungi_out, sep="\t")
arch_matrix.to_csv(arch_out,  sep="\t")
stats_df.to_csv(stats_out,    sep="\t", index=False)

print(f"\nOutputs saved:")
print(f"  Bacteria matrix : {bact_out}  ({len(all_bact_families)} families x {len(sample_ids)} samples)")
print(f"  Fungi matrix    : {fungi_out}  ({len(all_fungi_families)} families x {len(sample_ids)} samples)")
print(f"  Archaea matrix  : {arch_out}  ({len(all_arch_families)} families x {len(sample_ids)} samples)")
print(f"  Stats           : {stats_out}")
print("\nDone.")
