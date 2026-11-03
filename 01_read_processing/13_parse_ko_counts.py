#!/usr/bin/env python3
# =============================================================================
# Script:
# Parse eggNOG-mapper annotations to build KEGG Orthology (KO) count
# matrices for DESeq2 differential expression analysis and taxonomic
# breakdown of carbon metabolism pathways.
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
# Integrates three data sources via UniRef90_ID as common key:
#   - DIAMOND blastx outputs: read → UniRef90_ID (per sample)
#   - eggNOG-mapper annotations: UniRef90_ID → KO terms (col 12)
#   - taxonkit lineage: UniRef90_ID → TaxID → family
#
# Produces three output matrices:
#   1. counts_ko.tsv: full KO × samples matrix (all annotated KOs)
#   2. counts_ko_carbon.tsv: KO × samples for carbon metabolism KOs only
#   3. counts_ko_taxonomy.tsv: KO × family × sample long-format table
#      for taxonomic breakdown figures (as in Zhou et al. 2025)
#
# Carbon metabolism KO blocks (curated):
#   Block 1 - Glycolysis/Gluconeogenesis (map00010)
#   Block 2 - Pentose phosphate pathway (map00030)
#   Block 3 - TCA cycle (map00020)
#   Block 4 - Oxidative phosphorylation (map00190)
#   Block 5 - Microbial growth proxies (aminoacyl-tRNA, purine metabolism)
#
# KO handling for multi-KO annotations:
#   Sequences annotated with multiple KO terms (comma-separated) are
#   counted once per KO — i.e., the read count is distributed across
#   all annotated KOs. This is the most inclusive approach and is
#   appropriate for pathway-level analysis.
#
# Inputs:
#   EMAPPER_FILE : eggNOG-mapper annotations (uniref90_hits.emapper.annotations)
#   TAXID_FILE   : UniRef90_ID → TaxID (uniref90_hits_taxid.tsv)
#   LINEAGE_FILE : TaxID → lineage (taxid_lineage.tsv)
#   DIAMOND_DIR  : per-sample DIAMOND blastx outputs (*.uniref90.tsv.gz)
#   META_FILE    : sample metadata (sample_metadata.tsv)
#
# Outputs:
#   counts_ko.tsv          : KO × samples count matrix (DESeq2 input)
#   counts_ko_carbon.tsv   : KO × samples (carbon metabolism only)
#   counts_ko_taxonomy.tsv : long-format KO × family × sample counts
#   ko_parsing_stats.tsv   : per-sample annotation stats
#
# Usage:
#   Submit via SLURM (see wrapper script) — emapper file is 1.6GB
#   Estimated RAM: 32GB, Runtime: 1-2h
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

# ── Carbon metabolism KOs (curated) ──────────────────────────────────────────
# Block 1: Glycolysis / Gluconeogenesis (map00010)
GLYCOLYSIS_KOS = {
    "K00844", "K12407", "K00845", "K01810", "K06859", "K13810", "K15916",
    "K00850", "K16370", "K00927", "K11389", "K01834", "K15633", "K15634",
    "K15635", "K00873", "K12406", "K00016", "K00101"
}
# Block 2: Pentose phosphate pathway (map00030)
PPP_KOS = {
    "K00036", "K01057", "K07404", "K00033", "K01783", "K02082", "K01808",
    "K01807", "K00615", "K00616", "K01100", "K03841", "K04041"
}
# Block 3: TCA cycle (map00020)
TCA_KOS = {
    "K01647", "K01681", "K01682", "K00031", "K00030", "K01903", "K01902",
    "K00234", "K00235", "K00236", "K00237", "K01676", "K01677", "K01678",
    "K00239", "K00240", "K00241", "K00242", "K01902", "K00174", "K00175",
    "K00176", "K00177", "K01899", "K01900", "K00382", "K00627", "K00169",
    "K00170", "K00171", "K00172"
}
# Block 4: Oxidative phosphorylation (map00190)
OXPHOS_KOS = {
    "K00330", "K00331", "K00332", "K00333", "K00334", "K00335", "K00336",
    "K00337", "K00338", "K00339", "K00340", "K00341", "K00342", "K00343",
    "K00344", "K00345", "K02256", "K02257", "K02258", "K02259", "K02260",
    "K02261", "K02262", "K02263", "K02264", "K02265", "K02266", "K02267",
    "K02268", "K02274", "K02275", "K02276", "K02277", "K02278", "K02279",
    "K02110", "K02111", "K02112", "K02113", "K02114", "K02115", "K02116",
    "K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02123",
    "K02124", "K02125", "K02126", "K02127", "K02128", "K02129", "K02130",
    "K02131", "K02132", "K02133", "K02134", "K02135", "K02136"
}
# Block 5: Growth proxies
GROWTH_KOS = {
    # Aminoacyl-tRNA biosynthesis (map00970) — proxy for active growth
    "K01866", "K01867", "K01868", "K01869", "K01870", "K01872", "K01873",
    "K01874", "K01875", "K01876", "K01878", "K01879", "K01880", "K01881",
    "K01882", "K01883", "K01884", "K01885", "K01886", "K01887", "K01889",
    "K01890", "K01891", "K01892", "K01893", "K01894", "K01895", "K04566",
    # Purine metabolism (map00230)
    "K00627", "K01945", "K01939", "K00602", "K00764", "K01756", "K00759",
    "K00760", "K01196", "K00549"
}

CARBON_KOS = GLYCOLYSIS_KOS | PPP_KOS | TCA_KOS | OXPHOS_KOS | GROWTH_KOS

# KO block labels for output
KO_BLOCK = {}
for k in GLYCOLYSIS_KOS: KO_BLOCK[k] = "Glycolysis"
for k in PPP_KOS:        KO_BLOCK[k] = "Pentose_phosphate"
for k in TCA_KOS:        KO_BLOCK[k] = "TCA_cycle"
for k in OXPHOS_KOS:     KO_BLOCK[k] = "Oxidative_phosphorylation"
for k in GROWTH_KOS:     KO_BLOCK[k] = "Growth_proxies"

# ── Fungal phyla ──────────────────────────────────────────────────────────────
FUNGAL_PHYLA = {
    "Ascomycota", "Basidiomycota", "Chytridiomycota",
    "Mucoromycota", "Glomeromycota", "Mortierellomycota"
}

# ── Step 1: Load TaxID → lineage (domain + family) ───────────────────────────
print("Step 1: Loading taxonkit lineage...")
taxid_to_info = {}  # taxid -> (domain, phylum, family)
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
print("\nStep 2: Loading UniRef90_ID → TaxID mapping...")
uniref_to_taxinfo = {}  # uniref_id -> (domain, phylum, family)

with open(TAXID_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        uniref_id, taxid = parts
        if taxid in taxid_to_info:
            uniref_to_taxinfo[uniref_id] = taxid_to_info[taxid]

print(f"  UniRef90 IDs with taxonomy: {len(uniref_to_taxinfo):,}")

# ── Step 3: Load eggNOG-mapper annotations → UniRef90_ID → KOs ───────────────
print("\nStep 3: Loading eggNOG-mapper KO annotations...")
uniref_to_kos = {}  # uniref_id -> list of KO terms

with open(EMAPPER_FILE) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        uniref_id = parts[0]
        ko_field  = parts[11]
        if ko_field == "-" or ko_field == "":
            continue
        # Parse KO terms: "ko:K00844,ko:K01810" -> ["K00844", "K01810"]
        kos = [k.replace("ko:", "") for k in ko_field.split(",")]
        uniref_to_kos[uniref_id] = kos

print(f"  UniRef90 IDs with KO annotation: {len(uniref_to_kos):,}")

# ── Step 4: Parse DIAMOND outputs and build count matrices ───────────────────
print("\nStep 4: Parsing DIAMOND outputs and building count matrices...")
metadata = pd.read_csv(META_FILE, sep="\t")

# ko_counts[sample_id][ko] = count
ko_counts = defaultdict(lambda: defaultdict(int))
# ko_tax_counts[sample_id][ko][family] = count (for taxonomy breakdown)
ko_tax_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
stats = []

for _, row in metadata.iterrows():
    sample_id    = row["sample_id"]
    file_prefix  = row["file_prefix"]
    diamond_file = os.path.join(DIAMOND_DIR, f"{file_prefix}.uniref90.tsv.gz")

    if not os.path.exists(diamond_file):
        print(f"  WARNING: {diamond_file} not found, skipping.")
        continue

    total_hits  = 0
    ko_hits     = 0
    carbon_hits = 0

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
                ko_counts[sample_id][ko] += 1
                if ko in CARBON_KOS:
                    carbon_hits += 1
                    ko_tax_counts[sample_id][ko][family] += 1

    stats.append({
        "sample_id"   : sample_id,
        "treatment"   : row["treatment"],
        "timepoint"   : row["timepoint"],
        "total_hits"  : total_hits,
        "ko_hits"     : ko_hits,
        "carbon_hits" : carbon_hits,
        "pct_ko"      : round(ko_hits / total_hits * 100, 2) if total_hits > 0 else 0,
        "pct_carbon"  : round(carbon_hits / total_hits * 100, 2) if total_hits > 0 else 0,
    })
    print(f"  {sample_id}: {total_hits:,} hits | KO: {ko_hits:,} ({stats[-1]['pct_ko']}%) | carbon: {carbon_hits:,} ({stats[-1]['pct_carbon']}%)")

# ── Step 5: Build and save count matrices ─────────────────────────────────────
print("\nStep 5: Building count matrices...")

sample_ids = [row["sample_id"] for _, row in metadata.iterrows()
              if os.path.exists(os.path.join(DIAMOND_DIR, f"{row['file_prefix']}.uniref90.tsv.gz"))]

# Full KO matrix
all_kos = sorted(set(ko for s in sample_ids for ko in ko_counts[s].keys()))
ko_matrix = pd.DataFrame(
    {s: [ko_counts[s].get(k, 0) for k in all_kos] for s in sample_ids},
    index=all_kos
)
ko_matrix.index.name = "KO"

# Carbon KO matrix
carbon_kos = sorted(k for k in all_kos if k in CARBON_KOS)
carbon_matrix = ko_matrix.loc[carbon_kos] if carbon_kos else pd.DataFrame()

# Taxonomy breakdown (long format)
tax_rows = []
for sample_id in sample_ids:
    sample_row = metadata[metadata["sample_id"] == sample_id].iloc[0]
    for ko, family_counts in ko_tax_counts[sample_id].items():
        block = KO_BLOCK.get(ko, "Other")
        for family, count in family_counts.items():
            tax_rows.append({
                "sample_id" : sample_id,
                "treatment" : sample_row["treatment"],
                "timepoint" : sample_row["timepoint"],
                "KO"        : ko,
                "block"     : block,
                "family"    : family,
                "count"     : count
            })
tax_df = pd.DataFrame(tax_rows)

# Stats
stats_df = pd.DataFrame(stats)

# Save
ko_out     = os.path.join(OUT_DIR, "counts_ko.tsv")
carbon_out = os.path.join(OUT_DIR, "counts_ko_carbon.tsv")
tax_out    = os.path.join(OUT_DIR, "counts_ko_taxonomy.tsv")
stats_out  = os.path.join(OUT_DIR, "ko_parsing_stats.tsv")

ko_matrix.to_csv(ko_out,     sep="\t")
carbon_matrix.to_csv(carbon_out, sep="\t")
tax_df.to_csv(tax_out,       sep="\t", index=False)
stats_df.to_csv(stats_out,   sep="\t", index=False)

print(f"\nOutputs saved:")
print(f"  Full KO matrix    : {ko_out}  ({len(all_kos):,} KOs x {len(sample_ids)} samples)")
print(f"  Carbon KO matrix  : {carbon_out}  ({len(carbon_kos)} KOs x {len(sample_ids)} samples)")
print(f"  Taxonomy breakdown: {tax_out}  ({len(tax_df):,} rows)")
print(f"  Stats             : {stats_out}")
print("\nDone.")
