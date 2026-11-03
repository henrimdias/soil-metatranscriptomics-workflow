#!/bin/bash
# =============================================================================
# Script:
# Extract TaxIDs from UniRef90 FASTA headers for sequences with DIAMOND
# blastx hits, to enable taxonomic assignment of metatranscriptomic reads
# at the family level using taxonkit.
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
# UniRef90 FASTA headers contain the NCBI TaxID of the representative
# organism in the format "TaxID=XXXXXX". This script:
#   Step 1 — Scans the UniRef90 FASTA (45GB gz) for headers matching
#            the 3.4M unique sseqids that had DIAMOND blastx hits across
#            all 18 metatranscriptomic samples.
#   Step 2 — Extracts UniRef90_ID → TaxID pairs into a TSV file.
#
# The resulting mapping file is used downstream with taxonkit to obtain
# full taxonomic lineages (phylum → family) for each UniRef90 cluster,
# enabling taxonomic attribution of metatranscriptomic read counts.
#
# This enables construction of a scatter plot comparing:
#   X-axis: relative abundance of bacterial/fungal families in 16S/ITS amplicon
#   Y-axis: relative abundance of transcripts attributed to the same families
#   Color: treatment
#   Shape: timepoint (24h, 48h, 72h)
#
# Assumptions:
#   - UniRef90 FASTA headers contain TaxID in format "TaxID=XXXXXX"
#   - HITS_IDS file contains unique sseqids from DIAMOND output (one per line)
#   - No conda environment needed (uses system grep/sed/awk)
#
# Inputs:
#   UNIREF_FASTA : UniRef90 protein sequences (uniref90.fasta.gz)
#   HITS_IDS     : non-redundant list of UniRef90 sseqids with hits
#
# Outputs:
#   HITS_TAXID   : TSV file with columns: UniRef90_ID, TaxID
#
# Usage:
#   sbatch 10_extract_uniref90_taxids.sh
#
# Notes:
#   - grep -F uses fixed string matching (faster than regex for large files)
#   - Runtime estimate: 2-4h depending on I/O speed
#   - Output size: ~100MB for 3.4M entries
# =============================================================================
#SBATCH --job-name=extract_taxid
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/extract_taxid_%j.out
#SBATCH --error=logs/extract_taxid_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR=${PROJ}
UNIREF_FASTA=${BASE_DIR}/db/uniref90/uniref90.fasta.gz
HITS_IDS=${BASE_DIR}/db/uniref90/uniref90_hits_unique.txt
HITS_TAXID=${BASE_DIR}/db/uniref90/uniref90_hits_taxid.tsv

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "IDs to match: $(wc -l < ${HITS_IDS})"
echo "=========================================="

# ── Extract UniRef90_ID → TaxID from FASTA headers ───────────────────────────
echo "[$(date)] Scanning UniRef90 FASTA headers..."

zcat ${UNIREF_FASTA} | \
    grep '^>' | \
    grep -Ff ${HITS_IDS} | \
    sed 's/^>//; s/ .*TaxID=\([0-9]*\).*/\t\1/' \
    > ${HITS_TAXID}

N_MAPPED=$(wc -l < ${HITS_TAXID})
echo "[$(date)] Done: ${N_MAPPED} UniRef90 IDs mapped to TaxID"
echo "Output: ${HITS_TAXID}"

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
