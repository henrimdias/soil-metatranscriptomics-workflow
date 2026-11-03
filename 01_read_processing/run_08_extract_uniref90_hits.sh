#!/bin/bash
# =============================================================================
# Script:
# Extract UniRef90 sequences corresponding to DIAMOND blastx hits across
# all samples using seqtk, to generate a non-redundant FASTA input for
# eggNOG-mapper functional annotation.
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
# DIAMOND blastx against UniRef90 produces per-sample hit tables with
# sseqid (UniRef90 cluster IDs) in column 2. This script:
#   Step 1 — Aggregates all unique UniRef90 sseqids across 18 samples
#            into a single non-redundant ID list.
#   Step 2 — Extracts the corresponding protein sequences from the
#            UniRef90 FASTA using seqtk subseq in a single pass.
#
# The resulting FASTA (~3.4M sequences) is used as input for eggNOG-mapper
# (emapper.py -m diamond) to obtain KEGG Orthology (KO) annotations,
# which are then joined back to the per-sample read counts.
#
# Assumptions:
#   - Conda environment `metatx` contains seqtk v1.5+
#   - DIAMOND UniRef90 outputs exist at DIAMOND_DIR/*.uniref90.tsv.gz
#   - UniRef90 FASTA exists at UNIREF_FASTA
#   - Sufficient scratch space (~20GB) for the output FASTA
#
# Inputs:
#   DIAMOND_DIR : per-sample DIAMOND blastx outputs (*.uniref90.tsv.gz)
#   UNIREF_FASTA: UniRef90 protein sequences (uniref90.fasta.gz)
#
# Outputs:
#   HITS_IDS    : non-redundant list of UniRef90 IDs with hits
#   HITS_FASTA  : extracted protein sequences for eggNOG-mapper input
#
# Usage:
#   sbatch 08_extract_uniref90_hits.sh
#
# Notes:
#   - seqtk subseq performs a single pass through the FASTA (45GB gz)
#   - Runtime estimate: 1-3h depending on I/O speed
#   - TMPDIR set to scratch to avoid /tmp space issues during sort
# =============================================================================
#SBATCH --job-name=seqtk_extract
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/seqtk_extract_%j.out
#SBATCH --error=logs/seqtk_extract_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR=${PROJ}
DIAMOND_DIR=${BASE_DIR}/results/diamond/uniref90
UNIREF_FASTA=${BASE_DIR}/db/uniref90/uniref90.fasta.gz
HITS_IDS=${BASE_DIR}/db/uniref90/uniref90_hits_unique.txt
HITS_FASTA=${BASE_DIR}/db/uniref90/uniref90_hits.fasta
TMPDIR=${SCRATCH_DIR}/tmp
mkdir -p ${TMPDIR}

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "=========================================="

# ── Step 1: Aggregate unique UniRef90 IDs across all samples ─────────────────
echo "[$(date)] Step 1: Collecting unique UniRef90 hit IDs..."

for f in ${DIAMOND_DIR}/*.uniref90.tsv.gz; do
    zcat $f | awk '{print $2}'
done | sort -T ${TMPDIR} -u > ${HITS_IDS}

echo "  Unique UniRef90 IDs: $(wc -l < ${HITS_IDS})"

# ── Step 2: Extract sequences using seqtk ────────────────────────────────────
echo "[$(date)] Step 2: Extracting sequences with seqtk..."

seqtk subseq \
    ${UNIREF_FASTA} \
    ${HITS_IDS} \
    > ${HITS_FASTA}

N_SEQS=$(grep -c '>' ${HITS_FASTA})
echo "  Sequences extracted: ${N_SEQS}"
echo "  Output size: $(du -sh ${HITS_FASTA} | cut -f1)"

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
