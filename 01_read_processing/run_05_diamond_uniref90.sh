#!/bin/bash
# =============================================================================
# Script:
# DIAMOND blastx translated alignment of metatranscriptomic mRNA reads
# against UniRef90 for broad functional annotation and KEGG Orthology
# mapping via eggNOG-mapper or UniProt ID mapping.
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
# This SLURM array job runs DIAMOND blastx on 18 paired mRNA read files
# (forward reads only) against the UniRef90 DIAMOND index (January 2026
# release; 188,848,220 sequences). Hits are used for broad functional
# annotation and subsequent mapping to KEGG Orthology (KO) terms via
# the UniProt ID mapping file. Only the best hit per query read is retained.
#
# The script performs:
#   - DIAMOND blastx of mRNA fwd reads against uniref90.dmnd
#   - E-value threshold: 1e-5
#   - Sequence identity threshold: 60%
#   - Top 1 hit per read (--max-target-seqs 1)
#   - Extended tabular output (outfmt 6 + extra fields)
#   - Compressed output to save disk space
#
# Assumptions:
#   - Conda environment `metatx` contains DIAMOND v2.1.13
#   - mRNA reads have been quality-trimmed (fastp) and rRNA-depleted (SortMeRNA)
#   - UniRef90 DIAMOND index exists at DB_UNIREF path
#   - Only fwd (R1) reads are used as query
#
# Inputs:
#   MRNA_DIR  : rRNA-depleted mRNA reads (*_mRNA_fwd.fq.gz)
#   DB_UNIREF : UniRef90 DIAMOND binary index (.dmnd)
#
# Outputs:
#   OUT_DIR   : per-sample tabular DIAMOND output (.uniref90.tsv.gz)
#
# Usage:
#   sbatch 04_diamond_uniref90.sh
#
# Array:
#   one task per sample; task N reads row N of sample_metadata.tsv
#
# Notes:
#   - UniRef90 index is 86GB; DIAMOND loads it in chunks (--block-size)
#   - --block-size 4 uses ~64GB RAM per task; adjust if OOM occurs
#   - Runtime estimate: 4-12h per sample depending on read count
# =============================================================================
#SBATCH --job-name=diamond_uniref90
#SBATCH --partition=compute
#SBATCH --array=1-18
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/diamond_uniref90_%A_%a.out
#SBATCH --error=logs/diamond_uniref90_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

# ── Paths ─────────────────────────────────────────────────────────────────────
MRNA_DIR=${PROJ}/data/mRNA
DB_UNIREF=${PROJ}/db/uniref90/uniref90.dmnd
OUT_DIR=${PROJ}/results/diamond/uniref90
mkdir -p ${OUT_DIR}

# ── Sample for this array task: row N of sample_metadata.tsv, column file_prefix ─
#    (one array task per sample; set --array=1-<number of samples>)
SAMPLE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $2}' ${PROJ}/sample_metadata.tsv)
INPUT=${MRNA_DIR}/${SAMPLE}_mRNA_fwd.fq.gz
OUTPUT=${OUT_DIR}/${SAMPLE}.uniref90.tsv

echo "=========================================="
echo "Task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample  : $SAMPLE"
echo "Input   : $INPUT"
echo "Started : $(date)"
echo "=========================================="

# ── Run DIAMOND blastx ────────────────────────────────────────────────────────
diamond blastx \
    --query        ${INPUT} \
    --db           ${DB_UNIREF} \
    --out          ${OUTPUT} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --evalue       1e-5 \
    --id           80 \
    --max-target-seqs 1 \
    --block-size   4 \
    --threads      ${SLURM_CPUS_PER_TASK} \
    --compress     1 \
    --quiet

echo "=========================================="
echo "Finished: $(date)"
echo "Hits    : $(zcat ${OUTPUT}.gz | wc -l)"
echo "=========================================="
