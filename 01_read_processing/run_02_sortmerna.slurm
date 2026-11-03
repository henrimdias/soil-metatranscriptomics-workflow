#!/bin/bash
# =============================================================================
# Script:
# rRNA depletion of paired-end metatranscriptomic reads using SortMeRNA v4,
# separating rRNA and mRNA fractions for downstream functional analysis.
#
# Author (script):
# Henrique M. Dias
#
# Affiliation:
# South Dakota State University
#
# Date:
# 2025
#
# Description:
# This SLURM array job runs SortMeRNA on paired-end metatranscriptomic
# libraries, one array task per sample. Reads are separated into rRNA (retained for
# taxonomic validation) and mRNA (used for functional analysis via DIAMOND
# and taxon-aware mapping).
#
# The script performs:
#   - Activation of the conda environment `metatx`
#   - rRNA filtering against the smr_v4.3_default_db (SILVA-based, bacteria,
#     archaea, and eukarya rRNA)
#   - Output of rRNA reads (aligned) and mRNA reads (other) as gzipped FASTQ
#     pairs
#   - Cleanup of per-sample SortMeRNA work directories (kvdb) after each run
#
# Assumptions:
#   - Conda environment `metatx` contains SortMeRNA v4.3.7
#   - Trimmed reads follow the naming convention:
#       {PREFIX}_R1.trim.fq.gz / {PREFIX}_R2.trim.fq.gz
#   - The SortMeRNA database has been downloaded from:
#       https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
#   - Output directories (mRNA/, rRNA/) exist prior to submission
#
# Inputs:
#   TRIM_DIR : trimmed paired-end reads (.trim.fq.gz)
#   DB       : SortMeRNA default database (smr_v4.3_default_db.fasta)
#
# Outputs:
#   MRNA_DIR : mRNA-enriched reads (per sample, R1 + R2)
#   RRNA_DIR : rRNA reads (per sample, R1 + R2; retained for QC and taxonomy)
#
# Usage:
#   sbatch run_02_sortmerna.slurm
#
# Array:
#   one task per sample; task N reads row N of sample_metadata.tsv
# =============================================================================
#SBATCH --job-name=sortmerna
#SBATCH --partition=compute
#SBATCH --array=1-18
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sortmerna_%A_%a.out
#SBATCH --error=logs/sortmerna_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ~/miniconda3/etc/profile.d/conda.sh
conda activate metatx

# ── Paths ─────────────────────────────────────────────────────────────────────
TRIM_DIR=${PROJ}/data/trim
MRNA_DIR=${PROJ}/data/mRNA
RRNA_DIR=${PROJ}/data/rRNA
DB=${SCRATCH_DIR}/dbs/sortmerna/smr_v4.3_default_db.fasta
WORKDIR=${PROJ}/work

# ── Sample for this array task: row N of sample_metadata.tsv, column file_prefix ─
#    (one array task per sample; set --array=1-<number of samples>)
SAMPLE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $2}' ${PROJ}/sample_metadata.tsv)

echo "=========================================="
echo "Task ID : $SLURM_ARRAY_TASK_ID"
echo "Sample  : $SAMPLE"
echo "Started : $(date)"
echo "=========================================="

# ── SortMeRNA requires a clean workdir per sample ─────────────────────────────
# (kvdb locks and refuses to overwrite an existing run directory)
SAMPLE_WORKDIR=${WORKDIR}/sortmerna_${SAMPLE}
rm -rf ${SAMPLE_WORKDIR}
mkdir -p ${SAMPLE_WORKDIR}

# ── Run SortMeRNA ─────────────────────────────────────────────────────────────
sortmerna \
    --ref     ${DB} \
    --reads   ${TRIM_DIR}/${SAMPLE}_R1.trim.fq.gz \
    --reads   ${TRIM_DIR}/${SAMPLE}_R2.trim.fq.gz \
    --aligned ${RRNA_DIR}/${SAMPLE}_rRNA \
    --other   ${MRNA_DIR}/${SAMPLE}_mRNA \
    --fastx \
    --paired-out \
    --out2 \
    --threads ${SLURM_CPUS_PER_TASK} \
    --workdir ${SAMPLE_WORKDIR}

# ── Cleanup workdir (kvdb can be several GB per sample) ───────────────────────
rm -rf ${SAMPLE_WORKDIR}

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
