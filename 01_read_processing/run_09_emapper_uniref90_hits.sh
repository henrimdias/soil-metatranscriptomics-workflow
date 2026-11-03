#!/bin/bash
# =============================================================================
# Script:
# Run eggNOG-mapper v2.1.13 on UniRef90 hit sequences to obtain KEGG
# Orthology (KO) annotations for downstream DESeq2 analysis.
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
# Runs eggNOG-mapper in diamond mode on the non-redundant set of UniRef90
# protein sequences that had blastx hits across the 18 metatranscriptomic
# samples. The resulting annotations (KO terms, COG categories, GO terms)
# are used to map per-sample read counts to KEGG Orthology groups for
# differential expression analysis with DESeq2.
#
# Pipeline context:
#   DIAMOND blastx → UniRef90 hits → seqtk extract → eggNOG-mapper → KO terms
#   → join with per-sample DIAMOND counts → DESeq2
#
# Assumptions:
#   - Conda environment `eggnog2` contains eggNOG-mapper v2.1.13
#   - eggNOG database v5.0.2 exists at EGGNOG_DB_DIR
#   - Input FASTA (3.4M sequences, 1.8GB) exists at INPUT_FASTA
#   - Sufficient scratch space (~5GB) for output files
#
# Inputs:
#   INPUT_FASTA  : UniRef90 protein sequences with hits (uniref90_hits.fasta)
#   EGGNOG_DB_DIR: eggNOG database directory (eggnog.db, eggnog_proteins.dmnd)
#
# Outputs:
#   OUT_DIR/uniref90_hits.emapper.annotations : main annotation file
#   OUT_DIR/uniref90_hits.emapper.hits         : DIAMOND hits
#   OUT_DIR/uniref90_hits.emapper.seed_orthologs : seed orthologs
#
# Key output columns in .emapper.annotations:
#   Col 1  : query sequence ID (UniRef90 ID)
#   Col 12 : KEGG Orthology (KO) terms
#   Col 7  : GO terms
#   Col 5  : COG category
#
# Usage:
#   sbatch 09_emapper_uniref90_hits.sh
#
# Notes:
#   - -m diamond uses DIAMOND for the search step (faster than HMMER)
#   - --itype proteins: input is protein sequences
#   - --cpu 16: uses all allocated CPUs
#     (requires ~32GB RAM; remove if OOM)
#   - Runtime estimate: 24-48h for 3.4M sequences
# =============================================================================
#SBATCH --job-name=emapper
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/emapper_%j.out
#SBATCH --error=logs/emapper_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate eggnog2
export TMPDIR=${SCRATCH_DIR}/tmp
export HOME=${SCRATCH_DIR}/tmp

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR=${PROJ}
INPUT_FASTA=${BASE_DIR}/db/uniref90/uniref90_hits.fasta
EGGNOG_DB_DIR=${SCRATCH_DIR}/dbs/eggnog
OUT_DIR=${BASE_DIR}/results/emapper
mkdir -p ${OUT_DIR}

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "Input   : ${INPUT_FASTA}"
echo "Seqs    : $(grep -c '>' ${INPUT_FASTA})"
echo "=========================================="

# ── Run eggNOG-mapper ─────────────────────────────────────────────────────────
emapper.py \
    -i ${INPUT_FASTA} \
    --itype proteins \
    -m diamond \
    --data_dir ${EGGNOG_DB_DIR} \
    --output uniref90_hits \
    --output_dir ${OUT_DIR} \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --override \
    --scratch_dir ${SCRATCH_DIR}/tmp \
    --no_file_comments

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
