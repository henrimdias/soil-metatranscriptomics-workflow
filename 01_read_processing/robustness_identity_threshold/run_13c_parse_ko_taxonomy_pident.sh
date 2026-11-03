#!/bin/bash
# =============================================================================
# Script:
# SLURM wrapper to run 13c_parse_ko_taxonomy_pident.py -- rebuilds the
# KO x family x sample table at DIAMOND identity >= 80, 90 and 95 % (taxonomic-assignment robustness check).
# Three thresholds are held in memory at once, hence the larger allocation.
# inputs (this is a variant of script 13, not script 07).
#
# Usage:
#   sbatch run_13c_parse_ko_taxonomy_pident.sh
# =============================================================================
#SBATCH --job-name=ko_taxonomy_pident
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ko_taxonomy_pident_%j.out
#SBATCH --error=logs/ko_taxonomy_pident_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "=========================================="

python ${PROJ}/scripts/13c_parse_ko_taxonomy_pident.py

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
