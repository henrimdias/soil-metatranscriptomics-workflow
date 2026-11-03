#!/bin/bash
# =============================================================================
# Script:
# SLURM wrapper to run 13b_parse_ko_taxonomy_full.py -- builds the full
# (unfiltered) KO x family breakdown and the KO -> KEGG_Pathway lookup used
# for family-resolved analyses. Same resource profile as
# run_13_parse_ko_counts.sh, since it re-reads the same emapper/DIAMOND
# inputs (this is a variant of script 13, not script 07).
#
# Usage:
#   sbatch run_13b_parse_ko_taxonomy_full.sh
# =============================================================================
#SBATCH --job-name=ko_taxonomy_full
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=logs/ko_taxonomy_full_%j.out
#SBATCH --error=logs/ko_taxonomy_full_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "=========================================="

python ${PROJ}/scripts/13b_parse_ko_taxonomy_full.py

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
