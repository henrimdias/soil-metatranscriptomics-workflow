#!/bin/bash
# =============================================================================
# Script:
# SLURM wrapper to run 13_parse_ko_counts.py — parses eggNOG-mapper
# annotations to build KEGG Orthology (KO) count matrices for DESeq2
# and taxonomic breakdown of carbon metabolism pathways.
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
#   - eggNOG-mapper annotations: UniRef90_ID → KO terms
#   - taxonkit lineage: UniRef90_ID → TaxID → family
#
# Produces:
#   counts_ko.tsv          : full KO × samples matrix (all KOs, DESeq2 input)
#   counts_ko_carbon.tsv   : KO × samples (carbon metabolism KOs only)
#   counts_ko_taxonomy.tsv : long-format KO × family × sample (as in Zhou et al. 2025)
#   ko_parsing_stats.tsv   : per-sample annotation stats
#
# Assumptions:
#   - Conda environment `metatx` contains Python 3.x and pandas
#   - eggNOG-mapper output exists at results/emapper/
#   - taxonkit lineage file exists at db/uniref90/taxid_lineage.tsv
#   - UniRef90 TaxID mapping exists at db/uniref90/uniref90_hits_taxid.tsv
#   - DIAMOND UniRef90 outputs exist at results/diamond/uniref90/
#
# Usage:
#   sbatch run_13_parse_ko_counts.sh
#
# Notes:
#   - eggNOG-mapper annotations file is 1.6GB — loaded once into memory
#   - Estimated RAM: 32-64GB
#   - Estimated runtime: 1-2h
# =============================================================================
#SBATCH --job-name=ko_counts
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/ko_counts_%j.out
#SBATCH --error=logs/ko_counts_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "=========================================="

python ${PROJ}/scripts/13_parse_ko_counts.py

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
