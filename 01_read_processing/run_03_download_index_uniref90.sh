#!/bin/bash
# =============================================================================
# Script:
# Download and index UniRef90 database for DIAMOND blastx annotation of
# metatranscriptomic reads.
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
# Downloads the UniRef90 FASTA file from the UniProt FTP server and builds
# a DIAMOND binary database index (.dmnd) for use in translated alignment
# of metatranscriptomic mRNA reads. UniRef90 clusters UniProtKB sequences
# at 90% sequence identity, providing broad functional coverage while
# reducing redundancy relative to the full UniProtKB database.
#
# The script performs:
#   - Download of uniref90.fasta.gz from UniProt FTP (January 2026 release)
#   - DIAMOND makedb indexing of the downloaded FASTA
#   - Log of database statistics to makedb.log
#
# Assumptions:
#   - Conda environment `metatx` contains DIAMOND v2.1.13
#   - Compute nodes have outbound internet access to ftp.uniprot.org
#   - Sufficient scratch space (~150GB) for FASTA + DMND index
#
# Inputs:
#   UniProt FTP: https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/
#
# Outputs:
#   DB_DIR/uniref90.fasta.gz : downloaded FASTA (compressed)
#   DB_DIR/uniref90.dmnd     : DIAMOND binary index
#   DB_DIR/makedb_uniref90.log : indexing log
#
# Usage:
#   sbatch run_03_download_index_uniref90.sh
#
# Notes:
#   - Download alone may take 2-4h depending on network speed
#   - Indexing requires ~60-80GB RAM and ~1-2h with 16 threads
#   - Total disk usage: ~150GB (compressed FASTA ~50GB + DMND ~100GB)
# =============================================================================
#SBATCH --job-name=uniref90_db
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/uniref90_db_%j.out
#SBATCH --error=logs/uniref90_db_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# ── Activate conda environment ────────────────────────────────────────────────
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate metatx

# ── Paths ─────────────────────────────────────────────────────────────────────
DB_DIR=${PROJ}/db/uniref90
mkdir -p ${DB_DIR}

FASTA_URL=https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
FASTA_OUT=${DB_DIR}/uniref90.fasta.gz
DMND_OUT=${DB_DIR}/uniref90
LOG=${DB_DIR}/makedb_uniref90.log

echo "=========================================="
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo "=========================================="

# ── Step 1: Download UniRef90 ─────────────────────────────────────────────────
echo "[$(date)] Downloading UniRef90 from UniProt FTP..."
wget --no-verbose \
     --show-progress \
     --continue \
     -O ${FASTA_OUT} \
     ${FASTA_URL}

echo "[$(date)] Download complete. File size: $(du -sh ${FASTA_OUT} | cut -f1)"

# ── Step 2: Index with DIAMOND ────────────────────────────────────────────────
echo "[$(date)] Building DIAMOND database index..."
diamond makedb \
    --in ${FASTA_OUT} \
    --db ${DMND_OUT} \
    --threads ${SLURM_CPUS_PER_TASK} \
    2>&1 | tee ${LOG}

echo "=========================================="
echo "Finished: $(date)"
echo "DMND index: $(ls -lh ${DMND_OUT}.dmnd)"
echo "=========================================="
