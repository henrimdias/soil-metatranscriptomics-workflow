#!/bin/bash
#SBATCH --job-name=qc_trim
#SBATCH --partition=compute         # <-- e.g., short, normal, highmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu
#SBATCH --output=logs/qc_trim_%A_%a.out
#SBATCH --error=logs/qc_trim_%A_%a.err
#SBATCH --array=1-18

# ===============================================================
# Step 0: Environment setup
# ===============================================================
# Prevent /etc/bashrc issue
export BASHRCSOURCED=1
set -euo pipefail
PROJECT=${PROJ}
cd $PROJECT

source ~/.bashrc
conda activate metatx

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" samples.tsv)
sample=$(echo "$LINE" | awk '{print $1}')
R1=$(echo "$LINE" | awk '{print $2}')
R2=$(echo "$LINE" | awk '{print $3}')

echo "[$(date)] Starting QC/trim for $sample"

mkdir -p data/trim results/qc

# ===============================================================
# Step 1: FastQC before trimming
# ===============================================================
fastqc -t 8 -o results/qc "$R1" "$R2"

# ===============================================================
# Step 2: Trimming (adapters + low-quality reads)
# ===============================================================
fastp \
  -w 8 \
  -i "$R1" -I "$R2" \
  -o "data/trim/${sample}_R1.trim.fq.gz" \
  -O "data/trim/${sample}_R2.trim.fq.gz" \
  -h "results/qc/${sample}_fastp.html" \
  -j "results/qc/${sample}_fastp.json"

# ===============================================================
# Step 3: FastQC after trimming
# ===============================================================
fastqc -t 8 -o results/qc \
  "data/trim/${sample}_R1.trim.fq.gz" \
  "data/trim/${sample}_R2.trim.fq.gz"

echo "[$(date)] Finished QC/trim for $sample"

# ===============================================================
# Step 4: Run MultiQC automatically after the last job finishes
# ===============================================================
if [ "$SLURM_ARRAY_TASK_ID" -eq "$(($(wc -l < samples.tsv)-1))" ]; then
  echo "[$(date)] Running MultiQC summary..."
  multiqc results/qc -o results/qc
  echo "[$(date)] MultiQC report ready → results/qc/multiqc_report.html"
fi

echo "[$(date)] Done."
