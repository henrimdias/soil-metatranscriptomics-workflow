#!/bin/bash
#SBATCH --job-name=kofam_dl
#SBATCH --partition=quickq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/kofam_dl_%j.out
#SBATCH --error=logs/kofam_dl_%j.err
set -euo pipefail

KOFAM=${SCRATCH_DIR}/dbs/kofam
mkdir -p "$KOFAM"; cd "$KOFAM"

echo "[$(date)] ko_list"
[ -s ko_list ] || { wget -q https://www.genome.jp/ftp/db/kofam/ko_list.gz && gunzip -f ko_list.gz; }

echo "[$(date)] profiles"
[ -d profiles ] || { wget -q https://www.genome.jp/ftp/db/kofam/profiles.tar.gz \
                     && tar xzf profiles.tar.gz && rm -f profiles.tar.gz; }

echo "[$(date)] pronto"
echo "ko_list  : $(wc -l < ko_list) linhas"
echo "profiles : $(ls profiles | wc -l) arquivos"
du -sh "$KOFAM"
