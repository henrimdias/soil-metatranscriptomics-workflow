#!/bin/bash
#SBATCH --job-name=uniref90_taxdb
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --output=logs/uniref90_taxdb_%j.out
#SBATCH --error=logs/uniref90_taxdb_%j.err
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate metatx

DB=${PROJ}/db/uniref90
TAX=${SCRATCH_DIR}/dbs/taxdump
FASTA=$DB/uniref90.fasta.gz
MAP=$DB/uniref90.accession2taxid.gz

echo "[$(date)] passo 1: mapa accession->taxid a partir dos cabecalhos"
{ printf 'accession\taccession.version\ttaxid\tgi\n'
  zcat "$FASTA" | grep '^>' \
    | sed -n 's/^>\([^ ]*\).*TaxID=\([0-9]*\).*/\1\t\1\t\2\t0/p'
} | gzip > "$MAP.tmp.gz"
mv "$MAP.tmp.gz" "$MAP"

N=$(zcat "$MAP" | wc -l)
echo "[$(date)] mapa com $N linhas"
if [ "$N" -lt 50000000 ]; then
  echo "ERRO: mapa pequeno demais, o padrao do cabecalho nao casou" >&2
  exit 1
fi

echo "[$(date)] passo 2: diamond makedb com taxonomia"
diamond makedb --in "$FASTA" --db "$DB/uniref90_tax" \
  --taxonmap "$MAP" --taxonnodes "$TAX/nodes.dmp" --taxonnames "$TAX/names.dmp" \
  --threads "$SLURM_CPUS_PER_TASK"

echo "[$(date)] pronto"; ls -la "$DB/uniref90_tax.dmnd"
