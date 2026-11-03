#!/bin/bash
# =============================================================================
# reanálise fúngica por montagem de novo
# Etapa 4: nf-core/metatdenovo v1.4.0
#
#   co-assembly MEGAHIT das 18 amostras
#   + ORFs com TransDecoder (eucariótico)
#   + taxonomia DIAMOND/UniRef90 (índice com taxonomia, etapa 01)
#   + função com eggNOG-mapper e KofamScan
#
# Pré-requisitos, nesta ordem:
#   01_uniref90_taxdb.slurm   feito  -> uniref90_tax.dmnd
#   02_kofam_download.slurm   feito  -> dbs/kofam
#   03_make_samplesheet.sh    feito  -> samplesheet.csv
#   cluster.config e diamond_dbs.csv    no mesmo diretório
#
# Uso:      sbatch 04_metatdenovo.slurm
# Retomar:  o mesmo comando. O -resume já está ligado e o Nextflow continua de
#           onde parou. NUNCA apague o work dir antes de um resume.
#
# Autor: Henrique M. Dias / SDSU
# =============================================================================
#SBATCH --job-name=metatdenovo
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/metatdenovo_%j.out
#SBATCH --error=logs/metatdenovo_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

set -euo pipefail

# ── Ambiente ─────────────────────────────────────────────────────────────────
# Este job é só o DRIVER: ele fica vivo o tempo todo e submete cada etapa do
# pipeline como um job SLURM próprio (ver cluster.config). Por isso pede tempo
# longo mas quase nenhum CPU ou memória.
module load nextflow/25.04 apptainer/1.4.1

# O módulo é de maio/2025 e não enxerga o plugin nf-schema 2.7.2 que o
# metatdenovo 1.4.0 exige. NXF_VER faz o lançador baixar a versão certa para
# ~/.nextflow/framework sem mexer no módulo.
export NXF_VER=26.04.6
export NXF_OFFLINE=true

export NXF_APPTAINER_CACHEDIR=${SCRATCH_DIR}/singularity_cache
export NXF_OPTS='-Xms1g -Xmx8g'

# ── Caminhos ─────────────────────────────────────────────────────────────────
SCR=${SCRATCH_DIR}
PROJ=${SCR}/metatdenovo_crosscheck
DBS=${SCR}/dbs

OUTDIR=${PROJ}/results/megahit
WORKDIR=${SCR}/nxf_work_metatdenovo          # fora do projeto: chega a 1-2 TB e é descartável

mkdir -p logs "${OUTDIR}" "${WORKDIR}" "${NXF_APPTAINER_CACHEDIR}"

echo "=============================================="
echo "Job     : ${SLURM_JOB_ID}"
echo "Início  : $(date)"
echo "Saída   : ${OUTDIR}"
echo "Work    : ${WORKDIR}"
echo "=============================================="

# ── Corrida ──────────────────────────────────────────────────────────────────
# Justificativa de cada escolha, para ficar auditável depois:
#
#   --skip_trimming     os reads já passaram por fastp no workflow read-based (01_read_processing). Reaparar
#                       mudaria o denominador e quebraria a comparabilidade com
#                       a contagem de reads pós-fastp do workflow read-based.
#
#   --sequence_filter   remoção de rRNA por k-mer com BBDuk contra a referência
#                       do SortMeRNA. É ferramenta DIFERENTE da do workflow read-based (01_read_processing),
#                       então os números de rRNA não vão bater — e não devem.
#                       Reportar separado, sem tentar reconciliar.
#
#   --orf_caller        transdecoder é o chamador eucariótico. Esta é a maior
#                       mudança metodológica em relação ao workflow read-based (01_read_processing), que não
#                       tinha etapa eucariótica nenhuma.
#
#   --min_contig_length o mínimo padrão do TransDecoder é 100 aa = 300 nt, logo
#                       contig menor que isso não produz ORF completo.
#
#   --diamond_dbs       índice UniRef90 reconstruído COM taxonomia na etapa 01.
#                       O índice antigo do workflow read-based (01_read_processing) não tem taxonomia.
#
#   --skip_eukulele     o EUKulele é camada de validação, não pré-requisito.
#                       Ele roda sobre os ORFs, e os ORFs não mudam — dá para
#                       rodar depois sobre o mesmo .faa sem refazer a montagem.
#
#   --diamond_top       fica no padrão 10. O paper do pipeline mostra que 3
#                       empurra mais atribuições para nível de espécie; isso é
#                       análise de sensibilidade da fase 2, não a principal.

nextflow run "${PROJ}/pipeline_offline/1_4_0" \
    -profile apptainer \
    -c cluster.config \
    -work-dir "${WORKDIR}" \
    -resume \
    --input             samplesheet.csv \
    --outdir            "${OUTDIR}" \
    --skip_trimming \
    --sequence_filter   "${DBS}/sortmerna/smr_v4.3_default_db.fasta" \
    --assembler         megahit \
    --orf_caller        transdecoder \
    --diamond_dbs       diamond_dbs.csv \
    --skip_eukulele \
    --eggnog_dbpath     "${DBS}/eggnog" \
    --kofam_dir         "${DBS}/kofam" \
    -with-report        "${OUTDIR}/nf_report.html" \
    -with-timeline      "${OUTDIR}/nf_timeline.html" \
    -with-trace         "${OUTDIR}/nf_trace.txt"

echo "=============================================="
echo "Fim : $(date)"
echo "=============================================="
