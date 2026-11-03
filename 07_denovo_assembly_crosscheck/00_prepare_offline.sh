#!/bin/bash
# =============================================================================
# 00 — Preparar a corrida para rodar OFFLINE.
#
# Por que isto existe: os nós de computação deste cluster estão atrás de um
# appliance Cisco Umbrella cuja CA eles não confiam. Qualquer https a partir de
# um nó de computação morre com "SSL certificate problem: unable to get local
# issuer certificate". O nó de LOGIN confia na CA e tem rede normal.
#
# Consequência: tudo que a corrida precisa da internet — o código do pipeline,
# os configs institucionais do nf-core, o plugin nf-schema e as 33 imagens de
# container — tem que estar em disco ANTES de submeter.
#
# Rode este script no nó de LOGIN, dentro de um tmux (leva ~40 min).
#
#   tmux new -s prep
#   bash 00_prepare_offline.sh
# =============================================================================
set -euo pipefail

SCR=${SCRATCH_DIR}
PROJ=${SCR}/metatdenovo_crosscheck
CACHE=${SCR}/singularity_cache

# Diretório de blobs temporários do apptainer/singularity. Fora do cache final,
# senão os blobs OCI se misturam com os SIFs e o Nextflow se perde.
export SINGULARITY_CACHEDIR=${SCR}/sing_blobs
export APPTAINER_CACHEDIR=${SCR}/sing_blobs

mkdir -p "${CACHE}" "${SINGULARITY_CACHEDIR}" "${PROJ}/pipeline_offline"

module load apptainer/1.4.1

# ── 1. Código do pipeline + configs institucionais ───────────────────────────
# `nf-core pipelines download` traz o repositório e reescreve o
# custom_config_base para o diretório local ../configs, de modo que o
# nextflow.config não tente mais buscar nfcore_custom.config no GitHub.
#
# --container-system none: NÃO deixe o nf-core baixar as imagens. Ele usa o
# registro Docker (quay.io) e o que cai no cache são camadas OCI, não SIFs — o
# Nextflow não reconhece e tenta baixar tudo de novo do nó de computação.
# As imagens são puxadas no passo 2, uma a uma, dos URLs que os próprios
# módulos declaram.
if [ ! -d "${PROJ}/pipeline_offline/1_4_0" ]; then
    nf-core pipelines download metatdenovo \
        --revision 1.4.0 \
        --outdir "${PROJ}/pipeline_offline" \
        --compress none \
        --container-system none
fi

# ── 2. Imagens de container ──────────────────────────────────────────────────
# A lista real sai dos módulos do pipeline. São três origens:
#   https://depot.galaxyproject.org/singularity/...   SIF direto
#   https://community-cr-prod.seqera.io/.../data      SIF direto (Wave)
#   oras://community.wave.seqera.io/library/...       OCI via ORAS
#
# O nome do arquivo no cache segue a convenção do Nextflow: tira o protocolo,
# troca "/" e ":" por "-", acrescenta ".img". Se o nome não bater exatamente,
# o Nextflow baixa de novo — do nó de computação, onde não há rede.
P=${PROJ}/pipeline_offline/1_4_0
grep -rhoE "(https://[^'\" ]+|oras://[^'\" ]+)" "$P/modules" "$P/subworkflows" \
  | grep -E 'seqera|galaxyproject|depot' | sort -u > "${PROJ}/containers_needed.txt"

echo "Containers referenciados: $(wc -l < "${PROJ}/containers_needed.txt")"

while read -r url; do
    name=$(echo "$url" | sed -E 's|^https://||; s|^oras://||; s|[/:]|-|g').img
    if [ -s "${CACHE}/${name}" ]; then echo "JA TEM  $name"; continue; fi
    if apptainer pull --name "${CACHE}/${name}" "$url" >/dev/null 2>&1; then
        echo "OK      $name"
    else
        echo "FALHOU  $url"
    fi
done < "${PROJ}/containers_needed.txt"

# NOTA: duas linhas FALHOU são esperadas e inofensivas — as variantes arm64 de
# fastqc e multiqc (oras://). As amd64 correspondentes vêm de outra URL e são
# baixadas com sucesso. 33 de 35 é o resultado correto.

# ── 3. Plugin nf-schema ──────────────────────────────────────────────────────
# O metatdenovo 1.4.0 declara nf-schema@2.7.2 e o módulo nextflow/25.04 é velho
# demais para conhecê-lo. Puxar aqui, no login, deixa o plugin em
# ~/.nextflow/plugins e o NXF_OFFLINE do job não precisa da rede.
module load nextflow/25.04
export NXF_VER=26.04.6
nextflow plugin install nf-schema@2.7.2 || true
ls ~/.nextflow/plugins/

# ── 4. Conferência ───────────────────────────────────────────────────────────
echo
echo "=============================================="
echo "imagens no cache : $(find "${CACHE}" -maxdepth 1 -type f -name '*.img' | wc -l)"
echo "pipeline local   : $([ -d "$P" ] && echo OK || echo FALTA)"
echo "configs locais   : $([ -d "${PROJ}/pipeline_offline/configs" ] && echo OK || echo FALTA)"
echo "plugin nf-schema : $([ -d ~/.nextflow/plugins/nf-schema-2.7.2 ] && echo OK || echo FALTA)"
echo "=============================================="
