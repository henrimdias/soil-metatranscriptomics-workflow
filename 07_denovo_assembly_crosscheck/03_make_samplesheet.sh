#!/bin/bash
# =============================================================================
# 03 — Build the nf-core/metatdenovo samplesheet from sample_metadata.tsv.
#
# Uses the fastp-trimmed reads of 01_read_processing (data/trim/), which are
# correctly paired. It deliberately does NOT use data/mRNA/*_mRNA_{fwd,rev}.fq.gz:
# those are not synchronised pairs (see the SortMeRNA note in the main README).
# rRNA removal happens inside the pipeline instead (--sequence_filter).
#
# Sample IDs are taken from the sample_id column of sample_metadata.tsv, so the
# assembly-based tables join directly onto the read-based ones.
#
#   export SCRATCH_DIR=... PROJ=...   # PROJ = read-processing directory
#   bash 03_make_samplesheet.sh
# =============================================================================
set -euo pipefail

META=${PROJ}/sample_metadata.tsv     # columns: sample_id  file_prefix  ...
TRIM=${PROJ}/data/trim
OUT=samplesheet.csv

echo "sample,fastq_1,fastq_2" > "${OUT}"
missing=0
while IFS=$'\t' read -r sample prefix _; do
    r1=${TRIM}/${prefix}_R1.trim.fq.gz
    r2=${TRIM}/${prefix}_R2.trim.fq.gz
    for f in "$r1" "$r2"; do
        [ -s "$f" ] || { echo "MISSING: $f" >&2; missing=$((missing+1)); }
    done
    echo "${sample},${r1},${r2}" >> "${OUT}"
done < <(tail -n +2 "${META}")

if [ "${missing}" -gt 0 ]; then
    echo "ERROR: ${missing} read file(s) missing — samplesheet written but NOT usable." >&2
    exit 1
fi
echo "Wrote ${OUT}: $(( $(wc -l < "${OUT}") - 1 )) samples, all files present."
