# 01 — Read processing (HPC, SLURM)

Run in order from `$PROJ`, with `PROJ` and `SCRATCH_DIR` exported (see the main README). Array jobs run one task per sample, reading `sample_metadata.tsv`.

| Step | Script | Tool | Output used downstream |
|---|---|---|---|
| 1 | `run_01_qc_trim_multiqc.slurm` | FastQC, fastp (Q20, ≥ 50 bp), MultiQC | `data/trim/*_{R1,R2}.trim.fq.gz` |
| 2 | `run_02_sortmerna.slurm` | SortMeRNA, `smr_v4.3_default_db` | `data/mRNA/*_mRNA_fwd.fq.gz` |
| 3 | `run_03_download_index_uniref90.sh` | UniRef90 download, `diamond makedb` | `db/uniref90/uniref90.dmnd` |
| 5 | `run_05_diamond_uniref90.sh` | DIAMOND blastx, `--id 80 --evalue 1e-5`, best hit | `results/diamond/uniref90/*.tsv` |
| 8 | `run_08_extract_uniref90_hits.sh` | seqtk | FASTA of the UniRef90 proteins hit by ≥ 1 read |
| 9 | `run_09_emapper_uniref90_hits.sh` | eggNOG-mapper | protein → KO and KEGG_Pathway |
| 10 | `run_10_extract_uniref90_taxids.sh` | header parsing | UniRef90 ID → NCBI TaxID |
| 11 | `run_11_taxonkit_lineage.sh` | TaxonKit | TaxID → lineage |
| 12 | `12a_parse_taxonomy_counts.py` | | `counts_family_{bacteria,fungi,archaea}.tsv` |
| 13 | `run_13_parse_ko_counts.sh` → `13_parse_ko_counts.py` | | KO × sample counts |
| 13b | `run_13b_parse_ko_taxonomy_full.sh` → `13b_parse_ko_taxonomy_full.py` | | `counts_ko_taxonomy_full.tsv` (family × KO × sample, long format), `ko_taxonomy_full_stats.tsv`, `ko_to_pathway.tsv` |

Searching only the proteins hit by at least one read (steps 8–9) keeps eggNOG-mapper to a small fraction of UniRef90. The step numbers are not contiguous; the missing numbers belong to analyses that are not part of this workflow.

`robustness_identity_threshold/` re-parses the DIAMOND output at 80, 90 and 95 % identity in one pass (`counts_ko_taxonomy_pid{80,90,95}.tsv.gz`), so that family-resolved results can be checked at stricter assignment thresholds. `06_taxonomic_decomposition/decomposition_contrasts.py` accepts each table unchanged.
