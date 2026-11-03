# De novo assembly cross-check with nf-core/metatdenovo

The main workflow annotates individual reads: each 150 bp forward read is searched against UniRef90 and inherits the function and taxonomy of its best hit. That is sensitive, but a ~50 amino-acid fragment can be hard to place, and there is no eukaryote-aware gene model. As a double check, the same trimmed reads are processed independently by co-assembly with **nf-core/metatdenovo**, which assembles transcripts, calls open reading frames (ORFs), annotates them and counts reads per ORF. Agreement between the two routes on overall annotation recovery, domain composition and the number of testable KOs per domain indicates that conclusions do not depend on the read-based method; disagreement shows where they might.

**Pipeline:** nf-core/metatdenovo 1.4.0 — Di Leo D, Nilsson E, Krinos A, Pinhassi J, Lundin D. 2025. The Nextflow nf-core/metatdenovo pipeline for reproducible annotation of metatranscriptomes, and more. *PeerJ* 13:e20328. [doi:10.7717/peerj.20328](https://doi.org/10.7717/peerj.20328) · [nf-co.re/metatdenovo](https://nf-co.re/metatdenovo/1.4.0)

## Configuration used

| Choice | Setting | Reason |
|---|---|---|
| Input | trimmed pairs from `01_read_processing` (`data/trim/`), `--skip_trimming` | same starting reads as the read-based route; properly paired |
| rRNA removal | in-pipeline BBDuk (`--sequence_filter`) against the SortMeRNA default reference | independent of SortMeRNA, so rRNA fractions are compared, not assumed equal |
| Assembly | MEGAHIT co-assembly of all samples | pools depth for low-abundance taxa such as fungi |
| Minimum contig length | 300 nt (`min_contig_length` in `cluster.config`) | TransDecoder's default minimum ORF is 100 aa |
| ORF calling | TransDecoder | eukaryote-aware gene model |
| Taxonomy | DIAMOND against UniRef90 rebuilt with taxonomy (`db_preparation/01_uniref90_taxdb.slurm`) | the read-based index carries no taxonomy; the accession → TaxID map is taken from the `TaxID=` field of the UniRef90 headers |
| Function | eggNOG-mapper and KofamScan | KofamScan gives KOs, the unit of the read-based analysis |
| EUKulele | skipped (`--skip_eukulele`) | can be added later on the same ORFs |

## Running

```bash
export SCRATCH_DIR=/path/to/scratch
export PROJ=/path/to/read-processing/project       # holds data/trim/, db/uniref90/ and sample_metadata.tsv
sbatch db_preparation/01_uniref90_taxdb.slurm     # UniRef90 DIAMOND index with --taxonmap/--taxonnodes/--taxonnames
sbatch db_preparation/02_kofam_download.slurm     # KofamScan ko_list and profiles
bash   03_make_samplesheet.sh                     # samplesheet.csv from sample_metadata.tsv
bash   00_prepare_offline.sh                      # on a login node: pipeline, configs, plugin and container images
sbatch 04_metatdenovo.slurm                       # Nextflow driver job; reads cluster.config and diamond_dbs.csv
```

`00_prepare_offline.sh` is needed only where compute nodes cannot reach the internet: it downloads the pipeline with `nf-core pipelines download`, together with its configs, the `nf-schema` plugin and the container images, and the run then uses `NXF_OFFLINE=true`. `cluster.config` holds the SLURM executor, resource and Apptainer settings of the cluster where this was run and must be adapted elsewhere. Comments in these scripts are in Portuguese.

## What to compare

| Metric | Read-based route | Assembly route |
|---|---|---|
| rRNA fraction | SortMeRNA log | BBDuk logs |
| Reads reaching annotation | reads with a UniRef90 hit (`ko_taxonomy_full_stats.tsv`) | reads mapped back to contigs and assigned to ORFs (`summary_tables/*.overall_stats.tsv.gz`) |
| Domain composition | `counts_family_*.tsv` | `summary_tables/*.uniref90.diamond.taxonomy.tsv.gz` |
| Testable KOs per domain | `02_count_matrices` output | distinct KofamScan KOs per domain (`summary_tables/*.kofamscan-uniq.tsv.gz`) |

## Software run inside the pipeline

| Tool | Reference |
|---|---|
| nf-core/metatdenovo 1.4.0 | Di Leo D et al. 2025. *PeerJ* 13:e20328. [doi:10.7717/peerj.20328](https://doi.org/10.7717/peerj.20328) |
| nf-core | Ewels PA et al. 2020. *Nat Biotechnol* 38:276–278. [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x) |
| Nextflow 26.04.6 | Di Tommaso P et al. 2017. *Nat Biotechnol* 35:316–319. [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| Apptainer 1.4.1 (Singularity) | Kurtzer GM, Sochat V, Bauer MW. 2017. *PLoS One* 12:e0177459. [doi:10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459) |
| BBDuk, BBMap (BBTools) | Bushnell B. BBMap. [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/) |
| MEGAHIT | Li D, Liu CM, Luo R, Sadakane K, Lam TW. 2015. *Bioinformatics* 31:1674–1676. [doi:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033) |
| TransDecoder | Haas BJ et al. 2013. *Nat Protoc* 8:1494–1512. [doi:10.1038/nprot.2013.084](https://doi.org/10.1038/nprot.2013.084); [github.com/TransDecoder](https://github.com/TransDecoder/TransDecoder) |
| featureCounts (Subread) | Liao Y, Smyth GK, Shi W. 2014. *Bioinformatics* 30:923–930. [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656) |
| DIAMOND | Buchfink B, Reuter K, Drost HG. 2021. *Nat Methods* 18:366–368. [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x) |
| eggNOG-mapper | Cantalapiedra CP et al. 2021. *Mol Biol Evol* 38:5825–5829. [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293) |
| KofamScan | Aramaki T et al. 2020. *Bioinformatics* 36:2251–2252. [doi:10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859) |
