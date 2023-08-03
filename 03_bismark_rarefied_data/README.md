# `03_bismark_rarefied_data`: mapping rarefied Illumina reads #

Raw sequencing data are not available due to re-identifiability concerns.

## Rarefaction ##

In order to remove the effect of coverage on downstream analyses, we opted to rarefy the original reads so that all eight samples had equal depth (down to the shallowest depth of 166,282,895 paired end reads).

The rarefaction was carried out using a script https://github.com/lyijin/common/blob/master/subsample_fastq.py before trimming, producing the "ER" and "WR" files from the original files with "E" and "W" suffixes ("R" means "rarefied"). The subsampling process is deterministic, as the Python script incorporates the use of a set seed.

Sanity check on pre-trimmed FASTQ files (not available due to re-identifiability concerns) to check script is working as expected:

```shell
# this command uses a bash-ism (the echo part) which divides total lines by 4, might not work with other *NIX shells
$ for a in *_R1.fastq.gz; do printf "${a}:     "; b=`zcat ${a} | wc -l`; echo $((${b}/4)); done
WR025V1ER_R1.fastq.gz:     166282895
WR025V1WR_R1.fastq.gz:     166282895
WR025V9ER_R1.fastq.gz:     166282895
WR025V9WR_R1.fastq.gz:     166282895
WR069V1ER_R1.fastq.gz:     166282895
WR069V1WR_R1.fastq.gz:     166282895
WR069V9ER_R1.fastq.gz:     166282895
WR069V9WR_R1.fastq.gz:     166282895
```

## Pipeline ##

Post-rarefaction, downstream mapping and methylation tallying follows the standard `bismark` procedure. To facilitate working on our compute cluster, we wrote a pipeline nicknamed `bismsmark` (which combines `bismark` + `snakemake`; https://github.com/lyijin/bismsmark) that mapped reads against two genomes: the overall human genome "grch38p13_lambda_puc", and the human 45S rDNA consensus sequence. Both files are available to download from CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/data`).

We opted to use the GRCh38 patch 13 human genome because it serves as the baseline for GENCODE release 38. We retained only the main chromosomes (i.e. excluded unlocalised/alt sequences), and spiked in two control sequences: lambda and pUC19. Both of the latter sequences were from NEB.

As the intermediate `fastq` and `bam` files are too large to be uploaded (and not essential for downstream analyses), this folder only contains the log files produced by each step in the `bismark` pipeline. The script in `scripts/` summarises these log files into tsv files in this folder (`compiled_bismark_logs.*.tsv`).

The cov files, produced from `bismark_methylation_extractor`, are available from CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/03_bismark_rarefied_data/05_renamed_covs`).
