# `02_bismark_wgbs_emseq_data`: mapping raw Illumina reads to the human genome #

Raw data are not provided in this repo, but can be downloaded from (TODO: insert NCBI SRA/CSIRO DAP link here when data is uploaded).

## Pipeline ##

The `bismsmark` pipeline (which combines `bismark` + `snakemake`) was used to map the reads against the "grch38p13_lambda_puc" genome.

This genome, GRCh38 patch 13, was used because it serves as the baseline for GENCODE release 38. We retained only the main chromosomes (i.e. excluded unlocalised/alt sequences), and spiked in two control sequences: lambda and pUC19. Both of the latter sequences were from NEB.

As the intermediate `fastq` and `bam` files are too large to be uploaded, we have only included the final compressed `cov` files produced by `bismark`. These files are the inputs for the multi-headed comparisons.
