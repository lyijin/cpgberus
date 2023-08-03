# `14_methepic_vs_emgseq_wgbs/` folder #

## Data used for plots ##

Processing of EPIC data is described in `02_process_methepic_data/`, while the compiled betas for EM-seq and WGBS data sits in `04_parse_bismark_covs/`.

Genome annotations are from a processed RData in `data/`.

## Idea behind R code ##

Folder contains two mega R scripts.

The first one, `methepic_emseq_wgbs_per-pos_beta.R`, carries out a per-position, three-way comparisons between EM-seq, WGBS and EPIC readouts across common methylated positions across the human genome.

The second one, `methepic_emseq_wgbs_feature_beta.R`, builds upon the first, and specifically aggregates CpGs by their genomic annotations. Do, say, CpG islands have different readouts across different experimental method?

The HTML files in this folder is a quick way to view the R scripts + the plots that were made. Not all plots were included in the manuscript.
