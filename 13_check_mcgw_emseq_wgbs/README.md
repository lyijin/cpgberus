# `13_check_mcgw_emseq_wgbs/` folder #

## Idea behind R code ##

Folder contains two R scripts to answer two separate questions we had that were sparked mostly after coming across Adam et al., 2022 (https://www.nature.com/articles/s42003-022-03033-4) and Ravichandran et al. 2022 (https://www.science.org/doi/10.1126/sciadv.abm2427).

Note that Adam et al. does not mention MCGW explicitly, but when one zooms into TET2 and 5mC in particular (Fig. 1E), preference for MCGW is obvious.

1. Papers claim that MCGW is overrepresented in sequences that is preferentially demethylated by TET2. M = A/C, W = A/T. From the list of 124 CpGs with differential WGBS vs. EM-seq methylation, how many of them contain MCGW? How many do not? Do these numbers occur in a statistically skewed manner? Answered by `sigcpgs_mcgw_bias.R`

2. MCGW is a non-palindromic motif. Papers claims that strand with MCGW is preferentially demethylated over the non-MCGW strand. Potentially dovetailing with our observations that differentially methylated positions (EM-seq vs. WGBS) have uneven coverage and uneven methylation across strands. Is this observable in the overall dataset? Answered by `overall_mcgw_bias.R`

The HTML files in this folder is a quick way to see the R scripts in action + description of logic + the plots that were made. Not all plots were included in the manuscript.
