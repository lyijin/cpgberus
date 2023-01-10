# `13_check_mcgw_emseq_wgbs/` folder #

## Idea behind R code ##

Folder contains two R scripts to answer two separate questions we had that were sparked mostly after coming across Ravichandran et al. 2022, https://www.science.org/doi/10.1126/sciadv.abm2427.

1. Paper claims that MCGW is overrepresented in sequences that is preferentially demethylated by TET2. M = A/C, W = A/T. From the list of 124 CpGs with differential WGBS vs. EM-seq methylation, how many of them contain MCGW? how many do not? do these numbers occur in a statistically skewed manner? answered by `sigcpgs_mcgw_bias.R`

2. MCGW is a non-palindromic motif. paper claims that strand that has MCGW is preferentially demethylated over the non-MCGW strand. dovetails with our observation that diff meth positions have uneven coverage and uneven methylation across strands. my guess is that the strand with MCGW has lower methylation than the non-MCGW—and this observation should only happen for EM-seq, not for WGBS. answered by `overall_mcgw_bias.R`

The HTML file in this folder is a quick way to see the R script in action + the plots that were made. Not all plots were included in the manuscript.
