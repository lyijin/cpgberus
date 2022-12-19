#!/usr/bin/env Rscript

"> sigcpgs_mcgw_bias.R <

There are 124 CpGs where methylation readouts in EM-seq differs substantially
from WGBS. This inter-method difference could be driven by TET2's preference
for MCGW motifs (Ravichandran, Sci Adv, 2022). In that case, when we parse the
NCGN sequence context surrounding the 124 CpGs, could there be an unexpected
enrichment/depletion of MCGW amongst the 16 possible NCGN motifs?

When this script checks for 'MCGW', it actually checks for MCGW or WCGK. This
is because MCGW is not palindromic. The reverse complement of MCGW, WCGK, 
implies that MCGW is on the non-reference strand.

With this definition of MCGW, a Python script that parsed the MCGW/NCGN ratio
across the reference genome of GRCh38.p13 returned a baseline rate of 44.68%
(not 50%, there are only 7 MCGW motifs. 7/16 = 43.75%, so 44.68% is close-ish
to uninformed, naive expectations).
" -> doc

setwd('~/csiro/stopwatch/cpgberus/13_check_mcgw_emseq_wgbs/')

SIG_CPGS <- './significant_regions.csv'
MCGW_4MERS <- c('ACGA', 'ACGT', 'CCGA', 'CCGT')
WCGK_4MERS <- c('ACGG', 'ACGT', 'TCGG', 'TCGT')
MCGW_WCGK <- unique(c(MCGW_4MERS, WCGK_4MERS))
length(MCGW_WCGK)  # 7 is expected here

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# read CpG positions with significantly different meth readouts between
# EM-seq and WGBS
sig_df <- read.csv(SIG_CPGS, 
                   colClasses=c('NULL', NA, NA, NA, 'NULL', 'NULL', NA, NA))
head(sig_df)
# mu1 refers to the mean methylation beta in EM-seq samples; and likewise
# mu2 for WGBS samples

# this csv file does not have the NNNCGNNN context around it. use Biostrings
# and BSgenome.Hsapiens.UCSC.hg38 to extract the sequence contexts. but first,
# extract chrom:start-end and make sure ALL sequences are "CG"
sig_gr <- GRanges(seqnames=sig_df$seqnames,
                  ranges=IRanges(start=sig_df$start, end=sig_df$end))
stopifnot(all(getSeq(Hsapiens, sig_gr) == 'CG'))

# good. we can then proceed to extract the NCGN context (to check whether
# MCGW is over-represented
sig_gr <- GRanges(seqnames=sig_df$seqnames,
                  ranges=IRanges(start=sig_df$start-1, end=sig_df$end+1))
sig_df <- cbind(sig_df, ncgn=as.character(getSeq(Hsapiens, sig_gr)))
sig_df$bool_mcgw <- sig_df$ncgn %in% MCGW_WCGK
head(sig_df)
table(sig_df$bool_mcgw)

# calculate p value based on GRCh38.p13's expected rate of 44.68%. model
# observation against a binomial distribution.
# H0: proportion of MCGW is similar to genome's, i.e., no reason to suspect
#     that MCGW drives differences in methylation level observations in
#     EM-seq vs. WGBS
# H1: proportion of MCGW is greater/smaller than expected, motif does indeed
#     influence differences in inter-method readouts
pbinom(sum(sig_df$bool_mcgw), size=length(sig_df$bool_mcgw), prob=0.4468)
# p > 0.05, unable to reject H0. MCGW is present in expected quantities

# changing the probabilities to uninformed, naive (7 motifs of 16) does not
# change the p value much. H0 is unable to be rejected
pbinom(sum(sig_df$bool_mcgw), size=length(sig_df$bool_mcgw), prob=0.4375)

# list deps used in this script
sessionInfo()
