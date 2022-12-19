#!/usr/bin/env Rscript

"> overall_mcgw_bias.R <

As MCGW is not fully palindromic, the MCGW-containing strand could potentially
be preferentially demethylated over the non-MCGW-containing one (Ravichandran,
Sci Adv, 2022). Can we replicate that observation with our EM-seq vs. WGBS
samples, and partially explain our observed inter-strand coverage and
methylation betas differences?

This script parses through the intermediate files that calculated, on a per-
CpG dinucleotide basis, difference in coverage and betas across strands.

There are three categories:
  0. Neither strand are MCGW.                               9 possibilities
  1. One strand is MCGW, opposite is not.                   6 possibilities
  2. One strand is MCGW, opposite is ALSO MCGW (ACGT/ACGT). 1 possibility
                                                           ----------------
                                                     Total 16 possibilities (NCGN)

Null expectation is that (1) would look different from (0) and (2). Latter two
should have comparable amounts of inter-strand differences.

The kicker would be--if indeed MCGW-containing strands have higher demethylation
rates, is this specific to the EM-seq samples? As WGBS does not rely on TET2
in the conversion process, it should NOT have the same MCGW-related differences.
" -> doc

setwd('~/csiro/stopwatch/cpgberus/13_check_mcgw_emseq_wgbs/')

COMBINED_TSVS <- Sys.glob('../04_parse_bismark_covs/WR*R.grch38p13_lambda_puc.combined.tsv.gz')
CAT0_4MERS <- c('ACGC', 'CCGC', 'CCGG', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'TCGA', 'TCGC')
CAT1_4MERS <- c('ACGA', 'ACGG', 'CCGA', 'CCGT', 'TCGG', 'TCGT')
CAT2_4MERS <- 'ACGT'

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}

read_combined_tsv <- function(filename) {
  # use data.table to quickly read table
  dt <- fread(filename)
  diag_message('Pre-filtered table has ', nrow(dt), ' dinucleotides.')
  
  # filter for CpG dinucleotides with coverage >= 5 for computational
  # efficiency, and to focus analysis on better-covered dinucs
  dt <- dt[dt$meth_cov + dt$unmeth_cov >= 5]
  diag_message('Post-filtered table has ', nrow(dt), ' dinucleotides.')
  
  # create a temp array to sort out categorisation
  temp <- dt$cpg_context_nnncgnnn
  temp[substr(dt$cpg_context_nnncgnnn, 3, 6) %in% CAT0_4MERS] <- 'CAT0'
  temp[substr(dt$cpg_context_nnncgnnn, 3, 6) %in% CAT1_4MERS] <- 'CAT1'
  temp[substr(dt$cpg_context_nnncgnnn, 3, 6) %in% CAT2_4MERS] <- 'CAT2'
  
  # remove rows that do not fall into these three categories + sanity check
  dt <- dt[grepl('^CAT', temp)]
  temp <- temp[grepl('^CAT', temp)]
  temp <- as.factor(temp)
  stopifnot(nlevels(temp) <= 3)
  
  # split data table into cat1/2/3 based on NCGN sequence
  dt$category <- temp
  
  # categorisation shouldn't remove rows, double check
  diag_message('Post-categorised table has ', nrow(dt), ' dinucleotides.')
  
  # convert meth% into beta (divide by 100) for graph plotting
  dt$meth_pct <- dt$meth_pct / 100
  dt$abs_delta_meth_pct <- dt$abs_delta_meth_pct / 100
  colnames(dt)[4] <- 'beta'
  colnames(dt)[8] <- 'abs_delta_beta'
  
  # add sample ID into dt
  sample_id <- sub('.*_covs/', '', filename)
  sample_id <- sub('\\.grch38.*', '', sample_id)
  dt$sample_id <- sample_id
  
  dt
}

# hacky rbinding across 8 files
for (n in 1:length(COMBINED_TSVS)) {
  diag_message('Processing "', COMBINED_TSVS[n], '"...')
  
  if (n == 1) {
    dt <- read_combined_tsv(COMBINED_TSVS[n])
  } else {
    dt <- rbind(dt, read_combined_tsv(COMBINED_TSVS[n]))
  }
}

# plot separate box-and-whiskers plot for beta, abs_delta_beta and evenness
#+ fig.width=10, fig.height=5
ggplot(dt, aes(x=sample_id, y=beta, fill=category)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual('MCGW-containing strands', labels=c(0,1,2), 
                    values=c('#e7e1ef','#c994c7','#dd1c77')) +
  ggtitle('Effects of MCGW on beta') +
  xlab('Sample ID') +
  ylab('Beta') +
  theme_minimal(12) +
  theme(legend.position='top')
# interesting to note that beta for ACGT/ACGT is lower than the other two
# categories. not sure why; but differences are CONSISTENT ACROSS EM-SEQ
# AND WGBS. more likely the differences are related to genome properties, and
# not due to conversion method

ggplot(dt, aes(x=sample_id, y=abs_delta_beta, fill=category)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual('MCGW-containing strands', labels=c(0,1,2), 
                    values=c('#e7e1ef','#c994c7','#dd1c77')) +
  ggtitle('Effects of MCGW on absolute delta beta') +
  xlab('Sample ID') +
  ylab('Absolute delta beta') +
  theme_minimal(12) +
  theme(legend.position='top')
# similar observation as raw beta values, again mostly independent of method

ggplot(dt, aes(x=sample_id, y=evenness, fill=category)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual('MCGW-containing strands', labels=c(0,1,2), 
                    values=c('#e7e1ef','#c994c7','#dd1c77')) +
  ggtitle('Effects of MCGW on evenness',
          subtitle='Evenness: 0.5 = perfectly even; 1.0 = perfectly uneven') +
  xlab('Sample ID') +
  ylab('Evenness') +
  theme_minimal(12) +
  theme(legend.position='top')
# well, nice to see that at least for one measure, it is mostly independent
# of MCGW context...

# list deps used in this script
sessionInfo()
