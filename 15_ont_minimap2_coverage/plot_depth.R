#!/usr/bin/env Rscript

"> plot_depth.R <

After mapping FASTQ files against 45S reference with `minimap2`, and after
parsing the BAM files with `samtools depth` into plaintext TSV files, plot the
data out in those tables.
" -> doc

setwd('~/csiro/stopwatch/human_tests/05_ramac.ont_long_reads/02_minimap2_vs_45s/')
PARSED_MPILEUP_TSV <- c(Sys.glob('./W*.parsed_mp.tsv.gz'), Sys.glob('./non-*.parsed_mp.tsv.gz'))

suppressPackageStartupMessages({
  library(ggplot2)
})

combined_df <- data.frame()
for (pmt in PARSED_MPILEUP_TSV) {
  barcode <- strsplit(pmt, '.', fixed=TRUE)[[1]][2]
  barcode <- gsub('/', '', barcode, fixed=TRUE)
  
  # read df
  df <- read.delim(
    pmt, header=FALSE,
    col.names=c('scaf', 'pos', 'cov_total', 'cov_watson', 'cov_crick',
                'cov_a', 'cov_c', 'cov_g', 'cov_t'))
  
  # only keep relevant columns
  df <- df[, c(2, 4, 5)]
  
  # standardise pos wrt reference sequence used in WGBS/EM-seq work, and make
  # sure transcribed 45S region starts at pos +1. the nanopore refseq has 11 kb
  # of tail end sequence moved to the 5' end; while the short reads refseq only
  # has 1 kb tail end sequence moved to 5' end
  # TL;DR: subtract 11 kb from pos values from minimap2 results
  df$pos <- df$pos - 11000
  head(df)
  
  # reshape wide-to-long
  df <- reshape(
    df, direction='long',
    varying=c('cov_watson', 'cov_crick'), v.names='cov',
    idvar='pos', timevar='strand', times=c('watson', 'crick'))
  head(df)
  
  # add in provenance
  df$dataset <- barcode
  
  # rbind everything together
  combined_df <- rbind(combined_df, df)
}

#+ fig.width=8, fig.height=8
# plot coverage across loci
g <- ggplot(combined_df, aes(x=pos, y=cov)) +
  geom_bar(aes(fill=strand), stat='identity', width=1, alpha=0.8) +
  scale_fill_manual(values=c('#5ab4ac', '#d8b365')) +
  scale_x_continuous(breaks=c(-891, -290, 1, 13332, 16541, 18425), guide=guide_axis(n.dodge=2)) +
  # x-axis: go slightly upstream and downstream of where probes cut
  # y-axis: force 3,000 max across all plots to ease visual comparison
  coord_cartesian(xlim=c(-1000, 19000), ylim=c(0, 3000)) +
  # shaded region, 1:13332, is transcribed 45S region
  annotate('rect', xmin=1, xmax=13332, ymin=-Inf, ymax=Inf, alpha=0.2) +
  facet_grid(rows=vars(dataset)) +
  xlab('Position along 45S reference (bp)') +
  ylab('Coverage') +
  theme_minimal(14) +
  theme(legend.position='top',
        panel.grid.minor.x=element_blank(),
        panel.spacing.y = unit(1.5, 'lines'))
print(g)
ggsave('plot_depth.pdf', plot=g, width=8, height=8)
