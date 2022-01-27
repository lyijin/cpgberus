#!/usr/bin/env Rscript

"> three-way_per-pos_beta.R <

Reads in betas from EM-seq, WGBS and ONT Cas9 experiments, and compares them 
on a per-position basis.

Shares the same cutoffs used in the coverage script, so that outputs from both
scripts are comparable.
" -> doc

setwd('~/csiro/stopwatch/cpgberus/16_loci_specific_three_way/')

# define constants
RDNA_FASTA <- '../data/homo_sapiens.45s.fa.gz'
EMSEQ_COV <- Sys.glob('./WR*ER.hsap_45s.cov.gz')
WGBS_COV <- Sys.glob('./WR*WR.hsap_45s.cov.gz')
ONT_BED <- Sys.glob('./WR*.bed.gz')

EMSEQ_COLOR <- '#1b9e77'
ONT_COLOR <- '#d95f02'
WGBS_COLOR <- '#7570b3'

suppressPackageStartupMessages({
  library(cowplot)
  library(eulerr)
  library(ggplot2)
  library(RColorBrewer)
})

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}

# function to annotate line of best fit in scatterplot
lm_eqn <- function(df, x, y) {
  # modified from https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  model <- lm(as.formula(paste(y, '~', x)), df)
  eq <- substitute(
    italic(y) == c + m %.% italic(x)*','~italic(r)^2 == r2*','~italic(p) == pval,
    list(c = format(unname(coef(model)[1]), digits=3),
         m = format(unname(coef(model)[2]), digits=3),
         r2 = format(summary(model)$r.squared, digits=4),
         pval = format(summary(model)$coefficients[2,4], digits=1)))
  as.character(as.expression(eq))
}

# function to calculate sum/len of a binary string
pct_of_ones <- function(binary_string) {
  int_vector <- as.integer(unlist(strsplit(binary_string, split = "")))
  
  sum(int_vector) / length(int_vector) * 100
}

# read data
rdna_seq <- scan(RDNA_FASTA, what='character', skip=1)

# read EM-seq and WGBS data into separate dfs, but logic-wise they're fairly similar
emseq_df <- data.frame()
for (f in EMSEQ_COV) {
  sample_id <- strsplit(f, '.', fixed=TRUE)[[1]][2]
  sample_id <- gsub('/', '', sample_id, fixed=TRUE)
  
  temp_df <- read.delim(
    f, header=FALSE, 
    col.names=c('scaf', 'pos', 'pos2', 'beta', 'meth', 'unmeth'))
  temp_df$sample_id <- sample_id
  
  emseq_df <- rbind(emseq_df, temp_df)
}

wgbs_df <- data.frame()
for (f in WGBS_COV) {
  sample_id <- strsplit(f, '.', fixed=TRUE)[[1]][2]
  sample_id <- gsub('/', '', sample_id, fixed=TRUE)
  
  temp_df <- read.delim(
    f, header=FALSE, 
    col.names=c('scaf', 'pos', 'pos2', 'beta', 'meth', 'unmeth'))
  temp_df$sample_id <- sample_id
  
  wgbs_df <- rbind(wgbs_df, temp_df)
}

# tidy up EM-seq/WGBS df
# first 1 kb is promoter sequence; make sure transcription starts at +1
emseq_df$pos <- emseq_df$pos - 1000
wgbs_df$pos <- wgbs_df$pos - 1000
emseq_df$cov <- emseq_df$meth + emseq_df$unmeth
wgbs_df$cov <- wgbs_df$meth + wgbs_df$unmeth

# ONT BEDs are "bedMethyl" files
# https://www.encodeproject.org/data-standards/wgbs/
ont_df <- data.frame()
for (f in ONT_BED) {
  sample_id <- strsplit(f, '.', fixed=TRUE)[[1]][2]
  sample_id <- gsub('/', '', sample_id, fixed=TRUE)
  
  temp_df <- read.delim(
    f, header=FALSE, 
    col.names=c('scaf', 'pos0', 'pos', 'name', 'score', 'strand',
                'startcodon', 'stopcodon', 'color_rgb', 'cov', 'beta'))
  temp_df$sample_id <- sample_id
  
  ont_df <- rbind(ont_df, temp_df)
}
# tidy up ONT df
# first 11 kb is promoter sequence; make sure transcription starts at +1
ont_df$pos <- ont_df$pos - 11000

# set comparative window to [1, 13332], the transcribed 45S region. both df
# should have data in this range
emseq_df <- emseq_df[emseq_df$pos >= 1 & emseq_df$pos <= 13332, ]
wgbs_df <- wgbs_df[wgbs_df$pos >= 1 & wgbs_df$pos <= 13332, ]
ont_df <- ont_df[ont_df$pos >= 1 & ont_df$pos <= 13332, ]

# set a coverage threshold of 50 to increase confidence in beta values
diag_message('Filter `emseq_df` for positions with coverage >= 50. ',
             'Original nrow: ', nrow(emseq_df), '; ',
             'filtered nrow: ', nrow(emseq_df[emseq_df$cov >= 50, ]))
emseq_df <- emseq_df[emseq_df$cov >= 50, ]

diag_message('Filter `wgbs_df` for positions with coverage >= 50. ',
             'Original nrow: ', nrow(wgbs_df), '; ',
             'filtered nrow: ', nrow(wgbs_df[wgbs_df$cov >= 50, ]))
wgbs_df <- wgbs_df[wgbs_df$cov >= 50, ]

diag_message('Filter `ont_df` for positions with coverage >= 50. ',
             'Original nrow: ', nrow(ont_df), '; ',
             'filtered nrow: ', nrow(ont_df[ont_df$cov >= 50, ]))
ont_df <- ont_df[ont_df$cov >= 50, ]

# calculate mean coverage for each method at this point
diag_message('Mean coverage of remaining positions are: ',
             'EM-seq ', mean(emseq_df$cov), '; ',
             'WGBS ', mean(wgbs_df$cov), '; ',
             'ONT Cas9 ', mean(ont_df$cov))

# only keep relevant columns (keep betas from this point on)
emseq_df <- emseq_df[, c('pos', 'beta', 'sample_id')]
wgbs_df <- wgbs_df[, c('pos', 'beta', 'sample_id')]
ont_df <- ont_df[, c('pos', 'beta', 'sample_id')]

# convert long-to-wide
emseq_wide_df <- reshape(emseq_df, idvar='pos', timevar='sample_id', direction='wide')
colnames(emseq_wide_df) <- gsub('^beta.', '', colnames(emseq_wide_df))
emseq_wide_df <- emseq_wide_df[order(emseq_wide_df$pos), ]

wgbs_wide_df <- reshape(wgbs_df, idvar='pos', timevar='sample_id', direction='wide')
colnames(wgbs_wide_df) <- gsub('^beta.', '', colnames(wgbs_wide_df))
wgbs_wide_df <- wgbs_wide_df[order(wgbs_wide_df$pos), ]

ont_wide_df <- reshape(ont_df, idvar='pos', timevar='sample_id', direction='wide')
colnames(ont_wide_df) <- gsub('^beta.', '', colnames(ont_wide_df))
ont_wide_df <- ont_wide_df[order(ont_wide_df$pos), ]

# check how many rows have NA beta values in ANY of the four datasets (driven
# by insufficient coverage)
diag_message('Number of incomplete rows in `emseq_wide_df`: ',
             sum(!complete.cases(emseq_wide_df)))
diag_message('Number of incomplete rows in `wgbs_wide_df`: ',
             sum(!complete.cases(wgbs_wide_df)))
diag_message('Number of incomplete rows in `ont_wide_df`: ',
             sum(!complete.cases(ont_wide_df)))

# hmm, not a lot. drop positions where methylation levels were NA in ANY of the
# four datasets
emseq_wide_df <- emseq_wide_df[complete.cases(emseq_wide_df), ]
wgbs_wide_df <- wgbs_wide_df[complete.cases(wgbs_wide_df), ]
ont_wide_df <- ont_wide_df[complete.cases(ont_wide_df), ]

# do WGBS datasets have fewer positions with meth data than EM-seq/ONT? check
# overlap of positions with methylation data across three different techniques 
pos_list <- list(`EM-seq`=emseq_wide_df$pos,
                 WGBS=wgbs_wide_df$pos,
                 `ONT Cas9`=ont_wide_df$pos)
plot(venn(pos_list))
# ... yes.

# combine the wide dfs for more diagnostic plots
wide_df <- merge(emseq_wide_df, wgbs_wide_df, by='pos', all=TRUE,
                 suffixes=c('_emseq', '_wgbs'))
wide_df <- merge(wide_df, ont_wide_df, by='pos', all=TRUE)

# select for positions that are present in all 12 datasets
wide_df <- wide_df[complete.cases(wide_df), ]
rownames(wide_df) <- wide_df$pos
wide_df <- wide_df[, -1]

# calculate per-position MEAN methylation levels for every method separately
wide_df$meanER <- rowMeans(wide_df[, grepl('ER$', colnames(wide_df))])
wide_df$meanWR <- rowMeans(wide_df[, grepl('WR$', colnames(wide_df))])
wide_df$meanO <- rowMeans(wide_df[, grepl('O$', colnames(wide_df))])
head(wide_df)

# sanity check: should match the intersection of the three circles in venn diagram
nrow(wide_df)    # expected: 3336

# plot a PCA to visualise overall profile of datasets
set.seed(1337)
pca <- prcomp(t(wide_df), scale=TRUE)
pca_coords <- as.data.frame(pca$x)
eigs <- pca$sdev ^ 2

pca_df <- data.frame(PC1=pca_coords$PC1,
                     PC2=pca_coords$PC2,
                     row.names=rownames(pca_coords))
pca_df$method <- rownames(pca_coords)
pca_df$method <- gsub('^.*ER$', 'EM-seq', pca_df$method)
pca_df$method <- gsub('^.*WR$', 'WGBS', pca_df$method)
pca_df$method <- gsub('^.*O$', 'ONT Cas9', pca_df$method)
pca_df$sample <- gsub('ER$|WR$|O$', '', rownames(pca_df))

g <- ggplot(pca_df, aes(x=PC1, y=PC2, color=method, fill=method, shape=sample)) +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c(EMSEQ_COLOR, ONT_COLOR, WGBS_COLOR)) +
  scale_fill_manual(values=c(EMSEQ_COLOR, ONT_COLOR, WGBS_COLOR)) +
  scale_shape_manual(values=21:25) +
  xlab(paste0('PC1 (', round(eigs[1] / sum(eigs) * 100, 2), '%)')) +
  ylab(paste0('PC2 (', round(eigs[2] / sum(eigs) * 100, 2), '%)')) +
  theme_minimal(12) +
  theme(legend.position='top', legend.box='vertical', legend.margin=margin())

#+ fig.width=6, fig.height=6
print(g)
ggsave('three-way.beta.pca.pdf', width=6, height=6)

# PCA indicates variation is largest across methods, much less variation across
# samples. which is why per-method means were computed. subsequent plots are
# based off these mean values

#+ fig.width=8, fig.height=8
# plot EM-seq vs. WGBS
ggplot(wide_df, aes(x=meanER, y=meanWR)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanWR'), parse=TRUE) +
  ggtitle('Per-position methylation levels, EM-seq vs. WGBS') +
  xlab('EM-seq') +
  ylab('WGBS') +
  theme_minimal(12)

# plot EM-seq vs. ONT
ggplot(wide_df, aes(x=meanER, y=meanO)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanO'), parse=TRUE) +
  ggtitle('Per-position methylation levels, EM-seq vs. ONT Cas9') +
  xlab('EM-seq') +
  ylab('ONT Cas9') +
  theme_minimal(12)

# WGBS vs. ONT
ggplot(wide_df, aes(x=meanWR, y=meanO)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanWR', 'meanO'), parse=TRUE) +
  ggtitle('Per-position methylation levels, WGBS vs. ONT Cas9') +
  xlab('WGBS') +
  ylab('ONT Cas9') +
  theme_minimal(12)

# create a long df for plotting
long_df <- wide_df[, c('meanER', 'meanWR', 'meanO')]
long_df$pos <- as.numeric(rownames(long_df))
long_df <- reshape(long_df, direction='long',
                   idvar='pos',
                   varying=c('meanER', 'meanWR', 'meanO'), v.name='beta',
                   timevar='method', times=c('EM-seq', 'WGBS', 'ONT Cas9'))

# sort by pos, reset row numbering
long_df <- long_df[order(long_df$pos), ]
rownames(long_df) <- NULL
head(long_df)

# calculate GC% in a window of WINDOW_BP
WINDOW_BP <- 101  # this value should be odd, so that the base at midpoint
                  # is sandwiched by equal num of bases upstream and downstream
bp_per_side <- (WINDOW_BP - 1) / 2

# create a string where C/G == 1 and A/T == 0, then
#   GC% = sum(window) / len(window)
rdna_binary <- gsub('C', '1', rdna_seq, ignore.case=TRUE)
rdna_binary <- gsub('G', '1', rdna_binary, ignore.case=TRUE)
rdna_binary <- gsub('[^1]', '0', rdna_binary, ignore.case=TRUE)
gcpct_df <- data.frame(pos=(bp_per_side + 1):(nchar(rdna_seq) - bp_per_side),
                       gcpct=0)
gcpct_df$gcpct <- unlist(lapply(
  gcpct_df$pos, 
  function(x) pct_of_ones(substr(rdna_binary, x - bp_per_side, x + bp_per_side))))

# first 1 kb is promoter sequence; make sure transcription starts at +1
gcpct_df$pos <- gcpct_df$pos - 1000

#+ fig.width=10, fig.height=5
# plot beta across loci
g1 <- ggplot(long_df, aes(x=pos, y=beta, color=method)) +
  geom_point(alpha=0.1) +
  geom_smooth(method='loess', span=0.1) +
  scale_color_manual(values=c(EMSEQ_COLOR, ONT_COLOR, WGBS_COLOR)) +
  coord_cartesian(xlim=c(0, 13332)) +
  ggtitle('Methylation levels and GC% across 45S loci') +
  xlab('Position') +
  ylab('Methylation level (%)') +
  theme_minimal(12) +
  theme(legend.position='top',
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

g2 <- ggplot(gcpct_df, aes(x=pos, y=gcpct)) +
  geom_smooth(method='loess', span=0.1) +
  coord_cartesian(xlim=c(0, 13332)) +
  xlab('Position') +
  ylab('GC%') +
  theme_minimal(12)

plot_grid(g1, g2, ncol=1, align='v', rel_heights=c(0.75, 0.25))
ggsave('three-way.beta_across_loci.pdf', width=10, height=5)

# does GC% influence the discrepancies in beta? re-do scatterplot, but colour
# points by GC%
# 
# merge dataframes together (inner-join), reset row names
wide_df <- merge(wide_df, gcpct_df, by.x=0, by.y='pos', all=FALSE)
colnames(wide_df)[1] <- 'pos'
wide_df$pos <- as.numeric(wide_df$pos)
wide_df <- wide_df[order(wide_df$pos), ]
rownames(wide_df) <- NULL
head(wide_df)

# plot EM-seq vs. WGBS and colour points by GC%
g3 <- ggplot(wide_df, aes(x=meanER, y=meanWR, color=gcpct)) +
  geom_point(alpha=0.4) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(55, 95),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 55), ylim=c(0, 55)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanWR'), parse=TRUE) +
  ggtitle('Per-position methylation levels') +
  xlab('EM-seq (%)') +
  ylab('WGBS (%)') +
  theme_minimal(12) +
  theme(legend.position=c(0.75, 0.1))

# plot EM-seq vs. ONT Cas9 and colour points by GC%
g4 <- ggplot(wide_df, aes(x=meanER, y=meanO, color=gcpct)) +
  geom_point(alpha=0.4) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(55, 95),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 55), ylim=c(0, 55)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanO'), parse=TRUE) +
  ggtitle('') +
  xlab('EM-seq (%)') +
  ylab('ONT Cas9 (%)') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

# plot WGBS vs. ONT Cas9 and colour points by GC%
g5 <- ggplot(wide_df, aes(x=meanWR, y=meanO, color=gcpct)) +
  geom_point(alpha=0.4) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(55, 95),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 55), ylim=c(0, 55)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanWR', 'meanO'), parse=TRUE) +
  ggtitle('') +
  xlab('WGBS (%)') +
  ylab('ONT Cas9 (%)') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

# do, like, proper stats
wide_df$delta_wgbs_emseq <- wide_df$meanWR - wide_df$meanER
summary(lm(wide_df$delta_wgbs_emseq ~ wide_df$gcpct))
g6 <- ggplot(wide_df, aes(x=gcpct, y=delta_wgbs_emseq)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(40, 95), ylim=c(-20, 40)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_wgbs_emseq'), parse=TRUE) +
  ggtitle('Per-position residual methylation levels vs. GC%') +
  xlab('GC%') +
  ylab('Residual WGBS - EM-seq (%)') +
  theme_minimal(12)

wide_df$delta_ont_emseq <- wide_df$meanO - wide_df$meanER
g7 <- ggplot(wide_df, aes(x=gcpct, y=delta_ont_emseq)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(40, 95), ylim=c(-20, 40)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_ont_emseq'), parse=TRUE) +
  ggtitle('') +
  xlab('GC%') +
  ylab('Residual ONT Cas9 - EM-seq (%)') +
  theme_minimal(12)

wide_df$delta_ont_wgbs <- wide_df$meanO - wide_df$meanWR
g8 <- ggplot(wide_df, aes(x=gcpct, y=delta_ont_wgbs)) +
  geom_point(alpha=0.2) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(40, 95), ylim=c(-20, 40)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_ont_wgbs'), parse=TRUE) +
  ggtitle('') +
  xlab('GC%') +
  ylab('Residual ONT Cas9 - WGBS (%)') +
  theme_minimal(12)

#+ fig.width=10, fig.height=12
plot_grid(g3, g6, g4, g7, g5, g8, ncol=2, rel_heights=c(1, 1, 1),
          labels=c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave('three-way.beta_and_gc.pdf', width=10, height=12)

# list deps used in this script
sessionInfo()
