#!/usr/bin/env Rscript

"> methepic_emseq_wgbs_per-pos_beta.R <

Reads in betas from MethylationEPIC arrays, EM-seq and WGBS, and compares them 
on a per-position basis.

For short read tech (EM-seq/WGBS), coverage cutoffs are defined as having
at least 3 of 4 samples from the same method having a coverage of 5 and
above (i.e., EM-seq 3 of 4 with coverage >= 5 && WGBS 3 of 4 with coverage >= 5).
MethylationEPIC readouts are 'analogue' so no 'coverage' to deal with.
" -> doc

setwd('~/csiro/stopwatch/cpgberus/14_methepic_vs_emseq_wgbs/')

# colour codes are from Dark2 panel
EMSEQ_COLOR <- '#1b9e77'
EPIC_COLOR <- '#e7298a'
WGBS_COLOR <- '#7570b3'

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(cowplot)
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
    italic(y) == c + m %.% italic(x)*','~italic(r) == rval*','~italic(p) == pval*','~italic(n) == n_df,
    list(c = format(unname(coef(model)[1]), digits=3),
         m = format(unname(coef(model)[2]), digits=3),
         rval = format(cor(df[[x]], df[[y]]), digits=4),
         pval = format(summary(model)$coefficients[2,4], digits=1),
         n_df = format(nrow(df), big.mark=',')))
  as.character(as.expression(eq))
}

# function to calculate sum/len of a binary string
pct_of_ones <- function(binary_string) {
  int_vector <- as.integer(unlist(strsplit(binary_string, split = "")))
  
  sum(int_vector) / length(int_vector) * 100
}

# load processed EPIC data
load('../02_process_methepic_data/epic_betas_hg38.RData')

# fix EPIC GRanges--all of the current positions refer to the 'C' on the Watson
# strand. this is correct when strand is '+', but not correct when strand is
# '-' (methylated C on Crick would be the G on Watson)
head(epic_gr, 10)
getSeq(BSgenome.Hsapiens.UCSC.hg38, head(epic_gr, 10))
ranges(epic_gr)[strand(epic_gr) == '-'] <- shift(ranges(epic_gr)[strand(epic_gr) == '-'], 1)
head(epic_gr, 10)
getSeq(BSgenome.Hsapiens.UCSC.hg38, head(epic_gr, 10))

# load per-position data from EM-seq and WGBS
cov_df <- read.delim('../04_parse_bismark_covs/rarefied.coverages.tsv.gz')
beta_df <- read.delim('../04_parse_bismark_covs/rarefied.methpcts.tsv.gz')

# convert meth % values (0-100) in `beta_df` to betas (0-1)
beta_df[4:ncol(beta_df)] <- beta_df[4:ncol(beta_df)] / 100
head(beta_df)

# filtering `beta_df` for positions with min 5 coverage in ALL 8 rarefied
# samples was too stringent: initial union of all 55m positions (covered at
# least once in at least one sample) drops to < 1m post-filtering.
#
# requiring min 3 of 4 EM-seq samples && min 3 of 4 WGBS samples produces a
# more workable set of positions, ~7.7m.
#
# get numbers that summarises coverages across EM-seq and WGBS samples
diag_message('Number of positions with min 5 coverage in 4/4 EM-seq samples: ',
             sum(rowSums(cov_df[, grepl('ER$', colnames(cov_df))] >= 5) == 4))
diag_message('Number of positions with min 5 coverage in 4/4 WGBS samples: ',
             sum(rowSums(cov_df[, grepl('WR$', colnames(cov_df))] >= 5) == 4))
diag_message('Number of positions with min 5 coverage in 3/4 EM-seq samples: ',
             sum(rowSums(cov_df[, grepl('ER$', colnames(cov_df))] >= 5) >= 3))
diag_message('Number of positions with min 5 coverage in 3/4 EM-seq samples: ',
             sum(rowSums(cov_df[, grepl('WR$', colnames(cov_df))] >= 5) >= 3))
# ... WGBS not performing too hot here, tsk tsk

# note: for these lines to work, positions in `beta_df` and `cov_df` must have 
#       the same order (guaranteed by upstream Python script). else `beta_df` 
#       will not be sliced properly
diag_message('Filter `beta_df` for positions where 3/4 EM-seq samples AND 3/4 WGBS samples have coverage >= 5. ')
diag_message('Original nrow: ', nrow(beta_df), '; ')
beta_df <- beta_df[rowSums(cov_df[grepl('ER$', colnames(cov_df))] >= 5) >= 3 &
                     rowSums(cov_df[grepl('WR$', colnames(cov_df))] >= 5) >= 3, ]
cov_df <- cov_df[rowSums(cov_df[grepl('ER$', colnames(cov_df))] >= 5) >= 3 &
                   rowSums(cov_df[grepl('WR$', colnames(cov_df))] >= 5) >= 3, ]
diag_message('filtered nrow: ', nrow(beta_df))

diag_message('Mean coverage of remaining positions are: ',
             'EM-seq ', mean(as.matrix(cov_df[grepl('ER$', colnames(cov_df))])), '; ',
             'WGBS ', mean(as.matrix(cov_df[grepl('WR$', colnames(cov_df))])))
# again, WGBS not performing well in terms of coverage

# convert positions from `beta_df` into a GRanges object, so that intersecting
# EPIC data (already in GRanges object) is easy
beta_gr <- GRanges(paste(beta_df$scaffold, beta_df$start_pos, sep=':'))
mcols(beta_gr) <- beta_df[4:ncol(beta_df)]
head(beta_gr)

# aggressively subset both GRanges to the intersection of both
diag_message('# positions in `epic_gr` before intersection: ', length(epic_gr))
epic_gr <- subsetByOverlaps(epic_gr, beta_gr, type='equal', ignore.strand=TRUE)
beta_gr <- subsetByOverlaps(beta_gr, epic_gr, type='equal', ignore.strand=TRUE)
epic_gr <- sort(sortSeqlevels(epic_gr), ignore.strand=TRUE)
beta_gr <- sort(sortSeqlevels(beta_gr), ignore.strand=TRUE)
diag_message('# positions in `epic_gr` after intersection: ', length(epic_gr))

# sanity check after subsetting: lengths are identical, sequence names are in
# same order, as well as start/end values
stopifnot(length(epic_gr) == length(beta_gr))
stopifnot(identical(seqnames(epic_gr), seqnames(beta_gr)))
stopifnot(identical(start(epic_gr), start(beta_gr)))
stopifnot(identical(end(epic_gr), end(beta_gr)))

# this allows for the cbinding of metadata together. store the metadata
# in `epic_gr`, which is better annotated than `beta_gr`
values(epic_gr) <- cbind(values(beta_gr), values(epic_gr))

# transform the metadata into a proper data.frame (not the 'DataFrame' S4 object
# used in GRanges), then calculate per-method mean methylation levels
wide_df <- as.data.frame(values(epic_gr))
wide_df$meanER <- rowMeans(wide_df[, grepl('ER$', colnames(wide_df))], na.rm=TRUE)
wide_df$meanWR <- rowMeans(wide_df[, grepl('WR$', colnames(wide_df))], na.rm=TRUE)
wide_df$meanI <- rowMeans(wide_df[, grepl('I$', colnames(wide_df))])
head(wide_df)

# plot a PCA to visualise overall profile of datasets. note that PCA doesn't
# like dealing with NA, necessitating the use of complete.cases() to remove
# any rows with NA in them
set.seed(1337)
pca <- prcomp(t(wide_df[complete.cases(wide_df), ]), scale=TRUE)
pca_coords <- as.data.frame(pca$x)
eigs <- pca$sdev ^ 2

pca_df <- data.frame(PC1=pca_coords$PC1,
                     PC2=pca_coords$PC2,
                     PC3=pca_coords$PC3,
                     row.names=rownames(pca_coords))
pca_df$method <- rownames(pca_coords)
pca_df$method <- gsub('^.*ER$', 'EM-seq', pca_df$method)
pca_df$method <- gsub('^.*WR$', 'WGBS', pca_df$method)
pca_df$method <- gsub('^.*I$', 'EPIC', pca_df$method)
pca_df$sample <- gsub('ER$|WR$|I$', '', rownames(pca_df))

g1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=method, fill=method, shape=sample)) +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_fill_manual(values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_shape_manual(values=21:25) +
  xlab(paste0('PC1 (', round(eigs[1] / sum(eigs) * 100, 2), '%)')) +
  ylab(paste0('PC2 (', round(eigs[2] / sum(eigs) * 100, 2), '%)')) +
  theme_minimal(12) +
  theme(legend.position='top', legend.box='vertical', legend.margin=margin())
g2 <- ggplot(pca_df, aes(x=PC2, y=PC3, color=method, fill=method, shape=sample)) +
  geom_point(size=3, alpha=0.5) +
  scale_color_manual(values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_fill_manual(values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_shape_manual(values=21:25) +
  xlab(paste0('PC2 (', round(eigs[2] / sum(eigs) * 100, 2), '%)')) +
  ylab(paste0('PC3 (', round(eigs[3] / sum(eigs) * 100, 2), '%)')) +
  theme_minimal(12) +
  theme(legend.position='none', legend.box='vertical', legend.margin=margin())

#+ fig.width=6, fig.height=6
plot_grid(g1, g2, ncol=1, align='v', rel_heights=c(0.55, 0.45))
ggsave('three-way.beta.pca.pdf', width=6, height=6)
# PCA indicates variation is largest across methods (specifically EPIC vs.
# the other two short read methods). There is some variation across samples,
# seen in PC3 (approx PC2 in terms of % variation explained). In PC3, the
# four-sided shapes are on top (WR025) while the three-sided ones are bottom
# (WR069). However, most variation still originate from choice of method,
# hence the calculation of per-method means; subsequent plots are based off
# these mean values

#+ fig.width=8, fig.height=8
# plot EM-seq vs. WGBS
ggplot(wide_df, aes(x=meanER, y=meanWR)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanWR'), parse=TRUE) +
  ggtitle('Per-position beta values, EM-seq vs. WGBS') +
  xlab('EM-seq') +
  ylab('WGBS') +
  theme_minimal(12)

# plot EM-seq vs. EPIC
ggplot(wide_df, aes(x=meanER, y=meanI)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanI'), parse=TRUE) +
  ggtitle('Per-position beta values, EM-seq vs. EPIC') +
  xlab('EM-seq') +
  ylab('EPIC') +
  theme_minimal(12)

# WGBS vs. EPIC
ggplot(wide_df, aes(x=meanWR, y=meanI)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanWR', 'meanI'), parse=TRUE) +
  ggtitle('Per-position beta values, WGBS vs. EPIC') +
  xlab('WGBS') +
  ylab('EPIC') +
  theme_minimal(12)

# does GC% influence the discrepancies in beta? re-do scatterplot, but colour
# points by GC%
#
# calculate GC% in a window of WINDOW_BP
WINDOW_BP <- 101  # this value should be odd, so that the base at midpoint
                  # is sandwiched by equal num of bases upstream and downstream
bp_per_side <- (WINDOW_BP - 1) / 2

# letterFrequency(CG) / letterfrequency(ACGT) deals with edge cases where windows
# contain N or non-ACGT bases, which get ignored
window_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, epic_gr + bp_per_side)
wide_df$gcpct <- as.vector(
  letterFrequency(window_seqs, 'CG') / letterFrequency(window_seqs, 'ACGT') * 100)
head(wide_df)

# check distribution of GC% values, to select a colour scheme for the heatmap
# centered appropriately over the median
quantile(wide_df$gcpct, c(0.05, .1, .5, .9, .95))

# hmm. human genome GC% is 42%, but perhaps EPIC probes are preferentially
# chosen around CpG islands, bumping the GC% up slightly. choose a [30, 70]
# colour scale for subsequent plots
# 
# plot EM-seq vs. WGBS and colour points by GC%
g3 <- ggplot(wide_df, aes(x=meanER, y=meanWR, color=gcpct)) +
  geom_point(size=0.5, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanWR'), parse=TRUE) +
  ggtitle('Per-position beta values') +
  xlab('EM-seq') +
  ylab('WGBS') +
  theme_minimal(12) +
  theme(legend.position=c(0.75, 0.1))

# plot EM-seq vs. EPIC and colour points by GC%
g4 <- ggplot(wide_df, aes(x=meanER, y=meanI, color=gcpct)) +
  geom_point(size=0.5, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanER', 'meanI'), parse=TRUE) +
  ggtitle('') +
  xlab('EM-seq') +
  ylab('EPIC') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

# plot WGBS vs. EPIC and colour points by GC%
g5 <- ggplot(wide_df, aes(x=meanWR, y=meanI, color=gcpct)) +
  geom_point(size=0.5, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'meanWR', 'meanI'), parse=TRUE) +
  ggtitle('') +
  xlab('WGBS') +
  ylab('EPIC') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

# do, like, proper stats
wide_df$delta_wgbs_emseq <- wide_df$meanWR - wide_df$meanER
summary(lm(wide_df$delta_wgbs_emseq ~ wide_df$gcpct))
g6 <- ggplot(wide_df, aes(x=gcpct, y=delta_wgbs_emseq)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(20, 80), ylim=c(-0.5, 0.5)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_wgbs_emseq'), parse=TRUE) +
  ggtitle('Per-position residual beta values vs. GC%') +
  xlab('GC%') +
  ylab('Residual WGBS - EM-seq') +
  theme_minimal(12)

wide_df$delta_ont_emseq <- wide_df$meanI - wide_df$meanER
g7 <- ggplot(wide_df, aes(x=gcpct, y=delta_ont_emseq)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(20, 80), ylim=c(-0.5, 0.5)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_ont_emseq'), parse=TRUE) +
  ggtitle('') +
  xlab('GC%') +
  ylab('Residual EPIC - EM-seq') +
  theme_minimal(12)

wide_df$delta_ont_wgbs <- wide_df$meanI - wide_df$meanWR
g8 <- ggplot(wide_df, aes(x=gcpct, y=delta_ont_wgbs)) +
  geom_point(size=0.5, alpha=0.05) + 
  geom_smooth(method=lm, formula='y ~ x', alpha=0.5) +
  coord_cartesian(xlim=c(20, 80), ylim=c(-0.5, 0.5)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df, 'gcpct', 'delta_ont_wgbs'), parse=TRUE) +
  ggtitle('') +
  xlab('GC%') +
  ylab('Residual EPIC - WGBS') +
  theme_minimal(12)

#+ fig.width=10, fig.height=12
plot_grid(g3, g6, g4, g7, g5, g8, ncol=2, rel_heights=c(1, 1, 1),
          labels=c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave('three-way.beta_and_gc.pdf', width=10, height=12)

# sanity check: confirm that scatterplots, whilst slightly overplotted, are
# visually accurate--replots of these with fewer points (10k) should show
# similar trends
set.seed(42)
wide_df_10k <- wide_df[sample(nrow(wide_df), 10000), ]
g3_10k <- ggplot(wide_df_10k, aes(x=meanER, y=meanWR, color=gcpct)) +
  geom_point(size=1, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df_10k, 'meanER', 'meanWR'), parse=TRUE) +
  ggtitle('Per-position methylation levels') +
  xlab('EM-seq (%)') +
  ylab('WGBS (%)') +
  theme_minimal(12) +
  theme(legend.position=c(0.75, 0.1))

# plot EM-seq vs. EPIC and colour points by GC%
g4_10k <- ggplot(wide_df_10k, aes(x=meanER, y=meanI, color=gcpct)) +
  geom_point(size=1, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df_10k, 'meanER', 'meanI'), parse=TRUE) +
  ggtitle('') +
  xlab('EM-seq (%)') +
  ylab('EPIC (%)') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

# plot WGBS vs. EPIC and colour points by GC%
g5_10k <- ggplot(wide_df_10k, aes(x=meanWR, y=meanI, color=gcpct)) +
  geom_point(size=1, alpha=0.2) + 
  scale_color_distiller('GC%', palette='RdYlBu', limits=c(30, 70),
                        guide=guide_colorbar(direction='horizontal')) +
  geom_line(stat='smooth', method=lm, formula='y ~ x', alpha=0.3) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  annotate('text', x=-Inf, y=Inf, hjust=0, vjust=1,
           label=lm_eqn(wide_df_10k, 'meanWR', 'meanI'), parse=TRUE) +
  ggtitle('') +
  xlab('WGBS (%)') +
  ylab('EPIC (%)') +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.1))

plot_grid(g3, g3_10k, g4, g4_10k, g5, g5_10k, ncol=2, rel_heights=c(1, 1, 1),
          labels=c('A', 'B', 'C', 'D', 'E', 'F'))
# ... uh, looks ok to me?


# hmm. there's some interesting trends in the high GC% region. how many positions
# have high GC% in its immediate context?
diag_message('Positions with GC > 70% in immediate context: ', sum(wide_df$gcpct > 70))
diag_message('Positions with GC > 75% in immediate context: ', sum(wide_df$gcpct > 75))
# don't think analyses using these collection of points would be convincing...

# but could it be that there are more cross-reactive probes amongst the high
# GC (GC > 75%) probes than expected?
high_gc_probes <- rownames(wide_df[wide_df$gcpct > 75, ])
head(high_gc_probes)
cross_hyb_probes <- read.csv('13059_2016_1066_MOESM1_ESM.csv')$X
cross_hyb_probes <- intersect(cross_hyb_probes, rownames(wide_df))
head(cross_hyb_probes)

# numbers for fisher's exact / hypergeometric test
diag_message('There are ', nrow(wide_df), ' EPIC probes (not 850k).')
diag_message('There are ', length(high_gc_probes), ' high GC probes.')
diag_message('There are ', length(cross_hyb_probes),
             ' cross-reactive probes within the universe of EPIC probes.')
diag_message('There are ', sum(high_gc_probes %in% cross_hyb_probes),
             ' high GC & cross-reactive probes.')
diag_message("One-tailed Fisher's exact for P(X >= 20) = ", 
             fisher.test(matrix(c(20, 4785, 215, 98650), nrow=2), alternative='greater')$p.value)
# in plain English, it means that there are significantly MORE cross-reactive
# probes in the high-GC probes.

# list deps used in this script
sessionInfo()
