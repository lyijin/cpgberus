#!/usr/bin/env Rscript

"> methepic_emseq_wgbs_feature_beta.R <

Reads in betas from MethylationEPIC arrays, EM-seq and WGBS, and further
analyses the betas on a per-feature basis.

Studied features are:
1. CpG islands, CpG shores, neither islands/shores
2. Promoter, gene body, intergenic

Coverage cutoffs used to pick better-covered positions in EM-seq and WGBS
datasets are inherited from `methepic_emseq_wgbs_per-pos_beta.R`.
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
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
})

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}

# function to calculate coverage ratios
# courtesy Jason Ross
coverageRatio <- function(query, subject, ratio=TRUE) {
  strand(subject) <- "*"
  subject <- reduce(subject)
  names(query) <- 1:length(query)
  
  # Find the pairs of overlapping ranges
  qs_pairs <- findOverlapPairs(query, subject)
  # Map this back to the original query
  query_map <- as.integer(names(S4Vectors::first(qs_pairs)))
  
  if(ratio) {
    x <- numeric(length=length(query))
    x[query_map] <- width(pintersect(qs_pairs)) / width(S4Vectors::first(qs_pairs))
  } else {
    x <- integer(length=length(query))
    x[query_map] <- width(pintersect(qs_pairs))
  }
  x
}

# function provides annotations for CpGs in the GenomeRanges object
# courtesy Jason Ross
annotateCpgs <- function(meth_gr, tx, cpgIslands, cpgShores, tss_proximal=2000) {
  # want to do most analysis for "just" protein coding as well as all genes
  tx_pc <- tx[tx$transcript_type == "protein_coding"]
  
  numGeneTSS <- function(this_tx) {
    # number of different genes TSSs within the DMR
    tss_ovl <- as.matrix(findOverlaps(meth_gr, resize(this_tx, 1, fix="start")))
    tss_ovl <- tapply(this_tx$transcript_id[tss_ovl[, 2]], tss_ovl[, 1], function(x) length(unique(x)))
    nGeneTSS <- integer(length(meth_gr))
    nGeneTSS[as.integer(names(tss_ovl))] <- unname(tss_ovl)
    nGeneTSS
  }
  
  distTSS <- function(this_tx) {
    # distance to closest TSS
    reg.dist <- as.data.frame(distanceToNearest(meth_gr,
                                                resize(this_tx, 1,
                                                       fix="start"))
    )
    i <- reg.dist$subjectHits
    df <- data.frame(distanceTSS = reg.dist$distance,
                     tx_name = this_tx$tx_name[i],
                     transcript_type = this_tx$transcript_type[i],
                     gene_id = this_tx$gene_id[i],
                     gene_name = this_tx$gene_name[i],
                     stringsAsFactors = FALSE
    )
    df
  }
  
  meth_gr$nGeneTSS <- numGeneTSS(tx)
  meth_gr$nProtGeneTSS <- numGeneTSS(tx_pc)
  
  dist_df <- distTSS(tx)
  mcols(meth_gr) <- cbind(mcols(meth_gr), dist_df)
  
  distpc_df <- distTSS(tx_pc)
  names(distpc_df) <- paste(names(distpc_df), "prot", sep="_")
  mcols(meth_gr) <- cbind(mcols(meth_gr), distpc_df[, c(1:2, 4:5)])
  
  # distance to closest CpG island
  reg.dist <- distanceToNearest(meth_gr, cpgIslands)
  meth_gr$distanceCpGi <- NA
  meth_gr$CpGi <- NA
  meth_gr$distanceCpGi[queryHits(reg.dist)] <- mcols(reg.dist)$distance
  meth_gr$CpGi[queryHits(reg.dist)] <- subjectHits(reg.dist)
  
  # % promoter/genebody/intergenic
  tss_gr <- resize(tx, 1, fix="start")
  promoters_gr <- reduce(resize(tss_gr, tss_proximal * 2, fix="center"))
  
  genebody_gr <- setdiff(reduce(tx), promoters_gr)
  genome_gr <- GRanges(seqlevels(tx), IRanges(1, seqlengths(tx)))
  genic_gr <- union(promoters_gr, tx)
  strand(genic_gr) <- "*"
  intergenic_gr <- setdiff(genome_gr, genic_gr)
  
  meth_gr$promoter <- coverageRatio(meth_gr, promoters_gr)
  meth_gr$genebody <- coverageRatio(meth_gr, genebody_gr)
  meth_gr$intergenic <- coverageRatio(meth_gr, intergenic_gr)
  
  # % CpG island/CpG shore/nonCpG
  meth_gr$CpGisland <- coverageRatio(meth_gr, cpgIslands)
  meth_gr$CpGshores <- coverageRatio(meth_gr, cpgShores)
  cpgSet <- union(cpgIslands, cpgShores)
  meth_gr$nonCpG <- coverageRatio(meth_gr, setdiff(genome_gr, cpgSet))
  
  meth_gr
}

# function to wide-to-long convert a sliced df, and slap it with `feature_label`
genomic_feature_as_long_df <- function(sliced_wide_df, feature_label) {
  # last column stores GC% of positions in the sliced df. calculate mean and
  # then discard the column so it doesn't get reshaped
  mean_gcpct <- mean(sliced_wide_df$gcpct)
  sliced_wide_df <- sliced_wide_df[-ncol(sliced_wide_df)]
  
  long_df <- reshape(
    sliced_wide_df, varying=1:3, v.names='beta',
    timevar='library_type',
    times=colnames(sliced_wide_df[1:3]),
    direction='long')
  long_df <- long_df[-3]
  long_df$genomic_feature <- paste0(
    feature_label,
    '|(n = ', format(nrow(sliced_wide_df), big.mark=','), ',',
    '|GC% = ', format(mean_gcpct, digits=3), '%)')
  
  long_df
}

# function to calculate corr values within a sliced df
corr_values_as_long_df <- function(sliced_wide_df, feature_label) {
  long_df <- data.frame(
    genomic_feature=feature_label,
    corr_value=c(cor(sliced_wide_df[[1]], sliced_wide_df[[2]], method='pearson'),
                 cor(sliced_wide_df[[1]], sliced_wide_df[[3]], method='pearson'),
                 cor(sliced_wide_df[[2]], sliced_wide_df[[3]], method='pearson')),
    type=c(paste(colnames(sliced_wide_df)[1], 'v', colnames(sliced_wide_df)[2]),
           paste(colnames(sliced_wide_df)[1], 'v', colnames(sliced_wide_df)[3]),
           paste(colnames(sliced_wide_df)[2], 'v', colnames(sliced_wide_df)[3]))
  )
  
  long_df
}

# load pre-processed genomic annotations from `01_txdb/`, saved in `data`
load('../data/gencode.v38.annotation_2021-08-04.RData')

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

# from this point on, script diverges from the per-position analysis that
# "inspired" this script

# transform the metadata into a proper data.frame (not the 'DataFrame' S4 object
# used in GRanges), then calculate per-method mean methylation levels
wide_df <- as.data.frame(values(epic_gr))
wide_df$`EM-seq` <- rowMeans(wide_df[, grepl('ER$', colnames(wide_df))], na.rm=TRUE)
wide_df$WGBS <- rowMeans(wide_df[, grepl('WR$', colnames(wide_df))], na.rm=TRUE)
wide_df$EPIC <- rowMeans(wide_df[, grepl('I$', colnames(wide_df))])

# only keep means in wide_df
wide_df <- wide_df[, !grepl('^WR', colnames(wide_df))]
head(wide_df)

# calculate GC% in a window of WINDOW_BP for positions in wide_df, so that
# mean GC% can be calculated for each genomic context later
WINDOW_BP <- 101  # this value should be odd, so that the base at midpoint
# is sandwiched by equal num of bases upstream and downstream
bp_per_side <- (WINDOW_BP - 1) / 2

# letterFrequency(CG) / letterfrequency(ACGT) deals with edge cases where windows
# contain N or non-ACGT bases, which get ignored
window_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, epic_gr + bp_per_side)
wide_df$gcpct <- as.vector(
  letterFrequency(window_seqs, 'CG') / letterFrequency(window_seqs, 'ACGT') * 100)
head(wide_df)

# annotate GRanges object with genomic annotation
epic_anno_gr <- annotateCpgs(
  meth_gr=epic_gr, tx=gencode_tx, cpgIslands=cpgIslands, 
  cpgShores=cpgShores, tss_proximal=2000)

# create a df to plot beta trends across different genomic contexts
plot_df <- rbind(
  genomic_feature_as_long_df(wide_df[epic_anno_gr$CpGisland == 1, ], 'CpG island'),
  genomic_feature_as_long_df(wide_df[epic_anno_gr$CpGshores == 1, ], 'CpG shores'),
  genomic_feature_as_long_df(wide_df[epic_anno_gr$nonCpG == 1, ], 'Non-island/shore'),
  genomic_feature_as_long_df(wide_df[epic_anno_gr$promoter == 1, ], 'Promoter'),
  genomic_feature_as_long_df(wide_df[epic_anno_gr$genebody == 1, ], 'Gene body'),
  genomic_feature_as_long_df(wide_df[epic_anno_gr$intergenic == 1, ], 'Intergenic')
)
rownames(plot_df) <- NULL

corr_df <- rbind(
  corr_values_as_long_df(wide_df[epic_anno_gr$CpGisland == 1, ], 'CpG island'),
  corr_values_as_long_df(wide_df[epic_anno_gr$CpGshores == 1, ], 'CpG shores'),
  corr_values_as_long_df(wide_df[epic_anno_gr$nonCpG == 1, ], 'Non-island/shore'),
  corr_values_as_long_df(wide_df[epic_anno_gr$promoter == 1, ], 'Promoter'),
  corr_values_as_long_df(wide_df[epic_anno_gr$genebody == 1, ], 'Gene body'),
  corr_values_as_long_df(wide_df[epic_anno_gr$intergenic == 1, ], 'Intergenic')
)
corr_df

g1 <- ggplot(plot_df, aes(x=beta, y=fct_inorder(genomic_feature))) +
  geom_boxplot(aes(color=library_type, fill=library_type), alpha=0.5, width=0.85, outlier.shape=NA) +
  scale_color_manual('Library type', values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_fill_manual('Library type', values=c(EMSEQ_COLOR, EPIC_COLOR, WGBS_COLOR)) +
  scale_x_continuous(expression(beta), position='top') +
  scale_y_discrete('', labels=function(x) gsub('|', '\n', x, fixed=TRUE), limits=rev) +
  theme_classic(12) +
  theme(legend.position='top', 
        axis.line.y=element_blank(), axis.ticks.y=element_blank())

g2 <- ggplot(corr_df, aes(x=type, y=fct_inorder(genomic_feature))) +
  geom_text(data=corr_df, aes(label=round(corr_value, 3))) +
  scale_x_discrete(expression(paste('Pearson ', italic(r))), labels=function(x) sub(' v ', ' v\n', x, fixed=TRUE), position='top') +
  scale_y_discrete('', limits=rev) +
  theme_classic(12) +
  theme(legend.position='top', 
        axis.text.y=element_blank(), axis.line.y=element_blank(),
        axis.ticks.y=element_blank())

#+ fig.width=8, fig.height=5
plot_grid(g1, g2, ncol=2, align='h', axis='bt', rel_widths=c(0.6, 0.4))
ggsave('three-way.beta.genomic_features.pdf', width=8, height=5)

# list deps used in this script
sessionInfo()
