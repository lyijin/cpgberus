rm(list=ls())

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

path_base <- "/home/lie128/csiro/stopwatch/cpgberus/01_txdb"
path_code <- file.path(path_base)

source(file.path(path_code, "annotate_ranges.R"))


gla <- read.csv(file="/home/ros259/VirtualBox_share/GlaI_zscore.csv", stringsAsFactors = F)
gla_gr <- GRanges(seqnames = gla$Chromosome,
                  ranges=IRanges(start=gla$start + 1,
                                 end=gla$end),
                  seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
mcols(gla_gr) <- gla[, 4:ncol(gla)]
#dmrs <- gla_gr

gla_gr <- annotateDMRs(gla_gr, tx=gencode_tx, cpgIslands = cpgIslands, cpgShores = cpgShores)
gla_anno <- data.frame(gla_gr, stringsAsFactors = FALSE)

write.csv(gla_anno, file="/home/ros259/VirtualBox_share/GlaI_zscore_annotated_1based.csv")
