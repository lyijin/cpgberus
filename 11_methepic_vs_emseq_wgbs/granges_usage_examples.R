#####  Dependencies  #####

library(GenomicRanges)
library(data.table)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)


#####  Paths  #####

path_base = getwd()
path_cpgberus = "/media/hb-stopwatch/work/cpgberus"
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")
path_cpgberus_epic = file.path(path_cpgberus, "04_parse_methylationEPIC")

#####  Load data  #####

load(file=file.path(path_cpgberus_covs, "grch38p13_combined_covs_grl.RData"))
load(file=file.path(path_cpgberus_epic, "epic_betas_hg38.RData"))

#####  Find and report overlaps  #####


samples = names(covs_grl)

epic_gr



s = samples[1]
epic_s = substr(s, 1, 7)



hits = findOverlaps(covs_grl[[s]], epic_gr)

# query is seq

x = covs_grl[[s]][queryHits(hits), ]
y = epic_gr[subjectHits(hits), epic_s]


# subject is epic
epic_gr[, s]


sort(epic_gr[seqnames(epic_gr) == "chrY"])

my_gr = GRanges(seqnames="chr1", ranges=IRanges(start=10840, end=10855))
subsetByOverlaps(covs_grl[[1]], my_gr)




covs_grl[[1]][seqnames(covs_grl[[1]]) == "chrY"]
