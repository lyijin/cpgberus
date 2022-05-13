#####  Dependencies  #####

library(GenomicRanges)

#####  Paths  #####

path_base = getwd()
path_cpg_context = file.path(path_base, "12_cpg_context")
path_cpg_context_data = file.path(path_base, "data", "12_cpg_context")
path_txdb = file.path(path_base, "data", "data")
path_seq = file.path(path_base, "data", "04_parse_bismark_covs")

#####  Functions  #####

source(file.path(path_cpg_context, "annotate_ranges_functions.R"))


#####  Load data  #####

load(file=file.path(path_txdb, "gencode.v38.annotation_2021-08-04.RData"))
load(file=file.path(path_seq, "Rarefied_grch38p13_combined_covs_grl.RData"))

#####  Overlap txdb with data  #####


# Annotate once as CpGs appear in multiple lists.

# Make a set of all the unique CpG coords across the coverage files
cpgs_gr = unique(granges(unlist(Rarefied_covs_grl), use.names=FALSE))


cpgs_anno_gr = annotateCpgs(meth_gr=cpgs_gr, tx=gencode_tx,
                            cpgIslands=cpgIslands, cpgShores=cpgShores,
                            tss_proximal=2000)

save(cpgs_anno_gr, file=file.path(path_cpg_context_data, "Rarefied_covs_cpg_annotations.RData"))

















































x = GRangesList(sapply(Rarefied_covs_grl, head))
class(x)
x = x[1:2]
union(x)

class(x[[1]])
union(x[[1]], x[[2]])


x[1:2]
x[[1:2]]


# There is a union method for CompressedGRangesList
# x="CompressedGRangesList", y="missing"
# (inherited from: x="ANY", y="ANY")
showMethods(GenomicRanges::union)

## Construction with GRangesList():
gr1 <- GRanges("chr2", IRanges(3, 6),
               strand="+", score=5L, GC=0.45)
gr2 <- GRanges(c("chr1", "chr1"), IRanges(c(7,13), width=3),
               strand=c("+", "-"), score=3:4, GC=c(0.3, 0.5))
gr3 <- GRanges(c("chr1", "chr2"), IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)
class(grl)
# Error with CompressedGRangesList: Error in as.vector(x) : no method for coercing this S4 class to a vector
union(unlist(grl))
# But works fine if two GRanges are handed to union
union(grl[[1]], grl[[2]])
# Now convert grl to a 'SimpleGRangesList'
grl <- GRangesList(grl, compress = FALSE)
class(grl)
# Same error with SimpleGRangesList: Error in as.vector(x) : no method for coercing this S4 class to a vector
union(grl)

# Unlist works on a CompressedGRangesList object
x = reduce(unlist(grl[1:2]))
y = union(grl[[1]], grl[[2]])
identical(x, y)

# The same error if I unlist first to a GRanges object
x = union(unlist(grl))


union_cpgs_gr <- union(Rarefied_covs_grl)


