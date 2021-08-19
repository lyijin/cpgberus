library(GenomicRanges)


#####  Paths  #####

path_base = getwd()
# Assumes your ./data symlink points to hb-stopwatch/work/cpgberus
path_cpgberus = file.path(path_base, "data")
path_txdb = file.path(path_cpgberus, "01_txdb")
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")


#####  Load  #####

# Load nnncgnnn_universe

universe = read.table(file=file.path(path_cpgberus_covs, "grch38p13.nnncgnnn_universe.tsv"),
                      col.names = c("motif", "freq"))
# Discard motifs with N in them
universe = universe[!grepl("N", universe$motif), ]

# Load coverages
load(file=file.path(path_cpgberus_covs, "grch38p13_combined_covs_grl.RData"))

# Combine meth and unmeth coverage into total coverage
#for(s in names(covs_grl)) {
#  print(s)
#  covs_grl[[s]]$cov = covs_grl[[s]]$meth_cov + covs_grl[[s]]$unmeth_cov
#}


#####  Fill a coverage matrix  #####

# Fill a matrix of counts
mat = matrix(0,
             nrow = nrow(universe),
             ncol = length(covs_grl),
             dimnames = list(universe$motif, names(covs_grl)))

motifs = list(meth=mat,
              unmeth=mat)

for(s in names(covs_grl)) {
  print(s)
  meth_sums = tapply(covs_grl[[s]]$meth_cov, INDEX=factor(covs_grl[[s]]$cpg_context_nnncgnnn), FUN=sum)
  unmeth_sums = tapply(covs_grl[[s]]$unmeth_cov, INDEX=factor(covs_grl[[s]]$cpg_context_nnncgnnn), FUN=sum)
  
  # Discard motifs with N in them
  meth_sums = meth_sums[!grepl("N", names(meth_sums))]
  unmeth_sums = unmeth_sums[!grepl("N", names(unmeth_sums))]
  # Record data for this sample
  motifs$meth[names(meth_sums), s] = as.vector(meth_sums)
  motifs$unmeth[names(unmeth_sums), s] = as.vector(unmeth_sums)
}

motifs$all = motifs$meth + motifs$unmeth

motifs$gc_pct = mat
motifs$beta = mat

for(s in names(covs_grl)) {
  print(s)
  meth_pct = tapply(covs_grl[[s]]$meth_pct, INDEX=factor(covs_grl[[s]]$cpg_context_nnncgnnn), FUN=mean)
  gc_pct = tapply(covs_grl[[s]]$gc_pct, INDEX=factor(covs_grl[[s]]$cpg_context_nnncgnnn), FUN=mean)
  
  # Discard motifs with N in them
  meth_pct = meth_pct[!grepl("N", names(meth_pct))]
  gc_pct = gc_pct[!grepl("N", names(gc_pct))]
  # Record data for this sample
  motifs$beta[names(meth_pct), s] = as.vector(meth_pct) / 100
  motifs$gc_pct[names(gc_pct), s] = as.vector(gc_pct)
}


#####  Normalise matrix  #####

motifs$meth_enrich = mat
motifs$unmeth_enrich = mat
motifs$all_enrich = mat

# Normalise coverages to reflect the commonality of the motif in the genome
for(s in names(covs_grl)) {
  motifs$all_enrich[, s] = motifs$all[, s] / universe$freq
  motifs$meth_enrich[, s] = motifs$meth[, s] / universe$freq
  motifs$unmeth_enrich[, s] = motifs$unmeth[, s] / universe$freq 
}

# Divide the coverage per motif by the library size, expressed as per million reads
for(this in c("all_enrich", "meth_enrich", "unmeth_enrich")) {
  for(s in names(covs_grl)) {
    motifs[[this]][, s] = motifs[[this]][, s] / ( sum(motifs[[this]][, s]) / 1e6)
  }
}

save(motifs, file=file.path(path_cpgberus_covs, "motif_matrices.RData"))
