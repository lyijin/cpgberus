rm(list = ls())

#####  Dependencies  #####

library(GenomicRanges)
library(data.table)
library(GenomeInfoDb)


#####  Paths  #####

# path_base = getwd()
# Assumes your ./data symlink points to hb-stopwatch/work/cpgberus
# path_cpgberus = file.path(path_base, "data")
path_cpgberus = "/datasets/work/hb-stopwatch/work/cpgberus"
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")


#####  Load and process data  #####

# HPC clusters has issues with Seqinfo(genome="hg38") because: 
# 1) Connecting to the internet because of whitelisting
# 2) Outdated package. 
# Instead, generate a RData file with the Seqinfo required and load it.

# hg38 = Seqinfo(genome="hg38")
load(file.path(path_cpgberus, "04_parse_bismark_covs", "Seqinfo_hg38.RData"))

#####  Create matrix for rarefied data #####
############################################

rarefied_files = list.files(path_cpgberus_covs, pattern = "R.grch38p13_lambda_puc.combined.tsv.gz")
covs = list()

for(f in rarefied_files) {
  sample_name = sub("(\\w+)\\..*", "\\1", f)
  print(sample_name)
  this_cov = fread(file.path(path_cpgberus_covs, f), data.table = FALSE)
  # The coverage file contains lambda and pUC19_NEB genomes, discard them
  this_cov = this_cov[!this_cov$scaffold %in% c("lambda", "pUC19_NEB"), ]
  
  cov_gr = GRanges(seqnames = this_cov$scaffold,
                   ranges=IRanges(start = this_cov$dinuc_start, end = this_cov$dinuc_end),
                   strand = "*")
  mcols(cov_gr) = this_cov[, 4:10]
  seqinfo(cov_gr) = hg38[seqlevels(cov_gr)]
  covs[[sample_name]] = cov_gr
}
rm(cov_gr, this_cov, sample_name)

# Form a GRangesList object
Rarefied_covs_grl = do.call("GRangesList", covs)

save(Rarefied_covs_grl, file=file.path(path_cpgberus_covs, "Rarefied_grch38p13_combined_covs_grl.RData"))

rm(Rarefied_covs_grl)


#####  Create matrix for not rarefied data #####
################################################

not_rarefied_files = list.files(path_cpgberus_covs, pattern = "combined.tsv.gz")
not_rarefied_files = not_rarefied_files[!(not_rarefied_files %in% rarefied_files)]
covs = list()

for(f in not_rarefied_files) {
  sample_name = sub("(\\w+)\\..*", "\\1", f)
  print(sample_name)
  this_cov = fread(file.path(path_cpgberus_covs, f), data.table = FALSE)
  # The coverage file contains lambda and pUC19_NEB genomes, discard them
  this_cov = this_cov[!this_cov$scaffold %in% c("lambda", "pUC19_NEB"), ]
  
  cov_gr = GRanges(seqnames = this_cov$scaffold,
                   ranges=IRanges(start = this_cov$dinuc_start, end = this_cov$dinuc_end),
                   strand = "*")
  mcols(cov_gr) = this_cov[, 4:10]
  seqinfo(cov_gr) = hg38[seqlevels(cov_gr)]
  covs[[sample_name]] = cov_gr
}
rm(cov_gr, this_cov, sample_name)

# Form a GRangesList object
Not_rarefied_covs_grl = do.call("GRangesList", covs)

save(Not_rarefied_covs_grl, file=file.path(path_cpgberus_covs, "Not_rarefied_grch38p13_combined_covs_grl.RData"))

q("no")