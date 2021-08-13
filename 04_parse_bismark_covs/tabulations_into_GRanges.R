rm(list = ls())

#####  Dependencies  #####

library(GenomicRanges)
library(data.table)
library(GenomeInfoDb)


#####  Paths  #####

path_base = getwd()
path_cpgberus = "/media/hb-stopwatch/work/cpgberus"
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")


#####  Load and process data  #####

hg38 = Seqinfo(genome="hg38")

covs = list()
for(f in list.files(path_cpgberus_covs, pattern = "combined.tsv.gz")) {
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
covs_grl = do.call("GRangesList", covs)

save(covs_grl, file=file.path(path_cpgberus_covs, "grch38p13_combined_covs_grl.RData"))


q("no")

