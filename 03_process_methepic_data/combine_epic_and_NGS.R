#####  Dependencies  #####

library(GenomicRanges)
library(data.table)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)


#####  Paths  #####

path_base = getwd()
#path_cpgberus = "/media/hb-stopwatch/work/cpgberus"
#path_cpgberus_covs = file.path(path_cpgberus, "data", "04_parse_bismark_covs")
#path_cpgberus_epic = file.path(path_cpgberus, "data", "04_parse_methylationEPIC")
path_cpgberus_covs = file.path(path_base, "data", "04_parse_bismark_covs")
path_cpgberus_epic = file.path(path_base, "data", "04_parse_methylationEPIC")


#####  Load data  #####

#load(file=file.path(path_cpgberus_covs, "grch38p13_combined_covs_grl.RData"))
load(file=file.path(path_cpgberus_covs, "Rarefied_grch38p13_combined_covs_grl.RData"))
load(file=file.path(path_cpgberus_epic, "epic_betas_hg38.RData"))

#####  Find and report overlaps  #####

samples = names(Rarefied_covs_grl)

Rarefied_covs_epic = lapply(samples, function(s) {
  # For EPIC arrays, flip the ending E or W on the identifier to I
  # An optional R character might be there for rarified
  epic_s = sub("(\\w+)[EW]R?", "\\1I", s)
  # Find overlapping CpG sites
  # query is NGS, subject is EPIC
  hits = findOverlaps(Rarefied_covs_grl[[s]], epic_gr)
  
  x = Rarefied_covs_grl[[s]][queryHits(hits), ]
  #queryHits(hits)
  mcols(x)$EPIC_pct = mcols(epic_gr[subjectHits(hits)])[, epic_s] * 100
  mcols(x)$IlmnID = names(epic_gr[subjectHits(hits)])
  mcols(x)$Type = manifest[mcols(x)$IlmnID, "Infinium_Design_Type"]
  x
})
names(Rarefied_covs_epic) = samples

# Form a GRangesList object
Rarefied_covs_epic = do.call("GRangesList", Rarefied_covs_epic)

save(Rarefied_covs_epic, file=file.path(path_cpgberus_covs, "Rarefied_grch38p13_combined_covs_epic-overlap_grl.RData"))

