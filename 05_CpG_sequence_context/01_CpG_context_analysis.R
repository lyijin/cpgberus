#####  Dependencies  #####

library(GenomicRanges)
library(motifStack)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")

#####  Analysis  #####

# Extract first 100000 rows from each Granges for testing, seperate EM-Seq (E) and WGBS (W)
# Subset_covs_grl = endoapply(covs_grl,"[",1:100000)
# Subset_covs_grl = covs_grl[seqnames(covs_grl) == "chr1"]
Subset_covs_grl_E = covs_grl[grep("E$", names(covs_grl))]
Subset_covs_grl_W = covs_grl[grep("W$", names(covs_grl))]

# Find common regions for E and W in GrangesList and subset
Common_intersect_regions_E = Reduce(intersect, Subset_covs_grl_E)
Subset_covs_grl_common_E = endoapply(Subset_covs_grl_E, subsetByOverlaps, Common_intersect_regions_E)
rm(Subset_covs_grl_E)
rm(Common_intersect_regions_E)

Common_intersect_regions_W = Reduce(intersect, Subset_covs_grl_W)
Subset_covs_grl_common_W = endoapply(Subset_covs_grl_W, subsetByOverlaps, Common_intersect_regions_W)
rm(Subset_covs_grl_W)
rm(Common_intersect_regions_W)

# Function to tabulate nucleotide species for each postion of motif
Motif_frequency_table <- function(Granges_object) {

    # Create empty nucleotide frequency table
    cpg_context = Granges_object$cpg_context_nnncgnnn
    largest_motif = max(nchar(cpg_context))

    column_names = unique(unlist(strsplit(paste(cpg_context,collapse=""), "")))
    nucleotide_freq = data.frame(matrix(NA, nrow = 0, ncol = length(column_names)))
    colnames(nucleotide_freq) <- column_names

    # Count nucleotide population for first position of motif, loop + 1 until end of motif.
    for (i in seq(largest_motif)) {
        nucleotide_freq_temp = table(substr(cpg_context, i, i))
        nucleotide_freq_temp = data.frame(t(unclass(nucleotide_freq_temp)))
        
        # Fill in missing nucleotides in nucleotide_freq_temp with 0
        if (identical(sort(column_names), sort(colnames(nucleotide_freq_temp))) == FALSE) {
            columns_to_add = column_names[(!(column_names %in% colnames(nucleotide_freq_temp)))]
            temp_df = data.frame(matrix(0, nrow = 1, ncol = length(columns_to_add)))
            colnames(temp_df) <- columns_to_add
            nucleotide_freq_temp = cbind(nucleotide_freq_temp, temp_df)
        }
        nucleotide_freq = rbind(nucleotide_freq, nucleotide_freq_temp)
    }

    return(data.frame(t(nucleotide_freq)))
}

# Plot Motif
Plot_data <- lapply(Subset_covs_grl_common_W, Motif_frequency_table)

for (data_name in names(Plot_data)) {
    motif <- new("pcm", mat=as.matrix(Plot_data[[data_name]]), name=data_name)
    png(paste0(data_name, "_common.png"), width = 11.69, height = 8.3, units = "in", res = 600)
    plot(motif, ic.scale=FALSE, ylab="probability")
    dev.off()
}

Plot_data <- lapply(Subset_covs_grl_common_E, Motif_frequency_table)

for (data_name in names(Plot_data)) {
    motif <- new("pcm", mat=as.matrix(Plot_data[[data_name]]), name=data_name)
    png(paste0(data_name, "_common.png"), width = 11.69, height = 8.3, units = "in", res = 600)
    plot(motif, ic.scale=FALSE, ylab="probability")
    dev.off()
}

Plot_data <- lapply(covs_grl, Motif_frequency_table)

for (data_name in names(Plot_data)) {
    motif <- new("pcm", mat=as.matrix(Plot_data[[data_name]]), name=data_name)
    png(paste0(data_name, "_all.png"), width = 11.69, height = 8.3, units = "in", res = 600)
    plot(motif, ic.scale=FALSE, ylab="probability")
    dev.off()
}

rm(covs_grl)