#####  Dependencies  #####

library(GenomicRanges)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")
load("/scratch1/gua020/CpGberus_data/motif_matrices.RData")

#####  Analysis  #####

# # # Function to filter granges object in the following order: 1) Lower quartile discard; 2) meth percentage range; 3) Coverage cutoff's for meth / unmeth
Filter_meth_pct <- function(Granges_object, discard_lower_quartile, meth_pct_keep_range, total_coverage_cutoff) {
    
    temp_grange = Granges_object
    
    if (discard_lower_quartile == TRUE) {
        quartiles = quantile(temp_grange$meth_cov + temp_grange$unmeth_cov)
        lower_quartile = quartiles[["25%"]]
        temp_grange = temp_grange[(temp_grange$meth_cov + temp_grange$unmeth_cov) >= lower_quartile, ]
        print(paste0("Discarded < 25% lower quartile which equates to coverage < ", lower_quartile))
    }
    
    temp_grange = subset(temp_grange, meth_pct >= min(meth_pct_keep_range) & meth_pct <= max(meth_pct_keep_range))
    temp_grange = temp_grange[(temp_grange$meth_cov + temp_grange$unmeth_cov) >= total_coverage_cutoff, ]
    return(temp_grange)
}

# Filter CpGs for quartile, meth and coverage
Subset_covs_grl_all = endoapply(covs_grl, Filter_meth_pct, TRUE, c(0, 100), 0)

# # # Function to tabulate motifs
# # # Additional option to mutliply motif by total coverage
Motif_frequency_table <- function(Granges_object, number_of_granges, motif_coverage_multiplier) {

    counter <<- counter + 1
    print(paste0("Analysing Grange object: ", as.character(counter), " / ", number_of_granges))
    
    # Multiply motif cpg_context by total coverage if TRUE
    if (motif_coverage_multiplier == TRUE) {
        motif_multiplier = Granges_object$meth_cov + Granges_object$unmeth_cov
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = motif_multiplier)
    } else {
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = 1)
    }

    # Sum motifs and return dataframe
    motif_table_sum = aggregate(Motif_multiplier ~ Motif, cpg_context, FUN=sum)
    return(data.frame(motif_table_sum))
}

Plot_density <- function(Motif_list_df, file_name, sample_names, remove_N) {

    if (remove_N == TRUE) {       
        Motif_list_df = lapply(Motif_list_df, function(x) x[!grepl("N", x$Motif), ])
    }
    
    pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    
    for (samp_name in sample_names) {
        Motif_list_df_temp = Motif_list_df[names(Motif_list_df) == samp_name]
        plot(density(Motif_list_df_temp[[samp_name]]$Motif_multiplier))
        plot(density(log2(Motif_list_df_temp[[samp_name]]$Motif_multiplier)))
    }
    
    dev.off()    
}

# Consider normalising by frequency?
counter = 0
Motif_all <- lapply(Subset_covs_grl_all, Motif_frequency_table, length(Subset_covs_grl_all), FALSE)
Plot_density(Motif_all, "No_multiply_density", names(Motif_all), TRUE)

counter = 0
Motif_all_multiply <- lapply(Subset_covs_grl_all, Motif_frequency_table, length(Subset_covs_grl_all), TRUE)
Plot_density(Motif_all, "Multiply_density", names(Motif_all), TRUE)

counter = 0
Motif_all_no_filt <- lapply(covs_grl, Motif_frequency_table, length(covs_grl), FALSE)
Plot_density(Motif_all_no_filt, "No_filt_multiply_density", names(Motif_all_no_filt), TRUE)

# Note: I believe Jason doesn't filter on coverage before summing. 
# pdf("Test.pdf", width = 11.69, height = 8.3)
# plot(density(motifs$all_enrich[, 1]))
# plot(density(log2(motifs$all_enrich[, 1])))
# dev.off()