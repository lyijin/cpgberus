#####  Dependencies  #####

library(GenomicRanges)
library(motifStack)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")

#####  Analysis  #####

# Extract first 100000 rows from each Granges for testing, seperate EM-Seq (E) and WGBS (W)
# Subset_covs_grl = endoapply(covs_grl,"[",1:100000)

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

# # # Function to to filter granges list object for common CpGs
Filter_common_CpGs <- function(Granges_list_object) {
    Common_intersect_regions = Reduce(intersect, Granges_list_object)
    Subset_covs_grl_common = endoapply(Granges_list_object, subsetByOverlaps, Common_intersect_regions)
    return(Subset_covs_grl_common)
}

# Filter CpGs for quartile, meth and coverage
Subset_covs_grl_meth = endoapply(covs_grl, Filter_meth_pct, TRUE, c(100, 100), 0)
Subset_covs_grl_unmeth = endoapply(covs_grl, Filter_meth_pct, TRUE, c(0, 0), 0)
Subset_covs_grl_middle = endoapply(covs_grl, Filter_meth_pct, TRUE, c(40, 60), 0)
Subset_covs_grl_all = endoapply(covs_grl, Filter_meth_pct, TRUE, c(0, 100), 0)

# Subset_covs_grl_meth_common = Filter_common_CpGs(Subset_covs_grl_meth)
# Subset_covs_grl_unmeth_common = Filter_common_CpGs(Subset_covs_grl_unmeth)
# Subset_covs_grl_middle_common = Filter_common_CpGs(Subset_covs_grl_middle)
# Subset_covs_grl_all_common = Filter_common_CpGs(Subset_covs_grl_all)

# # # Function to tabulate nucleotide species for each postion of motif
# # # Additional option to mutliply motif by total coverage
Motif_frequency_table <- function(Granges_object, number_of_granges, motif_coverage_multiplier) {

    counter <<- counter + 1
    print(paste0("Analysing Grange object: ", as.character(counter), " / ", number_of_granges))
    
    # Create empty nucleotide frequency table
    cpg_context = unique(Granges_object$cpg_context_nnncgnnn)
    largest_motif = max(nchar(cpg_context))

    column_names = unique(unlist(strsplit(paste(cpg_context,collapse=""), "")))
    nucleotide_freq = data.frame(matrix(NA, nrow = 0, ncol = length(column_names)))
    colnames(nucleotide_freq) <- column_names
    
    # Multiply motif cpg_context by total coverage if TRUE
    if (motif_coverage_multiplier == TRUE) {
        motif_multiplier = Granges_object$meth_cov + Granges_object$unmeth_cov
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = motif_multiplier)
    } else {
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = 1)
    }

    # Count nucleotide population for first position of all motifs, loop + 1 until end of motif.
    motif_list = cpg_context$Motif
    
    for (i in seq(largest_motif)) {    
        print(paste0("Motif position: ", i))
        
        # Substitue motifs with n nucleotide position, then sum
        motif_n_elements = substr(motif_list, i, i)
        cpg_context$Motif <- motif_n_elements
        motif_table_sum = aggregate(Motif_multiplier ~ Motif, cpg_context, FUN=sum)
        
        # Clean up motif_table_sum
        nucleotide_freq_temp = as.data.frame(t(motif_table_sum))
        colnames(nucleotide_freq_temp) <- as.character(nucleotide_freq_temp["Motif", ])
        nucleotide_freq_temp = nucleotide_freq_temp[!(row.names(nucleotide_freq_temp) %in% "Motif"), , drop = FALSE]
        
        # Fill in missing nucleotides in nucleotide_freq_temp with 0
        if (identical(sort(column_names), sort(colnames(nucleotide_freq_temp))) == FALSE) {
            columns_to_add = column_names[(!(column_names %in% colnames(nucleotide_freq_temp)))]
            temp_df = data.frame(matrix(0, nrow = 1, ncol = length(columns_to_add)))
            colnames(temp_df) <- columns_to_add
            nucleotide_freq_temp = cbind(nucleotide_freq_temp, temp_df)
        }
        nucleotide_freq_temp = lapply(nucleotide_freq_temp, as.integer)
        nucleotide_freq = rbind(nucleotide_freq, nucleotide_freq_temp)
    }

    return(data.frame(t(nucleotide_freq)))
}

# # # Function to convert list of dataframes to motif class pcm
Motif_stack_conversion <- function(Motif_list_df, remove_N) {

    #Remove "N" nucleotide if true
    if (remove_N == TRUE) {
        Motif_list_df = lapply(Motif_list_df, function(x) x[rownames(x)!="N", ])
    }
    
    # Convert dataframe to motif class pcm
    for (data_name in names(Motif_list_df)) {
        Motif_list_df[[data_name]] <- new("pcm", mat=as.matrix(Motif_list_df[[data_name]]), name=data_name)
    }
    
    return(Motif_list_df)
}

# # # Funtion to plot by patient ID
Motif_plot <- function(Motif_list_df, file_name, sample_names, remove_N) {

    pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    
    for (samp_name in sample_names) {
        Motif_list_df_temp = Motif_list_df[names(Motif_list_df) == samp_name]
        plot(Motif_stack_conversion(Motif_list_df_temp, remove_N)[[samp_name]], ic.scale=FALSE, ylab="probability")
    }
    
    plot.new()
    motifStack(Motif_stack_conversion(Motif_list_df, remove_N), layout="tree", ic.scale=FALSE, ylab="probability")
    dev.off()
}

# Motif tables and plots where motif multiplied by coverage
counter = 0
Plot_data_meth <- lapply(Subset_covs_grl_meth, Motif_frequency_table, length(Subset_covs_grl_meth), TRUE)
Motif_plot(Plot_data_meth, "Multiply_Meth", names(Plot_data_meth), TRUE)

counter = 0
Plot_data_unmeth <- lapply(Subset_covs_grl_unmeth, Motif_frequency_table, length(Subset_covs_grl_unmeth), TRUE)
Motif_plot(Plot_data_unmeth, "Multiply_Unmeth", names(Plot_data_unmeth), TRUE)

counter = 0
Plot_data_middle <- lapply(Subset_covs_grl_middle, Motif_frequency_table, length(Subset_covs_grl_middle), TRUE)
Motif_plot(Plot_data_middle, "Multiply_Middle", names(Plot_data_middle), TRUE)

counter = 0
Plot_data_all <- lapply(Subset_covs_grl_all, Motif_frequency_table, length(Subset_covs_grl_all), TRUE)
Motif_plot(Plot_data_all, "Multiply_all", names(Plot_data_all), TRUE)

rm(covs_grl)

# Code to remove rows with 0 in list of dataframe process by Merge_CpGs function
Covs_grl_all_df_common_3 = sapply(names(Covs_grl_all_df_common), function(x) which(Covs_grl_all_df_common[[x]]$N == 0),
                        simplify = FALSE, USE.NAMES = TRUE)
                        
rows_to_remove = unique(unlist(Covs_grl_all_df_common_3))

Covs_grl_all_df_common_3 = sapply(names(Covs_grl_all_df_common), function(x) Covs_grl_all_df_common[[x]][-rows_to_remove, ],
                        simplify = FALSE, USE.NAMES = TRUE)
                        
                        
# counter = 0
# Plot_data_meth_common <- lapply(Subset_covs_grl_meth_common, Motif_frequency_table, length(Subset_covs_grl_meth_common), TRUE)
# Motif_plot(Plot_data_meth_common, "Multiply_Meth_common", names(Plot_data_meth_common), TRUE)

# counter = 0
# Plot_data_unmeth_common <- lapply(Subset_covs_grl_unmeth_common, Motif_frequency_table, length(Subset_covs_grl_unmeth_common), TRUE)
# Motif_plot(Plot_data_unmeth_common, "Multiply_Unmeth_common", names(Plot_data_unmeth_common), TRUE)

# counter = 0
# Plot_data_middle_common <- lapply(Subset_covs_grl_middle_common, Motif_frequency_table, length(Subset_covs_grl_middle_common), TRUE)
# Motif_plot(Plot_data_middle_common, "Multiply_Middle_common", names(Plot_data_middle_common), TRUE)


# # Motif tables and plots where motif not multiplied
# counter = 0
# Plot_data_meth_no_multiply <- lapply(Subset_covs_grl_meth, Motif_frequency_table, length(Subset_covs_grl_meth), FALSE)
# Motif_plot(Plot_data_meth_no_multiply, "Meth", names(Plot_data_meth_no_multiply), TRUE)

# counter = 0
# Plot_data_unmeth_no_multiply <- lapply(Subset_covs_grl_unmeth, Motif_frequency_table, length(Subset_covs_grl_unmeth), FALSE)
# Motif_plot(Plot_data_unmeth_no_multiply, "Unmeth", names(Plot_data_unmeth_no_multiply), TRUE)

# counter = 0
# Plot_data_middle_no_multiply <- lapply(Subset_covs_grl_middle, Motif_frequency_table, length(Subset_covs_grl_middle), FALSE)
# Motif_plot(Plot_data_middle_no_multiply, "Middle", names(Plot_data_middle_no_multiply), TRUE)

# counter = 0
# Plot_data_meth_common_no_multiply <- lapply(Subset_covs_grl_meth_common, Motif_frequency_table, length(Subset_covs_grl_meth_common), FALSE)
# Motif_plot(Plot_data_meth_common_no_multiply, "Meth_common", names(Plot_data_meth_common_no_multiply), TRUE)

# counter = 0
# Plot_data_unmeth_common_no_multiply <- lapply(Subset_covs_grl_unmeth_common, Motif_frequency_table, length(Subset_covs_grl_unmeth_common), FALSE)
# Motif_plot(Plot_data_unmeth_common_no_multiply, "Unmeth_common", names(Plot_data_unmeth_common_no_multiply), TRUE)

# counter = 0
# Plot_data_middle_common_no_multiply <- lapply(Subset_covs_grl_middle_common, Motif_frequency_table, length(Subset_covs_grl_middle_common), FALSE)
# Motif_plot(Plot_data_middle_common_no_multiply, "Middle_common", names(Plot_data_middle_common_no_multiply), TRUE)