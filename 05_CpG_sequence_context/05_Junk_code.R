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

###################################################################
# # Function to calculate the delta between samples and groups for:
# # 1) Coverage evenness
# # 2) Deviation in absolute methylation levels
# # Above are calculated in Yi Jin's python script.
###################################################################

Calculate_delta_evenness_meth_pct <- function(List_of_dataframes_common, biological_sample_names, reference_sample_names) {
    
    samp_id_list = names(List_of_dataframes_common)
    
    # QC to check if "chr" and "pos" columns are identical
    for (samp_id in samp_id_list) {
        if (identical(List_of_dataframes_common[[samp_id_list[1]]]$chr, List_of_dataframes_common[[samp_id]]$chr) == FALSE) {
            stop(paste0("Chromosomes not identical between samples: ", samp_id_list[1], " ", samp_id))
        } else if (identical(List_of_dataframes_common[[samp_id_list[1]]]$pos, List_of_dataframes_common[[samp_id]]$pos) == FALSE) {
            stop(paste0("CpG position not identical between samples: ", samp_id_list[1], " ", samp_id))
        }
    }
    
    # Calculate delta evenness and meth pct between groups
    mean_reference = List_of_dataframes_common[reference_sample_names]
    mean_reference = lapply(mean_reference, function(x) data.frame(group_delta_evenness = x$evenness))
    mean_reference = rowMeans(do.call(cbind, mean_reference))
    
    mean_condition = List_of_dataframes_common[samp_id_list[!(samp_id_list %in% reference_sample_names)]]
    mean_condition = lapply(mean_condition, function(x) data.frame(group_delta_evenness = x$evenness))
    mean_condition = rowMeans(do.call(cbind, mean_condition))
    
    evenness_diff = mean_reference - mean_condition
    
    mean_reference = List_of_dataframes_common[reference_sample_names]
    mean_reference = lapply(mean_reference, function(x) data.frame(group_abs_delta_meth_pct = x$abs_delta_meth_pct))
    mean_reference = rowMeans(do.call(cbind, mean_reference))
    
    mean_condition = List_of_dataframes_common[samp_id_list[!(samp_id_list %in% reference_sample_names)]]
    mean_condition = lapply(mean_condition, function(x) data.frame(group_abs_delta_meth_pct = x$abs_delta_meth_pct))
    mean_condition = rowMeans(do.call(cbind, mean_condition))
    
    meth_pct_diff = mean_reference - mean_condition
    
    # Calculate delta evenness and meth pct between samples
    for (samp_id in biological_sample_names) {
        samples_to_analyse = grep(samp_id, samp_id_list, value = TRUE)
        control_sample = grep(paste0(reference_sample_names, collapse = "|"), samples_to_analyse, value = TRUE)
        
        for (samp_id_2 in samples_to_analyse) {
            delta_evenness = List_of_dataframes_common[[control_sample]]$evenness - List_of_dataframes_common[[samp_id_2]]$evenness
            delta_meth_pct_sample = List_of_dataframes_common[[control_sample]]$abs_delta_meth_pct - List_of_dataframes_common[[samp_id_2]]$abs_delta_meth_pct
            List_of_dataframes_common[[samp_id_2]] = cbind(List_of_dataframes_common[[samp_id_2]], data.frame(sample_delta_evenness = delta_evenness, sample_abs_delta_meth_pct = delta_meth_pct_sample,
                                                            group_delta_evenness = evenness_diff, group_abs_delta_meth_pct = meth_pct_diff))
        }
    }
    
    return(List_of_dataframes_common)
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






# # Add mu1 and mu2 to metadata

# meta_data = endoapply(Covs_grl_common_bsseq_list_stats_no_smooth, subsetByOverlaps, covs_grl_sig_not_smoothed_common)



# Significant_regions = do.call("rbind", Covs_grl_common_bsseq_list_stats)
# Significant_regions = Significant_regions[(Significant_regions$fdr <= 0.05), ]
# Significant_regions = GRanges(seqnames = Significant_regions$chr,
    # IRanges(start = Significant_regions$pos,
    # end = (Significant_regions$pos + 1)))

# covs_grl_sig_smoothed_common = endoapply(covs_grl, subsetByOverlaps, Significant_regions)
# #rm(covs_grl)

# # Plot PCA graphs unsmoothed
# temp_df = sapply(names(covs_grl_sig_not_smoothed_common), function(x) data.frame(chr = seqnames(covs_grl_sig_not_smoothed_common[[x]]), pos = start(covs_grl_sig_not_smoothed_common[[x]]),
                        # N = (covs_grl_sig_not_smoothed_common[[x]]$meth_cov + covs_grl_sig_not_smoothed_common[[x]]$unmeth_cov), X = covs_grl_sig_not_smoothed_common[[x]]$meth_cov, 
                        # evenness = covs_grl_sig_not_smoothed_common[[x]]$evenness, abs_delta_meth_pct = covs_grl_sig_not_smoothed_common[[x]]$abs_delta_meth_pct, sample_name = x),
                        # simplify = FALSE, USE.NAMES = TRUE)

# # Find rows to discard based on coverage                        
# Rows_to_keep = List_df_to_matrix(temp_df, "N")
# Rows_to_keep = (rowSums(Rows_to_keep <= 15) == 0)
# Rows_to_keep = names(Rows_to_keep)[Rows_to_keep]
                        
# matrix = List_df_to_matrix(temp_df, "evenness")
# Abs_delta_meth_pct_matrix = List_df_to_matrix(temp_df, "abs_delta_meth_pct")
# Coverage_matrix = List_df_to_matrix(temp_df, "N")

# matrix = matrix[Rows_to_keep, ]
# Abs_delta_meth_pct_matrix = Abs_delta_meth_pct_matrix[Rows_to_keep, ]
# Coverage_matrix = Coverage_matrix[Rows_to_keep, ]

# colnames(matrix) = paste0("Evenness", "_", colnames(matrix))
# colnames(Abs_delta_meth_pct_matrix) = paste0("Abs_delta_meth", "_", colnames(Abs_delta_meth_pct_matrix))
# colnames(Coverage_matrix) = paste0("Coverage", "_", colnames(Coverage_matrix))

# Combined_matrix = cbind(matrix, Abs_delta_meth_pct_matrix, Coverage_matrix)
# List_of_matrices_plot = list(Evenness = matrix, Abs_delta_meth_pct = Abs_delta_meth_pct_matrix, Coverage = Coverage_matrix, Combined = Combined_matrix)









# res.pca <- prcomp(matrix, center = TRUE, scale = TRUE)
# res_2.pca <- prcomp(Abs_delta_meth_pct_matrix, center = TRUE, scale = TRUE)
# res_3.pca <- prcomp(Coverage_matrix, center = TRUE, scale = TRUE)
# res_4.pca <- prcomp(Combined_matrix, center = TRUE, scale = TRUE)

# pdf("Test_all.pdf", width = 11.69, height = 8.3)

# plot(fviz_eig(res.pca))
# plot(fviz_pca_biplot(res.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))

# plot(fviz_eig(res_2.pca))
# plot(fviz_pca_biplot(res_2.pca , axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res_2.pca , axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))

# plot(fviz_eig(res_3.pca))
# plot(fviz_pca_biplot(res_3.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res_3.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))

# plot(fviz_eig(res_4.pca))
# plot(fviz_pca_biplot(res_4.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res_4.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))

# dev.off()

# # Plot evenness and delta and coverage scatterplots, full and subset. Include PCA for CpG context.






# Covs_grl_all_df_common = sapply(names(Covs_grl_all_df_common), function(x) data.frame(chr = seqnames(Covs_grl_all_df_common[[x]]), pos = start(Covs_grl_all_df_common[[x]]),
                        # N = (Covs_grl_all_df_common[[x]]$meth_cov + Covs_grl_all_df_common[[x]]$unmeth_cov), X = Covs_grl_all_df_common[[x]]$meth_cov, 
                        # evenness = Covs_grl_all_df_common[[x]]$evenness, abs_delta_meth_pct = Covs_grl_all_df_common[[x]]$abs_delta_meth_pct, sample_name = x),
                        # simplify = FALSE, USE.NAMES = TRUE)

# Covs_grl_all_df_common = Calculate_delta_meth_pct(Covs_grl_all_df_common, c("WR025V1", "WR025V9", "WR069V1", "WR069V9"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"))

# save(Covs_grl_all_df_common, file = "Strand_bias.RData")


# # Create directroy to store plots
# dir.create("Strand_bias_plots_not_subsample")


# #############################

# # Convert to list of dataframe to matrix of variable of interest
# matrix = List_df_to_matrix(Covs_grl_all_df_common, "evenness")
# Abs_delta_meth_pct_matrix = List_df_to_matrix(Covs_grl_all_df_common, "abs_delta_meth_pct")
# Coverage_matrix = List_df_to_matrix(Covs_grl_all_df_common, "N")

# matrix = Filter_matrix_NA(matrix, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), 
                                    # c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), 2, 2)

# Abs_delta_meth_pct_matrix = Filter_matrix_NA(Abs_delta_meth_pct_matrix, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), 
                                    # c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), 2, 2)

# Coverage_matrix = head([Coverage_matrix == 0])
# Coverage_matrix = Filter_matrix_NA(Coverage_matrix, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), 
                                    # c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), 2, 2)
                                    
# Abs_delta_meth_pct_matrix = Filter_matrix_NA(Abs_delta_meth_pct_matrix, "abs_delta_meth_pct")
# Coverage_matrix = Filter_matrix_NA(Coverage_matrix, "N")

# # Extract top 3000 variable CpGs based on evenness 
# # ind = order(rowVars(matrix, na.rm = TRUE), decreasing = TRUE)[1:5000]
# # matrix = matrix[ind, ]
# # Abs_delta_meth_pct_matrix = Abs_delta_meth_pct_matrix[ind, ]
# # Coverage_matrix = Coverage_matrix[ind, ]

# colnames(matrix) = paste0("Evenness", "_", colnames(matrix))
# colnames(Abs_delta_meth_pct_matrix) = paste0("Abs_delta_meth", "_", colnames(Abs_delta_meth_pct_matrix))
# colnames(Coverage_matrix) = paste0("Coverage", "_", colnames(Coverage_matrix))

# Combined_matrix = cbind(matrix, Abs_delta_meth_pct_matrix, Coverage_matrix)

# # PCA
# # Create PCA of all values, show that the greatest seperation for evenness and Delta is sequecing type, while its patient for coverage
# # Create individual heatmaps (5000 most variable regions) with GC percentage + small coverage matrix?

# res.pca <- prcomp(matrix, center = TRUE, scale = TRUE)
# res_2.pca <- prcomp(Abs_delta_meth_pct_matrix, center = TRUE, scale = TRUE)
# res_3.pca <- prcomp(Coverage_matrix, center = TRUE, scale = TRUE)
# #res_4.pca <- prcomp(Combined_matrix, center = TRUE, scale = TRUE)

# # png("Test.png", width = 11.69, height = 11.69, units = "in", res = 600)
# # print(ggpairs(data.frame(Combined_matrix), title = "DEpiBurn clinical correlations", upper = list(continuous = wrap("cor", size = 1, method = "pearson")), lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.1))) + theme(strip.placement = "outside", text = element_text(size = 2)))
# # dev.off()

# pdf("Test_all.pdf", width = 11.69, height = 8.3)

# plot(fviz_eig(res.pca))
# plot(fviz_pca_biplot(res.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# # print(Heatmap(matrix, name = "Evenness", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE))

# plot(fviz_eig(res_2.pca))
# plot(fviz_pca_biplot(res_2.pca , axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res_2.pca , axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# # print(Heatmap(Abs_delta_meth_pct_matrix, name = "Abs_delta_meth_pct", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE))

# plot(fviz_eig(res_3.pca))
# plot(fviz_pca_biplot(res_3.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# #plot(fviz_pca_var(res_3.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# # print(Heatmap(Coverage_matrix, name = "Coverage", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE))

# # plot(fviz_eig(res_4.pca))
# # plot(fviz_pca_biplot(res_4.pca, axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# # #plot(fviz_pca_var(res_4.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
# # print(Heatmap(Combined_matrix, name = "Evenness and Abs_delta_meth_pct", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE))

# dev.off()


# # Make a heatmap here, test code.
# heatmap_1 = Heatmap(matrix, name = "Evenness", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE, col = colorRamp2(c(min(matrix), max(matrix)), c("white", "orange")))
# heatmap_2 = Heatmap(Abs_delta_meth_pct_matrix, name = "Abs_delta_meth_pct", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, col = colorRamp2(c(min(Abs_delta_meth_pct_matrix), max(Abs_delta_meth_pct_matrix)), c("white", "orange")))
# heatmap_3 = Heatmap(Coverage_matrix, name = "Coverage", column_title = "Samples", row_title = "CpG position", show_row_names = FALSE, show_row_dend = FALSE, col = colorRamp2(c(min(Coverage_matrix), max(Coverage_matrix)), c("white", "orange")))

# png("Test.png", width = 11.69, height = 8.3, units = "in", res = 300)
# print(heatmap_1 + heatmap_2 + heatmap_3)
# dev.off()






# #############################

# # Analyse and plot unmodified evenness and absolute meth delta
# Covs_grl_all_df_tally_abs_delta = Variable_distribution(Covs_grl_all_df_common, "Tally", "evenness", "density", "Strand_bias_plots_not_subsample",  c(0.5, 1.0), 1)
# Covs_grl_all_df_tally_evenness = Variable_distribution(Covs_grl_all_df_common, "Tally", "abs_delta_meth_pct", "Abs_delta_meth_pct_density", "Strand_bias_plots_not_subsample", c(0, 100), 1)


# # Analyse and plot group delta's, use first element in list because data is the same
# Covs_grl_all_df_tally_abs_delta_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_delta_evenness", "Group_density", "Strand_bias_plots_not_subsample",  c(-0.5, 0.5), 1)
# Covs_grl_all_df_tally_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_abs_delta_meth_pct", "Group_Abs_delta_meth_pct_density", "Strand_bias_plots_not_subsample", c(-100, 100), 1)

# # Extract WGBS samples (EM-Seq is reference, so all zero's) analyse and plot sample delta's
# Covs_grl_all_df_common = Covs_grl_all_df_common[grep("W$", names(Covs_grl_all_df_common))]
# Covs_grl_all_df_tally_abs_delta_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_delta_evenness", "Sample_density", "Strand_bias_plots_not_subsample", c(-0.5, 0.5), 1)
# Covs_grl_all_df_tally_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_abs_delta_meth_pct", "Sample_Abs_delta_meth_pct_density", "Strand_bias_plots_not_subsample", c(-100, 100), 1)


# #####  Plot graphs for subsampled data #####
# ############################################

# # Create directroy to store plots
# dir.create("Strand_bias_plots_subsample")

# # Filter for common CpGs and calculate deltas between samples and groups
# Covs_grl_all_df_common = Filter_common_CpGs(covs_grl_subsample)
# Covs_grl_all_df_common = sapply(names(Covs_grl_all_df_common), function(x) data.frame(chr = seqnames(Covs_grl_all_df_common[[x]]), pos = start(Covs_grl_all_df_common[[x]]),
                        # N = (Covs_grl_all_df_common[[x]]$meth_cov + Covs_grl_all_df_common[[x]]$unmeth_cov), X = Covs_grl_all_df_common[[x]]$meth_cov, 
                        # evenness = Covs_grl_all_df_common[[x]]$evenness, abs_delta_meth_pct = Covs_grl_all_df_common[[x]]$abs_delta_meth_pct, sample_name = x),
                        # simplify = FALSE, USE.NAMES = TRUE)
            
# Covs_grl_all_df_common = Calculate_delta_meth_pct(Covs_grl_all_df_common, c("WR025V1", "WR025V9", "WR069V1", "WR069V9"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"))

# # Analyse and plot group delta's, use first element in list because data is the same
# Covs_grl_all_df_tally_abs_delta_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_delta_evenness", "Group_density", "Strand_bias_plots_subsample",  c(-0.5, 0.5), 1)
# Covs_grl_all_df_tally_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_abs_delta_meth_pct", "Group_Abs_delta_meth_pct_density", "Strand_bias_plots_subsample", c(-100, 100), 1)

# # Extract WGBS samples (EM-Seq is reference, so all zero's) analyse and plot sample delta's
# Covs_grl_all_df_common = Covs_grl_all_df_common[grep("W$", names(Covs_grl_all_df_common))]
# Covs_grl_all_df_tally_abs_delta_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_delta_evenness", "Sample_density", "Strand_bias_plots_subsample", c(-0.5, 0.5), 1)
# Covs_grl_all_df_tally_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_abs_delta_meth_pct", "Sample_Abs_delta_meth_pct_density", "Strand_bias_plots_subsample", c(-100, 100), 1)

# #save(Covs_grl_all_df_tally_abs_delta_group, Covs_grl_all_df_tally_group, Covs_grl_all_df_tally_abs_delta_sample, Covs_grl_all_df_tally_sample, file ="Strand_bias.RData")
# quit(save = "no")


#############################################################
# # Filter matrix row to remove NA based on cutoff for each group.
# # Cutoff means count of NA > that you want to discard.
#############################################################

Filter_matrix_NA <- function(Matrix_input, group_1, group_2, group_1_cutoff, group_2_cutoff) {

    # Count NA per row for each group
    Row_NA_count_1 = rowSums(is.na((Matrix_input[, group_1])))
    Row_NA_count_2 = rowSums(is.na((Matrix_input[, group_2])))
    
    # Remove rows containing NA's > cutoff
    Matrix_input = Matrix_input[!(Row_NA_count_1 > group_1_cutoff | Row_NA_count_2 > group_2_cutoff), ]
    
    return(Matrix_input)
}


################################################################################
# # Function to to filter granges list object for common CpGs
################################################################################

Filter_common_CpGs <- function(Granges_list_object) {
    Common_intersect_regions = Reduce(intersect, Granges_list_object)
    Subset_covs_grl_common = endoapply(Granges_list_object, subsetByOverlaps, Common_intersect_regions)
    return(Subset_covs_grl_common)
}