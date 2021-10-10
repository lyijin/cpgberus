#####  Script information  #####
################################

# Explore evenness and delta meth between watson and crick.
# See if there are any differences between EM-Seq to WGBS.
# Coverage evenness and absolute delta are based on Yi Jin's python script.
# These values may have to be recalculated for subsampled data.

#####  Dependencies  #####
##########################

library(GenomicRanges)
# library(ggplot2)
# library(scales)
library(reshape2)
# library(ComplexHeatmap)
# library(matrixStats)
# library(circlize)
library(factoextra)
# library(GGally)

#####  Load data  #####
#######################

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")
load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_common.RData")
#load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Data_coverage_filtered.RData")
#load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Strand_bias.RData")


#####  Functions  #####
#######################

###########################################################
# # Funtion to plot either density or tally distribution
# # and return plotting values. Export as png because tally 
# # contains a lot of datapoints making a large PDF.
###########################################################

Variable_distribution <- function(List_of_dataframes, density_or_tally, variable_column_name, plot_name, folder_name, x_axis_cutoff, density_adjust) {

    n = names(List_of_dataframes)
    
    if (density_or_tally == "Density") {
        # Calculate density, then dataframe of x and y density coord + name
        print("Calculating densities")
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name], adjust = density_adjust))
        temp_df = lapply(n, function(n) data.frame(Sample_name = n, x_values = temp_df[[n]][["x"]], y_values = temp_df[[n]][["y"]]))
    } else if (density_or_tally == "Tally") {
        # Calculate tally density, convert factor to numeric, then dataframe of x and y density coord + name
        print("Calculating tallys")
        temp_df = lapply(List_of_dataframes, function(x) data.frame(table(x[, variable_column_name])))
        temp_df = lapply(temp_df, function(x) data.frame(Var1 = as.numeric(levels(x$Var1))[x$Var1], Freq = x$Freq))
        temp_df = lapply(n, function(n) data.frame(Sample_name = n, x_values = temp_df[[n]][["Var1"]], y_values = temp_df[[n]][["Freq"]]))
    }
    
    # Collapse list of dataframes
    temp_df = do.call(rbind, temp_df)
    
    # Add sample name without suffix and sequencing type column
    temp_df$Sample_name_2 = substr(temp_df$Sample_name, 1, nchar(temp_df$Sample_name) - 1)

    temp_df$Seq_type = temp_df$Sample_name
    temp_df$Seq_type[grep("E$", temp_df$Sample_name)] = "EM-Seq"
    temp_df$Seq_type[grep("W$", temp_df$Sample_name)] = "WGBS"

    # Plot
    png(paste0(folder_name, "/01_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0(folder_name, "/02_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0(folder_name, "/03_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0(folder_name, "/04_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + facet_wrap( ~ Sample_name_2, ncol = 2) + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0(folder_name, "/05_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=log2(y_values), color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + facet_wrap( ~ Sample_name_2, ncol = 2) + xlab(variable_column_name) + ylab(paste0("log2 ", density_or_tally)))
    dev.off()
    
    return(temp_df)
}

################################################################################
# # Function to to filter granges list object for common CpGs
################################################################################

Filter_common_CpGs <- function(Granges_list_object) {
    Common_intersect_regions = Reduce(intersect, Granges_list_object)
    Subset_covs_grl_common = endoapply(Granges_list_object, subsetByOverlaps, Common_intersect_regions)
    return(Subset_covs_grl_common)
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


#############################################################
# # Function to convert lists of dataframes to matrix
#############################################################

List_df_to_matrix <- function(List_of_dataframes_common, variable_column_name) {

    samp_id_list = names(List_of_dataframes_common)
    
    # QC to check if "seqnames" and "start" columns are identical
    for (samp_id in samp_id_list) {
        if (identical(List_of_dataframes_common[[samp_id_list[1]]]$seqnames, List_of_dataframes_common[[samp_id]]$seqnames) == FALSE) {
            stop(paste0("Chromosomes not identical between samples: ", samp_id_list[1], " ", samp_id))
        } else if (identical(List_of_dataframes_common[[samp_id_list[1]]]$start, List_of_dataframes_common[[samp_id]]$start) == FALSE) {
            stop(paste0("CpG position not identical between samples: ", samp_id_list[1], " ", samp_id))
        }
    }
    
    # Create matrix
    Final_matrix = matrix(nrow = nrow(List_of_dataframes_common[[1]]), ncol = length(samp_id_list))
    colnames(Final_matrix) = samp_id_list
    rownames(Final_matrix) = paste0(List_of_dataframes_common[[1]]$seqnames, "_", List_of_dataframes_common[[1]]$start)
    
    for (samp_id in samp_id_list) {
        Final_matrix[, samp_id] = List_of_dataframes_common[[samp_id]][, variable_column_name]
    }
    
    return(Final_matrix)
}

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

#############################################################
# # Plot PCA and scatter
#############################################################

Plot_PCA_cor <- function(List_of_dataframes, columns_of_interest, group_1, group_2, file_name) {

    # Create matrix list
    Matrix_list = vector(mode = "list", length = length(columns_of_interest))
    names(Matrix_list) = columns_of_interest

    for (col_var in columns_of_interest) {
        temp_matrix = List_df_to_matrix(List_of_dataframes, col_var)
        colnames(temp_matrix) = paste0(col_var, "_", colnames(temp_matrix))
        Matrix_list[[col_var]] = temp_matrix
    }

    Matrix_list$Combined = do.call("cbind", Matrix_list)
    
    # Calculate PCA for each matrix
    List_names = names(Matrix_list)
    PCA_results = vector(mode = "list", length = length(List_names))
    names(PCA_results) = List_names
    
    for (name_var in List_names) {
        PCA_results[[name_var]] = prcomp(Matrix_list[[name_var]], center = TRUE, scale = TRUE)
    }
    
    # Calculate mean dataframe for each name_var
    Dataframe_for_scatterplot = data.frame(mu1 = List_of_dataframes[[1]]$mu1, mu2 = List_of_dataframes[[1]]$mu2)
    
    for (name_var in List_names) {
        group_1_index = grep(paste0(group_1, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        group_2_index = grep(paste0(group_2, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        Dataframe_for_scatterplot[paste0(name_var, "_group_1")] <- rowMeans(Matrix_list[[name_var]][, group_1_index])
        Dataframe_for_scatterplot[paste0(name_var, "_group_2")] <- rowMeans(Matrix_list[[name_var]][, group_2_index])
    }
    
    # Plot Scree, PCA and scatter
    pdf(paste0(file_name, ".pdf"), width = 8.3, height = 8.3)
    
    for (name_var in List_names) {
        plot(fviz_eig(PCA_results[[name_var]]) + labs(title = name_var))
        
        plot(fviz_pca_biplot(PCA_results[[name_var]], axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) + labs(title = name_var))
        
        if (name_var == "Combined") {
            #pass
        } else if (name_var == "N") {
            name_var_column_names = grep(paste0("^", name_var), colnames(Dataframe_for_scatterplot), value = TRUE)
            temp_df = melt(Dataframe_for_scatterplot, id.vars = c("mu1", "mu2"), measure.vars = name_var_column_names)
            
            p = ggplot(temp_df, aes(x = mu1, y = mu2, color = log2(value)))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
              xlab("EM-Seq") + ylab("WGBS") + ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
        } else {
            name_var_column_names = grep(paste0("^", name_var), colnames(Dataframe_for_scatterplot), value = TRUE)
            temp_df = melt(Dataframe_for_scatterplot, id.vars = c("mu1", "mu2"), measure.vars = name_var_column_names)
            
            p = ggplot(temp_df, aes(x = mu1, y = mu2, color = value))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
              xlab("EM-Seq") + ylab("WGBS") + ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
        }
    }
    
    dev.off()
}


#####  Stand bias analysis for common and significant CpGs  #####
#####  for unsmoothed data.                                 #####
#################################################################

# Extract significant CpGs from DMLtest
Significant_regions = do.call("rbind", Covs_grl_common_bsseq_list_stats_no_smooth)
Significant_regions = Significant_regions[(Significant_regions$fdr <= 0.05), ]
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_not_smoothed_common = lapply(covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_not_smoothed_common = lapply(covs_grl_sig_not_smoothed_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})
                                            
covs_grl_sig_not_smoothed_common = sapply(names(covs_grl_sig_not_smoothed_common), function(x) data.frame(covs_grl_sig_not_smoothed_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_not_smoothed_common = sapply(names(covs_grl_sig_not_smoothed_common), function(x) cbind(covs_grl_sig_not_smoothed_common[[x]], data.frame(N = covs_grl_sig_not_smoothed_common[[x]]$meth_cov + covs_grl_sig_not_smoothed_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
Plot_PCA_cor(covs_grl_sig_not_smoothed_common, c("evenness", "abs_delta_meth_pct", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Not_smoothed_strand_bias")


#####  Stand bias analysis for common and significant CpGs  #####
#####  for smoothed data.                                   #####
#################################################################

# Extract significant CpGs from DMLtest
Significant_regions = do.call("rbind", Covs_grl_common_bsseq_list_stats)
Significant_regions = Significant_regions[(Significant_regions$fdr <= 0.05), ]
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})
                                            
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Smoothed_strand_bias")

quit(save = "no")

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