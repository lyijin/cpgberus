#####  Script information  #####
################################

# Explore evenness and delta meth between watson and crick.
# See if there are any differences between EM-Seq to WGBS.
# Coverage evenness and absolute delta are based on Yi Jin's python script.
# These values may have to be recalculated for subsampled data.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(ggplot2)
library(scales)
library(reshape2)
library(ComplexHeatmap)
# library(matrixStats)
# library(circlize)
library(factoextra)
# library(GGally)
library(BSgenome.Hsapiens.UCSC.hg38)

#####  Load data  #####
#######################

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")
load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_common.RData")
#load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Data_coverage_filtered.RData")
#load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Strand_bias.RData")


#####  Functions  #####
#######################

########################################################
# # Funtion to plot either density or tally distribution
# # and return plotting values.
########################################################

Variable_distribution <- function(List_of_dataframes, density_or_tally, variable_column_name, plot_name, x_axis_cutoff) {

    n = names(List_of_dataframes)
    
    if (density_or_tally == "Density") {
        # Calculate density, then dataframe of x and y density coord + name
        print("Calculating densities")
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name]))
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
    pdf(paste0(plot_name, ".pdf"), width = 11.69, height = 8.3)
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
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

Plot_PCA_cor <- function(List_of_dataframes, columns_of_interest, group_1, group_2, group_1_label, group_2_label, file_name) {

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
    Dataframe_for_scatterplot = data.frame(EM_Seq_beta = List_of_dataframes[[1]]$mu1, WGBS_beta = List_of_dataframes[[1]]$mu2)
    
    for (name_var in List_names) {
        group_1_index = grep(paste0(group_1, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        group_2_index = grep(paste0(group_2, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        Dataframe_for_scatterplot[paste0(name_var, "_", group_1_label)] <- rowMeans(Matrix_list[[name_var]][, group_1_index])
        Dataframe_for_scatterplot[paste0(name_var, "_", group_2_label)] <- rowMeans(Matrix_list[[name_var]][, group_2_index])
    }
    
    rownames(Dataframe_for_scatterplot) = rownames(Matrix_list[[1]])
    
    # Plot Scree, PCA and scatter
    pdf(paste0(file_name, ".pdf"), width = 8.3, height = 8.3)
    
    for (name_var in List_names) {
        plot(fviz_eig(PCA_results[[name_var]]) + labs(title = name_var))
        
        plot(fviz_pca_biplot(PCA_results[[name_var]], axes = c(1, 2), col.var = "contrib", label = "var", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) + labs(title = name_var))
        
        if (name_var == "Combined") {
            #pass
        } else if (name_var == "N") {
            name_var_column_names = grep(paste0("^", name_var), colnames(Dataframe_for_scatterplot), value = TRUE)
            temp_df = melt(Dataframe_for_scatterplot, id.vars = c("EM_Seq_beta", "WGBS_beta"), measure.vars = name_var_column_names)
            
            p = ggplot(temp_df, aes(x = EM_Seq_beta, y = WGBS_beta, color = log2(value)))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
                            ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
            
            temp_df = melt(Dataframe_for_scatterplot, id.vars = name_var_column_names, measure.vars = c("EM_Seq_beta", "WGBS_beta"))
            
            p = ggplot(temp_df, aes(x = get(name_var_column_names[1]), y = get(name_var_column_names[2]), color = value))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
                            xlab(paste0("EM-Seq ", name_var)) + ylab(paste0("WGBS ", name_var)) + ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
        } else {
            name_var_column_names = grep(paste0("^", name_var), colnames(Dataframe_for_scatterplot), value = TRUE)
            temp_df = melt(Dataframe_for_scatterplot, id.vars = c("EM_Seq_beta", "WGBS_beta"), measure.vars = name_var_column_names)
            
            p = ggplot(temp_df, aes(x = EM_Seq_beta, y = WGBS_beta, color = value))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
                            ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
            
            temp_df = melt(Dataframe_for_scatterplot, id.vars = name_var_column_names, measure.vars = c("EM_Seq_beta", "WGBS_beta"))

            p = ggplot(temp_df, aes(x = get(name_var_column_names[1]), y = get(name_var_column_names[2]), color = value))
            Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
                            xlab(paste0("EM-Seq ", name_var)) + ylab(paste0("WGBS ", name_var)) + ggtitle(name_var) + facet_wrap( ~ variable)
            plot(Final_graph)
        }
    }
    
    dev.off()
    
    return(list(Dataframe_for_scatterplot, Matrix_list))
}

#######################################################################
# # Funtion to plot extract motifs and calculcate CG percentage from
# # a dataframe of C positions
#######################################################################

Extract_motif_calculate_GC <- function(Input_dataframe, nucleotides_backward, nucleotides_forward) {

    Motif_grange = GRanges(seqnames = Input_dataframe$seqnames,
        IRanges(start = Input_dataframe$start - nucleotides_backward,
        end = (Input_dataframe$start + nucleotides_forward)))
    
    Final_dataframe_grange = GRanges(seqnames = Input_dataframe$seqnames,
        IRanges(start = Input_dataframe$start, end = (Input_dataframe$start + 1)))

    Final_dataframe_grange$cpg_context_nnncgnnn <- getSeq(BSgenome.Hsapiens.UCSC.hg38, Motif_grange)
    Final_dataframe_grange$GC_percentage <- as.numeric(letterFrequency(Final_dataframe_grange$cpg_context_nnncgnnn, letters = "GC", as.prob = TRUE))
    
    Final_dataframe_grange = data.frame(Final_dataframe_grange)
    rownames(Final_dataframe_grange) = paste0(Final_dataframe_grange$seqnames, "_", Final_dataframe_grange$start)
    
    return(Final_dataframe_grange)
}


#####  Stand bias analysis for common and significant CpGs  #####
#####  for unsmoothed data.                                 #####
#################################################################

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Covs_grl_common_bsseq_list_stats_no_smooth), function(x) Covs_grl_common_bsseq_list_stats_no_smooth[[x]][Covs_grl_common_bsseq_list_stats_no_smooth[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
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
covs_grl_sig_not_smoothed_common = sapply(names(covs_grl_sig_not_smoothed_common), function(x) cbind(covs_grl_sig_not_smoothed_common[[x]], data.frame(Beta = covs_grl_sig_not_smoothed_common[[x]]$meth_cov / covs_grl_sig_not_smoothed_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
Not_smoothed_common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_not_smoothed_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "EM-Seq", "WGBS", "Not_smoothed_strand_bias")

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_not_smoothed_common, "Tally", "N", "Significant_CpGs_common_tally_no_smooth", c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_not_smoothed_common, "Density", "N", "Significant_CpGs_common_density_no_smooth", c(-500, 4000))

# Make a heatmap here, test code.
Beta_matrix = Not_smoothed_common_means_and_matrices[[2]]$Beta
Evenness_matrix = Not_smoothed_common_means_and_matrices[[2]]$evenness
Abs_delta_matrix = Not_smoothed_common_means_and_matrices[[2]]$abs_delta_meth_pct
Coverage_matrix = Not_smoothed_common_means_and_matrices[[2]]$N
Coverage_matrix_delta = as.matrix(data.frame(WR025V1 = Coverage_matrix[,"N_WR025V1E"] - Coverage_matrix[,"N_WR025V1W"], WR025V9 = Coverage_matrix[,"N_WR025V9E"] - Coverage_matrix[,"N_WR025V9W"],
                                    WR069V1= Coverage_matrix[,"N_WR069V1E"] - Coverage_matrix[,"N_WR069V1W"], WR069V9 = Coverage_matrix[,"N_WR069V9E"] - Coverage_matrix[,"N_WR069V9W"]))
GC_matrix = Extract_motif_calculate_GC(covs_grl_sig_not_smoothed_common[[1]], nucleotides_backward = 3, nucleotides_forward = 4)
GC_matrix = as.matrix(GC_matrix["GC_percentage"])

colnames(Beta_matrix) = gsub("^.*_", "", colnames(Beta_matrix))
colnames(Evenness_matrix) = gsub("^.*_", "", colnames(Evenness_matrix))
colnames(Abs_delta_matrix) = gsub("^.*_", "", colnames(Abs_delta_matrix))
    
heatmap_1 = Heatmap(Beta_matrix, name = "Beta", column_title = "Beta", row_title = "CpG position", show_row_names = FALSE, show_row_dend = TRUE)
heatmap_2 = Heatmap(Evenness_matrix, name = "evenness", column_title = "Evenness", row_title = "CpG position", show_row_names = FALSE)
heatmap_3 = Heatmap(Abs_delta_matrix, name = "abs_delta_meth_pct", column_title = "Absolute delta meth %", row_title = "CpG position", show_row_names = FALSE)
heatmap_4 = Heatmap(Coverage_matrix_delta, name = "Coverage", column_title = "Coverage", row_title = "CpG position", show_row_names = FALSE)
heatmap_5 = Heatmap(GC_matrix, name = "GC_percentage", column_title = "GC_percentage", row_title = "CpG position", show_row_names = FALSE)

png("Test.png", width = 11.69, height = 8.3, units = "in", res = 300)
print(heatmap_1 + heatmap_2 + heatmap_3 + heatmap_4 + heatmap_5)
dev.off()

#####  Stand bias analysis for common and significant CpGs  #####
#####  for smoothed data.                                   #####
#################################################################

# # Extract significant CpGs from DMLtest
# Significant_regions = sapply(names(Covs_grl_common_bsseq_list_stats), function(x) Covs_grl_common_bsseq_list_stats[[x]][Covs_grl_common_bsseq_list_stats[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
# Significant_regions = do.call("rbind", Significant_regions)
# Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                # mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# # Extract significant regions from unfiltered data and append mu1 and mu2 metadata
# covs_grl_sig_common = lapply(covs_grl, mergeByOverlaps, Significant_regions)
# covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            # temp = x[, 1]
                                            # values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            # return(temp)})
                                            
# covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
# covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
# covs_grl_sig_common = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "EM-Seq", "WGBS", "Smoothed_strand_bias")

quit(save = "no")