#####  Script information  #####
################################

# Explore evenness and delta meth between watson and crick.
# See if there are any differences between EM-Seq to WGBS.
# Coverage evenness and absolute delta are based on Yi Jin's python script.
# These values may have to be recalculated for subsampled data, when subsampling
# the counts.
#
# Also, the mean beta calculated by the DSS package (mu1 and mu2) differs
# from manual mean beta calculations. Reading the raw code, there seems to be a small 
# constant added from documentation: "adding a small constant could bring trouble when there's no coverage!!!".
# This is why the graphs are different between manual and DSS calculated scatterplots for beta.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(ggplot2)
library(scales)
library(reshape2)
library(ComplexHeatmap)
library(factoextra)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gridExtra)
library(circlize)

#####  Load data  #####
#######################

path_to_cpgerus = "/scratch/user/uqdguanz/Projects/Meth/cpgberus"
source(paste0(path_to_cpgerus, "/05_CpG_sequence_context/07_Motif_functions.R"))

# Load data just before analysis below and remove objects after to save on memory

# load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Not_rarefied_grch38p13_combined_covs_grl.RData"))
# load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Rarefied_grch38p13_combined_covs_grl.RData"))
# load(file.path(path_to_cpgerus, "04_parse_bismark_covs/grch38p13_combined_covs_grl.RData"))

# load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Not_rarefied_CpG_stats_common.RData"))
# load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Rarefied_CpG_stats_common.RData"))
# load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Original_CpG_stats_common.RData"))
# load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Original_CpG_stats_existing.RData"))

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
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name], na.rm = TRUE))
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
# # Plot PCA and scatter
#
# Rows with NA will be removed for the PCA analysis.
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
        Matrix_list_no_NA = Matrix_list[[name_var]]
        Matrix_list_no_NA = Matrix_list_no_NA[complete.cases(Matrix_list_no_NA),]
        PCA_results[[name_var]] = prcomp(Matrix_list_no_NA, center = TRUE, scale = TRUE)
    }
    
    # Calculate mean dataframe for each name_var
    Dataframe_for_scatterplot = data.frame(EM_Seq_beta = List_of_dataframes[[1]]$mu1, WGBS_beta = List_of_dataframes[[1]]$mu2)
    
    for (name_var in List_names) {
        group_1_index = grep(paste0(group_1, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        group_2_index = grep(paste0(group_2, "$", collapse = "|"), colnames(Matrix_list[[name_var]]))
        Dataframe_for_scatterplot[paste0(name_var, "_", group_1_label)] <- rowMeans(Matrix_list[[name_var]][, group_1_index], na.rm = TRUE)
        Dataframe_for_scatterplot[paste0(name_var, "_", group_2_label)] <- rowMeans(Matrix_list[[name_var]][, group_2_index], na.rm = TRUE)
    }
    
    rownames(Dataframe_for_scatterplot) = rownames(Matrix_list[[1]])
    
    # Plot Scree, PCA and scatter
    pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    
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
#
# Needs seqnames and start columns
#######################################################################

Extract_motif_calculate_GC <- function(Input_dataframe, nucleotides_backward, nucleotides_forward) {

    Motif_grange = GRanges(seqnames = Input_dataframe$seqnames,
        IRanges(start = Input_dataframe$start - nucleotides_backward,
        end = (Input_dataframe$start + nucleotides_forward)))
    
    Final_dataframe_grange = GRanges(seqnames = Input_dataframe$seqnames,
        IRanges(start = Input_dataframe$start, end = (Input_dataframe$start + 1)))

    Final_dataframe_grange$cpg_context_nnncgnnn <- getSeq(BSgenome.Hsapiens.UCSC.hg38, Motif_grange)
    Final_dataframe_grange$Motif_GC_percentage <- as.numeric(letterFrequency(Final_dataframe_grange$cpg_context_nnncgnnn, letters = "GC", as.prob = TRUE))
    
    Final_dataframe_grange = data.frame(Final_dataframe_grange)
    rownames(Final_dataframe_grange) = paste0(Final_dataframe_grange$seqnames, "_", Final_dataframe_grange$start)
    
    return(Final_dataframe_grange)
}

#######################################################################
# # Complex heatmap wrapper.
#
# For column bar annotations, is based on grouping WGBS and EM-seq.
# Columns orders should follow Annotation_groups object below. 
# List_of_matrices should be the output from List_df_to_matrix function.
#
# Todo: A lot is hardcoded, probably better do some kind of loop
# through the List_of_matrices.
#######################################################################

Plot_complex_heatmap <- function(List_of_matrices, EM_Seq_names, WGBS_Seq_names, File_name) {

    Beta_matrix = List_of_matrices$Beta
    Evenness_matrix = List_of_matrices$evenness
    Abs_delta_matrix = List_of_matrices$abs_delta_meth_pct
    Coverage_matrix = List_of_matrices$N
    # Coverage_matrix_delta = as.matrix(data.frame(WR025V1 = Coverage_matrix[,"N_WR025V1E"] - Coverage_matrix[,"N_WR025V1W"], WR025V9 = Coverage_matrix[,"N_WR025V9E"] - Coverage_matrix[,"N_WR025V9W"],
                                        # WR069V1= Coverage_matrix[,"N_WR069V1E"] - Coverage_matrix[,"N_WR069V1W"], WR069V9 = Coverage_matrix[,"N_WR069V9E"] - Coverage_matrix[,"N_WR069V9W"]))
    
    seqnames_var = sapply(strsplit(rownames(List_of_matrices$Beta), "_"), "[", 1)
    start_var = as.numeric(sapply(strsplit(rownames(List_of_matrices$Beta), "_"), "[", -1))
    Coverage_dataframe = data.frame(seqnames = seqnames_var, start = start_var)
    
    GC_matrix = Extract_motif_calculate_GC(Coverage_dataframe, nucleotides_backward = 3, nucleotides_forward = 4)
    GC_matrix = as.matrix(GC_matrix["Motif_GC_percentage"])
	colnames(GC_matrix) = "Motif GC %"

    colnames(Beta_matrix) = gsub("^.*_", "", colnames(Beta_matrix))
    colnames(Evenness_matrix) = gsub("^.*_", "", colnames(Evenness_matrix))
    colnames(Abs_delta_matrix) = gsub("^.*_", "", colnames(Abs_delta_matrix))
    colnames(Coverage_matrix) = gsub("^.*_", "", colnames(Coverage_matrix))
    
	Annotation_groups = c("EM-seq", "WGBS", "EM-seq", "WGBS", "EM-seq", "WGBS", "EM-seq", "WGBS")
    ha = HeatmapAnnotation(Seq_method = Annotation_groups, col = list(Seq_method = c("EM-seq" = "#1b9e77", "WGBS" = "#7570b3")), show_annotation_name = FALSE, annotation_legend_param = list(title = "Library type"))
	
	# Calculate dendrogram. Fill Beta matrix NA with row (CpG) mean for clustering.
	Beta_matrix_mean_NA_EM_seq = t(apply(Beta_matrix[ , EM_Seq_names], 1, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))))
	Beta_matrix_mean_NA_WGBS_seq = t(apply(Beta_matrix[ , WGBS_Seq_names], 1, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))))
	Beta_matrix_mean_NA = cbind(Beta_matrix_mean_NA_EM_seq, Beta_matrix_mean_NA_WGBS_seq)
	
	Beta_matrix_mean_NA = Beta_matrix_mean_NA[ , colnames(Beta_matrix)]
	Clustering_NA_row = as.dendrogram(hclust(dist(Beta_matrix_mean_NA)))
	
    heatmap_1 = Heatmap(Beta_matrix, col = colorRamp2(seq(min(Beta_matrix, na.rm = TRUE), max(Beta_matrix, na.rm = TRUE), length = 3), c("#4575b4", "#fee090", "#d73027")), na_col = "white",
						name = "Beta", column_title = "Beta", row_title = "CpG dinucleotide", show_row_names = FALSE, cluster_rows = Clustering_NA_row, top_annotation = ha, row_split = 2, border = TRUE)
    heatmap_2 = Heatmap(Evenness_matrix, col = colorRamp2(seq(min(Evenness_matrix, na.rm = TRUE), max(Evenness_matrix, na.rm = TRUE), length = 3), c("#4575b4", "#fee090", "#d73027")), heatmap_legend_param = list(at = seq(min(Evenness_matrix, na.rm = TRUE), max(Evenness_matrix, na.rm = TRUE), length = 5)), na_col = "white",
						name = "Evenness", column_title = "Evenness", row_title = "CpG dinucleotide", show_row_names = FALSE, top_annotation = ha, border = TRUE)
    heatmap_3 = Heatmap(Abs_delta_matrix, col = colorRamp2(seq(min(Abs_delta_matrix, na.rm = TRUE), max(Abs_delta_matrix, na.rm = TRUE), length = 3), c("#4575b4", "#fee090", "#d73027")), na_col = "white",
						name = "Absolute delta meth %", column_title = "Absolute delta meth %", row_title = "CpG dinucleotide", show_row_names = FALSE, top_annotation = ha, border = TRUE)
    heatmap_4 = Heatmap(log2(Coverage_matrix), col = colorRamp2(seq(min(log2(Coverage_matrix), na.rm = TRUE), max(5, na.rm = TRUE), length = 3), c("#4575b4", "#fee090", "#d73027")), heatmap_legend_param = list(at = round(seq(min(log2(Coverage_matrix), na.rm = TRUE), max(log2(Coverage_matrix), na.rm = TRUE), length = 6), 1)), na_col = "white",
						name = "Log2 Coverage", column_title = "Log2 Coverage", row_title = "CpG dinucleotide", show_row_names = FALSE, top_annotation = ha, border = TRUE)
    heatmap_5 = Heatmap(GC_matrix, col = colorRamp2(seq(min(GC_matrix, na.rm = TRUE), max(GC_matrix, na.rm = TRUE), length = 3), c("#4575b4", "#fee090", "#d73027")), heatmap_legend_param = list(at = seq(min(GC_matrix, na.rm = TRUE), max(GC_matrix, na.rm = TRUE), length = 7)), na_col = "white",
						name = "Motif GC %", row_title = "CpG dinucleotide", show_row_names = FALSE, border = TRUE)

    png(paste0(File_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(heatmap_1 + heatmap_2 + heatmap_3 + heatmap_4 + heatmap_5)
    dev.off()
}

################################################################################
# # Function to merge CpGs which exist between all samples, fill empty with 0.
# 
# Needs chr, pos, N and X. N is total coverage, X is meth.
################################################################################

Merge_CpGs <- function(Dataframe_list_object) {

    # Find existing chromosomes
    number_of_chromosomes = list()
    
    for (data_frame in Dataframe_list_object) {
        temp_list = as.character(unique(data_frame$chr))
        number_of_chromosomes = unlist(unique(c(number_of_chromosomes, temp_list)))
    }
    
    # Find existing CpGs in all samples per chromosome
    existing_CpGs = data.frame(chr = factor(), pos = integer())
    
    for (chr_num in number_of_chromosomes) {
        temp_pos = list()
        
        for (data_frame in Dataframe_list_object) {
            chr_num_pos = data_frame[(data_frame$chr == chr_num), "pos"]
            temp_pos = unique(unlist(c(temp_pos, chr_num_pos)))
        }
        
        existing_CpGs = rbind(existing_CpGs, data.frame(chr = factor(chr_num), pos = temp_pos))
        print(paste0(chr_num, " find existing CpGs is done."))
    }
    
    # Compare sample CpGs to existing CpGs and fill differences with 0    
    final_dataframe_list = list()
    
    for (dataframe_index in seq(length(Dataframe_list_object))) {
        data_frame = Dataframe_list_object[[dataframe_index]]
        data_frame_name = names(Dataframe_list_object[dataframe_index])
        
        temp_df = data.frame(matrix(ncol = ncol(data_frame), nrow = 0))
        colnames(temp_df) = colnames(data_frame)

        for (chr_num in number_of_chromosomes) {
            chr_num_pos = data_frame[(data_frame$chr == chr_num), "pos"]
            chr_num_pos_ref = existing_CpGs[(existing_CpGs$chr == chr_num), "pos"]
            chr_num_pos_missing = chr_num_pos_ref[!(chr_num_pos_ref %in% chr_num_pos)]
            
            if (length(chr_num_pos_missing) != 0) {
                temp_df_2 = data.frame(matrix(ncol = ncol(data_frame), nrow = 0))
                colnames(temp_df_2) = colnames(data_frame)
                
                temp_df_2[1:length(chr_num_pos_missing), "chr"] = chr_num
                temp_df_2[1:length(chr_num_pos_missing), "pos"] = chr_num_pos_missing
                temp_df_2[1:length(chr_num_pos_missing), "N"] = 0
                temp_df_2[1:length(chr_num_pos_missing), "X"] = 0
                temp_df = rbind(temp_df, temp_df_2)
            }
        }
        
        temp_df = rbind(data_frame, temp_df)
        temp_df = temp_df[order(temp_df$chr, temp_df$pos), ]
        rownames(temp_df) <- seq(1:nrow(temp_df))
        temp_df$chr = factor(temp_df$chr)
        final_dataframe_list[[data_frame_name]] <- temp_df
        print(paste0(data_frame_name, " data is finished merging: ", dataframe_index, " / ", length(Dataframe_list_object)))
    }
    return(final_dataframe_list)
}

#############################################################
### Function to do split violins
### Code from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
#############################################################

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#############################################################
# # Plot correlation of coverage, motif GC percentage and 
# # absolute delta meth %
#############################################################

Plot_correlation <- function(Input_matrix, EM_Seq_names, WGBS_Seq_names, File_name) {
	
	seqnames_var = sapply(strsplit(rownames(Input_matrix$Beta), "_"), "[", 1)
	start_var = as.numeric(sapply(strsplit(rownames(Input_matrix$Beta), "_"), "[", -1))
	Coverage_dataframe = data.frame(seqnames = seqnames_var, start = start_var)
		
	GC_matrix = Extract_motif_calculate_GC(Coverage_dataframe, nucleotides_backward = 3, nucleotides_forward = 4)
	GC_matrix = as.matrix(GC_matrix["Motif_GC_percentage"])
							
	Data_for_correlation = data.frame(cbind(Input_matrix$Beta, Input_matrix$abs_delta_meth_pct, Input_matrix$N, GC_matrix))
	Data_for_correlation$Meth_pos = rownames(Data_for_correlation)

	Data_for_correlation = melt(Data_for_correlation, id.vars=colnames(Data_for_correlation)[grep("Motif_GC_percentage|Meth_pos", colnames(Data_for_correlation))], variable.name="Condition", value.name="Condition_value")

	Abs_data = Data_for_correlation[grep("^abs_delta_meth_pct_", Data_for_correlation$Condition), ]
	Beta_data = Data_for_correlation[grep("^Beta_", Data_for_correlation$Condition), ]
	Coverage_data = Data_for_correlation[grep("^N_", Data_for_correlation$Condition), ]

	colnames(Abs_data) = gsub("^Condition_value$", "Abs_delta_meth_pct", colnames(Abs_data))
	colnames(Beta_data) = gsub("^Condition_value$", "Beta", colnames(Beta_data))
	colnames(Coverage_data) = gsub("^Condition_value$", "N", colnames(Coverage_data))

	Abs_data$Condition = gsub("^abs_delta_meth_pct_", "", Abs_data$Condition)
	Beta_data$Condition = gsub("^Beta_", "", Beta_data$Condition)
	Coverage_data$Condition =  gsub("^N_", "", Coverage_data$Condition)

	Data_for_correlation = merge(Abs_data, Beta_data, by = c("Motif_GC_percentage", "Meth_pos", "Condition"))
	Data_for_correlation = merge(Data_for_correlation, Coverage_data, by = c("Motif_GC_percentage", "Meth_pos", "Condition"))
	Data_for_correlation$Library_type = "Unknown"
	Data_for_correlation[grep(paste0("^", EM_Seq_names, "$", collapse = "|"), Data_for_correlation$Condition), "Library_type"] = "EM-seq"
	Data_for_correlation[grep(paste0("^", WGBS_Seq_names, "$", collapse = "|"), Data_for_correlation$Condition), "Library_type"] = "WGBS"

	Data_for_correlation$Motif_GC_percentage = as.factor(Data_for_correlation$Motif_GC_percentage)

	stat_box_data <- function(y, upper_limit = max(log2(Data_for_correlation$N), na.rm = TRUE) * 1.15) {
	  return( 
		data.frame(
		  y = 0.95 * upper_limit,
		  label = paste('count =', length(y), '\n',
						'median =', round(median(y, na.rm = TRUE), 1), '\n')
		)
	  )
	}

	png(paste0(File_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
	print(ggplot(Data_for_correlation, aes(Motif_GC_percentage, log2(N))) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = Beta, size = Abs_delta_meth_pct), width = 0.2) +
		 stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=2.5) + labs(x = "Motif GC %", y = "Log2 (Coverage)") + scale_size_continuous(range = c(1, 4)) + theme_bw() + facet_wrap( ~ Library_type))
	dev.off()
	
	p1 = ggplot(Data_for_correlation, aes(Motif_GC_percentage, log2(N))) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05) +
			stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=2.5) + labs(y = "Log2 (Coverage)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_wrap( ~ Library_type)
		 
	p2 = ggplot(Data_for_correlation, aes(Motif_GC_percentage, Beta)) + geom_violin() + geom_jitter(width = 0.05) + labs(x = "Motif GC %", y = "Beta") + theme_bw() + facet_wrap( ~ Library_type)

	p3 = arrangeGrob(p1, p2, nrow = 2)
	
	ggsave(file=paste0(File_name, "_2.png"), p3, width = 11.69, height = 8.3, units = "in", dpi = 300)
	
	p1 = ggplot(Data_for_correlation, aes(Motif_GC_percentage, log2(N), color = Library_type)) + geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) + geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
			labs(y = "Log2 (Coverage)") + theme_minimal(13) + scale_color_manual(values=c("#1b9e77", "#7570b3")) + labs(color = "Library type") + theme(axis.title.x=element_blank())
	
	p2 = ggplot(Data_for_correlation, aes(Motif_GC_percentage, Beta, fill = Library_type)) + geom_split_violin() + labs(x = "Motif GC %", y = "Beta", fill = "Library type") + theme_minimal(13) + scale_fill_manual(values=c("#1b9e77", "#7570b3"))
	
	p3 = arrangeGrob(p1, p2, nrow = 2)
	
	ggsave(file=paste0(File_name, "_3.png"), p3, width = 11.69, height = 6, units = "in", dpi = 300)
	
}

#####  Stand bias analysis for *all existing* and significant CpGs  #####
#####  between samples for *unsmoothed data*. This differs from     #####
#####  common CpG because the CpG only has to exist in one sample   #####
#####  to be included in stats, while common CpGs has to exist in   #####
#####  all samples to be included in stats. Rarefied data is used.  #####
#########################################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/04_outputs/Rarefied_existing_CpGs")
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Rarefied_CpG_stats_existing.RData"))

dir.create(full_path, recursive = TRUE)

Rarefied_covs_grl = lapply(Rarefied_covs_grl, function(x) cbind(data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov), data.frame(mcols(x))))
Rarefied_covs_grl = Merge_CpGs(Rarefied_covs_grl)
Rarefied_covs_grl = lapply(Rarefied_covs_grl, function(x) { temp = GRanges(seqnames = x$chr, ranges=IRanges(start = x$pos, end = (x$pos + 1)))
                                                values(temp) = x[!(colnames(x) %in% c("chr", "pos", "N", "X"))]
                                                return(temp)})                                            
                                                
Rarefied_covs_grl = do.call("GRangesList", Rarefied_covs_grl)

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Rarefied_existing_stats_no_smooth), function(x) Rarefied_existing_stats_no_smooth[[x]][Rarefied_existing_stats_no_smooth[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(Rarefied_covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})
         
# Plot correlations         
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(Beta = covs_grl_sig_common[[x]]$meth_cov / covs_grl_sig_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), "EM-Seq", "WGBS", file.path(full_path, "Unsmoothed_significant_CpGs_PCA"))

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_common, "Tally", "N", file.path(full_path, "Unsmoothed_significant_CpGs_tally"), c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_common, "Density", "N", file.path(full_path, "Unsmoothed_significant_CpGs_density"), c(-500, 4000))

# Plot heatmap. For row clustering, missing values use the group CpG average.
common_means_and_matrices_subset = common_means_and_matrices[[2]][c("evenness", "abs_delta_meth_pct", "Beta", "N")]
Plot_complex_heatmap(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), WGBS_Seq_names = c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), file.path(full_path, "Unsmoothed_significant_CpGs_heatmap"))

# Analyse motifs of heatmap.
seqnames_var_2 = sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", 1)
start_var_2 = as.numeric(sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", -1))
Coverage_dataframe_2 = data.frame(seqnames = seqnames_var_2, start = start_var_2)

GC_matrix_2 = Extract_motif_calculate_GC(Coverage_dataframe_2, nucleotides_backward = 3, nucleotides_forward = 4)

# Convert to GC% bin list of motifs
GC_matrix_list = list()
for (bin_var in sort(unique(GC_matrix_2$Motif_GC_percentage))) {
	
	GC_matrix_subset = GC_matrix_2[GC_matrix_2$Motif_GC_percentage == bin_var, ]
	temp_df = list(GC_matrix_subset)
	names(temp_df) = bin_var
	GC_matrix_list = append(GC_matrix_list, temp_df)
}

# Plot motif
counter <<- 0
Motif_plot_data = lapply(GC_matrix_list, Motif_frequency_table, length(GC_matrix_list))
Motif_plot(Motif_plot_data, file.path(full_path, "Motif_plots"), names(Motif_plot_data), TRUE)

# Plot correlations
Plot_correlation(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), WGBS_Seq_names = c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), file.path(full_path, "Correlation_plot"))

# Remove objects to clear memory
rm(Rarefied_covs_grl, Significant_regions, covs_grl_sig_common, common_means_and_matrices, temp_df, common_means_and_matrices_subset, Rarefied_existing_bsseq, Rarefied_existing_stats_smooth, Rarefied_existing_stats_no_smooth)
gc()


#####  Stand bias analysis for *all existing* and significant CpGs      #####
#####  between samples for *unsmoothed data*. This differs from         #####
#####  common CpG because the CpG only has to exist in one sample       #####
#####  to be included in stats, while common CpGs has to exist in       #####
#####  all samples to be included in stats. Not rarefied data is used.  #####
#############################################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/04_outputs/Not_rarefied_existing_CpGs")
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Not_rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Not_rarefied_CpG_stats_existing.RData"))

dir.create(full_path, recursive = TRUE)

# Create covs_grl (a Granges list) which contains raw data for all existing CpGs which exist for each sample.
Not_rarefied_covs_grl = lapply(Not_rarefied_covs_grl, function(x) cbind(data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov), data.frame(mcols(x))))
Not_rarefied_covs_grl = Merge_CpGs(Not_rarefied_covs_grl)
Not_rarefied_covs_grl = lapply(Not_rarefied_covs_grl, function(x) { temp = GRanges(seqnames = x$chr, ranges=IRanges(start = x$pos, end = (x$pos + 1)))
                                                values(temp) = x[!(colnames(x) %in% c("chr", "pos", "N", "X"))]
                                                return(temp)})                                            
                                                
Not_rarefied_covs_grl = do.call("GRangesList", Not_rarefied_covs_grl)

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Not_rarefied_existing_stats_no_smooth), function(x) Not_rarefied_existing_stats_no_smooth[[x]][Not_rarefied_existing_stats_no_smooth[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(Not_rarefied_covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})
         
# Plot correlations         
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(Beta = covs_grl_sig_common[[x]]$meth_cov / covs_grl_sig_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "EM-Seq", "WGBS", file.path(full_path, "Unsmoothed_significant_CpGs_PCA"))

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_common, "Tally", "N", file.path(full_path, "Unsmoothed_significant_CpGs_tally"), c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_common, "Density", "N", file.path(full_path, "Unsmoothed_significant_CpGs_density"), c(-500, 4000))

# Plot heatmap. For row clustering, missing values use the group CpG average.
common_means_and_matrices_subset = common_means_and_matrices[[2]][c("evenness", "abs_delta_meth_pct", "Beta", "N")]
Plot_complex_heatmap(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), WGBS_Seq_names = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), file.path(full_path, "Unsmoothed_significant_CpGs_heatmap"))

# Analyse motifs of heatmap.
seqnames_var_2 = sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", 1)
start_var_2 = as.numeric(sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", -1))
Coverage_dataframe_2 = data.frame(seqnames = seqnames_var_2, start = start_var_2)

GC_matrix_2 = Extract_motif_calculate_GC(Coverage_dataframe_2, nucleotides_backward = 3, nucleotides_forward = 4)

# Convert to GC% bin list of motifs
GC_matrix_list = list()
for (bin_var in sort(unique(GC_matrix_2$Motif_GC_percentage))) {
	
	GC_matrix_subset = GC_matrix_2[GC_matrix_2$Motif_GC_percentage == bin_var, ]
	temp_df = list(GC_matrix_subset)
	names(temp_df) = bin_var
	GC_matrix_list = append(GC_matrix_list, temp_df)
}

# Plot motif
counter <<- 0
Motif_plot_data = lapply(GC_matrix_list, Motif_frequency_table, length(GC_matrix_list))
Motif_plot(Motif_plot_data, file.path(full_path, "Motif_plots"), names(Motif_plot_data), TRUE)

# Plot correlations
Plot_correlation(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), WGBS_Seq_names = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), file.path(full_path, "Correlation_plot"))

# Remove temporary covs_grl after all existing analysis.
rm(Not_rarefied_covs_grl, Significant_regions, covs_grl_sig_common, common_means_and_matrices, temp_df, common_means_and_matrices_subset, Not_rarefied_existing_bsseq, Not_rarefied_existing_stats_smooth, Not_rarefied_existing_stats_no_smooth)
gc()


#####  Stand bias analysis for *common* and significant CpGs  #####
#####  for *unsmoothed data*. Rarefied data is used.          #####
###################################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/04_outputs/Rarefied_common_CpGs")
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Rarefied_CpG_stats_common.RData"))

dir.create(full_path, recursive = TRUE)

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Rarefied_common_stats_no_smooth), function(x) Rarefied_common_stats_no_smooth[[x]][Rarefied_common_stats_no_smooth[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(Rarefied_covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})

# Plot correlations                                            
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(Beta = covs_grl_sig_common[[x]]$meth_cov / covs_grl_sig_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), "EM-Seq", "WGBS", file.path(full_path, "Unsmoothed_significant_CpGs_PCA"))

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_common, "Tally", "N", file.path(full_path, "Unsmoothed_significant_CpGs_tally"), c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_common, "Density", "N", file.path(full_path, "Unsmoothed_significant_CpGs_density"), c(-500, 4000))

# Plot heatmap. For row clustering, missing values use the group CpG average.
common_means_and_matrices_subset = common_means_and_matrices[[2]][c("evenness", "abs_delta_meth_pct", "Beta", "N")]
Plot_complex_heatmap(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), WGBS_Seq_names = c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), file.path(full_path, "Unsmoothed_significant_CpGs_heatmap"))

# Analyse motifs of heatmap.
seqnames_var_2 = sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", 1)
start_var_2 = as.numeric(sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", -1))
Coverage_dataframe_2 = data.frame(seqnames = seqnames_var_2, start = start_var_2)

GC_matrix_2 = Extract_motif_calculate_GC(Coverage_dataframe_2, nucleotides_backward = 3, nucleotides_forward = 4)

# Convert to GC% bin list of motifs
GC_matrix_list = list()
for (bin_var in sort(unique(GC_matrix_2$Motif_GC_percentage))) {
	
	GC_matrix_subset = GC_matrix_2[GC_matrix_2$Motif_GC_percentage == bin_var, ]
	temp_df = list(GC_matrix_subset)
	names(temp_df) = bin_var
	GC_matrix_list = append(GC_matrix_list, temp_df)
}

# Plot motif
counter <<- 0
Motif_plot_data = lapply(GC_matrix_list, Motif_frequency_table, length(GC_matrix_list))
Motif_plot(Motif_plot_data, file.path(full_path, "Motif_plots"), names(Motif_plot_data), TRUE)

# Plot correlations
Plot_correlation(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), WGBS_Seq_names = c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), file.path(full_path, "Correlation_plot"))

rm(Rarefied_covs_grl, Rarefied_common_bsseq, Significant_regions, covs_grl_sig_common, common_means_and_matrices, temp_df, common_means_and_matrices_subset, Rarefied_common_stats_smooth, Rarefied_common_stats_no_smooth)
gc()


#####  Stand bias analysis for *common* and significant CpGs  #####
#####  for *unsmoothed data*. Not rarefied data is used.      #####
###################################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/04_outputs/Not_rarefied_common_CpGs")
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Not_rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Not_rarefied_CpG_stats_common.RData"))

dir.create(full_path, recursive = TRUE)

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Not_rarefied_common_stats_no_smooth), function(x) Not_rarefied_common_stats_no_smooth[[x]][Not_rarefied_common_stats_no_smooth[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(Not_rarefied_covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})

# Plot correlations                                            
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(Beta = covs_grl_sig_common[[x]]$meth_cov / covs_grl_sig_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "EM-Seq", "WGBS", file.path(full_path, "Unsmoothed_significant_CpGs_PCA"))

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_common, "Tally", "N", file.path(full_path, "Unsmoothed_significant_CpGs_tally"), c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_common, "Density", "N", file.path(full_path, "Unsmoothed_significant_CpGs_density"), c(-500, 4000))

# Plot heatmap. For row clustering, missing values use the group CpG average.
common_means_and_matrices_subset = common_means_and_matrices[[2]][c("evenness", "abs_delta_meth_pct", "Beta", "N")]
Plot_complex_heatmap(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), WGBS_Seq_names = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), file.path(full_path, "Unsmoothed_significant_CpGs_heatmap"))

# Analyse motifs of heatmap.
seqnames_var_2 = sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", 1)
start_var_2 = as.numeric(sapply(strsplit(rownames(common_means_and_matrices_subset$Beta), "_"), "[", -1))
Coverage_dataframe_2 = data.frame(seqnames = seqnames_var_2, start = start_var_2)

GC_matrix_2 = Extract_motif_calculate_GC(Coverage_dataframe_2, nucleotides_backward = 3, nucleotides_forward = 4)

# Convert to GC% bin list of motifs
GC_matrix_list = list()
for (bin_var in sort(unique(GC_matrix_2$Motif_GC_percentage))) {
	
	GC_matrix_subset = GC_matrix_2[GC_matrix_2$Motif_GC_percentage == bin_var, ]
	temp_df = list(GC_matrix_subset)
	names(temp_df) = bin_var
	GC_matrix_list = append(GC_matrix_list, temp_df)
}

# Plot motif
counter <<- 0
Motif_plot_data = lapply(GC_matrix_list, Motif_frequency_table, length(GC_matrix_list))
Motif_plot(Motif_plot_data, file.path(full_path, "Motif_plots"), names(Motif_plot_data), TRUE)

# Plot correlations
Plot_correlation(common_means_and_matrices_subset, EM_Seq_names = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), WGBS_Seq_names = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), file.path(full_path, "Correlation_plot"))

rm(Not_rarefied_covs_grl, Not_rarefied_common_bsseq, Significant_regions, covs_grl_sig_common, common_means_and_matrices, temp_df, common_means_and_matrices_subset, Not_rarefied_common_stats_smooth, Not_rarefied_common_stats_no_smooth)
gc()

quit(save = "no")