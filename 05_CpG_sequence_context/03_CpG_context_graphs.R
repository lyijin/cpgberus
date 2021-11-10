#####  Script information  #####
################################

# This code graphs the statisical analysis produced from 02_CpG_context_stats.R
# Graphs include correlations and motifs of significant CpGs comparing 
# EM-Seq to WGBS (with and without smoothing). In addition, GC percentage and
# coverage are explored to see whether this affects significant CpG calls.
# GC percentage calculations are not based off Yi Jin's python script. By default
# the GC percentage is 3 nucleotides backward and 4 nucleotides forward from a C position,
# which is a 8 nucleotide motif region.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(motifStack)
library(bsseq)

#####  Load data  #####
#######################

# Load data just before analysis below and remove objects after to save on memory

# load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_not_subsampled.RData")
# load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_subsampled.RData")
# load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_common.RData")

#####  Functions  #####
#######################

#######################################################################
# # Funtion to plot extract motifs and calculcate CG percentage from
# # a dataframe of C positions
#
# The input dataframe needs to be the DMLtest results from 
# 02_CpG_context_stats.R script because column names are hardcoded.
#
# Todo: Line 50, make this cleaner, something like !(colnames %in%).
#######################################################################

Extract_motif_calculate_GC <- function(Input_dataframe, nucleotides_backward, nucleotides_forward) {

    Motif_grange = GRanges(seqnames = Input_dataframe$chr,
        IRanges(start = Input_dataframe$pos - nucleotides_backward,
        end = (Input_dataframe$pos + nucleotides_forward)))
    
    Final_dataframe_grange = GRanges(seqnames = Input_dataframe$chr,
        IRanges(start = Input_dataframe$pos, end = (Input_dataframe$pos + 1)))
        
    mcols(Final_dataframe_grange) <- Input_dataframe[c("mu1", "mu2", "diff", "diff.se", "stat", "phi1", "phi2", "pval", "fdr")]
    Final_dataframe_grange$cpg_context_nnncgnnn <- getSeq(BSgenome.Hsapiens.UCSC.hg38, Motif_grange)
    Final_dataframe_grange$GC_percentage <- as.numeric(letterFrequency(Final_dataframe_grange$cpg_context_nnncgnnn, letters = "GC", as.prob = TRUE))
    Final_dataframe_grange$High_low <- Final_dataframe_grange$diff
    Final_dataframe_grange$High_low[Final_dataframe_grange$High_low >= 0] <- "Low"
    Final_dataframe_grange$High_low[Final_dataframe_grange$High_low <= 0] <- "High"
    
    return(Final_dataframe_grange)
}

#######################################################################
# # Function to tabulate nucleotide species for each postion of motif
#######################################################################

Motif_frequency_table <- function(Granges_object, number_of_granges) {

    counter <<- counter + 1
    print(paste0("Analysing Grange object: ", as.character(counter), " / ", number_of_granges))
    
    # Create empty nucleotide frequency table
    cpg_context = unique(Granges_object$cpg_context_nnncgnnn)
    largest_motif = max(nchar(cpg_context))

    column_names = unique(unlist(strsplit(paste(cpg_context,collapse=""), "")))
    nucleotide_freq = data.frame(matrix(NA, nrow = 0, ncol = length(column_names)))
    colnames(nucleotide_freq) <- column_names
    
    cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = 1)

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

#######################################################################
# # Function to convert list of dataframes to motif class pcm
#######################################################################

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

#######################################################################
# # Funtion to motif plot by patient ID
#######################################################################

Motif_plot <- function(Motif_list_df, file_name, sample_names, remove_N) {

    pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    
    for (samp_name in sample_names) {
        Motif_list_df_temp = Motif_list_df[names(Motif_list_df) == samp_name]
        plot(Motif_stack_conversion(Motif_list_df_temp, remove_N)[[samp_name]], ic.scale=FALSE, ylab="probability")
    }
    
    dev.off()
}

#######################################################################
# # Plotting wrapper function for correlation
#######################################################################

Plot_cor_motif <- function(Granges_object, file_name) {

    pdf(paste0(file_name, ".pdf"), width = 8.3, height = 8.3)
    
    p = ggplot(data.frame(Granges_object), aes(x = mu1, y = mu2, color = GC_percentage))
    Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
      xlab("EM-Seq") + ylab("WGBS") + ggtitle("GC percentage")
    plot(Final_graph)

    p = ggplot(data.frame(Granges_object), aes(x = mu1, y = mu2, color = log2(cov_raw_average_m1)))
    Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
      xlab("EM-Seq") + ylab("WGBS") + ggtitle("EM-Seq raw coverage")
    plot(Final_graph)

    p = ggplot(data.frame(Granges_object), aes(x = mu1, y = mu2, color = log2(cov_raw_average_m2)))
    Final_graph = p + geom_point() + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
      xlab("EM-Seq") + ylab("WGBS") + ggtitle("WGBS raw coverage")
    plot(Final_graph)
    
    dev.off()
}

#######################################################################
# # Convert flattened dataframe of DSS stats to a granges with metadata
#######################################################################

Convert_df_to_grange <- function(Dataframe_object, bsseq_object, Sample_group_1, Sample_group_2) {

    # Extract significant regions and generate granges object with sequences / GC percentage
    Covs_grl_stats_flatten_sig = Dataframe_object[(Dataframe_object$fdr <= 0.05), ]
    Covs_grl_stats_flatten_sig = Extract_motif_calculate_GC(Covs_grl_stats_flatten_sig, nucleotides_backward = 3, nucleotides_forward = 4)

    # Extract coverages from bbseq object and add to granges object metadata
    raw_coverages_sig = getCoverage(bsseq_object, regions = Covs_grl_stats_flatten_sig)
    raw_coverages_sig = data.frame(do.call("rbind", raw_coverages_sig))
    colnames(raw_coverages_sig) <- sampleNames(bsseq_object)

    Covs_grl_stats_flatten_sig$cov_raw_average_m1 <- as.integer(ceiling(rowMeans(raw_coverages_sig[Sample_group_1], na.rm = TRUE)))
    Covs_grl_stats_flatten_sig$cov_raw_average_m2 <- as.integer(ceiling(rowMeans(raw_coverages_sig[Sample_group_2], na.rm = TRUE)))

    Covs_grl_stats_flatten_sig$cov_raw_average_m1[(Covs_grl_stats_flatten_sig$cov_raw_average_m1 == 0)] <- NA
    Covs_grl_stats_flatten_sig$cov_raw_average_m2[(Covs_grl_stats_flatten_sig$cov_raw_average_m2 == 0)] <- NA

    return(Covs_grl_stats_flatten_sig)
}

#######################################################
# # Wrapper function to pre-process and plot stats list 
# # object (DMLtest) from 02_CpG_context_stats.
#######################################################

Process_and_plot <- function(Stats_list_DMLtest, Bsseq_object_DMLtest, Sample_group_1, Sample_group_2, file_name_correlation, file_name_motif) {

    # Flatten list of stats dataframes, covert to grange with metadata
    Covs_grl_stats = do.call("rbind", Stats_list_DMLtest)
    Covs_grl_stats = Convert_df_to_grange(Covs_grl_stats, Bsseq_object_DMLtest, Sample_group_1, Sample_group_2)

    # Plot signficant scatter of high / low clusters
    print("Plotting correlations")
    Plot_cor_motif(Covs_grl_stats, file_name_correlation)

    # Plot signficant motifs of high / low clusters
    High_clust = Covs_grl_stats[Covs_grl_stats$High_low == "High", ]
    Low_clust = Covs_grl_stats[Covs_grl_stats$High_low == "Low", ]

    High_clust_no_NA = High_clust[!is.na(High_clust$cov_raw_average_m1), ]
    High_clust_no_NA = High_clust_no_NA[!is.na(High_clust_no_NA$cov_raw_average_m2), ]

    Low_clust_no_NA = Low_clust[!is.na(Low_clust$cov_raw_average_m1), ]
    Low_clust_no_NA = Low_clust_no_NA[!is.na(Low_clust_no_NA$cov_raw_average_m2), ]

    temp_grange_list = list(High = High_clust, Low = Low_clust, High_no_NA = High_clust_no_NA, Low_no_NA = Low_clust_no_NA)
    
    print("Plotting motifs")
    counter <<- 0
    Plot_data_meth <- lapply(temp_grange_list, Motif_frequency_table, length(temp_grange_list))
    Motif_plot(Plot_data_meth, file_name_motif, names(Plot_data_meth), TRUE)
    
    return(Covs_grl_stats)
}

#####  Analysis for significant CpGs between EM-Seq and WGBS #####
##################################################################

dir.create("All_existing_CpGs_GC_context")
load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_not_subsampled.RData")
# Pre process and plot smoothed significant CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_all_bsseq_list_stats, Covs_grl_all_bsseq, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "All_existing_CpGs_GC_context/Smoothed_significant_CpGs", "All_existing_CpGs_GC_context/Smoothed_motifs")
rm(temp)
rm(Covs_grl_all_bsseq_list_stats)

# Pre process and plot not smoothed significant CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_all_bsseq_list_stats_no_smooth, Covs_grl_all_bsseq, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "All_existing_CpGs_GC_context/Unsmoothed_significant_CpGs", "All_existing_CpGs_GC_context/Unsmoothed_motifs")
rm(temp)
rm(Covs_grl_all_bsseq_list_stats_no_smooth)
rm(Covs_grl_all_bsseq)


dir.create("Subsample_CpGs_GC_context")
load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_subsampled.RData")
# Pre process and plot smoothed significant subsampled CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_all_bsseq_list_stats_subsample, Covs_grl_all_bsseq_subsample, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Subsample_CpGs_GC_context/Smoothed_significant_CpGs_subsample", "Subsample_CpGs_GC_context/Smoothed_motifs_subsample")
rm(temp)
rm(Covs_grl_all_bsseq_list_stats_subsample)

# Pre process and plot not smoothed significant subsampled CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_all_bsseq_list_stats_no_smooth_subsample, Covs_grl_all_bsseq_subsample, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Subsample_CpGs_GC_context/Unsmoothed_significant_CpGs_subsample", "Subsample_CpGs_GC_context/Unsmoothed_motifs_subsample")
rm(temp)
rm(Covs_grl_all_bsseq_list_stats_no_smooth_subsample)
rm(Covs_grl_all_bsseq_subsample)


dir.create("Common_CpGs_GC_context")
load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats_common.RData")
# Pre process and plot smoothed significant common CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_common_bsseq_list_stats, Covs_grl_common_bsseq, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Common_CpGs_GC_context/Smoothed_significant_CpGs_common", "Common_CpGs_GC_context/Smoothed_motifs_common")
rm(temp)
rm(Covs_grl_common_bsseq_list_stats)

# Pre process and plot not smoothed significant common CpGs between EM-Seq and WGBS
temp = Process_and_plot(Covs_grl_common_bsseq_list_stats_no_smooth, Covs_grl_common_bsseq, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "Common_CpGs_GC_context/Unsmoothed_significant_CpGs_common", "Common_CpGs_GC_context/Unsmoothed_motifs_common")
rm(temp)
rm(Covs_grl_common_bsseq_list_stats_no_smooth)
rm(Covs_grl_common_bsseq)

quit(save = "no")