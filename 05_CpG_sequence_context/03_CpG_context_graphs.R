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
library(bsseq)
library(showtext)

font_add_google("Source Sans Pro", "source_sans")
showtext_auto()
showtext_opts(dpi = 300)

#####  Load data  #####
#######################

path_to_cpgerus = "/scratch/user/uqdguanz/Projects/Meth/cpgberus"
source(paste0(path_to_cpgerus, "/05_CpG_sequence_context/07_Motif_functions.R"))

# Load data just before analysis below and remove objects after to save on memory

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
    Final_dataframe_grange$GC_percentage <- as.numeric(letterFrequency(Final_dataframe_grange$cpg_context_nnncgnnn, letters = "GC", as.prob = TRUE)) * 100
    Final_dataframe_grange$High_low <- Final_dataframe_grange$diff
    Final_dataframe_grange$High_low[Final_dataframe_grange$High_low >= 0] <- "Low"
    Final_dataframe_grange$High_low[Final_dataframe_grange$High_low <= 0] <- "High"
    
    return(Final_dataframe_grange)
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


#####  Analysis for *existing* significant CpGs between EM-Seq and WGBS (Rarefied) #####
########################################################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/03_outputs/Rarefied_existing_CpGs_GC_context")
dir.create(full_path, recursive = TRUE)

load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Rarefied_CpG_stats_existing.RData"))
# Pre process and plot smoothed significant common CpGs between EM-Seq and WGBS
temp = Process_and_plot(Rarefied_existing_stats_smooth, Rarefied_existing_bsseq, c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), 
                        file.path(full_path, "Smoothed_significant_CpGs_existing"), file.path(full_path, "Smoothed_motifs_existing"))

rm(temp, Rarefied_existing_stats_smooth)
gc()

# Pre process and plot not smoothed significant common CpGs between EM-Seq and WGBS
temp = Process_and_plot(Rarefied_existing_stats_no_smooth, Rarefied_existing_bsseq, c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER"), c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR"), 
                        file.path(full_path, "Unsmoothed_significant_CpGs_existing"), file.path(full_path, "Unsmoothed_motifs_existing"))

rm(temp, Rarefied_existing_stats_no_smooth, Rarefied_existing_bsseq)
gc()


quit(save = "no")
