#####  Script information  #####
################################

# This code graphs the statisical analysis produced from 01_CpG_context_stats.R
# Graphs include correlations and motifs of significant CpGs comparing 
# EM-Seq to WGBS (with and without smoothing). In addition, GC percentage and
# coverage are explored to see whether this affects significant CpG calls.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(motifStack)
library(bsseq)

#####  Load data  #####
#######################

load("/home/gua020/Development/CPGberus/cpgberus/05_CpG_sequence_context/Motif_stats.RData")

#####  Functions  #####
#######################

#######################################################################
# # Funtion to plot extract motifs and calculcate CG percentage from
# # a dataframe of C positions
#######################################################################

Extract_motif_calculate_GC <- function(Input_dataframe, nucleotides_backward, nucleotides_forward) {

    Motif_grange = GRanges(seqnames = Input_dataframe$chr,
        IRanges(start = Input_dataframe$pos - nucleotides_backward,
        end = (Input_dataframe$pos + nucleotides_forward)))
    
    Final_dataframe_grange = GRanges(seqnames = Input_dataframe$chr,
        IRanges(start = Input_dataframe$pos, end = (Input_dataframe$pos + 1)))
        
    mcols(Final_dataframe_grange) <- Input_dataframe[c("mu1", "mu2", "diff", "diff.se","stat", "phi1", "phi2", "pval", "fdr")]
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

Convert_df_to_grange <- function(Granges_object) {

    # Extract significant regions and generate granges object with sequences / GC percentage
    Covs_grl_stats_flatten_sig = Granges_object[(Granges_object$fdr <= 0.05), ]
    Covs_grl_stats_flatten_sig = Extract_motif_calculate_GC(Covs_grl_stats_flatten_sig, nucleotides_backward = 3, nucleotides_forward = 4)

    # Extract coverages from bbseq object and add to granges object metadata
                    #### Note: Have to QC to check if coverage agrees with Jason's original granges object ####
    raw_coverages_sig = getCoverage(Covs_grl_all_bsseq, regions = Covs_grl_stats_flatten_sig)
    raw_coverages_sig = data.frame(do.call("rbind", raw_coverages_sig))
    colnames(raw_coverages_sig) <- sampleNames(Covs_grl_all_bsseq)

    Covs_grl_stats_flatten_sig$cov_raw_average_m1 <- as.integer(round(rowMeans(raw_coverages_sig[c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E")])))
    Covs_grl_stats_flatten_sig$cov_raw_average_m2 <- as.integer(round(rowMeans(raw_coverages_sig[c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W")])))

    Covs_grl_stats_flatten_sig$cov_raw_average_m1[(Covs_grl_stats_flatten_sig$cov_raw_average_m1 == 0)] <- NA
    Covs_grl_stats_flatten_sig$cov_raw_average_m2[(Covs_grl_stats_flatten_sig$cov_raw_average_m2 == 0)] <- NA

    return(Covs_grl_stats_flatten_sig)
}


#####  Analysis for smoothed significant CpGs between EM-Seq and WGBS #####
###########################################################################

# Flatten list of stats dataframes, covert to grange with metadata
Covs_grl_stats_grange_smooth = do.call("rbind", Covs_grl_all_bsseq_list_stats)
Covs_grl_stats_grange_smooth = Convert_df_to_grange(Covs_grl_stats_grange_smooth)

# Plot signficant scatter of high / low clusters
Plot_cor_motif(Covs_grl_stats_grange_smooth, "Smoothed_significant_CpGs")

# Plot signficant motifs of high / low clusters
temp_grange_list = list(High = Covs_grl_stats_grange_smooth[Covs_grl_stats_grange_smooth$High_low == "High", ], 
    Low = Covs_grl_stats_grange_smooth[Covs_grl_stats_grange_smooth$High_low == "Low", ],
    EM_seq_NA = Covs_grl_stats_grange_smooth[is.na(Covs_grl_stats_grange_smooth$cov_raw_average_m1), ],
    WGBS_NA = Covs_grl_stats_grange_smooth[is.na(Covs_grl_stats_grange_smooth$cov_raw_average_m2), ])
    
counter = 0
Plot_data_meth <- lapply(temp_grange_list, Motif_frequency_table, length(temp_grange_list))
Motif_plot(Plot_data_meth, "Smoothed_motifs", names(Plot_data_meth), TRUE)


#####  Analysis for not smoothed significant CpGs between EM-Seq and WGBS #####
###############################################################################

# Flatten list of stats dataframes
Covs_grl_stats_grange_no_smooth = do.call("rbind", Covs_grl_all_bsseq_list_stats_no_smooth)
Covs_grl_stats_grange_no_smooth = Convert_df_to_grange(Covs_grl_stats_grange_no_smooth)

# Plot signficant scatter of high / low clusters
Plot_cor_motif(Covs_grl_stats_grange_no_smooth, "Not_smoothed_significant_CpGs")

# Plot signficant motifs of high / low clusters
temp_grange_list = list(High = Covs_grl_stats_grange_no_smooth[Covs_grl_stats_grange_no_smooth$High_low == "High", ], 
    Low = Covs_grl_stats_grange_no_smooth[Covs_grl_stats_grange_no_smooth$High_low == "Low", ],
    EM_seq_NA = Covs_grl_stats_grange_no_smooth[is.na(Covs_grl_stats_grange_no_smooth$cov_raw_average_m1), ],
    WGBS_NA = Covs_grl_stats_grange_no_smooth[is.na(Covs_grl_stats_grange_no_smooth$cov_raw_average_m2), ])
    
counter = 0
Plot_data_meth <- lapply(temp_grange_list, Motif_frequency_table, length(temp_grange_list))
Motif_plot(Plot_data_meth, "Not_smoothed_motifs", names(Plot_data_meth), TRUE)

quit(save = "no")