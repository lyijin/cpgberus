#####  Dependencies  #####

library(GenomicRanges)
library(bsseq)
library(DSS)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")
# load("/scratch1/gua020/CpGberus_data/motif_matrices.RData")

################################################################################
# # Function to merge CpGs which exist between all samples, fill empty with 0
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
        temp_df = data.frame(chr = factor(), pos = integer(), N = integer(), X = integer())

        for (chr_num in number_of_chromosomes) {
            chr_num_pos = data_frame[(data_frame$chr == chr_num), "pos"]
            chr_num_pos_ref = existing_CpGs[(existing_CpGs$chr == chr_num), "pos"]
            chr_num_pos_missing = chr_num_pos_ref[!(chr_num_pos_ref %in% chr_num_pos)]
            
            if (length(chr_num_pos_missing) != 0) {
                temp_df = rbind(temp_df, data.frame(chr = factor(chr_num), pos = chr_num_pos_missing, N = 0, X = 0))
            }
        }
        
        temp_df = rbind(data_frame, temp_df)
        temp_df = temp_df[order(temp_df$chr, temp_df$pos), ]
        rownames(temp_df) <- seq(1:nrow(temp_df))
        final_dataframe_list[[data_frame_name]] <- temp_df
        print(paste0(data_frame_name, " data is finished merging: ", dataframe_index, " / ", length(Dataframe_list_object)))
    }
    return(final_dataframe_list)
}


####################################################################
# # Function here to create BS object
####################################################################

Create_BS_object <- function(Dataframe_list_object) {

    # QC to check if chr and pos are identical for all samples, then cbind
    chr_ref = Dataframe_list_object[[1]]$chr
    pos_ref = Dataframe_list_object[[1]]$pos
    row_ref = rownames(Dataframe_list_object[[1]])
    temp_mat_N = matrix( , nrow = length(row_ref), ncol = 0)
    temp_mat_X = matrix( , nrow = length(row_ref), ncol = 0)
    
    for (dataframe_index in seq(length(Dataframe_list_object))) {

        if (identical(Dataframe_list_object[[dataframe_index]]$chr, chr_ref) == FALSE) {
            stop(paste0("Column chr (Sample: ", names(Dataframe_list_object[dataframe_index]), ") doesn't match reference column chr (Sample: ", names(Dataframe_list_object[1]), ")"))
        } 
        
        if (identical(Dataframe_list_object[[dataframe_index]]$pos, pos_ref) == FALSE) {
            stop(paste0("Column pos (Sample: ", names(Dataframe_list_object[dataframe_index]), ") doesn't match reference column pos (Sample: ", names(Dataframe_list_object[1]), ")"))
        }
        
        temp_mat_N = cbind(temp_mat_N, matrix(Dataframe_list_object[[dataframe_index]]$N))
        temp_mat_X = cbind(temp_mat_X, matrix(Dataframe_list_object[[dataframe_index]]$X))
    }
    
    BStmp <- BSseq(chr = chr_ref, pos = pos_ref, M = temp_mat_X, Cov = temp_mat_N, sampleNames = names(Dataframe_list_object))
}

##############################################################
# # Function here to create smooth fit by chromsome

# This was created so that the BS smoothed values can be put into
# DMLtest statistics function. DMLtest itself uses a smoothing
# average algorithim before statistics, where smoothing can be 
# turned off. However, I don't think there is a way to extract the
# smoothed values, hence the creation of this function. I haven't 
# tested whether substitution with BS smoothed values into DMLtest 
# is correct.
##############################################################

Create_smooth_fit_list_per_chr <- function(bsseq_object, smoothing_span, number_of_workers) {

    final_list = list()
    
    for (chr_num in levels(seqnames(bsseq_object))) {
        print(paste0("Smoothing: ", chr_num))
        BSobj_temp_chr = chrSelectBSseq(bsseq_object, seqnames = chr_num, order = TRUE)
        
        # For parallel processing (MulticoreParam), it seems to be producing zombie processes causing memory leaks
        # I don't think MulticoreParam (fork) can be used on the HPC clusters
        # Use SnowParam instead (sock)
        final_list[[chr_num]] <- BSmooth(BSseq = BSobj_temp_chr, h = smoothing_span, BPPARAM = SnowParam(workers = number_of_workers, progressbar = TRUE))
    }
    
    return(final_list)
}

##############################################################
# # Function here to perform stats by chromsome

# Approximately 1.2 gb per chromosome / worker (sample).
# On one HPC node, could run theoretically 20 workers (20 cpu's).
# Should be fine to run 20 samples per node (20gb needed approx).
# Create option to parallelize chromosome loop.
##############################################################

Create_stat_list_per_chr <- function(bsseq_object, smoothing_bool, smoothing_span, number_of_workers, group_1, group_2) {

    final_list = list()

    for (chr_num in levels(seqnames(bsseq_object))) {
        print(paste0("Smoothing: ", chr_num))
        BSobj_temp_chr = chrSelectBSseq(bsseq_object, seqnames = chr_num, order = TRUE)
        
        # For parallel processing (MulticoreParam), it seems to be producing zombie processes causing memory leaks
        # I don't think MulticoreParam (fork) can be used on the HPC clusters
        # Use SnowParam instead (sock)
        if (smoothing_bool == TRUE) {
            final_list[[chr_num]] <- DMLtest(BSobj_temp_chr, group1 = group_1, group2 = group_2, BPPARAM = SnowParam(workers = number_of_workers, progressbar = TRUE), smoothing.span = smoothing_span, smoothing = smoothing_bool)
        }
        if (smoothing_bool == FALSE) {
            final_list[[chr_num]] <- DMLtest(BSobj_temp_chr, group1 = group_1, group2 = group_2, BPPARAM = SnowParam(workers = number_of_workers, progressbar = TRUE), smoothing = smoothing_bool)
        }
    }
    
    return(final_list)
}

#######################################################################
# # Function to tabulate nucleotide species for each postion of motif
# # Additional option to mutliply motif by total coverage
#######################################################################

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
    
    # plot.new()
    # motifStack(Motif_stack_conversion(Motif_list_df, remove_N), layout="tree", ic.scale=FALSE, ylab="probability")
    dev.off()
}


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


#####  Analysis  #####

# Find all CpGs which exist in all samples, merge and convert to bbseq object
Covs_grl_all_df = lapply(covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov))
Covs_grl_all_df = Merge_CpGs(Covs_grl_all_df)
Covs_grl_all_bsseq = Create_BS_object(Covs_grl_all_df)
rm(covs_grl)
rm(Covs_grl_all_df)

# Smooth bbseq object per chromosome, return a list of bbseq smoothed objects.
# Covs_grl_all_bsseq_list_smooth = Create_smooth_fit_list_per_chr(Covs_grl_all_bsseq, 500, 5)

# Perform stats using DMLtest from DSS package
sample_group_1 = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E")
sample_group_2 = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W")
Covs_grl_all_bsseq_list_stats = Create_stat_list_per_chr(Covs_grl_all_bsseq, smoothing_bool = TRUE, smoothing_span = 500, number_of_workers = 8, group_1 = sample_group_1, group_2 = sample_group_2)



#####  Graphing, move to motif script  #####

library(GenomicRanges)
library(motifStack)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)

# Flatten list of stats dataframes
Covs_grl_stats_flatten = do.call("rbind", Covs_grl_all_bsseq_list_stats)

# Extract significant regions and generate granges object with sequences / GC percentage
Covs_grl_stats_flatten_sig = Covs_grl_stats_flatten[(Covs_grl_stats_flatten$fdr <= 0.05), ]
Covs_grl_stats_flatten_sig = Extract_motif_calculate_GC(Covs_grl_stats_flatten_sig, nucleotides_backward = 3, nucleotides_forward = 4)

# Plot signficant scatter and motifs of high / low clusters
p <- ggplot(data.frame(Covs_grl_stats_flatten_sig), aes(x = mu1, y = mu2, color = GC_percentage))
Final_graph = p + geom_point() + theme_bw() + theme(legend.position="top") +
  xlab("EM-Seq") + ylab("WGBS")
ggsave("Significant_correlation_scatter.pdf", Final_graph)

temp_grange_list = list(High = Covs_grl_stats_flatten_sig[Covs_grl_stats_flatten_sig$High_low == "High", ], 
    Low = Covs_grl_stats_flatten_sig[Covs_grl_stats_flatten_sig$High_low == "Low", ])
    
counter = 0
Plot_data_meth <- lapply(temp_grange_list, Motif_frequency_table, length(temp_grange_list), FALSE)
Motif_plot(Plot_data_meth, "High_low_motif", names(Plot_data_meth), TRUE)

# Plot raw coverages of significant CpGs
# QC check this with Jason meth matrix
raw_coverages_sig = getCoverage(Covs_grl_all_bsseq, regions = Covs_grl_stats_flatten_sig)
raw_coverages_sig = data.frame(do.call("rbind", raw_coverages_sig))
colnames(raw_coverages_sig) <- sampleNames(Covs_grl_all_bsseq)

Covs_grl_stats_flatten_sig$cov_raw_average_m1 <- as.integer(round(rowMeans(raw_coverages_sig[c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E")])))
Covs_grl_stats_flatten_sig$cov_raw_average_m2 <- as.integer(round(rowMeans(raw_coverages_sig[c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W")])))

Covs_grl_stats_flatten_sig$cov_raw_average_m1[(Covs_grl_stats_flatten_sig$cov_raw_average_m1 == 0)] <- NA
Covs_grl_stats_flatten_sig$cov_raw_average_m2[(Covs_grl_stats_flatten_sig$cov_raw_average_m2 == 0)] <- NA

p <- ggplot(data.frame(Covs_grl_stats_flatten_sig), aes(x = mu1, y = mu2, color = log2(cov_raw_average_m1)))
Final_graph = p + geom_point() + theme_bw() + theme(legend.position="top") +
  xlab("EM-Seq") + ylab("WGBS")
ggsave("Significant_correlation_scatter_EM_seq.pdf", Final_graph)

p <- ggplot(data.frame(Covs_grl_stats_flatten_sig), aes(x = mu1, y = mu2, color = log2(cov_raw_average_m2)))
Final_graph = p + geom_point() + theme_bw() + theme(legend.position="top") +
  xlab("EM-Seq") + ylab("WGBS")
ggsave("Significant_correlation_scatter_WGBS.pdf", Final_graph)

p <- ggplot(data.frame(Covs_grl_stats_flatten_sig), aes(x = log2(cov_raw_average_m1), y = log2(cov_raw_average_m2), color = GC_percentage))
Final_graph = p + geom_point() + theme_bw() + theme(legend.position="top") +
  xlab("EM-Seq") + ylab("WGBS")
ggsave("Significant_correlation_scatter_coverage.pdf", Final_graph)