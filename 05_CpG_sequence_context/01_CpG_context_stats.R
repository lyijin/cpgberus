#####  Script information  #####
################################

# Performs statistical analysis comparing EM-Seq to WGBS, using the DSS R package
# with smoothing and no smoothing parameters.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(bsseq)
library(DSS)

#####  Load data  #####
#######################

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")

#####  Functions  #####
#######################

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
        print(paste0("Analysing: ", chr_num))
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

#####  Analysis  #####
######################

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
Covs_grl_all_bsseq_list_stats_no_smooth = Create_stat_list_per_chr(Covs_grl_all_bsseq, smoothing_bool = FALSE, number_of_workers = 8, group_1 = sample_group_1, group_2 = sample_group_2)

# Save workspace
rm(sample_group_1)
rm(sample_group_2)
rm(Create_BS_object)
rm(Create_smooth_fit_list_per_chr)
rm(Create_stat_list_per_chr)
rm(Merge_CpGs)

save.image("Motif_stats.RData")
quit(save = "no")