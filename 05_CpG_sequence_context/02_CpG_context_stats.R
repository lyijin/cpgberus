#####  Dependencies  #####

library(GenomicRanges)
library(motifStack)
library(bsseq)
library(DSS)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")
# load("/scratch1/gua020/CpGberus_data/motif_matrices.RData")

#####  Analysis  #####

# # # # Function to filter granges object in the following order: 1) Lower quartile discard; 2) meth percentage range; 3) Coverage cutoff's for meth / unmeth
# Filter_meth_pct <- function(Granges_object, discard_lower_quartile, meth_pct_keep_range, total_coverage_cutoff) {
    
    # temp_grange = Granges_object
    
    # if (discard_lower_quartile == TRUE) {
        # quartiles = quantile(temp_grange$meth_cov + temp_grange$unmeth_cov)
        # lower_quartile = quartiles[["25%"]]
        # temp_grange = temp_grange[(temp_grange$meth_cov + temp_grange$unmeth_cov) >= lower_quartile, ]
        # print(paste0("Discarded < 25% lower quartile which equates to coverage < ", lower_quartile))
    # }
    
    # temp_grange = subset(temp_grange, meth_pct >= min(meth_pct_keep_range) & meth_pct <= max(meth_pct_keep_range))
    # temp_grange = temp_grange[(temp_grange$meth_cov + temp_grange$unmeth_cov) >= total_coverage_cutoff, ]
    # return(temp_grange)
# }

# # Filter CpGs for quartile, meth and coverage
# # Subset_covs_grl_all = endoapply(covs_grl, Filter_meth_pct, TRUE, c(0, 100), 0)

# # # # Function to tabulate motifs
# # # # Additional option to mutliply motif by total coverage
# Motif_frequency_table <- function(Granges_object, number_of_granges, motif_coverage_multiplier) {

    # counter <<- counter + 1
    # print(paste0("Analysing Grange object: ", as.character(counter), " / ", number_of_granges))
    
    # # Multiply motif cpg_context by total coverage if TRUE
    # if (motif_coverage_multiplier == TRUE) {
        # motif_multiplier = Granges_object$meth_cov + Granges_object$unmeth_cov
        # cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = motif_multiplier)
    # } else {
        # cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = 1)
    # }

    # # Sum motifs and return dataframe
    # motif_table_sum = aggregate(Motif_multiplier ~ Motif, cpg_context, FUN=sum)
    # return(data.frame(motif_table_sum))
# }

# Plot_density <- function(Motif_list_df, file_name, sample_names, remove_N) {

    # if (remove_N == TRUE) {       
        # Motif_list_df = lapply(Motif_list_df, function(x) x[!grepl("N", x$Motif), ])
    # }
    
    # pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    
    # for (samp_name in sample_names) {
        # Motif_list_df_temp = Motif_list_df[names(Motif_list_df) == samp_name]
        # plot(density(Motif_list_df_temp[[samp_name]]$Motif_multiplier))
        # plot(density(log2(Motif_list_df_temp[[samp_name]]$Motif_multiplier)))
    # }
    
    # dev.off()    
# }

# # # Function to merge CpGs which exist between all samples, fill empty with 0
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

# # # Function here to create BS object
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

# # # Function here to create smooth fit by chromsome
Create_smooth_fit_list_per_chr <- function(bsseq_object, smoothing_span, number_of_workers) {

    final_list = list()
    
    for (chr_num in levels(seqnames(bsseq_object))) {
        print(paste0("Smoothing: ", chr_num))
        BSobj_temp_chr = chrSelectBSseq(bsseq_object, seqnames = chr_num, order = TRUE)
        
        # For parallel processing (MulticoreParam), it seems to be producing zombie processes causing memory leaks
        # Use SnowParam instead
        final_list[[chr_num]] <- BSmooth(BSseq = BSobj_temp_chr, h = smoothing_span, BPPARAM = SnowParam(workers = number_of_workers, progressbar = TRUE))
    }
    
    return(final_list)
}
    
# Find all CpGs which exist in all samples, merge and convert to bbseq object
Covs_grl_all_df = lapply(covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov))
Covs_grl_all_df = Merge_CpGs(Covs_grl_all_df)
Covs_grl_all_bsseq = Create_BS_object(Covs_grl_all_df)
rm(covs_grl)
rm(Covs_grl_all_df)

# Smooth bbseq object
Covs_grl_all_bsseq_list_smooth = Create_smooth_fit_list_per_chr(Covs_grl_all_bsseq, 500, 5)


# Perform stats with DDS package
# dmlTest = DMLtest(final_list$chrB, group1=c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), group2=c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), smoothing = FALSE)


# Covs_grl_all_bsseq.fit_per_chr = chrSelectBSseq(Covs_grl_all_bsseq, seqnames = "chr22", order = TRUE)

# Subset_covs_grl_all_bsseq = Create_BS_object(Subset_covs_grl_all_df)
# Subset_covs_grl_all_bsseq.fit <- BSmooth(BSseq = Subset_covs_grl_all_bsseq, h = 500, BPPARAM = MulticoreParam(workers = 4, progressbar = TRUE), verbose = TRUE)

# BSobj = makeBSseqData(Subset_covs_grl_all_df, names(Subset_covs_grl_all_df)) # This is too slow, do own function.



# BSobj.fit <- BSmooth(BSseq = BSobj, h = 500, BPPARAM = MulticoreParam(workers = 4, progressbar = TRUE), verbose = TRUE)

# Perform stats with DDS package
# dmlTest = DMLtest(BSobj.fit, group1=c("C1", "C2"), group2=c("N1", "N2"))


# regions <- GRanges(seqnames = c("chr18"), ranges = IRanges(start = 3014740, end = 3014769))
# getMeth(BSobj, regions, type = "smooth")

# DESeq2 Test
# cts = motifs$all
# coldata = data.frame(condition = factor(c("Enzymatic", "Bisulfite", "Enzymatic", "Bisulfite", "Enzymatic", "Bisulfite", "Enzymatic", "Bisulfite"), levels = c("Bisulfite","Enzymatic")))
# rownames(coldata) <- colnames(cts)

# dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
# dds = DESeq(dds)
# res = results(dds)
# res_order = res[order(res$padj), ]
# res_significant = rownames(res_order[res_order$padj <= 0.01, ])


# Data_all_sig <- endoapply(covs_grl, function(x) x[x$cpg_context_nnncgnnn %in% res_significant, ])
# counter = 0
# Plot_data_all_sig <- lapply(Data_all_sig, Motif_frequency_table, length(Data_all_sig), TRUE)
# Motif_plot(Plot_data_all_sig, "Multiply_all_significant", names(Plot_data_all_sig), TRUE)




# Consider normalising by frequency?
# counter = 0
# Motif_all <- lapply(Subset_covs_grl_all, Motif_frequency_table, length(Subset_covs_grl_all), FALSE)
# Plot_density(Motif_all, "No_multiply_density", names(Motif_all), TRUE)

# counter = 0
# Motif_all_multiply <- lapply(Subset_covs_grl_all, Motif_frequency_table, length(Subset_covs_grl_all), TRUE)
# Plot_density(Motif_all, "Multiply_density", names(Motif_all), TRUE)

# counter = 0
# Motif_all_no_filt <- lapply(covs_grl, Motif_frequency_table, length(covs_grl), FALSE)
# Plot_density(Motif_all_no_filt, "No_filt_multiply_density", names(Motif_all_no_filt), TRUE)

# # Note: I believe Jason doesn't filter on coverage before summing. 
# pdf("Test.pdf", width = 11.69, height = 8.3)
# plot(density(motifs$all[, 1]))
# plot(density(log2(motifs$all[, 1])))
# dev.off()