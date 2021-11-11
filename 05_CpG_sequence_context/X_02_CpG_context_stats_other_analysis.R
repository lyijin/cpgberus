#####  Analysis where samples with high coverage was subsampled #####
#####################################################################

# Find all CpGs which exist in all samples, merge and convert to bbseq object
Covs_grl_all_df = lapply(covs_grl_subsample, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov))
Covs_grl_all_df = Merge_CpGs(Covs_grl_all_df)
Covs_grl_all_bsseq_subsample = Create_BS_object(Covs_grl_all_df)
rm(covs_grl_subsample)
rm(Covs_grl_all_df)

# Smooth bbseq object per chromosome, return a list of bbseq smoothed objects.
# Covs_grl_all_bsseq_list_smooth = Create_smooth_fit_list_per_chr(Covs_grl_all_bsseq_subsample, 500, 5)

# Perform stats using DMLtest from DSS package
sample_group_1 = c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E")
sample_group_2 = c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W")
Covs_grl_all_bsseq_list_stats_subsample = Create_stat_list_per_chr(Covs_grl_all_bsseq_subsample, smoothing_bool = TRUE, smoothing_span = 500, number_of_workers = 8, group_1 = sample_group_1, group_2 = sample_group_2)
Covs_grl_all_bsseq_list_stats_no_smooth_subsample = Create_stat_list_per_chr(Covs_grl_all_bsseq_subsample, smoothing_bool = FALSE, number_of_workers = 8, group_1 = sample_group_1, group_2 = sample_group_2)

# Save workspace
save(Covs_grl_all_bsseq_subsample, Covs_grl_all_bsseq_list_stats_subsample, Covs_grl_all_bsseq_list_stats_no_smooth_subsample, file = "Motif_stats_subsampled.RData")
rm(Covs_grl_all_bsseq_subsample)
rm(Covs_grl_all_bsseq_list_stats_subsample)
rm(Covs_grl_all_bsseq_list_stats_no_smooth_subsample)


#####  Analysis for unfiltered coverage #####
#############################################

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
save(Covs_grl_all_bsseq, Covs_grl_all_bsseq_list_stats, Covs_grl_all_bsseq_list_stats_no_smooth, file = "Motif_stats_not_subsampled.RData")
rm(Covs_grl_all_bsseq)
rm(Covs_grl_all_bsseq_list_stats)
rm(Covs_grl_all_bsseq_list_stats_no_smooth)