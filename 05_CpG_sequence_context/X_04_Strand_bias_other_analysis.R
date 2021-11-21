




#####  Stand bias analysis for *all existing* and significant CpGs  #####
#####  between samples for *smoothed data*. This differs from       #####
#####  common CpG because the CpG only has to exist in one sample   #####
#####  to be included in stats, while common CpGs has to exist in   #####
#####  all samples to be included in stats.                         #####
#########################################################################

# Create covs_grl (a Granges list) which contains raw data for all existing CpGs which exist for each sample.
covs_grl = lapply(covs_grl, function(x) cbind(data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov), data.frame(mcols(x))))
covs_grl = Merge_CpGs(covs_grl)
covs_grl = lapply(covs_grl, function(x) { temp = GRanges(seqnames = x$chr, ranges=IRanges(start = x$pos, end = (x$pos + 1)))
                                                values(temp) = x[!(colnames(x) %in% c("chr", "pos", "N", "X"))]
                                                return(temp)})                                            
                                                
covs_grl = do.call("GRangesList", covs_grl)

# Extract significant CpGs from DMLtest
Significant_regions = sapply(names(Covs_grl_all_bsseq_list_stats), function(x) Covs_grl_all_bsseq_list_stats[[x]][Covs_grl_all_bsseq_list_stats[[x]]$fdr <= 0.05, ], simplify = FALSE, USE.NAMES = TRUE)
Significant_regions = do.call("rbind", Significant_regions)
Significant_regions = GRanges(seqnames = Significant_regions$chr, IRanges(start = Significant_regions$pos, end = (Significant_regions$pos + 1)), 
                                mu1 = Significant_regions$mu1, mu2 = Significant_regions$mu2)

# Extract significant regions from unfiltered data and append mu1 and mu2 metadata
covs_grl_sig_common = lapply(covs_grl, mergeByOverlaps, Significant_regions)
covs_grl_sig_common = lapply(covs_grl_sig_common, function(x) {
                                            temp = x[, 1]
                                            values(temp) <- c(values(temp), x[c("mu1", "mu2")])
                                            return(temp)})
         
# Plot correlations         
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) data.frame(covs_grl_sig_common[[x]]), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(N = covs_grl_sig_common[[x]]$meth_cov + covs_grl_sig_common[[x]]$unmeth_cov)), simplify = FALSE, USE.NAMES = TRUE)
covs_grl_sig_common = sapply(names(covs_grl_sig_common), function(x) cbind(covs_grl_sig_common[[x]], data.frame(Beta = covs_grl_sig_common[[x]]$meth_cov / covs_grl_sig_common[[x]]$N)), simplify = FALSE, USE.NAMES = TRUE)
common_means_and_matrices = Plot_PCA_cor(covs_grl_sig_common, c("evenness", "abs_delta_meth_pct", "Beta", "N"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), "EM-Seq", "WGBS", "All_existing_CpGs_strand_bias/Smoothed_significant_CpGs_existing_strand_bias")

# Plot coverage distribution
temp_df = Variable_distribution(covs_grl_sig_common, "Tally", "N", "All_existing_CpGs_strand_bias/Smoothed_significant_CpGs_existing_tally", c(-5, 3000))
temp_df = Variable_distribution(covs_grl_sig_common, "Density", "N", "All_existing_CpGs_strand_bias/Smoothed_significant_CpGs_existing_density", c(-500, 4000))

# Filter matrix CpGs where only one sample contains NA for EM-Seq and WGBS groups, then plot heatmaps.
# Anything higher than this threshold will result in a clustering error when plotting heatmaps.
common_means_and_matrices_subset = common_means_and_matrices[[2]][c("evenness", "abs_delta_meth_pct", "Beta", "N")]
common_means_and_matrices_subset = lapply(common_means_and_matrices_subset, Filter_matrix_NA, c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"), c("WR025V1W", "WR025V9W", "WR069V1W", "WR069V9W"), 1, 1)
Plot_complex_heatmap(common_means_and_matrices_subset, "All_existing_CpGs_strand_bias/Smoothed_signficant_CpGs_existing")