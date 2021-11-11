#####  Subsample and create granges list for stat analysis script  #####
########################################################################

# Subsample (using rbinom) high coverage data to be more equivalent to low coverage data.
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) cbind(x, data.frame(Subsample_coverage = x[, "N"])))

set.seed(1)
Covs_grl_all_df_subset[["WR069V1E"]]["Subsample_coverage"] = rbinom(nrow(Covs_grl_all_df_subset[["WR069V1E"]]), Covs_grl_all_df_subset[["WR069V1E"]][, "N"], 0.30)
set.seed(1)
Covs_grl_all_df_subset[["WR025V1E"]]["Subsample_coverage"] = rbinom(nrow(Covs_grl_all_df_subset[["WR025V1E"]]), Covs_grl_all_df_subset[["WR025V1E"]][, "N"], 0.55)
set.seed(1)
Covs_grl_all_df_subset[["WR025V9E"]]["Subsample_coverage"] = rbinom(nrow(Covs_grl_all_df_subset[["WR025V9E"]]), Covs_grl_all_df_subset[["WR025V9E"]][, "N"], 0.75)

# Remove subsample_coverage = 0 and add metadata required for stats script.
# After talking to Jason, it is probably more appropriate to subsample meth and unmeth. Implement this option in the next version.
# In addition, the gold standard would be to subsample raw FASTQ files, then process these.
# Most likely, I will just use the common CpGs as the final matrix for analysis.
Covs_grl_all_df_subset = lapply(Covs_grl_all_df_subset, function(x) x[x[, "Subsample_coverage"] > 0, ])
Covs_grl_all_df_subset = lapply(Covs_grl_all_df_subset, function(x) cbind(x, data.frame(meth_cov = x[, "X"], unmeth_cov = (x[, "N"] - x[, "X"]))))
                                    
# Correct meth / unmeth for subsampled data.
Samples_to_correct = c("WR069V1E", "WR025V1E", "WR025V9E")

for (data_name in Samples_to_correct) {
    print(paste0("Correcting sample: ", data_name))
    temp_data = Covs_grl_all_df_subset[[data_name]]
    meth_cov = ceiling(temp_data$X / temp_data$N * temp_data$Subsample_coverage)
    
    # Randomise CpG meth_cov if decimal is x.5, indicating that could be meth or unmeth
    pos_to_randomise = (meth_cov - (temp_data$X / temp_data$N * temp_data$Subsample_coverage)) == 0.5
    pos_to_randomise[pos_to_randomise == TRUE] = 1
    set.seed(1)
    pos_to_randomise = rbinom(length(pos_to_randomise), pos_to_randomise, 0.5)
    
    # Subtract randomise CpG from meth_cov; 1 = unmeth, 0 = meth. Also calculate unmeth
    meth_cov = meth_cov - pos_to_randomise
    unmeth_cov = temp_data$Subsample_coverage - meth_cov
    
    # Replace list of dataframes
    temp_data$meth_cov = as.integer(meth_cov)
    temp_data$unmeth_cov = as.integer(unmeth_cov)
    Covs_grl_all_df_subset[[data_name]] = temp_data
}

# Calculate and plot coverage tally without filtering for coverage.
Covs_grl_all_df_tally_subsample = Variable_distribution(Covs_grl_all_df_subset, "Tally", "Subsample_coverage", "Coverage_analysis_tally_subsample", c(-5, 70))

# Convert to granges list for stats analysis.
covs_grl_subsample = lapply(Covs_grl_all_df_subset, function(x) GRanges(seqnames = x$chr, ranges=IRanges(start = x$pos, end = (x$pos + 1)), 
                            meth_cov = x$meth_cov, unmeth_cov = x$unmeth_cov, cpg_context_nnncgnnn = x$cpg_context_nnncgnnn,
                            evenness = x$evenness, abs_delta_meth_pct = x$abs_delta_meth_pct))

covs_grl_subsample = do.call("GRangesList", covs_grl_subsample)