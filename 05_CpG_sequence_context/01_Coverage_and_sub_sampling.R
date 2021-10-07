#####  Script information  #####
################################

# This code calculates coverage distribution and can also subsample using rbinom.
# This information can be used to make a coverage cutoff before statistical analysis.
# In addition, subsampling of high coverage samples can be performed before statistical analysis.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(ggplot2)
library(scales)

#####  Load data  #####
#######################

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")

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
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name], adjust = 5))
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

#####  Analysis using density and tally  #####
##############################################

# Total coverage calculation, convert to dataframe
Covs_grl_all_df = lapply(covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov, cpg_context_nnncgnnn = x$cpg_context_nnncgnnn,
                                                            evenness = x$evenness, abs_delta_meth_pct = x$abs_delta_meth_pct))
rm(covs_grl)

# Subset N (coverage) <= 100
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) x[x[, "N"] <= 100, ])

# Calculate and plot coverage density, coverage filtered < 100 because density calculations are affected.
Covs_grl_all_df_density = Variable_distribution(Covs_grl_all_df_subset, "Density", "N", "Coverage_analysis_density", c(-5, 70))

# Calculate and plot coverage tally without filtering for coverage
Covs_grl_all_df_tally = Variable_distribution(Covs_grl_all_df, "Tally", "N", "Coverage_analysis_tally", c(-5, 70))


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

# Save workspace
save(Covs_grl_all_df_density, Covs_grl_all_df_tally, Covs_grl_all_df_tally_subsample, covs_grl_subsample, file ="Data_coverage_filtered.RData")
quit(save = "no")