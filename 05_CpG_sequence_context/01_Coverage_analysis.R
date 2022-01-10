#####  Script information  #####
################################

# This code calculates coverage distribution.
# This information can be used to make a coverage cutoff before statistical analysis.

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(ggplot2)
library(scales)

#####  Load data  #####
#######################

path_to_cpgerus = "/datasets/work/hb-stopwatch/work/cpgberus"

load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Not_rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/Rarefied_grch38p13_combined_covs_grl.RData"))
load(file.path(path_to_cpgerus, "04_parse_bismark_covs/grch38p13_combined_covs_grl.RData"))

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
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name], na.rm = TRUE))
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
    Rarefaction_position = grepl("R$", temp_df$Sample_name)
    Not_rarefaction_position = !Rarefaction_position
    
    temp_df$Sample_name_2 = "temp_string"
    temp_df$Sample_name_2[Rarefaction_position] = substr(temp_df$Sample_name[Rarefaction_position], 1, nchar(temp_df$Sample_name[Rarefaction_position]) - 2)
    temp_df$Sample_name_2[Not_rarefaction_position] = substr(temp_df$Sample_name[Not_rarefaction_position], 1, nchar(temp_df$Sample_name[Not_rarefaction_position]) - 1)

    temp_df$Seq_type = temp_df$Sample_name
    temp_df$Seq_type[grep("E$|ER$", temp_df$Sample_name)] = "EM-Seq"
    temp_df$Seq_type[grep("W$|WR$", temp_df$Sample_name)] = "WGBS"

    # Plot
    pdf(paste0(plot_name, ".pdf"), width = 11.69, height = 8.3)
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    plot(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    return(temp_df)
}


#####  Analysis using density and tally (not rarefied)  #####
#############################################################

full_path = file.path(path_to_cpgerus, "05_CpG_sequence_context/01_outputs")
dir.create(full_path)

# Total coverage calculation, convert to dataframe
Covs_grl_all_df = lapply(Not_rarefied_covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov, cpg_context_nnncgnnn = x$cpg_context_nnncgnnn,
                                                            evenness = x$evenness, abs_delta_meth_pct = x$abs_delta_meth_pct))
rm(Not_rarefied_covs_grl)

# Subset N (coverage) <= 100
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) x[x[, "N"] <= 100, ])

# Calculate and plot coverage density, coverage filtered < 100 because density calculations are affected.
Not_rarefied_covs_grl_all_df_density = Variable_distribution(Covs_grl_all_df_subset, "Density", "N", file.path(full_path, "Not_rarefied_coverage_analysis_density"), c(-5, 70))

# Calculate and plot coverage tally without filtering for coverage
Not_rarefied_covs_grl_all_df_tally = Variable_distribution(Covs_grl_all_df, "Tally", "N", file.path(full_path, "Not_rarefied_coverage_analysis_tally"), c(-5, 70))

rm(Covs_grl_all_df)


#####  Analysis using density and tally (rarefied)  #####
#########################################################

# Total coverage calculation, convert to dataframe
Covs_grl_all_df = lapply(Rarefied_covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov, cpg_context_nnncgnnn = x$cpg_context_nnncgnnn,
                                                            evenness = x$evenness, abs_delta_meth_pct = x$abs_delta_meth_pct))
rm(Rarefied_covs_grl)

# Subset N (coverage) <= 100
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) x[x[, "N"] <= 100, ])

# Calculate and plot coverage density, coverage filtered < 100 because density calculations are affected.
Rarefied_covs_grl_all_df_density = Variable_distribution(Covs_grl_all_df_subset, "Density", "N", file.path(full_path, "Rarefied_coverage_analysis_density"), c(-5, 70))

# Calculate and plot coverage tally without filtering for coverage
Rarefied_covs_grl_all_df_tally = Variable_distribution(Covs_grl_all_df, "Tally", "N", file.path(full_path, "Rarefied_coverage_analysis_tally"), c(-5, 70))

rm(Covs_grl_all_df)


#####  Analysis using density and tally (original)  #####
#########################################################

# Total coverage calculation, convert to dataframe
Covs_grl_all_df = lapply(covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov, cpg_context_nnncgnnn = x$cpg_context_nnncgnnn,
                                                            evenness = x$evenness, abs_delta_meth_pct = x$abs_delta_meth_pct))
rm(covs_grl)

# Subset N (coverage) <= 100
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) x[x[, "N"] <= 100, ])

# Calculate and plot coverage density, coverage filtered < 100 because density calculations are affected.
Original_covs_grl_all_df_density = Variable_distribution(Covs_grl_all_df_subset, "Density", "N", file.path(full_path, "Original_coverage_analysis_density"), c(-5, 70))

# Calculate and plot coverage tally without filtering for coverage
Original_covs_grl_all_df_tally = Variable_distribution(Covs_grl_all_df, "Tally", "N", file.path(full_path, "Original_coverage_analysis_tally"), c(-5, 70))

# Save workspace
save(Not_rarefied_covs_grl_all_df_density, Not_rarefied_covs_grl_all_df_tally, Rarefied_covs_grl_all_df_density, Rarefied_covs_grl_all_df_tally, Original_covs_grl_all_df_density, Original_covs_grl_all_df_tally,
    file = file.path(full_path, "Coverage_analysis.RData"))

quit(save = "no")