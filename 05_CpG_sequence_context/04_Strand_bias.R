#####  Script information  #####
################################

# Explore eveness and delta meth between watson and crick.
# See if there is any differences between EM-Seq to WGBS

#####  Dependencies  #####
##########################

library(GenomicRanges)
library(ggplot2)
library(scales)
library(reshape2)

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
        temp_df = lapply(List_of_dataframes, function(x) density(x[, variable_column_name], adjust = 2))
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
    png(paste0("01_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0("02_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0("03_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0("04_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=y_values, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + facet_wrap( ~ Sample_name_2, ncol = 2) + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    png(paste0("05_", plot_name, ".png"), width = 11.69, height = 8.3, units = "in", res = 300)
    print(ggplot(temp_df, aes(x=x_values, y=log2(y_values), color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = x_axis_cutoff) + theme_bw() + facet_wrap( ~ Sample_name_2, ncol = 2) + xlab(variable_column_name) + ylab(density_or_tally))
    dev.off()
    
    return(temp_df)
}

#############################################################
# # Function to filter granges list object for common CpGs
#############################################################

Filter_common_CpGs <- function(Granges_list_object) {

    Common_intersect_regions = Reduce(intersect, Granges_list_object)
    Subset_covs_grl_common = endoapply(Granges_list_object, subsetByOverlaps, Common_intersect_regions)
    return(Subset_covs_grl_common)
}

###################################################################
# # Function to calculate the delta between samples and groups for:
# # 1) Coverage evenness
# # 2) Deviation in absolute methylation levels
# # Above are calculated in Yi Jin's python script.
###################################################################

Calculate_delta_evenness_meth_pct <- function(List_of_dataframes_common, biological_sample_names, reference_sample_names) {
    
    samp_id_list = names(List_of_dataframes_common)
    
    # QC to check if "chr" and "pos" columns are identical
    for (samp_id in samp_id_list) {
        if (identical(List_of_dataframes_common[samp_id_list[1]]$chr, List_of_dataframes_common[samp_id]$chr) == FALSE) {
            stop(paste0("Chromosomes not identical between samples: ", samp_id_list[1], " ", samp_id))
        } else if (identical(List_of_dataframes_common[samp_id_list[1]]$pos, List_of_dataframes_common[samp_id]$pos) == FALSE) {
            stop(paste0("CpG position not identical between samples: ", samp_id_list[1], " ", samp_id))
        }
    }
    
    # Calculate delta evenness and meth pct between groups
    mean_reference = List_of_dataframes_common[reference_sample_names]
    mean_reference = lapply(mean_reference, function(x) data.frame(group_delta_evenness = x$evenness))
    mean_reference = rowMeans(do.call(cbind, mean_reference))
    
    mean_condition = List_of_dataframes_common[samp_id_list[!(samp_id_list %in% reference_sample_names)]]
    mean_condition = lapply(mean_condition, function(x) data.frame(group_delta_evenness = x$evenness))
    mean_condition = rowMeans(do.call(cbind, mean_condition))
    
    evenness_diff = mean_reference - mean_condition
    
    mean_reference = List_of_dataframes_common[reference_sample_names]
    mean_reference = lapply(mean_reference, function(x) data.frame(group_abs_delta_meth_pct = x$abs_delta_meth_pct))
    mean_reference = rowMeans(do.call(cbind, mean_reference))
    
    mean_condition = List_of_dataframes_common[samp_id_list[!(samp_id_list %in% reference_sample_names)]]
    mean_condition = lapply(mean_condition, function(x) data.frame(group_abs_delta_meth_pct = x$abs_delta_meth_pct))
    mean_condition = rowMeans(do.call(cbind, mean_condition))
    
    meth_pct_diff = mean_reference - mean_condition
    
    # Calculate delta evenness and meth pct between samples
    for (samp_id in biological_sample_names) {
        samples_to_analyse = grep(samp_id, samp_id_list, value = TRUE)
        control_sample = grep(paste0(reference_sample_names, collapse = "|"), samples_to_analyse, value = TRUE)
        
        for (samp_id_2 in samples_to_analyse) {
            delta_evenness = List_of_dataframes_common[[control_sample]]$evenness - List_of_dataframes_common[[samp_id_2]]$evenness
            delta_meth_pct_sample = List_of_dataframes_common[[control_sample]]$abs_delta_meth_pct - List_of_dataframes_common[[samp_id_2]]$abs_delta_meth_pct
            List_of_dataframes_common[[samp_id_2]] = cbind(List_of_dataframes_common[[samp_id_2]], data.frame(sample_delta_evenness = delta_evenness, sample_abs_delta_meth_pct = delta_meth_pct_sample,
                                                            group_delta_evenness = evenness_diff, group_abs_delta_meth_pct = meth_pct_diff))
        }
    }
    
    return(List_of_dataframes_common)
}

#####  Plot graphs  #####
#########################

# Filter for common CpGs and calculate deltas between samples and groups
Covs_grl_all_df_common = Filter_common_CpGs(covs_grl)
#Covs_grl_all_df_common_test = endoapply(Covs_grl_all_df_common,"[",1:1000000)
Covs_grl_all_df_common = sapply(names(Covs_grl_all_df_common), function(x) data.frame(chr = seqnames(Covs_grl_all_df_common[[x]]), pos = start(Covs_grl_all_df_common[[x]]),
                        N = (Covs_grl_all_df_common[[x]]$meth_cov + Covs_grl_all_df_common[[x]]$unmeth_cov), X = Covs_grl_all_df_common[[x]]$meth_cov, 
                        evenness = Covs_grl_all_df_common[[x]]$evenness, abs_delta_meth_pct = Covs_grl_all_df_common[[x]]$abs_delta_meth_pct, sample_name = x),
                        simplify = FALSE, USE.NAMES = TRUE)
            
Covs_grl_all_df_common = Calculate_delta_evenness_meth_pct(Covs_grl_all_df_common, c("WR025V1", "WR025V9", "WR069V1", "WR069V9"), c("WR025V1E", "WR025V9E", "WR069V1E", "WR069V9E"))

# Analyse and plot group delta's, use first element in list because data is the same
Covs_grl_all_df_tally_abs_delta_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_delta_evenness", "Group_Evenness_tally", c(-0.5, 0.5))
Covs_grl_all_df_tally_evenness_group = Variable_distribution(Covs_grl_all_df_common[1], "Density", "group_abs_delta_meth_pct", "Group_Abs_delta_meth_pct_tally", c(-100, 100))

# Extract WGBS samples (EM-Seq is reference, so all zero's) analyse and plot sample delta's
Covs_grl_all_df_common = Covs_grl_all_df_common[grep("W$", names(Covs_grl_all_df_common))]
Covs_grl_all_df_tally_abs_delta_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_delta_evenness", "Sample_Evenness_tally", c(-0.5, 0.5))
Covs_grl_all_df_tally_evenness_sample = Variable_distribution(Covs_grl_all_df_common, "Density", "sample_abs_delta_meth_pct", "Sample_Abs_delta_meth_pct_tally", c(-100, 100))

save(Covs_grl_all_df_tally_abs_delta_group, Covs_grl_all_df_tally_evenness_group, Covs_grl_all_df_tally_abs_delta_sample, Covs_grl_all_df_tally_evenness_sample, file ="Strand_bias.RData")
quit(save = "no")