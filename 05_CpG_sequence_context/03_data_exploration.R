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

#####  Analysis using density  #####
####################################

# Total coverage calculation, convert to dataframe
Covs_grl_all_df = lapply(covs_grl, function(x) data.frame(chr = seqnames(x), pos = start(x), N = (x$meth_cov + x$unmeth_cov), X = x$meth_cov))

# Subset first million lines for testing
#Covs_grl_all_df_subset = lapply(Covs_grl_all_df,"[", 1:1500000, )

# Subset N (coverage) <= 100
Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) x[x[, "N"] <= 100, ])

# Calculate density
Covs_grl_all_df_density = lapply(Covs_grl_all_df_subset, function(x) density(x[, "N"], adjust = 5))

# Create dataframe x and y density coord + name, then collapse.
n = names(Covs_grl_all_df_density)
Covs_grl_all_df_density_xy = lapply(n, function(n) data.frame(Sample_name = n, x_values = Covs_grl_all_df_density[[n]][["x"]], y_values = Covs_grl_all_df_density[[n]][["y"]]))
Covs_grl_all_df_density_xy_collapse = do.call(rbind, Covs_grl_all_df_density_xy)

# Add sample name without suffix and sequencing type column
Covs_grl_all_df_density_xy_collapse$Sample_name_2 = substr(Covs_grl_all_df_density_xy_collapse$Sample_name, 1, nchar(Covs_grl_all_df_density_xy_collapse$Sample_name) - 1)

Covs_grl_all_df_density_xy_collapse$Seq_type = Covs_grl_all_df_density_xy_collapse$Sample_name
Covs_grl_all_df_density_xy_collapse$Seq_type[grep("E$", Covs_grl_all_df_density_xy_collapse$Sample_name)] <- "EM-Seq"
Covs_grl_all_df_density_xy_collapse$Seq_type[grep("W$", Covs_grl_all_df_density_xy_collapse$Sample_name)] <- "WGBS"

# Plot
pdf("Analysis_density.pdf", width = 11.69, height = 8.3)
ggplot(Covs_grl_all_df_density_xy_collapse, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw()
ggplot(Covs_grl_all_df_density_xy_collapse, aes(x=x_values, y=y_values, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_density_xy_collapse, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_density_xy_collapse, aes(x=x_values, y=y_values, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
dev.off()

#####  Analysis using tally  #####
##################################

# Tally data
n = names(Covs_grl_all_df_subset)
Covs_grl_all_df_subset_table = lapply(n, function(n) data.frame(table(Covs_grl_all_df[[n]][["N"]]), Sample_name = n))
Covs_grl_all_df_subset_table = do.call(rbind, Covs_grl_all_df_subset_table)

# Add sample name without suffix and sequencing type column
Covs_grl_all_df_subset_table$Sample_name_2 = substr(Covs_grl_all_df_subset_table$Sample_name, 1, nchar(Covs_grl_all_df_subset_table$Sample_name) - 1)

Covs_grl_all_df_subset_table$Seq_type = Covs_grl_all_df_subset_table$Sample_name
Covs_grl_all_df_subset_table$Seq_type[grep("E$", Covs_grl_all_df_subset_table$Sample_name)] <- "EM-Seq"
Covs_grl_all_df_subset_table$Seq_type[grep("W$", Covs_grl_all_df_subset_table$Sample_name)] <- "WGBS"

# Plot
pdf("Analysis_tally.pdf", width = 11.69, height = 8.3)
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
dev.off()

#####  Subsample  #####
#######################

#Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) cbind(x, data.frame(Subsample_coverage = rbinom(nrow(x), x[, "N"], 0.5))))

Covs_grl_all_df_subset = lapply(Covs_grl_all_df, function(x) cbind(x, data.frame(Subsample_coverage = x[, "N"])))
Covs_grl_all_df_subset[["WR069V1E"]]["Subsample_coverage"] <- rbinom(nrow(Covs_grl_all_df_subset[["WR069V1E"]]), Covs_grl_all_df_subset[["WR069V1E"]][, "N"], 0.25)
Covs_grl_all_df_subset[["WR025V1E"]]["Subsample_coverage"] <- rbinom(nrow(Covs_grl_all_df_subset[["WR025V1E"]]), Covs_grl_all_df_subset[["WR025V1E"]][, "N"], 0.6)

# Tally data
n = names(Covs_grl_all_df_subset)
Covs_grl_all_df_subset_table = lapply(n, function(n) data.frame(table(Covs_grl_all_df_subset[[n]][["Subsample_coverage"]]), Sample_name = n))
Covs_grl_all_df_subset_table = do.call(rbind, Covs_grl_all_df_subset_table)

# Add sample name without suffix and sequencing type column
Covs_grl_all_df_subset_table$Sample_name_2 = substr(Covs_grl_all_df_subset_table$Sample_name, 1, nchar(Covs_grl_all_df_subset_table$Sample_name) - 1)

Covs_grl_all_df_subset_table$Seq_type = Covs_grl_all_df_subset_table$Sample_name
Covs_grl_all_df_subset_table$Seq_type[grep("E$", Covs_grl_all_df_subset_table$Sample_name)] <- "EM-Seq"
Covs_grl_all_df_subset_table$Seq_type[grep("W$", Covs_grl_all_df_subset_table$Sample_name)] <- "WGBS"

# Plot
pdf("Analysis_tally_subset.pdf", width = 11.69, height = 8.3)
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Seq_type)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Sample_name)) + geom_line() + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
ggplot(Covs_grl_all_df_subset_table, aes(x=as.integer(Var1), y=Freq, group = Sample_name, color = Sample_name_2)) + geom_line(aes(linetype = Seq_type)) + scale_x_continuous(breaks = pretty_breaks(10), limits = c(-5, 70)) + theme_bw()
dev.off()