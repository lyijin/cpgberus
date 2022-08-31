#####  Script information  #####
################################

# Plot GC_bias curves, quality and bin counts produced by Picard (CollectGcBiasMetrics).

#####  Dependencies  #####
##########################

library(ggplot2)
library(reshape2)
library(scales)

#####  Load data  #####
#######################

path_to_cpgerus = "/scratch/user/uqdguanz/Projects/Meth/cpgberus"
full_path = paste0(path_to_cpgerus, "/05_CpG_sequence_context/05_outputs/")
dir.create(full_path, recursive = TRUE)

# Combine dataframes

Files_to_analyse = list.files(paste0(path_to_cpgerus, "/05_CpG_sequence_context/collectgcbiasmetrics"), "gc_bias_metrics.txt", full.names = TRUE)
Merged_dataframe = data.frame()

for (file_name in Files_to_analyse) {
  Input_file = read.csv(file_name, skip = 6, sep = "\t", row.names = NULL)
  Input_file$Sample_name = strsplit(basename(file_name), split = "\\.")[[1]][1]
  Merged_dataframe = rbind(Merged_dataframe, Input_file)
}

# Extract and assign library type
Merged_dataframe_filt = Merged_dataframe[Merged_dataframe$READS_USED == "ALL", ]
Merged_dataframe_filt$Library_type = "Temp"
Merged_dataframe_filt[grep("WR$", Merged_dataframe_filt$Sample_name), "Library_type"] = "WGBS"
Merged_dataframe_filt[grep("ER$", Merged_dataframe_filt$Sample_name), "Library_type"] = "EM-seq"

# Convert to long and rename groups
Merged_dataframe_filt_long = melt(Merged_dataframe_filt, id.vars=c("Sample_name", "Library_type", "GC"), measure.vars = c("NORMALIZED_COVERAGE", "MEAN_BASE_QUALITY", "READ_STARTS"))
Merged_dataframe_filt_long$variable = as.character(Merged_dataframe_filt_long$variable)
Merged_dataframe_filt_long[Merged_dataframe_filt_long$variable == "NORMALIZED_COVERAGE", "variable"] = "Normalised Coverage"
Merged_dataframe_filt_long[Merged_dataframe_filt_long$variable == "MEAN_BASE_QUALITY", "variable"] = "Average base quality"
Merged_dataframe_filt_long[Merged_dataframe_filt_long$variable == "READ_STARTS", "variable"] = "Read counts"

Merged_dataframe_filt_long$variable = factor(Merged_dataframe_filt_long$variable, levels = c("Normalised Coverage", "Average base quality", "Read counts"))
#Merged_dataframe_filt_long$GC = as.factor(Merged_dataframe_filt_long$GC)

# Plot graph

p = ggplot(Merged_dataframe_filt_long, aes(x = GC, y = value, color = Library_type, group = Library_type)) + geom_line(aes(color = Library_type, group = Sample_name), alpha = 0.5) + 
    stat_summary(geom = "line", fun = mean, size = 1.5, alpha = 0.7) + scale_y_continuous(breaks = breaks_pretty(10), labels = comma) + scale_x_continuous(breaks = breaks_pretty(10), labels = comma) + 
    labs(x = "GC %", y = "Value") + facet_wrap( ~ variable,  scales="free_y") + theme_minimal(13) + scale_color_manual(values=c("#1b9e77", "#7570b3")) + labs(color = "Library type")

png(paste0(full_path, "GC_bias.png"), units="in", width=11.7, height=4.15, res=300)
print(p)
dev.off()


Merged_dataframe_filt_long_2 = Merged_dataframe_filt_long
Merged_dataframe_filt_long_2 = Merged_dataframe_filt_long_2[!(Merged_dataframe_filt_long_2$variable == "Average base quality"), ]

p = ggplot(Merged_dataframe_filt_long_2, aes(x = GC, y = value, color = Library_type, group = Library_type)) + geom_line(aes(color = Library_type, group = Sample_name), alpha = 0.5) + 
    stat_summary(geom = "line", fun = mean, size = 1.5, alpha = 0.7) + scale_y_continuous(breaks = breaks_pretty(10), labels = comma) + scale_x_continuous(breaks = breaks_pretty(10), labels = comma) + 
    labs(x = "GC %", y = "Value") + facet_wrap( ~ variable,  scales="free_y") + theme_minimal(13) + scale_color_manual(values=c("#1b9e77", "#7570b3")) + labs(color = "Library type")

png(paste0(full_path, "GC_bias_2.png"), units="in", width=11.7, height=4.15, res=300)
print(p)
dev.off()

quit(save = "no")