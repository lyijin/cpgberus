#####  Script information  #####
################################

# Correlation plots of global CpGs.

#####  Dependencies  #####
##########################

library(ggplot2)
library(bsseq)
library(GGally)
library(reshape2)
library(scales)

#####  Load data  #####
#######################

path_to_cpgerus = "/scratch/user/uqdguanz/Projects/Meth/cpgberus"
full_path = paste0(path_to_cpgerus, "/05_CpG_sequence_context/06_outputs/")
dir.create(full_path, recursive = TRUE)

load(file.path(path_to_cpgerus, "05_CpG_sequence_context/02_outputs/Rarefied_CpG_stats_existing.RData"))

# Extract meth and CpG position
Beta_meth = getMeth(Rarefied_existing_bsseq, type = "raw")
rownames_var = mapply(rep, attributes(seqnames(Rarefied_existing_bsseq))$values, attributes(seqnames(Rarefied_existing_bsseq))$lengths)
rownames_var = paste0(unlist(rownames_var), "_", start(Rarefied_existing_bsseq))
rownames(Beta_meth) = rownames_var

# Select random rows and remove rows with NA
Beta_meth = Beta_meth[sample(nrow(Beta_meth), 4000000), ]
Beta_meth = na.omit(Beta_meth)
nrow(Beta_meth)

# Function for using geom_density_2d_filled in ggpairs (GGally)
ggpair_density_custom <- function(data, mapping, ...) {
    
    p <- ggplot(data = data, mapping = mapping) + geom_density_2d_filled(...) + geom_smooth(method = lm, color = "red", size = 0.5)
	p
	
}

png(paste0(full_path, "Correlation_plot.png"), units="in", width=11.7, height=8.3, res=300)
ggpairs(data.frame(Beta_meth), upper = list(continuous = wrap("cor", method = "pearson")), 
		lower = list(continuous = wrap(ggpair_density_custom, contour_var = "ndensity", adjust = 18, bins = 75), combo = "box_no_facet")) + ggplot2::theme_minimal(12) + ggplot2::theme(axis.text = element_text(size = 6))
dev.off()

# Calculate correlation matrix for comparisions
corr_matrix = cor(data.frame(Beta_meth), method = "pearson")

# Reshape dataframe for plotting
df_corr = data.frame(corr_matrix)
df_corr$comparision = rownames(df_corr)
df_corr = melt(df_corr, id.vars = "comparision")

# Find identical patients for plotting colour groups
df_corr$patient_name = sapply(strsplit(df_corr$comparision, "V1|V9"), `[`, 1)
df_corr$patient_name_2 = sapply(strsplit(as.character(df_corr$variable), "V1|V9"), `[`, 1)
df_corr$patient_name_3 = sapply(strsplit(as.character(df_corr$variable), "WR$|ER$"), `[`, 1)

ind = df_corr[ , "patient_name"] == df_corr[ , "patient_name_2"]
df_corr$Sample = "Different"
df_corr[ind, "Sample"] = "Same"

# Factor x-axis
df_corr$patient_name_3 = factor(df_corr$patient_name_3, levels = c("WR025V1", "WR025V9", "WR069V1", "WR069V9"))

# Find identical library for filtering
df_corr$library_name = sapply(strsplit(df_corr$comparision, "V1|V9"), `[`, 2)
df_corr$library_name_2 = sapply(strsplit(as.character(df_corr$variable), "V1|V9"), `[`, 2)

ind = (df_corr[ , "library_name"] == df_corr[ , "library_name_2"])
df_corr$Library = "Different"
df_corr[ind, "Library"] = "Same"

df_corr$library_name_3 = df_corr$library_name_2
df_corr[df_corr$library_name_3 == "ER", "library_name_3"] = "EM-seq"
df_corr[df_corr$library_name_3 == "WR", "library_name_3"] = "WGBS"

# Filter dataframe
df_corr_filt = df_corr[df_corr["Library"] == "Same", ]
df_corr_filt = df_corr_filt[!(df_corr_filt["value"] == 1), ]

# Median text for plot
stat_box_data <- function(y, upper_limit = max(df_corr_filt$value, na.rm = TRUE) * 1.075) {
  return( 
	data.frame(
	  y = 0.95 * upper_limit,
	  label = paste('median =', round(median(y, na.rm = TRUE), 3), '\n')
	)
  )
}
	
# Plot dataframe
p = ggplot(df_corr_filt, aes(x = patient_name_3, y = value))
Final_graph = p + geom_point(aes(color = Sample)) + theme_bw() + theme(legend.position="bottom", plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = pretty_breaks(10), labels = label_number(accuracy = 0.001)) +
				stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
				theme_minimal(15) + xlab("Sample") + theme(axis.title.x=element_blank()) + ylab("Pearson correlation") + facet_wrap( ~ library_name_3)
  
png(paste0(full_path, "Correlation_plot_2.png"), units="in", width=11.7, height=4.15, res=300)
plot(Final_graph)  
dev.off()

quit(save = "no")