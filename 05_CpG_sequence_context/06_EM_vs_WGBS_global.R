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

# # Extract WGBS and EM samples, then convert from wide to long format
# Subset_test_WGBS = Beta_meth[ , c("WR025V1WR", "WR025V9WR", "WR069V1WR", "WR069V9WR")]
# Subset_test_EM = Beta_meth[ , c("WR025V1ER", "WR025V9ER", "WR069V1ER", "WR069V9ER")]

# Subset_test_WGBS = melt(Subset_test_WGBS, value.name = "WGBS_Beta", varnames = "WGBS_CpG_position")
# colnames(Subset_test_WGBS)[is.na(colnames(Subset_test_WGBS))] = "WGBS_Sample"
# Subset_test_WGBS$WGBS_Sample_ID = gsub("WR$", "R", Subset_test_WGBS$WGBS_Sample)

# Subset_test_EM = melt(Subset_test_EM, value.name = "EM_Beta", varnames = "EM_CpG_position")
# colnames(Subset_test_EM)[is.na(colnames(Subset_test_EM))] = "EM_Sample"
# Subset_test_EM$EM_Sample_ID = gsub("ER$", "R", Subset_test_EM$EM_Sample)

# # Check if columns are identical. Stop if false
# if (identical(Subset_test_WGBS$WGBS_CpG_position, Subset_test_EM$EM_CpG_position) != TRUE) {
	# stop("CpG position not identical")
# } 

# if (identical(Subset_test_WGBS$WGBS_Sample_ID, Subset_test_EM$EM_Sample_ID) != TRUE) {
	# stop("Sample ID not identical")
# }

# # Merge and plot
# Subset_test_merged = cbind(Subset_test_WGBS, Subset_test_EM)
# Subset_test_merged_test = Subset_test_merged[Subset_test_merged$WGBS_Sample_ID == "WR025V1R", ]

# png(paste0("/QRISdata/Q4967/Working/Temp/", "Test_2.png"), units="in", width=11.7, height=8.3, res=300)
# test = ggplot(Subset_test_merged, aes(WGBS_Beta, EM_Beta)) + geom_density_2d_filled(contour_var = "ndensity", adjust = 18, bins = 75, show.legend = FALSE) + theme_bw() + facet_wrap( ~ WGBS_Sample_ID)
# test
# dev.off()

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