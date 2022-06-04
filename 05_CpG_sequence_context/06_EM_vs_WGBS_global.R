#####  Script information  #####
################################

# Correlation plots of global CpGs.

#####  Dependencies  #####
##########################

library(ggplot2)
library(bsseq)
library(GGally)

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

# Select random rows.
Beta_meth = Beta_meth[sample(nrow(Beta_meth), 4000000), ]

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
		lower = list(continuous = wrap(ggpair_density_custom, contour_var = "ndensity", adjust = 18, bins = 75), combo = "box_no_facet")) + ggplot2::theme_bw() + ggplot2::theme(axis.text = element_text(size = 6))
dev.off()

quit(save = "no")