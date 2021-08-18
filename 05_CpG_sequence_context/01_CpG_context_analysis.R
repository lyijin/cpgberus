#####  Dependencies  #####

library(GenomicRanges)
library(motifStack)

#####  Load data  #####

load("/scratch1/gua020/CpGberus_data/grch38p13_combined_covs_grl.RData")

#####  Analysis  #####

# Extract first 100000 rows from each Granges for testing, seperate EM-Seq (E) and WGBS (W)
# Subset_covs_grl = endoapply(covs_grl,"[",1:100000)
# Subset_covs_grl_E = covs_grl[grep("E$", names(covs_grl))]
# Subset_covs_grl_W = covs_grl[grep("W$", names(covs_grl))]

# Function to filter granges object for meth percentage range, and coverage cutoff's for meth / unmeth
Filter_meth_pct <- function(Granges_object, meth_pct_keep_range, meth_unmeth_cutoff) {
    temp_grange = subset(Granges_object, meth_pct >= min(meth_pct_keep_range) & meth_pct <= max(meth_pct_keep_range))
    temp_grange = subset(temp_grange, meth_cov >= meth_unmeth_cutoff[1] & unmeth_cov >= meth_unmeth_cutoff[2])
    return(temp_grange)
}

# Function to to filter granges list object for common CpGs
Filter_common_CpGs <- function(Granges_list_object) {
    Common_intersect_regions = Reduce(intersect, Granges_list_object)
    Subset_covs_grl_common = endoapply(Granges_list_object, subsetByOverlaps, Common_intersect_regions)
    return(Subset_covs_grl_common)
}

# Filter CpGs for meth and coverage
Subset_covs_grl_meth = endoapply(covs_grl, Filter_meth_pct, c(100, 100), c(5, 0))
Subset_covs_grl_unmeth = endoapply(covs_grl, Filter_meth_pct, c(0, 0), c(0, 5))

Subset_covs_grl_meth_common = Filter_common_CpGs(Subset_covs_grl_meth)
Subset_covs_grl_unmeth_common = Filter_common_CpGs(Subset_covs_grl_unmeth)

# Function to tabulate nucleotide species for each postion of motif
# Additional option to mutliply motif by total coverage
Motif_frequency_table <- function(Granges_object, number_of_granges, motif_coverage_multiplier) {

    counter <<- counter + 1
    print(paste0("Analysing Grange object: ", as.character(counter), " / ", number_of_granges))
    
    # Create empty nucleotide frequency table
    cpg_context = unique(Granges_object$cpg_context_nnncgnnn)
    largest_motif = max(nchar(cpg_context))

    column_names = unique(unlist(strsplit(paste(cpg_context,collapse=""), "")))
    nucleotide_freq = data.frame(matrix(NA, nrow = 0, ncol = length(column_names)))
    colnames(nucleotide_freq) <- column_names
    
    # Multiply motif cpg_context by total coverage if TRUE
    if (motif_coverage_multiplier == TRUE) {
        motif_multiplier = Granges_object$meth_cov + Granges_object$unmeth_cov
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = motif_multiplier)
    } else {
        cpg_context = data.frame(Motif = Granges_object$cpg_context_nnncgnnn, Motif_multiplier = 1)
    }

    # Count nucleotide population for first position of all motifs, loop + 1 until end of motif.
    motif_list = cpg_context$Motif
    
    for (i in seq(largest_motif)) {    
        print(paste0("Motif position: ", i))
        
        # Substitue motifs with n nucleotide position, then sum
        motif_n_elements = substr(motif_list, i, i)
        cpg_context$Motif <- motif_n_elements
        motif_table_sum = aggregate(Motif_multiplier ~ Motif, cpg_context, FUN=sum)
        
        # Clean up motif_table_sum
        nucleotide_freq_temp = as.data.frame(t(motif_table_sum))
        colnames(nucleotide_freq_temp) <- as.character(nucleotide_freq_temp["Motif", ])
        nucleotide_freq_temp = nucleotide_freq_temp[!(row.names(nucleotide_freq_temp) %in% "Motif"), , drop = FALSE]
        
        # Fill in missing nucleotides in nucleotide_freq_temp with 0
        if (identical(sort(column_names), sort(colnames(nucleotide_freq_temp))) == FALSE) {
            columns_to_add = column_names[(!(column_names %in% colnames(nucleotide_freq_temp)))]
            temp_df = data.frame(matrix(0, nrow = 1, ncol = length(columns_to_add)))
            colnames(temp_df) <- columns_to_add
            nucleotide_freq_temp = cbind(nucleotide_freq_temp, temp_df)
        }
        nucleotide_freq_temp = lapply(nucleotide_freq_temp, as.integer)
        nucleotide_freq = rbind(nucleotide_freq, nucleotide_freq_temp)
    }

    return(data.frame(t(nucleotide_freq)))
}

# Function to convert list of dataframes to motif class pcm
Motif_stack_conversion <- function(Motif_list_df, remove_N) {

    #Remove "N" nucleotide if true
    if (remove_N == TRUE) {
        Motif_list_df = lapply(Motif_list_df, function(x) x[rownames(x)!="N", ])
    }
    
    # Convert dataframe to motif class pcm
    for (data_name in names(Motif_list_df)) {
        Motif_list_df[[data_name]] <- new("pcm", mat=as.matrix(Motif_list_df[[data_name]]), name=data_name)
    }
    
    return(Motif_list_df)
}

# Funtion to plot by patient ID
Motif_plot <- function(Motif_list_df, file_name, sample_names, remove_N) {

    pdf(paste0(file_name, ".pdf"), width = 11.69, height = 8.3)
    for (samp_name in sample_names) {
        Motif_list_df_temp = Motif_list_df[names(Motif_list_df) == samp_name]
        plot(Motif_stack_conversion(Motif_list_df_temp, TRUE)[[samp_name]], ic.scale=FALSE, ylab="probability")
    }
    dev.off()
}

# Motif tables and plots where motif multiplied by coverage
counter = 0
Plot_data_meth <- lapply(Subset_covs_grl_meth, Motif_frequency_table, length(Subset_covs_grl_meth), TRUE)
Motif_plot(Plot_data_meth, "Meth_multiply", names(Plot_data_meth), TRUE)

counter = 0
Plot_data_unmeth <- lapply(Subset_covs_grl_unmeth, Motif_frequency_table, length(Subset_covs_grl_unmeth), TRUE)
Motif_plot(Plot_data_unmeth, "Unmeth_multiply", names(Plot_data_unmeth), TRUE)

counter = 0
Plot_data_meth_common <- lapply(Subset_covs_grl_meth_common, Motif_frequency_table, length(Subset_covs_grl_meth_common), TRUE)
Motif_plot(Plot_data_meth_common, "Meth_common_multiply", names(Plot_data_meth_common), TRUE)

counter = 0
Plot_data_unmeth_common <- lapply(Subset_covs_grl_unmeth_common, Motif_frequency_table, length(Subset_covs_grl_unmeth_common), TRUE)
Motif_plot(Plot_data_unmeth_common, "Unmeth_common_multiply", names(Plot_data_unmeth_common), TRUE)


# Motif tables and plots where motif not multiplied
counter = 0
Plot_data_meth_no_multiply <- lapply(Subset_covs_grl_meth, Motif_frequency_table, length(Subset_covs_grl_meth), FALSE)
Motif_plot(Plot_data_meth_no_multiply, "Meth", names(Plot_data_meth_no_multiply), TRUE)

counter = 0
Plot_data_unmeth_no_multiply <- lapply(Subset_covs_grl_unmeth, Motif_frequency_table, length(Subset_covs_grl_unmeth), FALSE)
Motif_plot(Plot_data_unmeth_no_multiply, "Unmeth", names(Plot_data_unmeth_no_multiply), TRUE)

counter = 0
Plot_data_meth_common_no_multiply <- lapply(Subset_covs_grl_meth_common, Motif_frequency_table, length(Subset_covs_grl_meth_common), FALSE)
Motif_plot(Plot_data_meth_common_no_multiply, "Meth_common", names(Plot_data_meth_common_no_multiply), TRUE)

counter = 0
Plot_data_unmeth_common_no_multiply <- lapply(Subset_covs_grl_unmeth_common, Motif_frequency_table, length(Subset_covs_grl_unmeth_common), FALSE)
Motif_plot(Plot_data_unmeth_common_no_multiply, "Unmeth_common", names(Plot_data_unmeth_common_no_multiply), TRUE)

rm(covs_grl)