# Comparative analysis of WGBS and EM-Seq data #

This folder contains comparative analysis of WGBS and EM-Seq data, organised into subfolders, as follows:
1. 01_outputs: Coverage analysis
   1. "Combined_tally.png" used for manuscript.
2. 02_outputs: Statistical analysis - Rarefied_CpG_stats_existing.RData containing the statistical analysis results using the DSS R package.
   1. The R object "Rarefied_existing_stats_no_smooth" was used for further downstream analysis.
3. 03_outputs: Motif analysis for statistically significant CpGs used for the heatmap figure.
   1. "Unsmoothed_motifs_existing_1.png" and "Unsmoothed_motifs_existing_2.png" used for manuscript.
4. 04_outputs: Coverage, beta and motif analysis for statistically significant CpGs.
   1. "Unsmoothed_significant_CpGs_heatmap.png" used for manuscript.
   2. "Correlation_plot_3.png" used for manuscript.
   3. "Motif_plots_X.png" used for manuscript, with X signifying a integer.
5. 05_outputs: Plot GC bias curves, quality and bin counts produced by Picard (CollectGcBiasMetrics)
   1. "GC_bias_2.png" used for manuscript.
6. 06_outputs: Pearson correlation analysis for 4 million random CpGs.
   1. "Correlation_plot.png" and "Correlation_plot_2.png" used for the manuscript.

Before running scripts in this folder, you need to generate and access to data in:
1. "cpgberus/04_parse_bismark_covs/Not_rarefied_grch38p13_combined_covs_grl.RData"
2. "cpgberus/04_parse_bismark_covs/Rarefied_grch38p13_combined_covs_grl.RData"
3. "cpgberus/05_CpG_sequence_context/collectgcbiasmetrics/" containing output from Picard (CollectGcBiasMetrics)

## Setting up the enviroment ##

Installation of miniconda.
```
Download miniconda source from: https://docs.conda.io/en/latest/miniconda.html#linux-installers
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
Specify path option: /home/uqdguanz/Programs/miniconda3
Specify conda init option: yes
```

Clone Pipeline_small_ncRNA pipline from github.
```
git clone https://github.com/lyijin/cpgberus.git
```

Create conda enviroment to install dependencies using yaml file.
```
cd cpgberus/05_CpG_sequence_context
conda install -n base conda-forge::mamba
mamba env create -n CpGBerus_env -f CpGBerus_env.yaml
cd ..
```

## Running the code ##

Activate the conda enviroment.
```
conda activate CpGBerus_env
```

Generate required Rdata files.
```
cd 04_parse_bismark_covs
R CMD BATCH tabulations_into_GRanges.R
cd ..
```

Run the data analysis scripts.
```
cd 05_CpG_sequence_context

R CMD BATCH 01_Coverage_analysis.R
sleep 5
R CMD BATCH 02_CpG_context_stats.R
sleep 5
R CMD BATCH 03_CpG_context_graphs.R
sleep 5
R CMD BATCH 04_Strand_bias.R
sleep 5
R CMD BATCH 05_GC_bias_analysis.R
sleep 5
R CMD BATCH 06_EM_vs_WGBS_global.R
```

Deactivate the conda enviroment.
```
conda deactivate
```