#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --job-name=CpGBerus
#SBATCH --time=6:00:00
#SBATCH --partition=general
#SBATCH --account=a_salomon
#SBATCH -o slurm.output
#SBATCH -e slurm.error

conda activate CpGBerus_env

#R CMD BATCH 01_Coverage_analysis.R

#sleep 5

#R CMD BATCH 02_CpG_context_stats.R

#sleep 5

R CMD BATCH 03_CpG_context_graphs.R

sleep 5

R CMD BATCH 04_Strand_bias.R

sleep 5

R CMD BATCH 05_GC_bias_analysis.R

sleep 5

R CMD BATCH 06_EM_vs_WGBS_global.R

conda deactivate
