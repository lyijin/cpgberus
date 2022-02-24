# `cpgberus`: comparison of whole-genome DNA methylation techniques #

A four-headed beast stands guard at the gates of Methylation. All four heads bark differently--which one's the best at revealing the secrets behind the guarded gates?

![Gustav Dore's illustration of Cerberus in Dante's *Inferno*; public domain](cerberus.jpg)

(Yes, pedantically speaking, Cerberuses tend to have a max of three heads. We, uh, had had a surprise inclusion of a fourth technique later during the project.)

## Common abbreviations ##

WGBS: Whole genome bisulphite sequencing, the current gold standard.

EM-seq: Enzymatic Methyl-seq, uses two enzymes to convert DNA instead of sodium bisulphite.

EPIC: Infinium MethylationEPIC arrays, checks methylation status of ~850,000 positions. Provides "analogue" methylation readouts.

ONT: Oxford Nanopore Technology, allows for long reads and detection of base modifications directly with fancy machine learning.

## Placeholder notes ##

Naming of folders in this project is mainly determined by task.

`data`: put reference data here e.g. genomes, annots, array manifests, ...

`00_series_folders`: these folders contain scripts/pipelines that process raw data into intermediary tables.

`10_series_folders`: these folders contain scripts that parse intermediary tables into plots, which become (supplementary) figures in the manuscript.\
Examples: `10_emseq_vs_wgbs`, `11_methepic_vs_emseq_wgbs`, ...

Standardisation of ggplot2 aesthetics: use `theme_minimal()` if you don't mind?

Standardisation of colour scheme: no opinion atm. Use default ggplot2 colour schemes for now, can be tweaked easily for publication later.

## Analysis notes ##

WR025 and WR069 are both females. The tables still do contain chrY and chrM (but not lambda or pUC19 though) so you might want to filter these chrs out downstream.

Rarefied files are denoted with an "R" in the suffix of the sample name, e.g. WR069V1W --> WR069V1WR. All rarefied files contain 166,282,895 reads at raw FASTQ level (i.e. can differ post-mapping/-dedup because different files had had different mapping/duplication rates).

### Location of `bam` files ###

These files are huge, won't ever be in the repo; they're accessible on the mounted `hb-stopwatch/` folder.

```bash
lie128@blackpuma-ri:/datasets/work/hb-stopwatch/work/key_intermediate_files/human_tests/02_agrf_ramac.emseq_wgbs_methepic.varsamp/03_dedup_grch38p13_lambda_puc$ ls -l WR025* WR069*
-rw-r--r-- 1 lie128 lie128  35G Oct 20 15:15 WR025V1ER_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  73G Sep 12 08:41 WR025V1E_S2_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  32G Oct 20 03:38 WR025V1WR_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  47G Sep 12 17:23 WR025V1W_S21_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  36G Oct 20 05:16 WR025V9ER_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  44G Sep 11 20:18 WR025V9E_S12_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  31G Oct 20 03:41 WR025V9WR_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  34G Sep 12 02:46 WR025V9W_S23_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  33G Oct 20 09:02 WR069V1ER_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125 122G Sep 17 00:22 WR069V1E_S4_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  32G Oct 20 03:49 WR069V1WR_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  53G Sep 12 01:59 WR069V1W_S22_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  35G Oct 20 05:48 WR069V9ER_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  37G Sep 11 18:52 WR069V9E_S14_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 lie128  31G Oct 20 04:19 WR069V9WR_R1_val_1_bismark_bt2_pe.deduplicated.bam
-rw-r--r-- 1 lie128 319125  31G Sep 12 19:55 WR069V9W_S24_R1_val_1_bismark_bt2_pe.deduplicated.bam
```
(there are more files in the folder mentioned above, but for this project, ignore files that are not WR025/WR069)

### Location of `bismark bam2nuc` output ###

These files are tiny, they're uploaded in this repo.

Original files: `cpgberus/02_bismark_wgbs_emseq_data/03_dedup_grch38p13_lambda_puc/*.nucleotide_stats.txt`

Rarified files: `cpgberus/06_rarefaction/03_dedup_grch38p13_lambda_puc/*.nucleotide_stats.txt`
