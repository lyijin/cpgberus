# `15_ont_minimap2_coverage/` folder #

Folder documents method used to check coverages of ONT reads (following the amplification-free, Cas9-directed protocol) mapping to the targeted 45S rDNA loci in the human genome. The loci was chosen because it had a high GC% (mean of 72%) relative to the human genome (mean of 42%), and it is known that WGBS readouts get a bit more biased at high GC% regions.

Let's go through how things were done in reverse.

## "plot_depth" script and visualisation ##

An R script reads the parsed `samtools mpileup` files for each of the separately barcoded samples, and plots coverage for both the Watson and Crick strands on a stacked bar chart. This was done to illustrate the "U-shaped" nature of the resulting coverages, but strangely this shape was not present in the non-barcoded reads. Perhaps non-barcoded reads did not have double-strand breaks produced by the Cas9 machinery, and the reads come from stochastic fracturing of the DNA ends in the loci of interest? The broken ends thread through the pores and get sequenced. Statistically, there's an equal chance of the DNA getting fractured along the entire loci--hence the flat coverage?

## Parsed mpileup files ##

These files are the `.parsed_mp.tsv.gz` files in this folder. They were parsed from raw tabular outputs from `samtools mpileup`.

## `samtools mpileup` ##

This command was run on the BAM files generated from the mapping of ONT FASTQ reads (those categorised as PASS) against the 45S reference sequence. These files are just under the file size upload limit but honestly, they aren't very interesting.

## Mapping ONT FASTQ against 45S reference ##

The command used to do this was

```shell
for a in ../00_raw_reads/*.fastq.gz; do b=`echo ${a:16:100} | sed 's/fastq.gz/bam/'` && minimap2 -a ~/csiro/stopwatch/raw_data/45s_sequences/homo_sapiens_ont/homo_sapiens.KY962518.last_11kb_first.ont.mmi ${a} | samtools sort -O BAM - > ${b}; done
```

Erm.

In essence, it's `minimap2 -a <reference genome mmi file> <FASTQ file>`. `*.mmi` files are minimap2 indices, check out the instructions at https://github.com/lh3/minimap2 about this.

## What the files look like ##

```shell
$ ls -l
-rw-r--r--  1 lie128 lie128   42M Sep 20 11:32 FAQ88026_pass_5ad1247b.barcode17.bam
-rw-r--r--  1 lie128 lie128   24M Sep 20 13:09 FAQ88026_pass_5ad1247b.barcode17.mpileup.tsv.gz
-rw-r--r--  1 lie128 lie128  362K Sep 20 13:23 FAQ88026_pass_5ad1247b.barcode17.parsed_mp.tsv.gz
-rw-r--r--  1 lie128 lie128   40M Sep 20 11:32 FAQ88026_pass_5ad1247b.barcode18.bam
-rw-r--r--  1 lie128 lie128   23M Sep 20 13:09 FAQ88026_pass_5ad1247b.barcode18.mpileup.tsv.gz
-rw-r--r--  1 lie128 lie128  362K Sep 20 13:23 FAQ88026_pass_5ad1247b.barcode18.parsed_mp.tsv.gz
-rw-r--r--  1 lie128 lie128   31M Sep 20 11:32 FAQ88026_pass_5ad1247b.barcode19.bam
-rw-r--r--  1 lie128 lie128   19M Sep 20 13:09 FAQ88026_pass_5ad1247b.barcode19.mpileup.tsv.gz
-rw-r--r--  1 lie128 lie128  346K Sep 20 13:23 FAQ88026_pass_5ad1247b.barcode19.parsed_mp.tsv.gz
-rw-r--r--  1 lie128 lie128   51M Sep 20 11:33 FAQ88026_pass_5ad1247b.barcode20.bam
-rw-r--r--  1 lie128 lie128   30M Sep 20 13:09 FAQ88026_pass_5ad1247b.barcode20.mpileup.tsv.gz
-rw-r--r--  1 lie128 lie128  374K Sep 20 13:23 FAQ88026_pass_5ad1247b.barcode20.parsed_mp.tsv.gz
-rw-r--r--  1 lie128 lie128  606M Sep 20 11:36 FAQ88026_pass_5ad1247b.unclassified.bam
-rw-r--r--  1 lie128 lie128   55M Sep 20 13:09 FAQ88026_pass_5ad1247b.unclassified.mpileup.tsv.gz
-rw-r--r--  1 lie128 lie128  467K Sep 20 13:24 FAQ88026_pass_5ad1247b.unclassified.parsed_mp.tsv.gz
```
