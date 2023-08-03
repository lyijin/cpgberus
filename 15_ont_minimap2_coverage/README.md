# `15_ont_minimap2_coverage/` folder #

This folder documents our method to check the coverages of ONT reads (following the amplification-free, Cas9-directed protocol) that mapped to the targeted 45S rDNA loci in the human genome. The loci was chosen because it had a high GC% (mean of 72%) relative to the human genome (mean of 42%), and it is known that WGBS readouts are biased at high GC% regions.

The process of producing the `*.parsed_mp.tsv` intermediate tables can be a bit convoluted. Due to upload size constraints, we only make `*.parsed_mp.tsv` files available in a compressed format (available in this folder).

## Mapping ONT FASTQ against 45S reference ##

The command used to do this was

```shell
for a in ../00_raw_reads/*.fastq.gz; do b=`echo ${a:16:100} | sed 's/fastq.gz/bam/'` && minimap2 -a ~/csiro/stopwatch/raw_data/45s_sequences/homo_sapiens_ont/homo_sapiens.KY962518.last_11kb_first.ont.mmi ${a} | samtools sort -O BAM - > ${b}; done
```

Erm.

In essence, it's `minimap2 -a <reference genome mmi file> <FASTQ file>`. `*.mmi` files are minimap2 indices, check out the instructions at https://github.com/lh3/minimap2 about this.

## `samtools mpileup` ##

This command was run on the BAM files generated from the mapping of ONT FASTQ reads (those categorised as PASS) against the 45S reference sequence. These files were generated with the generic command `samtools mpileup -A -B -d 0 -Q 0 -x ${input_bam} | gzip > ${output_file}`.

In our pipeline, these files have a `*.mpileup.tsv.gz` suffix.

## Parsed mpileup files ##

The compressed mpileup tables from the previous step were then parsed using `parse_mpileup_tsvs.py`, resulting in the `.parsed_mp.tsv.gz` that we make available in this folder.

## `plot_depth.R` script ##

An R script reads the `*.parsed_mp.tsv` files for each of the separately barcoded samples, and plots coverage for both the Watson and Crick strands on a stacked bar chart.

This was done to illustrate the "U-shaped" nature of the resulting coverages, but strangely this shape was not present in the non-barcoded reads. Perhaps non-barcoded reads did not have double-strand breaks produced by the Cas9 machinery, and the reads come from stochastic fracturing of the DNA ends in the loci of interest? The broken ends thread through the pores and get sequenced. Statistically, there's an equal chance of the DNA getting fractured along the entire loci--hence the flat coverage?
