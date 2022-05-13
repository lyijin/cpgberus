# `07_parse_bismark_bams/` folder #

Folder contains files obtained from parsing the deduplicated BAMs produced by the `bismark` pipeline, specifically the `deduplicated_bismark` step.

We wanted to mimic one of the figs in ONT's promotional material--the correlation of relative coverage of reads with the GC% content of the underlying genome (100 bp windows). They used Picard's CollectGcBiasMetrics, with no further explanation, so this was how we implemented ours from the files that we obtained.

## Exact commands ##

Generate sorted bams as required by Picard's tools (`bismark` tends to output unsorted BAMs) \
`for a in *.deduplicated.bam; do samtools sort ${a} -@ 10 > ${a:0:9}.sorted.bam; done`

Running the `CollectGcBiasMetrics` tool \
`for a in *.sorted.bam; do PicardCommandLine CollectGcBiasMetrics -I ${a} -O ${a:0:9}.gc_bias_metrics.txt -CHART ${a:0:9}.gc_bias_metrics.pdf -S ${a:0:9}.summary_metrics.txt -R grch38p13_lambda_puc.fa -IS_BISULFITE_SEQUENCED true -ALSO_IGNORE_DUPLICATES true & done; wait`

Outputs the `*.gc_bias_metrics.txt` files used in our eventual plot. The sorted BAM files are too large to be uploaded, hence their exclusion.
