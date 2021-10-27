# `04_parse_bismark_covs/` folder #

Folder contains scripts that parses bismark covs, and produces intermediate files that eases downstream analysis/plotting.

First of all, you need to mount the bowen storage drive\
`sudo mount -t nfs -o vers=3 fs1-cbr.nexus.csiro.au:/data/{hb-stopwatch}/work/ /<destination_mount_path>`

There should be a couple of files here when symlinked to the proper folder at\
`/<destination_mount_path>/cpgberus/04_parse_bismark_covs`

```bash
lie128@blackpuma-ri:/datasets/work/hb-stopwatch/work/cpgberus/04_parse_bismark_covs$ ls -l
total 23G
-rw-r--r-- 1 lie128 lie128  2.0K Aug  5 14:43 README.md
-rw-r--r-- 1 lie128 lie128  436M Sep 24 12:53 WR025V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  380M Oct 25 10:30 WR025V1ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  389M Sep 24 13:02 WR025V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  364M Oct 25 10:39 WR025V1WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  394M Sep 24 13:11 WR025V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  380M Oct 25 10:48 WR025V9ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  366M Sep 24 13:20 WR025V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  361M Oct 25 10:56 WR025V9WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  482M Sep 24 13:29 WR069V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  375M Oct 25 11:05 WR069V1ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  400M Sep 24 13:38 WR069V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  367M Oct 25 11:14 WR069V1WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  380M Sep 24 13:47 WR069V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  378M Oct 25 11:23 WR069V9ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  359M Sep 24 13:56 WR069V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  359M Oct 25 11:31 WR069V9WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128  560M Sep 24 13:09 all.coverages.tsv.gz
-rw-r--r-- 1 lie128 lie128  626M Sep 24 13:12 all.methpcts.tsv.gz
-rwxr-xr-x 1 lie128 lie128   11K Aug 11 15:37 combine_cov_cpgs.py*
-rw-r--r-- 1 lie128 lie128   58K Aug 16 10:42 grch38p13.nnncgnnn_universe.tsv
-rw-r--r-- 1 229072 1074248 653M Aug 16 11:59 grch38p13_combined_covs_epic-overlap_grl.RData
-rw-r--r-- 1 229072 1074248  14G Aug 16 11:44 grch38p13_combined_covs_grl.RData
-rw-r--r-- 1 229072 1074248 2.6M Aug 19 12:45 motif_matrices.RData
-rw-r--r-- 1 lie128 lie128  504M Oct 25 10:43 rarefied.coverages.tsv.gz
-rw-r--r-- 1 lie128 lie128  530M Oct 25 10:46 rarefied.methpcts.tsv.gz
-rwxr-xr-x 1 lie128 lie128  4.3K Aug 11 15:38 tabulate_cov_coverages_methpct.py*
-rwxr-xr-x 1 lie128 lie128  2.0K Aug 16 10:40 tabulate_nnncgnnn_universe.py*
```

all.*.gz files contain merged coverages/methylation percentages across the eight original input files. Make sure you set meaningful filters before you plot stuff out.

rarefied.*.gz files contain merged coverages/methylation percentages across the eight rarefied input files (all of them have 166,282,895 at raw FASTQ level). Make sure you set meaningful filters before you plot stuff out.

*.combined.tsv.gz merges coverages of C-on-Watson and C-on-Crick, and calculate some stats for it (e.g. unevenness, CpG context...). Check out the Python script for more info. Note that sample names with "R" suffix (e.g. "WR025V1ER" are from the rarefied inputs).
