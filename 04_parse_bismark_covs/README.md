# `04_parse_bismark_covs/` folder #

Folder contains scripts that parses bismark covs, and produces intermediate files that eases downstream analysis/plotting.

First of all, you need to mount the bowen storage drive\
`sudo mount -t nfs -o vers=3 fs1-cbr.nexus.csiro.au:/data/{hb-stopwatch}/work/ /<destination_mount_path>`

There should be a couple of files here when symlinked to the proper folder at\
`/<destination_mount_path>/cpgberus/04_parse_bismark_covs`

```
lie128@blackpuma-ri:~/csiro/stopwatch/cpgberus/04_parse_bismark_covs$ ls -l
total 4.3G
-rw-r--r-- 1 lie128 lie128 2.0K Aug  5 14:43 README.md
-rw-r--r-- 1 lie128 lie128 436M Aug 12 19:24 WR025V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 389M Aug 12 19:33 WR025V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 394M Aug 12 19:42 WR025V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 366M Aug 12 19:51 WR025V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 482M Aug 12 20:00 WR069V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 400M Aug 12 20:08 WR069V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 380M Aug 12 20:17 WR069V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 359M Aug 12 20:26 WR069V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 560M Aug 12 19:38 all.coverages.tsv.gz
-rw-r--r-- 1 lie128 lie128 626M Aug 12 19:41 all.methpcts.tsv.gz
-rwxr-xr-x 1 lie128 lie128  11K Aug 11 15:37 combine_cov_cpgs.py*
-rwxr-xr-x 1 lie128 lie128 4.3K Aug 11 15:38 tabulate_cov_coverages_methpct.py*
```

all.*.gz files contain merged coverages/methylation percentages across the eight input files. Make sure you set meaningful filters before you plot stuff out.

*.combined.tsv.gz merges coverages of C-on-Watson and C-on-Crick, and calculate some stats for it (e.g. unevenness, CpG context...). Check out the Python script for more info.
