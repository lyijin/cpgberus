# `04_parse_bismark_covs/` folder #

Folder contains scripts that parses bismark covs, and produces intermediate files that eases downstream analysis/plotting.

First of all, you need to mount the bowen storage drive\
`sudo mount -t nfs -o vers=3 fs1-cbr.nexus.csiro.au:/data/{hb-stopwatch}/work/ /<destination_mount_path>`

There should be a couple of files here when symlinked to the proper folder at\
`/<destination_mount_path>/cpgberus/04_parse_bismark_covs`

```
lie128@blackpuma-ri:~/csiro/stopwatch/cpgberus/04_parse_bismark_covs$ ls -l
total 4.3G
-rw-r--r-- 1 lie128 lie128  746 Aug  5 14:36 README.md
-rw-r--r-- 1 lie128 lie128 436M Aug  5 12:56 WR025V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 384M Aug  5 13:05 WR025V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 394M Aug  5 13:14 WR025V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 362M Aug  5 13:23 WR025V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 482M Aug  5 13:33 WR069V1E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 396M Aug  5 13:41 WR069V1W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 380M Aug  5 13:50 WR069V9E.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 354M Aug  5 13:59 WR069V9W.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r-- 1 lie128 lie128 558M Aug  5 13:10 all.coverages.tsv.gz
-rw-r--r-- 1 lie128 lie128 605M Aug  5 13:26 all.methpcts.tsv.gz
-rwxr-xr-x 1 lie128 lie128  11K Aug  5 12:40 combine_cov_cpgs.py*
-rwxr-xr-x 1 lie128 lie128 3.4K Aug  5 12:38 tabulate_cov_coverages_methpct.py*
```

all.*.gz files contain merged coverages/methylation percentages across the eight input files. Make sure you set meaningful filters before you plot stuff out.

*.combined.tsv.gz merges coverages of C-on-Watson and C-on-Crick, and calculate some stats for it (e.g. unevenness, CpG context...). Check out the Python script for more info.
