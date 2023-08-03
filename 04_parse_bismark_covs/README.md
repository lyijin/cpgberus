# `04_parse_bismark_covs/` folder #

Folder contains scripts that parses bismark covs, and produces intermediate files that eases downstream analysis/plotting.
CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/04_parse_bismark_covs`)

```bash
# file sizes are in bytes, total size across all files here are ~6.3 GB
$ ls -l
-rw-r--r--   2415850006 Nov 10  2021 Not_rarefied_grch38p13_combined_covs_grl.RData
-rw-r--r--   2209069326 Nov 10  2021 Rarefied_grch38p13_combined_covs_grl.RData
-rw-r--r--         3586 Nov 10  2021 Seqinfo_hg38.RData
-rw-r--r--    398044480 Oct 24  2021 WR025V1ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    380944825 Oct 24  2021 WR025V1WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    397910209 Oct 24  2021 WR025V9ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    377673836 Oct 24  2021 WR025V9WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    392705205 Oct 25  2021 WR069V1ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    383888635 Oct 25  2021 WR069V1WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    395509560 Oct 25  2021 WR069V9ER.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--    375838572 Oct 25  2021 WR069V9WR.grch38p13_lambda_puc.combined.tsv.gz
-rw-r--r--      2623251 Aug 19  2021 motif_matrices.RData
-rw-r--r--    528015724 Oct 24  2021 rarefied.coverages.tsv.gz
-rw-r--r--    554923448 Oct 24  2021 rarefied.methpcts.tsv.gz
```

`rarefied.*.gz` files contain merged coverages/methylation percentages across the eight rarefied input files (all of them have 166,282,895 at raw FASTQ level).

`*.combined.tsv.gz` merges coverages of C-on-Watson and C-on-Crick (i.e., tabulates per-CpG stats from per-position input files), and calculates some additional stats for it (e.g. unevenness, CpG context...). Check out the Python script for more info.
