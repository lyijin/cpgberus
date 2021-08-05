# `data/` folder #

First of all, you need to mount the bowen storage drive\
`sudo mount -t nfs -o vers=3 fs1-cbr.nexus.csiro.au:/data/{hb-stopwatch}/work/ /<destination_mount_path>`

There should be a couple of files here when symlinked to the proper folder at\
`/<destination_mount_path>/cpgberus/data`

```
lie128@blackpuma-ri:~/csiro/stopwatch/cpgberus/data$ ls -l
total 1016M
-rw-r--r-- 1 lie128 lie128  305 Aug  5 14:34 README.md
-rw-r--r-- 1 lie128 lie128  45M Aug  3 13:25 gencode.v38.annotation.gtf.gz
-rw-r--r-- 1 lie128 lie128 123M Aug  4 15:48 gencode.v38.annotation.sqlite
-rw-r--r-- 1 lie128 lie128  11M Aug  4 15:49 gencode.v38.annotation_2021-08-04.RData
-rw-r--r-- 1 lie128 lie128 839M Aug  3 13:39 grch38p13_lambda_puc.fa.gz
```