# `cpgberus`: comparison of whole-genome DNA methylation techniques #

A four-headed beast stands guard at the gates of Methylation. All four heads bark differently--which one's the best at revealing the secrets behind the guarded gates?

## Common abbreviations ##

WGBS: Whole genome bisulphite sequencing, the gold standard.

EM-seq: Enzymatic Methyl-seq, uses two enzymes to convert DNA instead of sodium bisulphite.

MethylationEPIC: Infinium MethylationEPIC arrays, checks methylation status of 850k positions. Provides "analogue" output.

ONT: Oxford Nanopore Technology, allows for long reads and detection of base modifications directly with fancy machine learning.

## Placeholder notes ##

Naming of folders in this project is mainly determined by task.

`data`: put reference data here e.g. genomes, annots, array manifests, ...

`00_series_folders`: these folders contain scripts/pipelines that process raw data into intermediary tables.

`10_series_folders`: these folders contain scripts that parse intermediary tables into plots, which become (supplementary) figures in the manuscript.\
Examples: `10_emseq_vs_wgbs`, `11_methepic_vs_emseq_wgbs`, ...

Standardisation of ggplot2 aesthetics: use `theme_minimal()` if you don't mind?

Standardisation of colour scheme: no opinion atm. Use default ggplot2 colour schemes for now, can be tweaked easily for publication later.

## TODO ##

WGBS datasets were treated as EM-seq datasets during processing (which is not correct). I am currently reprocessing the data, ETA 1 week. This means that all the intermediate files you're seeing now will change as well.

However, the filenames and column structure of those files are correct though, so go ahead and write R code to parse stuff and create plots. I'll send an email to confirm everything is correct next week, so you can simply re-run your scripts to generate new correct plots.
