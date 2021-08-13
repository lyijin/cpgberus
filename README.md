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

## Analysis notes ##

WR025 and WR069 are both females. The tables still do contain chrY and chrM (but not lambda or pUC19 though) so you might want to filter these chrs out downstream.
