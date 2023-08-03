# `cpgberus`: comparison of whole-genome DNA methylation techniques #

A four-headed beast stands guard at the gates of Methylation. All four heads bark differently--which one's the best at revealing the secrets behind the guarded gates?

![random cute dog picture sourced from the internet](cerberus.jpg)

(Yes, Cerberuses tend to have three heads. We, uh, had had a surprise inclusion of a fourth technique later during the project.)

## The four heads ##

**WGBS:** Whole genome bisulphite sequencing, the current gold standard.

**EM-seq:** Enzymatic Methyl-seq, uses two enzymes to convert DNA instead of sodium bisulphite.

**EPIC:** Infinium MethylationEPIC arrays, checks methylation status of > 850,000 cytosines. Provides more "analogue" methylation readouts. We used v1 arrays, v2 came out later after we finished our work.

**ONT:** Oxford Nanopore Technology, allows for long reads and detection of base modifications directly with ~~black magic~~ fancy machine learning.

## Folder notes ##

Naming of folders in this project is mainly determined by task. Sequential naming is intentional, and the later folders might depend on code from earlier folders--never the other way around.

`data`: contain reference data used by scripts in other folders e.g. genomes, annots, array manifests, etc. These files are too large to be uploaded to GitHub, but available via CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/data`).

`01_txdb`: prepare human annotations from GENCODE v38.

`02_process_methepic_data`: prepare MethylationEPIC data. Some larger files are housed on CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/02_process_methepic_data`).

`03_bismark_rarefied_data`: prepare WGBS and EM-seq data. Some larger files are at CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/03_bismark_rarefied_data`).

`04_parse_bismark_covs`: compile the per-sample methylation coverages and betas from `03` into giant tables used in later scripts. Some larger files can be found on CSIRO DAP (https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/04_parse_bismark_covs`).

`05_CpG_sequence_context`: perform EM-seq vs. WGBS analyses.

`06_process_ont_data`: prepare ONT data.

(07-12 is intentionally missing, were overflow folder numbers in case we needed to do more processing/analyses before the next batch of folders.)

`13_check_mcgw_emseq_wgbs`: do EM-seq reads have biased conversion around MCGW?

`14_methepic_vs_emseq_wgbs`: perform EPIC vs. EM-seq vs. WGBS analyses.

`15_ont_minimap2_coverage`: check coverage for ONT reads, estimate enrichment of 45S rDNA relative to non-targeted regions.

`16_loci_specific_three_way`: perform ONT vs. EM-seq vs. WGBS analyses centred around the 45S rDNA loci.

Detailed READMEs can be found in every subfolder!

## Note on viewing HTML output from R scripts ##

Annoyingly, HTML files are treated as plain-text when viewed in GitHub repos. We do use `knitr` to generate HTML output for our R scripts, as we understand that re-running our scripts to view script outputs can be difficult to set up, and might require large files. To view these HTML files in their entire glory, you can either:

1. `git clone` the entire repo, view the folder contents on your desktop operating system, and double-click to view the downloaded HTML files in your browser. Best if you want to view many/all of the hosted HTML files in this repo?

2. Append "https://htmlpreview.github.io/?" in front of the HTML file. e.g., "https://github.com/lyijin/cpgberus/blob/master/14_methepic_vs_emseq_wgbs/methepic_emseq_wgbs_per-pos_beta.html" becomes "https://htmlpreview.github.io/?https://github.com/lyijin/cpgberus/blob/master/14_methepic_vs_emseq_wgbs/methepic_emseq_wgbs_per-pos_beta.html". Better to view the one-off HTML file that piques your interest. To understand how this hack works, check out their README at https://github.com/htmlpreview/htmlpreview.github.com.
