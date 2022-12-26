# `16_loci_specific_three_way/` folder #

## Data used for plots ##

The ONT bed files containing per-position methylation levels in this folder were produced by `megalodon`. See `06_process_ont_data/` for a fuller description--methylation calling was a finicky process, best just use the bed files.

EM-seq and WGBS cov files with per-position methylation levels were called using the `bismark` pipeline, detailed in `03_bismark_rarefied_data/`.

## Idea behind R code ##

Folder contains two mega R script that carries out three-way comparisons between EM-seq, WGBS and ONT Cas9 readouts from the 45S rDNA loci in the human genome. One script looks at methylation levels (betas), another looks at coverages across the loci. The loci was chosen because it had a high GC% (mean of 72%) relative to the human genome (mean of 42%), and it is known that WGBS readouts get a bit more biased at high GC% regions.

Probably easier to view the HTML file in this folder to see the R script in action + the plots that were made. Not all plots were included in the manuscript.

## Note about the `rerio_beds/` folder ##

Over the course of our work, as ONT-associated bioinformatics is undergoing rapid development, the models used in calculating the confidence of whether a C in the CpG context is methylated/unmethylated has also changed rapidly. We used to use models from `rerio` (the only available models ~late 2021), and we observed that the ONT methylation calls were closer to WGBS than EM-seq. The most surprising finding was that `rerio` and WGBS calls were independent of the GC% context ("three-way_per-pos_beta.html", final plot, panel F; residuals were almost flat along origin). Something felt wrong, and we reached out to ONT to confirm how the `rerio` models were trained.

Through direct correspondence, we found out that `rerio` was trained on a mixture of WGBS and synthetic datasets ("synthetic" == unmethylated DNA was produced using PCR-ed reads, while methylated DNA was from M.SssI-treated reads). We interpreted this to mean that the ONT calls might tend to recapitulate properties of WGBS reads. This meant that in a three-way fight, `rerio` would probably take WGBS' side more than EM-seq's.

In the same email, ONT alerted us to the existence of another model that was purely trained on synthetic data: `remora`. We modified our methylation calling pipeline early 2022 to use pre-trained models from `remora` instead, to produce plots/conclusions in the manuscript.

For transparency, we include the previous analysis using `rerio` in the `rerio_beds` folder--the HTML files contain equivalent plots using the old model. The compressed BED files in this folder was from `rerio`, so one can regenerate the old plots simply by overwriting the compressed BEDs in the current folder; one can also do a head-to-head comparison of `rerio` and `remora` by comparing similarly named BED files (we did this as a sanity check, they look fairly similar). Feel free to peruse and draw your own conclusions!
