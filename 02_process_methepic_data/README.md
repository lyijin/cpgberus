# `02_process_methepic_data`: getting MethylationEPIC betas from four samples #

The MethylationEPIC data from four samples (WR025V1, WR025V9, WR069V1, WR069V9) were amongst a group of 24 other samples meant for another publication.

As (light-touch) batch correction is often done in pre-processing EPIC data, generation of this subset of data is not possible without providing raw data for all samples in the same batch. Batch correction cannot be done on just the four samples.

For this study, the exact genomic coordinates of the > 800k positions interrogated by MethylationEPIC and the beta values for the four samples is provided as a RData file (for use with R's `load()` function).
