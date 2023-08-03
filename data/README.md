# `data/` folder #

This folder contains genomes (for mapping) and gene annotation files that are too large for GitHub, they are available at https://data.csiro.au/collection/csiro:58492, navigate to "Files" > `/data`.

## Genomes ##

`grch38p13_lambda_puc.fa.gz`: GRCh38 human genome, patch 13 + control sequences from NEB's EM-seq library prep kit (E7120) to compute non-conversion rates.

`homo_sapiens.45s.fa.gz`: 1 kb promoter + transcribed region of human 45S, accession KY962518. For short reads to map against.

`homo_sapiens.KY962518.last_11kb_first.fa.gz`: entire human 45S from accession KY962518, but the last 11 kb has been moved to the front of the sequence. For ONT long reads to map against.

## Annotation files ##

`gencode.v38.annotation.gtf.gz`: original GENCODE v38 GTF file.

`gencode.v38.annotation.sqlite` and `gencode.v38.annotation_2021-08-04.RData`: produced from R scripts in `01_txdb/`, for use in other R scripts.
