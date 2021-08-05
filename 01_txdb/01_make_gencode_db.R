#!/usr/bin/env Rscript

"> 01_make_gencode_db.R <

When given a converted genome and GENCODE annotations, create an annotated
genome in the form sqlitedb and RData files.

Original repo (that processed hg19 stuff)
https://bitbucket.csiro.au/users/ros259/repos/txdb/browse
" -> doc

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(rtracklayer)
    library(BSgenome.Hsapiens.UCSC.hg38)
})


#####  Paths  #####
path_base <- "/home/lie128/csiro/stopwatch/cpgberus/01_txdb"
path_data <- file.path(path_base, "../data")


#####  Functions  #####
extractAttributes <- function(df, threads=10) {
    
    attributesForRow <- function(x) {
        # Should always have a key-value pair split by space
        y <- strsplit(x, " ")
        z <- sapply(y, '[[', 2)
        names(z) <- sapply(y, '[[', 1)
        z
    }
    
    # Sometimes there is a space after the semi-colon
    x <- strsplit(df$attribute, "; ?")
    
    attributes <- mclapply(x, attributesForRow, mc.preschedule = TRUE,
                           mc.cores = threads)
    attribute_set <- unique(unlist(lapply(attributes, names)))
    
    m <- matrix(NA, nrow=length(attributes), ncol=length(attribute_set),
                dimnames=list(1:length(attributes), attr=attribute_set))
    
    for(i in 1:length(attributes)) {
        a <- attributes[[i]]
        m[i, names(a)] <- as.vector(a)
    }
    data.frame(m, stringsAsFactors = FALSE)
}


rangesTxdb <- function(txdb=gencode, gtf, feature="gene") {
    
    stopifnot(feature %in% c("gene", "transcript"))
    
    if(feature == "gene") {
        gr <- genes(txdb)
    }
    if(feature == "transcript") {
        gr <- transcripts(txdb)
    }
    
    df <- gtf[gtf$feature == feature, ]
    df_attr <- extractAttributes(df)
    
    stopifnot(identical(nrow(df), length(gr)))
    
    merge_by <- NA
    last_mcol <- mcols(gr)[, ncol(mcols(gr))]
    for(i in 1:ncol(df_attr)) {
        is_same <- last_mcol %in% df_attr[, i]
        if(sum(is_same) == nrow(mcols(gr))) {
            merge_by <- i
            break
        }
    }
    
    attr_index <- match(last_mcol, df_attr[, i])
    mcols(gr) <- cbind(mcols(gr), df_attr[attr_index, ])
    gr
}


#####  Load files and make TxDb  #####

# Use release 38 for hg38; https://www.gencodegenes.org/human/release_38.html
build <- 38
source_url <- paste("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                    build, "/gencode.v", build, ".annotation.gtf.gz", sep="")

gtf_suffix <- paste("gencode.v", build, ".annotation", sep="")
gtf_file <- file.path(path_data, paste(gtf_suffix, 'gtf.gz', sep="."))
rda_file <- file.path(path_data, paste(gtf_suffix, "_", Sys.Date(), '.RData', sep=""))
sqlite_file <- file.path(path_data, paste(gtf_suffix, 'sqlite', sep="."))

if(!file.exists(gtf_file)) {
    download.file(source_url, destfile = gtf_file, method="wget")    
}

gencode <- makeTxDbFromGFF(gtf_file, 
                           format = "gtf", dataSource = source_url,
                           organism = "Homo sapiens",
                           chrominfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38),
                           miRBaseBuild = "Homo sapiens")

saveDb(gencode, file=sqlite_file)


#####  Annotate  #####
gtf <- read.delim(gtf_file,
                  stringsAsFactors=FALSE, comment.char="#", header = FALSE,
                  colClasses = c("character", "character", "character",
                                 "integer", "integer", "character",
                                 "character", "character", "character"),
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attribute"))
attr(gtf, "source_url") <- source_url
attr(gtf, "gtf_file") <- gtf_file

gencode_gene <- rangesTxdb(txdb=gencode, gtf, feature="gene")
gencode_tx <- rangesTxdb(txdb=gencode, gtf, feature="transcript")

#gencode_cds <- cds(gencode)
#gencode_exon <- exons(gencode)
#gencode_intergenic <- gaps(gencode_gene)
#gencode_introns <- unlist(intronsByTranscript(gencode))


#####  CpGislands  #####

session <- browserSession("UCSC")
genome(session) <- "hg38"
cpgIslands <- track(session, "cpgIslandExt")
cpgIslands <- as(cpgIslands, "GRanges")

#CpG island shores & 5kb
shore_size <- 2000

cpgShores <- setdiff(resize(cpgIslands, width(cpgIslands) + shore_size * 2,
                            fix="center"), cpgIslands)
cpgShores <- trim(cpgShores)

cpg5kb <- resize(cpgIslands, 5000)
cpg5kb <- trim(cpg5kb)

save(gencode, gencode_gene, gencode_tx, cpgIslands, cpgShores, cpg5kb,
     file=rda_file)  # gencode_cds, gencode_exon, gencode_intergenic, gencode_introns
