stopifnot(require(GenomicRanges))

loadLatest <- function(mypath, suffix) {
    mypattern <- paste(suffix, "_", ".*.RData", sep="")
    latest_file <- tail(sort(list.files(mypath, pattern = mypattern)), n = 1)
    file.path(mypath, latest_file)
}


nearestAnnotation <- function(query, subject, keep_hits = FALSE) {
    
    gr_dists <- as.data.frame(distanceToNearest(query, subject))

    subject_ovl <- as.data.frame(mcols(subject)[gr_dists$subjectHits, ])
    names(subject_ovl) <- names(mcols(subject))

    if(keep_hits == FALSE) {
        gr_dists <- data.frame(distance = gr_dists$distance)
    }
    mcols(query) <- cbind(mcols(query), gr_dists, subject_ovl)
    query
}


coverageRatio <- function(query, subject, ratio=TRUE) {
    
    strand(subject) <- "*"
    subject <- reduce(subject)
    names(query) <- 1:length(query)
    
    # Find the pairs of overlapping ranges
    qs_pairs <- findOverlapPairs(query, subject)
    # Map this back to the original query
    query_map <- as.integer(names(S4Vectors::first(qs_pairs)))

    if(ratio) {
        x <- numeric(length=length(query))
        x[query_map] <- width(pintersect(qs_pairs)) / width(S4Vectors::first(qs_pairs))
    } else {
        x <- integer(length=length(query))
        x[query_map] <- width(pintersect(qs_pairs))
    }
    x
}


annotateDMRs <- function(dmrs, tx, cpgIslands, cpgShores, tss_proximal=2000) {
    
    # Want to do most analysis for "just" protein coding as well as all genes
    tx_pc <- tx[tx$transcript_type == "protein_coding"]
    
    numGeneTSS <- function(this_tx) {
    # Number of different genes TSSs within the DMR
        tss_ovl <- as.matrix(findOverlaps(dmrs, resize(this_tx, 1, fix="start")))
        tss_ovl <- tapply(this_tx$transcript_id[tss_ovl[, 2]], tss_ovl[, 1], function(x) length(unique(x)))
        nGeneTSS <- integer(length(dmrs))
        nGeneTSS[as.integer(names(tss_ovl))] <- unname(tss_ovl)
        nGeneTSS
    }

    distTSS <- function(this_tx) {
        # Distance to closest TSS
        reg.dist <- as.data.frame(distanceToNearest(dmrs,
                                                    resize(this_tx, 1,
                                                           fix="start"))
                                  )
        i <- reg.dist$subjectHits
        df <- data.frame(distanceTSS = reg.dist$distance,
                         tx_name = this_tx$tx_name[i],
                         transcript_type = this_tx$transcript_type[i],
                         gene_id = this_tx$gene_id[i],
                         gene_name = this_tx$gene_name[i],
                         stringsAsFactors = FALSE
                         )
        df
    }
    
    
    dmrs$nGeneTSS <- numGeneTSS(tx)
    dmrs$nProtGeneTSS <- numGeneTSS(tx_pc)
    
    dist_df <- distTSS(tx)
    mcols(dmrs) <- cbind(mcols(dmrs), dist_df)
    
    distpc_df <- distTSS(tx_pc)
    names(distpc_df) <- paste(names(distpc_df), "prot", sep="_")
    mcols(dmrs) <- cbind(mcols(dmrs), distpc_df[, c(1:2, 4:5)])
    
    # Distance to closest CpG island
    reg.dist <- distanceToNearest(dmrs, cpgIslands)
    dmrs$distanceCpGi <- NA
    dmrs$CpGi <- NA
    dmrs$distanceCpGi[queryHits(reg.dist)] <- mcols(reg.dist)$distance
    dmrs$CpGi[queryHits(reg.dist)] <- subjectHits(reg.dist)
    
    # % promoter/genebody/intergenic
    
    tss_gr <- resize(tx, 1, fix="start")
    promoters_gr <- reduce(resize(tss_gr, tss_proximal * 2, fix="center"))

    genebody_gr <- setdiff(reduce(tx), promoters_gr)
    genome_gr <- GRanges(seqlevels(tx), IRanges(1, seqlengths(tx)))
    genic_gr <- union(promoters_gr, tx)
    strand(genic_gr) <- "*"
    intergenic_gr <- setdiff(genome_gr, genic_gr)
    
    dmrs$promoter <- coverageRatio(dmrs, promoters_gr)
    dmrs$genebody <- coverageRatio(dmrs, genebody_gr)
    dmrs$intergenic <- coverageRatio(dmrs, intergenic_gr)
    
    # % CpG island/CpG shore/nonCpG
    dmrs$CpGisland <- coverageRatio(dmrs, cpgIslands)
    dmrs$CpGshores <- coverageRatio(dmrs, cpgShores)
    cpgSet <- union(cpgIslands, cpgShores)
    dmrs$nonCpG <- coverageRatio(dmrs, setdiff(genome_gr, cpgSet))
    
    dmrs
}

build <- 38
path_base <- "/home/lie128/csiro/stopwatch/cpgberus/01_txdb"
path_data <- file.path(path_base, "../data")
gtf_suffix <- paste("gencode.v", build, ".annotation", sep="")
load(loadLatest(mypath = path_data, suffix = gtf_suffix))
