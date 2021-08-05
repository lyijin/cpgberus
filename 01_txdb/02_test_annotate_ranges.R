#!/usr/bin/env Rscript

"> 02_test_annotate_ranges.R <

Check that the compiled files work as expected.

Sources `annotate_ranges.R`, which doesn't do much when run on its own.

Original repo (that processed hg19 stuff)
https://bitbucket.csiro.au/users/ros259/repos/txdb/browse
" -> doc

suppressPackageStartupMessages({
    library(testthat)
    library(GenomicRanges)
})

build <- 38
path_base <- "/home/lie128/csiro/stopwatch/cpgberus/01_txdb"
path_code <- file.path(path_base)

source(file.path(path_code, "annotate_ranges.R"))


query <- GRanges(seqnames = c("chr2", rep("chr1", 9)),
                 ranges = IRanges(start=c(1, 1, 1, 1, 81, 1, 201, 50, 42, 50),
                                  width=c(100, 100, 100, 100, 20, 20, 100, 100, 1, 20)),
                 strand = c("+", "+", "-", "*", "+", "+", "+", "-", "+", "-"))

subject <- GRanges(seqnames = c(rep("chr1", 4)),
                   ranges = IRanges(start=c(60, 60, 50, 40),
                                    width=c(100, 50, 50, 30)),
                   strand = c("+", "+", "+", "-"))


test_that("coverage", {
    
    # Case for non-matching chromosome, should have 0 overlap
    expect_identical(coverageRatio(query[1], subject, ratio = FALSE), 0)
    
    # Cases for same region but different strand information in the subject
    expect_identical(coverageRatio(query[2:4], subject, ratio = FALSE), rep(61, 3))
    expect_identical(coverageRatio(query[2:4], subject[c(1, 3, 4)], ratio = FALSE), rep(61, 3))
    
    # Test for correct calcution of ratio versus width output
    expect_identical(coverageRatio(query[5:10], subject, ratio = FALSE), c(20, 0, 0, 100, 1, 20))
    expect_identical(coverageRatio(query[5:10], subject, ratio = TRUE), c(1, 0, 0, 1, 1, 1))

})
