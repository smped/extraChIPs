## These are extremeley rudimentary tests, just probing the input classes
## More detailed tests checking the output consistency are still required
## Detailed tests for bam files are also required

test_that("Erroneous input specs are caught", {

    rse <- SummarizedExperiment(rowRanges = GRanges())
    rse2 <- SummarizedExperiment(rowRanges = GRanges("chr1:1"))
    expect_error(dualFilter(""))
    expect_error(dualFilter(rse, ""))
    expect_error(dualFilter(rse, rse, ""))
    expect_error(dualFilter(rse, rse, GRanges(), ""))
    expect_error(dualFilter(rse, rse, GRanges(), 0.5, ""))
    expect_error(dualFilter(rse, rse2, GRanges(), 0.5, TRUE))
    expect_error(dualFilter(rse, rse, GRanges(), 0.5, TRUE))

})

## Now the sections from the vignette
library(tidyverse)
library(Rsamtools)
library(csaw)
library(BiocParallel)
library(rtracklayer)
bfl <- system.file(
    "extdata", "bam", c("ex1.bam", "ex2.bam", "input.bam"), package = "extraChIPs"
) %>%
    BamFileList()
names(bfl) <- c("ex1", "ex2", "input")
rp <- readParam(
    pe = "none",
    dedup = TRUE,
    restrict = "chr10"
)
wincounts <- windowCounts(
    bam.files = bfl,
    spacing = 60,
    width = 180,
    ext = 200,
    filter = 1,
    param = rp
)
wincounts$totals <- c(964076L, 989543L, 1172179L)
wincounts$sample <- colnames(wincounts)
wincounts$treat <- as.factor(c("ctrl", "treat", NA))
peaks <- import.bed(
    system.file("extdata", "peaks.bed.gz", package = "extraChIPs")
)
peaks <- granges(peaks)

test_that("Missing bam files error", {

    path <- unique(dirname(colData(wincounts)$bam.files))
    colData(wincounts)$bam.files <- basename(colData(wincounts)$bam.files)
    expect_error(
        dualFilter(
            x = wincounts[, !is.na(wincounts$treat)],
            bg = wincounts[, is.na(wincounts$treat)],
            ref = peaks,
            q = 0.8 # Better to use q = 0.5 on real data
        ),
        'all\\(file.exists.+is not TRUE'
    )

    colData(wincounts)$bam.files[1:2] <- file.path(
        path, colData(wincounts)$bam.files[1:2]
    )
    expect_error(
        dualFilter(
            x = wincounts[, !is.na(wincounts$treat)],
            bg = wincounts[, is.na(wincounts$treat)],
            ref = peaks,
            q = 0.8 # Better to use q = 0.5 on real data
        ),
        'all\\(file.exists.+is not TRUE'
    )

})

test_that("No peak overlaps errors", {
    expect_error(
        dualFilter(
            x = wincounts[, !is.na(wincounts$treat)],
            bg = wincounts[, is.na(wincounts$treat)],
            ref = GRanges(),
            q = 0.8 # Better to use q = 0.5 on real data
        ),
        'any.+is not TRUE'
    )
})

test_that("dualFilter runs as expected", {

    filtcounts <- dualFilter(
        x = wincounts[, !is.na(wincounts$treat)],
        bg = wincounts[, is.na(wincounts$treat)],
        ref = peaks,
        q = 0.8 # Better to use q = 0.5 on real data
    )
    expect_equal(dim(filtcounts), c(108, 2))
    expect_equal(assayNames(filtcounts), c("counts", "logCPM"))
    expect_equal(colnames(rowData(filtcounts)), "overlaps_ref")

})


