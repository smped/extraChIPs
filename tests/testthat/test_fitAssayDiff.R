nrows <- 200; ncols <- 6
counts <- matrix(as.integer(runif(nrows * ncols, 0, 1e4)), nrows)
counts[1, 1] <- 0
colnames(counts) <- paste0("Sample_", seq_len(ncols))
df <- DataFrame(treat = c("A", "A", "A", "B", "B", "B"))
df$treat <- as.factor(df$treat)
df$totals <- colSums(counts)
se <- SummarizedExperiment(
    assays = SimpleList(counts = counts), colData = df
)
X <- model.matrix(~treat, colData(se))
assay(se, "logCPM") <- edgeR::cpm(counts, log = TRUE, prior.count = 0)

test_that("Incorrect assay types error correctly", {

    expect_error(fitAssayDiff(se, "counts", method = "lt", design = X))
    expect_error(fitAssayDiff(se, "logCPM", method = "qlf", design = X))
    expect_error(fitAssayDiff(se, lib.size = "totals"))

})

test_that("Results appear correct", {

    ## glmFits
    new_se <- suppressMessages(
        fitAssayDiff(se, design = X, fc = 1.2)
    )
    row_data <- rowData(new_se)
    expect_equal(colnames(row_data), c("logFC", "logCPM", "PValue", "FDR"))
    new_se <- suppressMessages(
        fitAssayDiff(se, design = X, lib.size = NULL, fc = 1)
    )
    row_data <- rowData(new_se)
    expect_equal(colnames(row_data), c("logFC", "logCPM", "F", "PValue", "FDR"))
    # Limma Fits
    new_se <- fitAssayDiff(se, assay = "logCPM", method = "lt", design = X)
    row_data <- rowData(new_se)
    expect_equal(colnames(row_data), c("logFC", "logCPM", "t", "PValue", "FDR"))

})

test_that("setting asRanges = TRUE works", {
    expect_warning(fitAssayDiff(se, design = X, asRanges = TRUE))
    rowRanges(se) <- GRanges(paste("chr1:", seq_len(nrows)))
    gr <- fitAssayDiff(se, design = X, asRanges = TRUE)
    expect_true(is(gr, "GenomicRanges"))
})

test_that("Coefs & the design are checked", {

    expect_error(fitAssayDiff(se, lib.size = NULL))
    expect_error(
        fitAssayDiff(se, lib.size = NULL, design = X, coef = ncol(X) + 1)
    )
    expect_error(fitAssayDiff(se, lib.size = NULL, design = X, coef = ""))

})
