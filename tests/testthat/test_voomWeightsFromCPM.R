bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
wc <- csaw::windowCounts(bamFiles, filter=1)
cpm <- edgeR::cpm(wc, log = TRUE)
el <- voomWeightsFromCPM(cpm, lib.size = wc$totals)

test_that("Correct structure parses", {
    expect_s4_class(voomWeightsFromCPM(cpm, lib.size = wc$totals), "EList")
    colnames(cpm) <- c()
    expect_equal(
        rownames(voomWeightsFromCPM(cpm, lib.size = wc$totals)$design),
        as.character(1:2)
    )
})

test_that("lib.size must be provided", {
    expect_error(voomWeightsFromCPM(cpm))
})

test_that("log transforms are consistent", {
    expect_equal(
        voomWeightsFromCPM(cpm, lib.size = wc$totals),
        voomWeightsFromCPM(2^cpm, lib.size = wc$totals, isLogCPM = FALSE)
    )
})

test_that("voomChecks error correctly", {
    expect_error(
        .voomChecks(matrix(cpm[,1], nrow = 1)),
        "Need at least two genes to fit a mean-variance trend"
    )
    expect_error(
        .voomChecks(cpm,isLogCPM = TRUE, w0 = 1),
        "Supplied weights do not match the data"
    )
    expect_error(
        .voomChecks(cpm, w0 = c(1, 1), wc$totals[1], TRUE),
        "Library sizes do not match the data"
    )
    expect_error(
        .voomChecks(cpm, w0 = c(1, 1), rep(0, ncol(cpm)), TRUE),
        "Library sizes must be > 0"
    )
    expect_error(
        .voomChecks(matrix(rnorm(100), ncol = 2), c(1,1), c(1, 1), isLogCPM = FALSE),
        "Negative CPM values not allowed"
    )
    expect_true(
        .voomChecks(cpm, NULL, wc$totals, TRUE)
    )
    cpm_na <- cpm
    cpm_na[1] <- NA
    expect_error(.voomChecks(cpm_na), "NA values not allowed")
    cpm_0 <- 2^cpm
    cpm_0[1, 1] <- 0
    expect_error(
        .voomChecks(cpm_0, c(1, 1), isLogCPM = FALSE),
        "Please ensure an offset is used for estimation"
    )
}
)

test_that("VoomWeights are returned", {
    elw <- voomWeightsFromCPM(cpm, lib.size = wc$totals, w0 = rep(1, ncol(cpm)))
    expect_true("sample.weights" %in% colnames(elw$targets))
})
