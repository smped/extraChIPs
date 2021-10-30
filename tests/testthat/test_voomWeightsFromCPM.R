bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
wc <- csaw::windowCounts(bamFiles, filter=1)
cpm <- edgeR::cpm(wc, log = TRUE)
el <- voomWeightsFromCPM(cpm, lib.size = wc$totals)

test_that("Correct structure parses", {
  expect_s4_class(voomWeightsFromCPM(cpm, lib.size = wc$totals), "EList")
})

test_that("lib.size must be provided", {
  expect_error(voomWeightsFromCPM(cpm))
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
    .voomChecks(matrix(rnorm(10), ncol = 2), c(1,1), c(1, 1), isLogCPM = FALSE),
    "Negative CPM values not allowed"
  )
  expect_true(
    .voomChecks(cpm, NULL, wc$totals, TRUE)
  )
  cpm_na <- cpm
  cpm_na[1] <- NA
  expect_error(.voomChecks(cpm_na), "NA values not allowed")
}
)

test_that("VoomWeights are returned", {
  elw <- voomWeightsFromCPM(cpm, lib.size = wc$totals, w0 = rep(1, ncol(cpm)))
  expect_true("sample.weights" %in% colnames(elw$targets))
})
