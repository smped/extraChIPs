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
  expect_error(.voomChecks(cpm[,1]))
  expect_error(.voomChecks(cpm, w0 = 1))
  expect_error(.voomChecks(cpm, w0 = c(1, 1), wc$totals[1]))
  expect_error(.voomChecks(cpm, w0 = c(1, 1), wc$totals))
  expect_error(
    .voomChecks(matrix(rnorm(10), ncol = 2), c(1,1), c(1, 1), isLogCPM = FALSE)
  )
  expect_true(
    .voomChecks(cpm, NULL, wc$totals, TRUE)
  )
}
)
