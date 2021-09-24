cpm <- matrix(rnorm(10), ncol = 2)
lib_sizes <- rep(5, 2)

test_that("isLogCPM behaves correctly", {
  expect_error(voomWeightsFromCPM(cpm, isLogCPM = FALSE, lib.size = lib_sizes))
  expect_s4_class(voomWeightsFromCPM(cpm, lib.size = lib_sizes), "EList")
})

test_that("lib.size must be provided", {
  expect_error(voomWeightsFromCPM(cpm))
})
