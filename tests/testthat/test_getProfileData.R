bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:500")

test_that("Check structure of output", {
  bwfl <- BigWigFileList(c(a = bw, b = bw))
  out <- getProfileData(bwfl, gr, upstream = 100, bins = 10)
  expect_equal(names(out), c("a", "b"))
  expect_s4_class(out, "GRangesList")
  expect_s4_class(out$a$profile_data, "DataFrameList")
  expect_equal(
    colnames(out$a$profile_data[[1]]), c("score", "position", "bp")
  )
  expect_equal(dim(out$a$profile_data[[1]]), c(10, 3))
})

test_that("Paths behave correctly", {
  expect_s4_class(getProfileData(bw, gr, upstream = 10, bins = 10), "GRanges")
  expect_s4_class(getProfileData(c(bw, bw), gr, upstream = 10, bins = 10), "GRangesList")
})

test_that("Errors", {
  expect_error(getProfileData("", gr))
  expect_error(getProfileData(bw))
})
