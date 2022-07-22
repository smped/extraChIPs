bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:1000")

test_that("Check structure of output", {
  bwfl <- BigWigFileList(c(a = bw, b = bw))
  out <- getProfileData(bwfl, gr, upstream = 500, bins = 10)
  expect_equal(names(out), c("a", "b"))
  expect_s4_class(out, "GRangesList")
  expect_s4_class(out$a$profile_data, "SplitDataFrameList")
  expect_equal(
    colnames(out$a$profile_data[[1]]), c("score", "position", "bp")
  )
  expect_equal(dim(out$a$profile_data[[1]]), c(10, 3))
})

test_that("log transformation works", {
    pd <- getProfileData(bw, gr, upstream = 500, bins = 10, log = FALSE)
    expect_equal(range(pd$profile_data[[1]]$score), c(-0.75, 0))
    pd <- getProfileData(bw, gr, upstream = 500, bins = 10, log = TRUE)
    expect_equal(range(pd$profile_data[[1]]$score), c(-2, 0))
})

test_that("Errors", {
  expect_error(getProfileData("", gr))
  expect_error(getProfileData(bw))
  expect_error(getProfileData(bw, gr, offset = "a"))
})
