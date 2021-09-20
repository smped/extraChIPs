x <- GRanges(c("chr1:1-10", "chr1:101-110", "chr1:1001-1010", "chr2:1-10"))

test_that("Exclude works correctly", {

  expect_length(stitchRanges(x, exclude = GRanges("chr1:200:+")), 3)
  expect_length(stitchRanges(x), 2)

})
