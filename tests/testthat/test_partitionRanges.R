x <- GRanges(c("chr1:1-10", "chr1:6-15"))
x$id <- paste0("range", seq_along(x))
y <- GRanges(c("chr1:2-5", "chr1:5-12"))
y$id <- paste0("range", seq_along(y))

test_that("Partitions are correct", {
  gr <- partitionRanges(x, y)
  expect_equal(
    gr$id.x, paste0("range", c(1, 1, 1, 2, 2))
  )
  expect_equal(
    gr$id.y, c(NA, paste0("range", c(1, 2, 2)), NA)
  )
})

test_that("Empty ranges return x unchanged", {
  expect_equal(
    partitionRanges(x, GRanges()), x
  )
})
