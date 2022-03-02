test_that("overlapsProp returns the correct values",{
  x <- GRanges("chr1:1-10")
  y <- GRanges("chr1:1-5")
  expect_equal(overlapsProp(x, y), 0.5)
  expect_equal(overlapsProp(x, GRanges()), 0)
}
)

test_that("overlapsProp returns NULL values correctly", {
  expect_equal(overlapsProp(GRanges(), GRanges()), numeric())
})

