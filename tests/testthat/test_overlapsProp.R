test_that("propOverlap returns the correct values",{
  x <- GRanges("chr1:1-10")
  y <- GRanges("chr1:1-5")
  expect_equal(propOverlap(x, y), 0.5)
  expect_equal(propOverlap(x, GRanges()), 0)
}
)

test_that("propOverlap returns NULL values correctly", {
  expect_equal(propOverlap(GRanges(), GRanges()), numeric())
})

