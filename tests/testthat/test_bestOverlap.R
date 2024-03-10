gr <- GRanges("chr1:1-10")
gr_cat <- GRanges(c("chr1:2-10", "chr1:5-10"))
gr_cat$category <- c("a", "b")
grl <- splitAsList(gr_cat, gr_cat$category)

test_that("bestOverlaps errors when expected", {

  expect_error(bestOverlap(gr, gr_cat, var = ""))
  expect_error(
    bestOverlap(gr, setNames(grl, c())), "'y' must be a named GRangesList"
  )

})

test_that("Correct values are returned by bestOverlap", {

  x <- bestOverlap(gr, gr_cat, var = "category")
  expect_equal(x, "a")
  expect_equal(bestOverlap(GRanges(), grl), character())
  expect_equal(
    bestOverlap(gr, GRangesList(a = GRanges()), missing = "none"), "none"
  )

})
