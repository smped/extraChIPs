genes <- c("A10", "A2")

test_that("collapseGenes returns the correct output", {
  x <- collapseGenes(genes)
  expect_true(is(x, "glue"))
  expect_equal(as.character(x), "_A2_ and _A10_")
})
