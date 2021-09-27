nrows <- 200; ncols <- 4
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
df <- DataFrame(treat = c("A", "A", "B", "B"))
se <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  colData = df
)

test_that("Assay Density plots behave correctly", {

  expect_error(plotAssayDensities(se, colour = "col"))
  expect_error(plotAssayDensities(se, group = "col"))
  p <- plotAssayDensities(se)
  expect_equal(dim(p$data), c(512*4, 4))
  expect_equal(colnames(p$data), c("sample", "x", "y", "treat"))

})

test_that("Arguments are passed to density correctly", {
  p <- plotAssayDensities(se, n = 500)
  expect_equal(dim(p$data), c(500*4, 4))
})

test_that("Assay PCA plots error correctly", {
  expect_error(plotAssayPCA(se, colour = "col"))
  expect_error(plotAssayPCA(se, shape = "col"))
  expect_error(plotAssayPCA(se, label = "col"))
})
