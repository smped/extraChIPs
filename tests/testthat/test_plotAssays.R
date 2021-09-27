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

test_that("show_points behaves as expected", {
  p <- plotAssayPCA(se)
  expect_equal(length(p$layers), 1)
  expect_true(is(p$layers[[1]]$geom, "GeomPoint"))
  expect_null(p$mapping$colour)
  p <- plotAssayPCA(se, show_points = FALSE)
  expect_equal(length(p$layers), 0)
})

test_that("colours are added correctly", {
  p <- plotAssayPCA(se, colour = "treat")
  expect_equal(rlang::as_label(p$mapping$colour), "treat")
  expect_equal(
    grepl("PC", unlist(p$labels)), c(TRUE, TRUE, FALSE)
  )
  expect_equal(p$labels$colour, "treat")
})

test_that("labels repel correctly", {
  p <- plotAssayPCA(se, label = "sample")
  expect_equal(length(p$layers), 2)
  expect_s3_class(p$layers[[2]]$geom, "GeomTextRepel")
  p <- plotAssayPCA(se, label = "sample", show_points = FALSE)
  expect_equal(length(p$layers), 1)
  expect_s3_class(p$layers[[1]]$geom, "GeomText")
})

test_that("data is transformed correctly", {
  expect_error(plotAssayPCA(se, trans = ""))
  ## Still need to test this
})
