nrows <- 200
ncols <- 4
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
df <- DataFrame(treat = c("A", "A", "B", "B"))
se <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  colData = df
)

test_that("Assay Density plots behave correctly", {

  expect_error(plotAssayDensities(se, colour = "col"))
  p <- plotAssayDensities(se)
  expect_equal(dim(p$data), c(512*4, 4))
  expect_equal(colnames(p$data), c("colnames", "x", "y", "treat"))
  expect_equal(
    unlist(p$labels),
    c(x = "counts", y = "Density", group = "colnames")
  )
  p <- plotAssayDensities(se, colour = "treat", linetype = "treat")
  expect_equal(
    unlist(p$labels),
    c(
      x = "counts", y = "Density", colour = "treat", linetype = "treat",
      group = "colnames"
    )
  )

})

test_that("Assay Density transformations error", {
  expect_error(plotAssayDensities(se, trans = ""))
  p <- plotAssayDensities(se, trans = "log2")
  expect_true(median(p$data$x) < log2(1e4))
  expect_equal(p$labels$x, "log2 counts")
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

test_that("plotAssayRle errors correctly", {
  err <- "'arg' should be one of \"sample\", \"treat\""
  expect_error(plotAssayRle(se, "counts", x_col = ""), err)
  expect_error(plotAssayRle(se, "counts", colour = ""), err)
  expect_error(plotAssayRle(se, "counts", fill = ""), err)
  expect_error(plotAssayRle(se, "counts", rle_group = ""), err)
  expect_error(
    plotAssayRle(se, "counts", trans = "a"),
    "object 'a' of mode 'function' was not found"
  )
  expect_error(
    plotAssayRle(se, "counts", trans = "mean"),
    "This transformation is not applicable"
  )

})

test_that("plotAssayRle creates a plot", {
  p <- plotAssayRle(se, "counts", fill = "treat", n_max = 100)
  expect_true(is(p, "gg"))
  expect_equal(dim(p$data), c(100*ncols, 6))
  expect_equal(
    unlist(p$labels),
    c(y = "RLE", fill = "treat", colour = "colour", x = "sample")
  )
})
