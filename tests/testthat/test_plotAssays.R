nrows <- 200
ncols <- 4
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
df <- DataFrame(treat = c("A", "A", "B", "B"))
se <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  colData = df
)
se$totals <- colSums(counts)

test_that("Assay Density plots behave correctly", {

  expect_error(plotAssayDensities(se, colour = "col"))
  p <- plotAssayDensities(se)
  expect_equal(dim(p$data), c(800, 4))
  expect_equal(colnames(p$data), c("colnames", "vals", "treat", "totals"))
  lab_vec <- c(
      x = "counts", y = "Density", group = "colnames", colour = "colour",
      fill = "fill", alpha = "alpha", linetype = "linetype",
      linewidth = "linewidth", weight = "weight"
  )
  expect_equal(unlist(lapply(p$labels, as.character)), lab_vec)
  p <- plotAssayDensities(
      se, colour = "treat", group = "treat", linetype = "treat", fill = "treat"
  )
  expect_true(
      unique(unlist(p$labels[c("group", "linetype", "fill")])) == "treat"
  )

  se$vals <- runif(ncol(se))
  expect_warning(
      plotAssayDensities(se),
      "Any columns named 'colnames' or 'vals' will be overwritten"
  )

})

test_that("Assay Density transformations error", {
  expect_error(plotAssayDensities(se, trans = ""))
  p <- plotAssayDensities(se, trans = "log2")
  expect_true(median(p$data$vals) < log2(1e4))
  expect_equal(p$labels$x, "log2 counts")
  expect_error(plotAssayDensities(se, trans = "max"), "This transformation")
})

test_that("fixed params are handled as expected", {
    p <- plotAssayDensities(se, colour = "#809050", linewidth = 2, linetype = 2)
    expect_equal(
        unlist(p$layers[[1]]$aes_params),
        c(colour = "#809050", linetype = "2", linewidth = "2")
    )
    p <- plotAssayDensities(se, fill = "#809050", alpha = 0.2)
    expect_equal(
        unlist(p$layers[[1]]$aes_params),
        c(fill = "#809050", alpha = "0.2")
    )
})


test_that("Assay PCA plots error correctly", {
  expect_error(plotAssayPCA(se, colour = "a"))
  expect_error(plotAssayPCA(se, shape = "a"))
  expect_error(plotAssayPCA(se, label = "a"))
})

test_that("show_points behaves as expected", {
  p <- plotAssayPCA(se)
  expect_equal(length(p$layers), 2)
  expect_true(is(p$layers[[1]]$geom, "GeomPoint"))
  expect_null(p$mapping$colour)
  p <- plotAssayPCA(se, show_points = FALSE)
  expect_equal(length(p$layers), 1)
})

test_that("colours/size are added correctly", {

  p <- plotAssayPCA(se, colour = "treat", shape = 4, size = "totals")
  mappings <- c(
      x = "PC1", y = "PC2", colour = "treat", size = "totals", label = "colnames"
  )
  expect_equal(vapply(p$mapping, rlang::as_label, character(1)), mappings)
  expect_equal(rlang::as_label(p$mapping$colour), "treat")
  expect_equal(rlang::as_label(p$mapping$size), "totals")
  expect_equal(
    grepl("PC", unlist(p$labels))[1:4], c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_equal(p$labels$colour, "treat")
  expect_true(p$layers[[1]]$aes_params$shape == 4)
  # expect_equal(p$labels$size, "log10(totals)")

  p <- plotAssayRle(se, colour = "blue", fill = "white", trans = "log2")
  expect_equal(p$mapping$fill, "white")
  expect_equal(p$mapping$colour, "blue")
  expect_true(grepl("log2", p$labels$y))

})



test_that("labels repel correctly", {
  p <- plotAssayPCA(se, label = "treat")
  expect_equal(length(p$layers), 2)
  expect_s3_class(p$layers[[2]]$geom, "GeomTextRepel")
  p <- plotAssayPCA(se, label = "treat", show_points = FALSE)
  expect_equal(length(p$layers), 1)
  expect_s3_class(p$layers[[1]]$geom, "GeomText")
})

test_that("data is transformed correctly", {
    expect_error(plotAssayPCA(se, trans = ""))
    expect_error(plotAssayPCA(se, trans = "max"), "This transformation is not")
    expect_true(is(plotAssayPCA(se, trans = "log2"), "gg"))
})

test_that("plotAssayRle errors correctly", {
  err <- "'arg' should be one of .+"
  expect_error(plotAssayRle(se, "counts", colour = "a"), err)
  expect_error(plotAssayRle(se, "counts", fill = "a"), err)
  expect_error(plotAssayRle(se, "counts", rle_group = "a"), err)
  expect_error(plotAssayRle(se, "counts", by_x = "a"), err)
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
  expect_equal(dim(p$data), c(100*ncols, 5))
  expect_equal(
    unlist(p$labels),
    c(x = "Sample", y = "RLE (counts)", fill = "treat", colour = "colour")
  )
  p <- plotAssayRle(se, "counts", by_x = "treat")
  expect_equal(unlist(p$labels)[["x"]], "treat")
})

test_that("plotAssayHeatmap creates a plot", {

    p <- plotAssayHeatmap(se[1:10,], trans = "log10")
    expect_equal(c("index", "colnames", "value", "treat", "totals"), colnames(p$data))
    expect_equal(dim(p$data), c(40, 5))
    rowRanges(se) <- GRanges(paste0("chr:", seq_len(nrow(se))))
    p <- plotAssayHeatmap(se[1:10,], trans = "log10")
    expect_equal(c("range", "colnames", "value", "treat", "totals"), colnames(p$data))
    expect_true(is(p, "gg"))
    p <- plotAssayHeatmap(se[1:10,], trans = "log10", ysideline = TRUE)
    expect_true(is(p, "ggside"))
    expect_error(
        plotAssayHeatmap(se, n_max = 1),
        "Only 1 ranges can be drawn. Please change the n_max parameter if you wish to draw more."
    )

})
