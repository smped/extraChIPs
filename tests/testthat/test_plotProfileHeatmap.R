bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:500")
pd <- getProfileData(bw, gr, upstream = 100, bins = 10)

test_that("Correct plot is returned", {
  p <- plotProfileHeatmap(pd$profile_data, ids = granges(pd))
  expect_equal(is(p), "ggside")
  expect_equal(nrow(p$data), 10)
  expect_equal(colnames(p$data), c("id", "score", "position", "bp"))
  expect_null(p$labels$y)
  expect_equal(levels(p$data$id), "chr2:400-599")
  expect_equal(is(p$facet), "FacetSideNull")
})

test_that("Histogram is omitted when requested", {
  p <- plotProfileHeatmap(pd$profile_data, ids = granges(pd), addHist = FALSE)
  expect_equal(is(p), "gg")
  expect_equal(is(p$facet), "FacetNull")
})

test_that("Samples are facetted when provided as a list", {
  p <- plotProfileHeatmap(
    list(a = pd$profile_data, b = pd$profile_data), ids = granges(pd)
  )
  expect_equal(is(p$facet), "FacetSideWrap")
  expect_equal(levels(p$data$name), c("a", "b"))
})

test_that("IDs are specified", {
  expect_error(plotProfileHeatmap(pd$profile_data))
})

test_that("Missing names are added", {
  expect_message(
    plotProfileHeatmap(
      list(pd$profile_data, pd$profile_data), ids = granges(pd)
    ),
    "No names provided"
  )
})

test_that("xCol is correct", {
  expect_error(plotProfileHeatmap(pd$profile_data, ids = granges(pd), xCol = ""))
})

test_that("Each sample/ranges are the same size", {
  dodge <- pd$profile_data
  dodge[[1]] <- dodge[[1]][1,]
  expect_error(
    plotProfileHeatmap(
      list(a = pd$profile_data, b = dodge), ids = granges(pd)
    )
  )
  dodge <- pd$profile_data
  dodge[[2]] <- dodge[[1]]
  expect_error(
    plotProfileHeatmap(
      list(a = pd$profile_data, b = dodge), ids = granges(pd)
    )
  )
})
