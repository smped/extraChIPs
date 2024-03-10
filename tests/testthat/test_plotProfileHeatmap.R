bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:500")
pd <- getProfileData(bw, gr, upstream = 100, bins = 10)

test_that("Correct plot is returned", {
  p <- plotProfileHeatmap(pd, profileCol = "profile_data")
  expect_equal(is(p), "ggside")
  expect_equal(nrow(p$data), 10)
  expect_equal(colnames(p$data), c("range", "score", "position", "bp"))
  expect_null(p$labels$y)
  expect_equal(levels(p$data$range), "chr2:400-599")
  expect_equal(is(p$facet), "FacetSideNull")
  pd$facet <- "a"
  p <- plotProfileHeatmap(pd, profileCol = "profile_data", facetY = "facet")
})

test_that("Histogram is omitted when requested", {
  p <- plotProfileHeatmap(pd, profileCol = "profile_data", summariseBy = "none")
  expect_equal(is(p), "gg")
  expect_equal(is(p$facet), "FacetNull")
})

test_that("Samples are facetted when provided as a list", {
  p <- plotProfileHeatmap(
    GRangesList(a = pd, b = pd), profileCol = "profile_data"
  )
  expect_equal(is(p$facet), "FacetSideGrid")
  expect_equal(levels(p$data$name), c("a", "b"))
  p <- plotProfileHeatmap(
      GRangesList(pd, pd), profileCol = "profile_data"
  )
  expect_equal(levels(p$data$name), c("X.1", "X.2"))
})

test_that("checkProfileDataFrames behaves as expected", {
  expect_true(.checkProfileDataFrames(pd$profile_data, "bp", "score"))
  expect_message(
    .checkProfileDataFrames(pd$profile_data, "b", "scor"),
    "Column b is missingColumn scor is missing"
    )
  expect_false(
    suppressMessages(.checkProfileDataFrames(pd$profile_data, "bp", "scor"))
  )

  dfl <- list(pd$profile_data[[1]], pd$profile_data[[1]])
  expect_true(.checkProfileDataFrames(dfl, "bp", "score"))

  dfl[[2]] <- dfl[[2]][1,]
  expect_message(
    .checkProfileDataFrames(dfl, "bp", "score"),
    "Each data element must have the same number of rows"
  )
  dfl[[2]] <- "a"
  expect_message(
    .checkProfileDataFrames(dfl, "bp", "score"),
    ".+ list element must contain DataFrame or data.frame"
  )
  dfl[[2]] <- dfl[[1]][,3:1]
  expect_message(
    .checkProfileDataFrames(dfl, "bp", "score"),
    "All elements must have the same column names"
  )

})

