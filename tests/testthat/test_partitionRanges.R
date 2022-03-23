x <- GRanges(c("chr1:1-10:+", "chr1:5-10:+", "chr1:5-6:+"))
x$id <- "gene1"
x$tx_id <- paste0("transcript", seq_along(x))
y <- GRanges(c("chr1:1-5", "chr1:9-12"))
y$id <- paste0("range", seq_along(y))
gr <- partitionRanges(x, y)

test_that("Partitions are correct", {
  expect_equal(gr$id.x, rep("gene1", 5))
  expect_equal(gr$id.y, c("range1", "range1", NA, NA, "range2"))
  expect_equal(
    gr$tx_id,
    new(
      "CompressedCharacterList",
      elementType = "character", elementMetadata = NULL,
      metadata = list(),
      unlistData = c(
        "transcript1", "transcript2", "transcript3", "transcript3",
        "transcript1", "transcript2", "transcript1", "transcript2"
        ),
      partitioning = new(
        "PartitioningByEnd", end = c(1L, 3L, 4L, 6L, 8L), NAMES = NULL,
        elementType = "ANY", elementMetadata = NULL, metadata = list())
    )
  )
})

test_that("Arguments behave correctly", {
  gr_un <- partitionRanges(x, y, ignore.strand = TRUE)
  expect_equal(
    granges(partitionRanges(x, y, y_as_both = FALSE)),
    sort(granges(x))
  )
  expect_equal(mcols(gr), mcols(gr_un))
  expect_equal(as.character(strand(gr_un)), rep("*", 5))
  expect_equal(
    partitionRanges(x, y, simplify = FALSE)$id.x, rep("gene1", 8)
  )
})

test_that("Empty ranges return x unchanged", {
  expect_equal(
    partitionRanges(x, GRanges()), x
  )
})
