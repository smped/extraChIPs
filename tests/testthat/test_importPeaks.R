test_that("narrowPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )
  expect_true(.isValidNarrow(fl))
  expect_false(.isValidBroad(fl))

  grl <- importPeaks(fl, type = "narrow")
  gr <- unlist(grl)

  expect_equal(length(grl), 1L)
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 5L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))
})

test_that("broadPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.broadPeak", package = "extraChIPs"
  )
  expect_false(.isValidNarrow(fl))
  expect_true(.isValidBroad(fl))

  grl <- importPeaks(fl, type = "broad")
  gr <- unlist(grl)

  expect_equal(length(grl), 1L)
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 4L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))
})

test_that("seqinfo objects behave correctly for narrowPeak files", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )
  ## Succeed
  expect_type(
    importPeaks(
      fl, type = "narrow",
      seqinfo = Seqinfo(seqnames = "chr1", seqlengths = 8100000)
    ),
    "S4"
  )
  ## Empty GRanges
  expect_message(
    importPeaks(
      fl, type = "narrow",
      seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "coarse"
    ),
    "No ranges match the supplied seqinfo object"
  )
  ## Error
  expect_error(
    importPeaks(
      fl, type = "narrow",
      seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "error"
    )
  )
  ## Error
  expect_error(
    importPeaks(fl, "narrow", seqinfo = NULL)
  )
})

test_that("blacklists behave correctly",{
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )
  expect_error(
    importPeaks(gr, "narrow", blacklist = NULL)
  )
  expect_equal(
    length(importPeaks(fl, blacklist = GRanges("chr1:1299700"))),
    1L
  )
})

test_that("seqinfo objects behave correctly for broadPeak files", {
  fl <- system.file(
    "extdata/testFiles", "test.broadPeak", package = "extraChIPs"
  )
  ## Succeed
  expect_type(
    importPeaks(
      fl, type = "broad",
      seqinfo = Seqinfo(seqnames = "chr1", seqlengths = 8100000)
    ),
    "S4"
  )
  ## Empty GRanges
  expect_message(
    importPeaks(
      fl, type = "broad",
      seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "coarse"
    ),
    "No ranges match the supplied seqinfo object"
  )
  ## Error
  expect_error(
    importPeaks(
      fl, type = "broad",
      seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "error"
    )
  )
  ## Error
  expect_error(
    importPeaks(fl, "broad", seqinfo = NULL)
  )
})

test_that("blacklists behave correctly for broadPeak files",{
  fl <- system.file(
    "extdata/testFiles", "test.broadPeak", package = "extraChIPs"
  )
  expect_error(
    importPeaks(gr, "broad", blacklist = NULL)
  )
  expect_equal(
    length(importPeaks(fl, "broad", blacklist = GRanges("chr1:1299700"))),
    1L
  )
})

test_that("Empty files parse empty GRanges", {
  fl <- file.path(tempdir(), "empty.txt")
  file.create(fl)
  np1 <- importPeaks(fl, type = "narrow")
  bp1 <- importPeaks(fl, type = "broad")
  expect_true(is(c(np1, bp1), "GRangesList"))
  expect_equal(vapply(c(np1, bp1), length, integer(1)), c(0, 0))

  # And with a seqinfo
  sq <- Seqinfo(seqnames = "chr1")
  np2 <- importPeaks(fl, type = "narrow", seqinfo = sq)
  bp2 <- importPeaks(fl, type = "broad", seqinfo = sq)
  expect_true(is(c(np1, bp1), "GRangesList"))
  expect_equal(vapply(c(np1, bp1), length, integer(1)), c(0, 0))
})
