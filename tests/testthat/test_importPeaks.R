test_that("narrowPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "chipExtra"
  )
  expect_true(.isValidNarrow(fl))
  expect_false(.isValidBroad(fl))

  gr <- importPeaks(fl, type = "narrow")
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 5L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))
})

test_that("broadPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.broadPeak", package = "chipExtra"
  )
  expect_false(.isValidNarrow(fl))
  expect_true(.isValidBroad(fl))
  gr <- importPeaks(fl, type = "broad")
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 4L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))
})

test_that("seqinfo objects behave correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "chipExtra"
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
    "extdata/testFiles", "test.narrowPeak", package = "chipExtra"
  )
  expect_error(
    importPeaks(gr, "narrow", blacklist = NULL)
  )
  expect_equal(
    length(importPeaks(fl, blacklist = GRanges("chr1:1299700"))),
    1L
  )
})
