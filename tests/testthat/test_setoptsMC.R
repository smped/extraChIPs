i <- c(1, 1:9)
j <- 1:10

test_that("setdiffMC produces a List as expected", {
  ## Return a CharacterList
  cl <- .returnListColumn(letters[j], i, j, .simplify = FALSE)
  expect_equal(length(cl), max(i))
  expect_equal(length(cl[[1]]), 2)
  expect_true(is(cl, "CompressedCharacterList"))
})

test_that("setdiffMC handles List input", {
  x <- as(as.list(letters[j]), "CompressedCharacterList")
  x[[1]] <- c(x[[1]], "z")
  cl <- .returnListColumn(x, i, j, .simplify = TRUE)
  expect_equal(length(cl[[1]]), 3)
  expect_true(is(cl, "CompressedCharacterList"))
})

test_that("setdiffMC doesn't return a list when not expected", {
  cv <- .returnListColumn(letters[1:10], 1:10, 1:10, .simplify = TRUE)
  expect_true(is(cv, "character"))
  expect_equal(length(cv), 10)
})

test_that(".mapMcols2Ranges returns empty mcols when required", {
  x <- GRanges("chr1:1-10")
  y <- GRanges("chr1:5")
  sd <- .mapMcols2Ranges(x, y, TRUE, TRUE)
  expect_equal(
    mcols(sd),
    new(
      "DFrame", rownames = NULL, nrows = 1L,
      listData = structure(list(), .Names = character(0)),
      elementType = "ANY", elementMetadata = NULL, metadata = list()
    )
  )
})

test_that("reduceMC returns a list when expected", {

  x <- GRanges(c("chr1:1-10:+", "chr1:6-12:-"))
  x$id <- c("range1", "range2")
  id <- reduceMC(x, ignore.strand = TRUE)$id
  expect_true(is(id, "CompressedCharacterList"))
  expect_equal(id[[1]], c("range1", "range2"))

  id <- reduceMC(x, ignore.strand = FALSE)$id
  expect_equal(id, c("range1", "range2"))

})
