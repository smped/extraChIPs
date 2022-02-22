i <- c(1, 1:9)
j <- 1:10

test_that("setdiffWithMCols produces a List as expected", {
  ## Return a CharacterList
  cl <- .returnListColumn(letters[j], i, j)
  expect_equal(length(cl), max(i))
  expect_equal(length(cl[[1]]), 2)
  expect_true(is(cl, "CompressedCharacterList"))
})

test_that("setdiffWithMcols handles List input", {
  x <- as(as.list(letters[j]), "CompressedCharacterList")
  x[[1]] <- c(x[[1]], "z")
  cl <- .returnListColumn(x, i, j)
  expect_equal(length(cl[[1]]), 3)
  expect_true(is(cl, "CompressedCharacterList"))
})

test_that("setdiffWithMCols doesn't return a list when not expected", {
  cv <- .returnListColumn(letters[1:10], 1:10, 1:10)
  expect_true(is(cv, "character"))
  expect_equal(length(cv), 10)
})

test_that("setdiffWithMcols passes to default when no mcols provided", {
  x <- GRanges("chr1:1-10")
  y <- GRanges("chr1:5")
  sd <- setdiffWithMcols(x, y)
  expect_equal(
    mcols(sd),
    new(
      "DFrame", rownames = NULL, nrows = 2L,
      listData = structure(list(), .Names = character(0)),
      elementType = "ANY", elementMetadata = NULL, metadata = list()
    )
  )
})
