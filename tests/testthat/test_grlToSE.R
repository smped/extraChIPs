sq <- Seqinfo("chr1", 100)
a <- GRanges("chr1:1-10", seqinfo = sq)
a$feature <- "Gene"
a$logFC <- 1
b <- GRanges("chr1:6-15", seqinfo = sq)
b$feature <- "Promoter"
b$logFC <- -1
grl <- GRangesList(a = a, b = b)
merged <- GenomicRanges::reduce(unlist(grl))

test_that("grlToSE returns correct object", {
  se <- grlToSE(
    grl, assayCols = "logFC", metaCols = "feature", keyvals = "logFC"
  )
  expect_equal(
    assay(se, "logFC"),
    matrix(c(1, -1), ncol = 2, dimnames = list(c(), c("a", "b")))
  )
  expect_equal(
    rowRanges(se)$feature, as(list(c("Gene", "Promoter")),"CompressedList")
  )
})

test_that(".cols2Assay returned expected values", {
  expect_equal(.cols2Assays(merged, grl, c()), SimpleList())
  test_merge <- .cols2Assays(merged, grl, "logFC", "logFC", "min", TRUE)
  expect_true(is(test_merge, "SimpleList"))
  expect_equal(
    test_merge[[1]],
    matrix(c(1, -1), ncol = 2, dimnames = list(c(), c("a", "b")))
  )
})

test_that(".addMcols returnes expected values", {
  expect_equal(merged, .addMcols(merged, c()))
  test_merge <- .addMcols(merged, "feature", grl, TRUE)
  expect_true(is(test_merge$feature, "CompressedCharacterList"))
  expect_equal(test_merge$feature[[1]], c("Gene", "Promoter"))
})

test_that(".checkArgsGrlToSe catches everything", {
  expect_true(
    .checkArgsGrlToSe(
      .x = grl, .assayCols = "logFC", .metaCols = "feature",
      .keyvals = "logFC", .all_mcols = c("logFC", "feature")
    )
  )
  expect_warning(
    .checkArgsGrlToSe(
      .x = setNames(grl, c()), .assayCols = "logFC", .metaCols = "feature",
      .keyvals = "logFC", .all_mcols = c("logFC", "feature")
    ), "All elements of the x should be given informative names"
  )
  expect_warning(
    .checkArgsGrlToSe(
      .x = grl, .assayCols = "", .metaCols = "feature",
      .keyvals = "logFC", .all_mcols = c("logFC", "feature")
    ), "Some columns requested as assays.+"
  )
  expect_warning(
    .checkArgsGrlToSe(
      .x = grl, .assayCols = "logFC", .metaCols = "",
      .keyvals = "logFC", .all_mcols = c("logFC", "feature")
    ), "Some columns requested as mcols.+"
  )
  expect_warning(
    .checkArgsGrlToSe(
      .x = grl, .assayCols = "logFC", .metaCols = "feature",
      .keyvals = c(), .all_mcols = c("logFC", "feature")
    ), "Missing keyvals. Values will be chosen at random"
  )
  expect_warning(
    fail_no_keys <- .checkArgsGrlToSe(
      .x = grl, .assayCols = "logFC", .metaCols = "feature",
      .keyvals = "", .all_mcols = c("logFC", "feature")
    ), "No specified 'keyvals'were found"
  )
  expect_false(fail_no_keys)
  expect_warning(
    fail <- .checkArgsGrlToSe(
      .x = grl, .assayCols = "logFC", .metaCols = "logFC",
      .keyvals = "logFC", .all_mcols = c("logFC", "feature")
    ), "Columns for assays and mcols must be distinct."
  )
  expect_false(fail)
})

test_that("Empty SE is correctly structured", {
  expect_warning(se <- .emptySE(GRanges(seqinfo = sq)), "Returned object.+")
  expect_equal(length(se), 0)
  expect_equal(sq, seqinfo(se))
})
