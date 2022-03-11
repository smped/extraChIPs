test_that("rangeAsChar behaves correctly", {
    gr <- GRanges("chr1:1-10")
    gr$p <- 0.2
    df <- as_tibble(gr)
    expect_true(is(df, "tbl_df"))
    expect_equal(c("range", "p"), colnames(df))
    expect_equal(ncol(as_tibble(gr, rangeAsChar = FALSE)), 6)
    expect_error(as_tibble(gr, name = "p"), "A column named")
})

test_that("DFrame conversion is successful", {
    df <- DataFrame(
        id = seq_along(letters),
        listCol = as(as.list(letters), "CharacterList"),
        gr = GRanges(paste0("chr1:", seq_along(letters))),
        row.names = LETTERS
    )
    tbl <- as_tibble(df)
    expect_equal(df$id, tbl$id)
    expect_true(is(tbl$listCol, "list"))
    expect_true(is(tbl$gr, "character"))
    expect_equal(
        colnames(as_tibble(df, rangeAsChar = FALSE)),
        c("id", "listCol", "gr.seqnames", "gr.start", "gr.end", "gr.width",
          "gr.strand")
    )
    expect_equal(rownames(tbl), as.character(seq_along(letters)))
    tbl2 <- as_tibble(df, rownames = "rownames")
    expect_equal(tbl2$rownames, LETTERS)

})

test_that("Seqinfo objects are successfully coerced", {
    sq <- Seqinfo("chr1", 10, FALSE, "test")
    tbl <- tibble(
        seqnames = "chr1", seqlengths = 10,
        is_circular = FALSE, genome = "test"
    )
    expect_equal(as_tibble(sq), tbl)
})
