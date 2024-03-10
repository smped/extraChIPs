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

test_that("GInteractions are converted correctly", {
    hic <- InteractionSet::GInteractions(
        GRanges("chr1:1-10"), GRanges("chr1:201-210")
    )
    # Test without mcols
    tbl <- as_tibble(hic)
    expect_equal(names(tbl), paste0("anchor", 1:2))
    expect_true(all(vapply(tbl, is, logical(1), class2 = "character")))
    tbl <- as_tibble(hic, rangeAsChar = FALSE)
    expect_equal(
        names(tbl),
        paste0(
            c("seqnames", "start", "end", "width", "strand"),
            rep(c(".x", ".y"), each = 5)
        )
    )
    hic$id <- "interaction1"
    tbl <- as_tibble(hic)
    expect_equal(names(tbl), c("anchor1", "anchor2", "id"))

})

test_that("Empty objects return empty tibble objects", {
    expect_equal(as_tibble(DataFrame()), tibble())
    expect_equal(as_tibble(GRanges()), tibble(range = character()))
    expect_equal(
        as_tibble(InteractionSet::GInteractions()),
        tibble(anchor1 = character(), anchor2 = character())
    )
})
