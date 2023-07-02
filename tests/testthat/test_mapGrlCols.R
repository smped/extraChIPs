a <- GRanges(paste0("chr1:", seq(1, 51, length.out = 10)))
width(a) <- 4
a$logFC <- rnorm(length(a))
g <- as.list(paste0("Gene", seq_along(a)))
g[[1]] <- c("Gene0", g[[1]])
a$genes <- as(g, "CompressedList")
b <- GRanges(paste0("chr1:", floor(seq(15, 50, length.out = 5))))
width(b) <- 8
b$logFC <- rnorm(length(b))
g <- vector("list", length(b))
g[[length(b)]] <- a$genes[[length(a)]]
b$genes <- as(g, "CompressedList")
grl <- GRangesList(A = a, B = b)

test_that("mapGrlCols, works as expected", {
    gr <- mapGrlCols(grl, var = "logFC")
    expect_equal(.mcolnames(gr), c("A_logFC", "B_logFC"))

    gr <- mapGrlCols(grl, var = c("logFC", "genes"))
    expect_equal(.mcolnames(gr), c("A_logFC", "B_logFC", "A_genes", "B_genes"))
    expect_true(is(gr$A_genes, "List"))
    expect_true(is(gr$A_logFC, "numeric"))

    gr <- mapGrlCols(grl, var = c("logFC", "genes"), collapse = "genes")
    expect_equal(.mcolnames(gr), c("A_logFC", "B_logFC", "genes"))
    expect_true(is(gr$genes, "character"))

    expect_equal(mapGrlCols(grl), granges(makeConsensus(grl)))

})

test_that("collapse_sep works for list objects", {
    gr <- mapGrlCols(grl, var = c("logFC", "genes"), collapse = "genes", collapse_sep = list(genes = "; "))
    expect_true(any(grepl("; ", gr$genes)))
    expect_error(
        mapGrlCols(grl, collapse = "genes", collapse_sep = list(g = "; ")),
        "All columns being collapsed need a separator provided"
    )
})

test_that(".coerceList works", {
    x <- sample(500, 20)
    y0 <- splitAsList(x, x %% 4)
    y1 <- .coerceList(y0)
    expect_true(is(y1, "list"))
})
