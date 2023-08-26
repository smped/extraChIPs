test_that("GRCh38 works correctly", {
    sq <- defineSeqinfo()
    grch38 <- new(
        "Seqinfo",
        seqnames = c(
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
            "chrX", "chrY"
        ),
        seqlengths = c(
            248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
            170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
            135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
            90338345L, 83257441L, 80373285L, 58617616L, 64444167L, 46709983L,
            50818468L, 156040895L, 57227415L
        ),
        is_circular = rep(FALSE,  24), genome = rep("GRCh38", 24)
    )
    expect_equal(sq, grch38)
})

test_that("chr & mito work correctly", {
    sq <- defineSeqinfo("GRCh38", chr = FALSE, "MT")
    grch38 <- new(
        "Seqinfo",
        seqnames = c(seq_len(22), "X", "Y", "MT"),
        seqlengths = c(
            248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
            170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
            135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
            90338345L, 83257441L, 80373285L, 58617616L, 64444167L, 46709983L,
            50818468L, 156040895L, 57227415L, 16569L
        ),
        is_circular = c(rep(FALSE,  24), TRUE), genome = rep("GRCh38", 25)
    )
    expect_equal(sq, grch38)
})

test_that("Errors correctly", {
    expect_error(defineSeqinfo("hg19"))
    expect_error(defineSeqinfo(chr = NULL))
    expect_error(defineSeqinfo(mito = ""))
})
