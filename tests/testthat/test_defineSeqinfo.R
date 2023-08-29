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

test_that("GRCh37 works correctly", {
    sq <- defineSeqinfo("GRCh37")
    grch37 <- new(
        "Seqinfo",
        seqnames = c(
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
            "chrX", "chrY"
        ),
        seqlengths = c(
            249250621L, 243199373L, 198022430L, 191154276L, 180915260L,
            171115067L, 159138663L, 146364022L, 141213431L, 135534747L,
            135006516L, 133851895L, 115169878L, 107349540L, 102531392L,
            90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 48129895L,
            51304566L, 155270560L, 59373566L
        ),
        is_circular = rep(FALSE,  24), genome = rep("GRCh37", 24)
    )
    expect_equal(sq, grch37)
})

test_that("GRCm38 works correctly", {
    sq <- defineSeqinfo("GRCm38")
    grcm38 <- new(
        "Seqinfo",
        seqnames = c(
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"
        ),
        seqlengths = c(
            197195432L, 181748087L, 159599783L, 155630120L, 152537259L,
            149517037L, 152524553L, 131738871L, 124076172L, 129993255L, 121843856L,
            121257530L, 120284312L, 125194864L, 103494974L, 98319150L, 95272651L,
            90772031L, 61342430L, 166650296L, 15902555L
        ),
        is_circular = rep(FALSE,  21), genome = rep("GRCm38", 21)
    )
    expect_equal(sq, grcm38)
})

test_that("GRCm39 works correctly", {
    sq <- defineSeqinfo("GRCm39")
    grcm39 <- new(
        "Seqinfo",
        seqnames = c(
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"
        ),
        seqlengths = c(
            195471971L, 182113224L, 160039680L, 156508116L, 151834684L,
            149736546L, 145441459L, 129401213L, 124595110L, 130694993L, 122082543L,
            120129022L, 120421639L, 124902244L, 104043685L, 98207768L, 94987271L,
            90702639L, 61431566L, 171031299L, 91744698L
        ),
        is_circular = rep(FALSE,  21), genome = rep("GRCm39", 21)
    )
    expect_equal(sq, grcm39)
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
