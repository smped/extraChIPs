gr1 <- GRanges("chr1:1-10")
gr2 <- GRanges("chr1:21-30")
hic <- InteractionSet::GInteractions(gr1, gr2)
hiccol <- list(anchors = "lightblue", interactions = "red")
feat <- GRangesList(a = GRanges("chr1:8"))
genes <- gr1
mcols(genes) <- DataFrame(
    feature = "exon", gene = "ENSG", exon = "ENSE", transcript = "ENST",
    symbol = "ID"
)
cyto_df <- data.frame(
    chrom = "chr1", chromStart = 1, chromEnd = 100, name = "p1", gieStain = "gneg"
)

test_path <- system.file("tests", "test.bw", package = "rtracklayer")
test_bw <- rtracklayer::BigWigFileList(test_path)
names(test_bw) <- "a"

test_that("Correct plotting works", {
    data(grch37.cytobands)
    gr <- GRanges("chr1:11869-12227")
    feat_gr <- GRangesList(
        Promoter = GRanges("chr1:11800-12000"),
        Enhancer = GRanges("chr1:13000-13200")
    )
    hic <- InteractionSet::GInteractions(feat_gr$Promoter, feat_gr$Enhancer)
    genes <- c("chr1:11869-12227:+", "chr1:12613-12721:+", "chr1:13221-14409:+")
    genes <- GRanges(genes)
    mcols(genes) <- DataFrame(
        feature = "exon", gene = "ENSG00000223972", exon = 1:3,
        transcript = "ENST00000456328", symbol = "DDX11L1"
    )
    p <- plotHFGC(
        gr, hic = hic, features = feat_gr, genes = genes,
        zoom = 2, cytobands = grch37.cytobands, rotation.title = 90,
        featcol = c(Promoter = "red", Enhancer = "yellow")
    )
    expect_true(is(p, "list"))
    expect_equal(length(p), 6)
    expect_true(is(p[[1]], "IdeogramTrack"))
    expect_true(is(p[[2]], "GenomeAxisTrack"))
    expect_true(is(p[[3]], "InteractionTrack"))
    expect_true(is(p[[4]], "AnnotationTrack"))
    expect_true(is(p[[5]], "GeneRegionTrack"))
    expect_true(is(p[[6]], "ImageMap"))

    p <- plotHFGC(
        GRanges("chr2:1-1000"), coverage = test_bw, cytobands = grch37.cytobands,
        axistrack = FALSE, annotation = GRangesList(up = GRanges("chr2:501-505")),
        annotcol = list(up = "red"), ylim = c(-1, 5)
    )
    expect_equal(length(p), 4)
    expect_true(is(p[[2]], "AnnotationTrack"))
    expect_equal(p[[2]]@dp@pars$fill, c(up = "red"))
    expect_equal(p[[2]]@name, "")
    expect_true(is(p[[3]], "DataTrack"))
    expect_equal(p[[3]]@dp@pars$ylim, c(-1, 5))

})

test_that(".makeIdeoTrack behaves as expected", {
    expect_message(.makeIdeoTrack(gr1), "Could not find cytogenetic bands")
    genome(gr1) <- "GRCh37"
    tr <- .makeIdeoTrack(gr1, .fontsize = 12)
    expect_equal(tr@genome, setNames("hg19", "GRCh37"))
    expect_equal(dim(tr@bandTable), c(1136, 5))
    expect_equal(tr@dp@pars$fontsize, 12)
})
