# This script contains test moved to longtest
# They can be temporarily activated by uncommenting here and pushing without
# changing the version number. That will update the code_cov metrics to double
# check that key lines are hit, and can then be de-activate again along with
# a version bump to initiate a BioC build.


####################
## getProfileData ##
####################

bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:500")

# Moved to longtest
test_that("Paths behave correctly", {
  expect_s4_class(getProfileData(bw, gr, upstream = 10, bins = 10), "GRanges")
  expect_s4_class(getProfileData(c(bw, bw), gr, upstream = 10, bins = 10), "GRangesList")
})

##################
## mapByFeature ##
##################

## We need a combination of ranges which do & don't map
genes <- GRanges(c("chr1:2-10:*", "chr1:25-30:-", "chr1:31-40:+"))
genes$gene_id <- paste0("gene", seq_along(genes))

## Add a promoter for each gene
prom <- promoters(genes, upstream = 1, downstream = 1)

## Add an enhancer in gene2 & to the right of gene3
enh <- GRanges(c("chr1:20", "chr1:50"))
enh$which <- c("gene2", "gene1 (gi)")

## Map the last enhancer to gene1's promoter
gi <- InteractionSet::GInteractions(prom[1], enh[2])
mcols(gi) <- c()

## Now define key ranges
gr <- GRanges(paste0("chr1:", seq(0, 50, by = 10)))
gr$via_ol <- CharacterList(
    list(NULL, "gene1", NULL, c("gene2", "gene3"), "gene3", NULL)
)
gr$via_prom <- CharacterList(
    list("gene1", NULL, NULL, c("gene2", "gene3"), NULL, NULL)
)
gr$via_enh_gap5 <- CharacterList(
    list(NULL, NULL, "gene2", NULL, NULL, NULL)
)
gr$via_gi <- CharacterList(
    list("gene1", NULL, NULL, NULL, NULL, "gene1")
)

test_that("mapByFeature produces the correct output", {
    ol_only <- mapByFeature(gr, genes, cols = "gene_id", gr2gene = 0)$gene_id
    expected <- new(
        "CompressedCharacterList", elementType = "character",
        elementMetadata = NULL, metadata = list(),
        unlistData = c("gene1", "gene2", "gene3"),
        partitioning = new(
            "PartitioningByEnd", end = c(0L, 1L, 1L, 2L, 3L, 3L),
            NAMES = NULL, elementType = "ANY", elementMetadata = NULL, metadata = list()
        )
    )
    expect_equal(ol_only, expected)
    ## Update promoter only mappings
    expected[[1]] <- "gene1"; expected[[4]] <- c("gene2", "gene3")
    expect_equal(
        mapByFeature(gr, genes, prom = prom, cols = "gene_id", gr2gene = 0)$gene_id,
        expected
    )
    ## Update enhancer only mappings
    expected[[3]] <- "gene2"
    expect_equal(
        mapByFeature(
            gr, genes, prom = prom, enh = enh, cols = "gene_id",
            gr2gene = 0, enh2gene = 5
        )$gene_id,
        expected
    )
    ## Update interaction mappings
    expected[[6]] <- "gene1"
    expect_equal(
        mapByFeature(
            gr, genes, prom = prom, enh = enh, gi = gi,
            cols = "gene_id", gr2gene = 0, enh2gene = 5
        )$gene_id,
        expected
    )

})

##############
## plotHFGC ##
##############

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

