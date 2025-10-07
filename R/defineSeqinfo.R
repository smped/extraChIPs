#' @title Use package data to define a Seqinfo object
#'
#' @description
#' Use package data to define a Seqinfo object
#'
#' @details
#' This function will create a Seqinfo object from pre-defined data from the
#' Genome Reference Consortium.
#' Returned objects will always be restricted to assembled molecules only.
#' Currently implemented genome builds represent the four most common builds
#' for ChIP-Seq analysis
#'
#'
#' @param build The Genome build used
#' @param chr logical(1) Include the prefix "chr"
#' @param mito Specify M or MT to include the mitochondrial chromosome. Omitted
#' by default
#' @param ... Not used
#'
#' @return A Seqinfo object
#'
#' @examples
#' defineSeqinfo("GRCh37", TRUE)
#' defineSeqinfo("GRCh37", FALSE, "MT")
#'
#'
#' @export
defineSeqinfo <- function(
        build = c(
            "GRCh38", "GRCh37", "GRCm39", "GRCm38", "hg19", "hg38",
            "T2T-CHM13v2.0", "mm39","mm10"
        ),
        chr = TRUE, mito, ...
){
    build <- match.arg(build)
    if (build == "T2T-CHM13v2.0") build <- "T2T"
    stopifnot(is.logical(chr))
    mito_arg <- NULL
    if (!missing(mito)) {
        stopifnot(mito %in% c("M", "MT"))
        mito_arg <- paste0(", '", mito, "'")
    }
    fn <- paste0(".make_", build, "(", chr, mito_arg, ")")
    eval(parse(text = fn))
}

#' @importFrom GenomeInfoDb Seqinfo
.make_T2T <- function(chr) {
    seqnames <- c(seq_len(22), "X", "Y")
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- as.integer(
        c(
        248387328, 242696752, 201105948, 193574945, 182045439, 172126628,
        160567428, 146259331, 150617247, 134758134, 135127769, 133324548,
        113566686, 101161492, 99753195, 96330374, 84276897, 80542538,
        61707364, 66210255, 45090682, 51324926, 154259566, 62460029
    )[seq_along(seqnames)])
    circ <- c(rep(FALSE, 24))[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "T2T-CHM13v2.0")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_GRCh38 <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(22), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
        170805979L, 159345973L, 145138636L, 138394717L, 133797422L, 135086622L,
        133275309L, 114364328L, 107043718L, 101991189L, 90338345L, 83257441L,
        80373285L, 58617616L, 64444167L, 46709983L, 50818468L, 156040895L,
        57227415L, 16569L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 24), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "GRCh38")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_hg38 <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(22), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
        170805979L, 159345973L, 145138636L, 138394717L, 133797422L, 135086622L,
        133275309L, 114364328L, 107043718L, 101991189L, 90338345L, 83257441L,
        80373285L, 58617616L, 64444167L, 46709983L, 50818468L, 156040895L,
        57227415L, 16569L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 24), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "hg38")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_GRCh37 <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(22), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        249250621L, 243199373L, 198022430L, 191154276L, 180915260L,
        171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L,
        133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L,
        78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L,
        59373566L, 16571L, 16569L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 24), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "GRCh37")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_hg19 <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(22), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        249250621L, 243199373L, 198022430L, 191154276L, 180915260L,
        171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L,
        133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L,
        78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L,
        59373566L, 16571L, 16569L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 24), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "hg19")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_GRCm39  <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(19), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        195471971L, 182113224L, 160039680L, 156508116L, 151834684L,
        149736546L, 145441459L, 129401213L, 124595110L, 130694993L, 122082543L,
        120129022L, 120421639L, 124902244L, 104043685L, 98207768L, 94987271L,
        90702639L, 61431566L, 171031299L, 91744698L, 16299L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 21), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "GRCm39")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_mm39  <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(19), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        195471971L, 182113224L, 160039680L, 156508116L, 151834684L,
        149736546L, 145441459L, 129401213L, 124595110L, 130694993L, 122082543L,
        120129022L, 120421639L, 124902244L, 104043685L, 98207768L, 94987271L,
        90702639L, 61431566L, 171031299L, 91744698L, 16299L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 21), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "mm39")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_GRCm38  <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(19), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        197195432L, 181748087L, 159599783L, 155630120L, 152537259L,
        149517037L, 152524553L, 131738871L, 124076172L, 129993255L, 121843856L,
        121257530L, 120284312L, 125194864L, 103494974L, 98319150L, 95272651L,
        90772031L, 61342430L, 166650296L, 15902555L, 16299L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 21), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "GRCm38")
}

#' @importFrom GenomeInfoDb Seqinfo
.make_mm10  <- function(chr, mito = NULL) {
    seqnames <- c(seq_len(19), "X", "Y", mito)
    if (chr) seqnames <- paste0("chr", seqnames)
    seqlen <- c(
        197195432L, 181748087L, 159599783L, 155630120L, 152537259L,
        149517037L, 152524553L, 131738871L, 124076172L, 129993255L, 121843856L,
        121257530L, 120284312L, 125194864L, 103494974L, 98319150L, 95272651L,
        90772031L, 61342430L, 166650296L, 15902555L, 16299L
    )[seq_along(seqnames)]
    circ <- c(rep(FALSE, 21), TRUE)[seq_along(seqnames)]
    Seqinfo(seqnames, seqlen, circ, "mm10")
}



