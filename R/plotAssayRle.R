#' @title Plot RLE for a given assay within a SummarizedExperiment
#'
#' @description Plot RLE for a given assay within a SummarizedExperiment
#'
#' @details
#' Uses ggplot2 to create an RLE plot for the selected assay. Any numerical
#' transformation prior to performing the RLE can be specified using the
#' `trans` argument
#'
#' @return
#' A ggplot2 object
#'
#' @param x A SummarizedExperiment object
#' @param assay The assay to plot
#' @param colour Column from `colData(x)` to outline the boxplots
#' @param fill Column from `colData(x)` to fill the boxplots
#' @param rle_group Column from `colData(x)` to calculate RLE within groups
#' Commonly an alternative sample label.
#' @param n_max Maximum number of points to plot
#' @param trans character(1). NUmerical transformation to apply to the data
#' prior to RLE calculation
#' @param ... Not used
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayRle(se, "counts", fill = "treat")
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr left_join group_by mutate
#' @importFrom rlang syms '!!!' sym
#' @importFrom stats median
#' @import ggplot2
#'
#' @rdname plotAssayRle-methods
#' @aliases plotAssayRle
#' @export
setMethod(
    "plotAssayRle",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour = NULL, fill = NULL, rle_group = NULL,
        n_max = Inf, trans = c(), ...
    ) {

        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        df$sample <- rownames(df)
        args <- colnames(df)
        if (!is.null(colour)) colour <- sym(match.arg(colour, args))
        if (!is.null(fill)) fill <- sym(match.arg(fill, args))
        if (!is.null(rle_group)) rle_group <- match.arg(rle_group, args)
        grps <- c("ind", rle_group)

        n_max <- min(nrow(x), n_max)
        ind <- seq_len(n_max)
        if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)

        mat <- assay(x, assay)[ind,]
        if (!is.null(trans)) {
            mat <- match.fun(trans)(mat)
            trans_ok <- all(
                is.matrix(mat), nrow(mat) == length(ind),
                colnames(mat) == colnames(x)
            )
            if (!trans_ok) stop("This transformation is not applicable")
        }

        mat_df <- as.data.frame(mat)
        mat_df$ind <- ind
        mat_df <- pivot_longer(
            mat_df,
            all_of(colnames(x)), names_to = "sample", values_to = "assay"
        )
        mat_df <- left_join(mat_df, df, by = "sample")
        mat_df <- group_by(mat_df, !!!syms(grps))
        mat_df <- mutate(mat_df, rle = assay - median(assay))

        sample <- rle <- NULL ## R CMD check
        ggplot(
            mat_df,
            aes(sample, rle, fill = {{ fill }}, colour = {{ colour }})
        ) +
            geom_boxplot() +
            labs(y = "RLE")

    }
)
