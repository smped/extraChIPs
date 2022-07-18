
plotOverlaps <- function(data, col = hcl.colors(3, "Zissou")) {


    ## Collect areas & totals, defining which intersections we have
    areas <- .calculateAreas(data)
    reps <- intersect(c("x", "y", "z"), names(areas))
    grps <- names(areas)
    pairs <- intersect(c("xy", "xz", "yz"), grps)
    three_way <- "xyz" %in% grps
    totals <- sapply(
        grps,
        function(x){
            i <- grepl(x, names(areas))
            sum(unlist(areas[i]))
        },
        simplify = FALSE
    )
    ## xz will be missed in the above sapply
    if ("xz" %in% grps) totals$xz <- areas$xz + areas$xyz

    ## Define Colours & Labels
    nm <- names(data)
    if (length(nm) != length(data)) nm <- c("x", "y", "z")[seq_along(data)]
    col <- rep(col, length(data))[seq_along(data)]
    names(col) <- reps

    ## From here is taken largely from Biovenn::draw.venn. By Tim Hulsen
    ## https://cran.r-project.org/web/packages/BioVenn/index.html

    ## Plotting area
    p <- list(w = 1, h = 1)
    ## Amplification factor
    amp <- 1e5 / sum(unlist(areas))

    ## Radius calculations
    r <- lapply(reps, function(x) .sqr(totals[[x]] * amp / pi))
    names(r) <- reps

    ## Distances
    d <- sapply(
        pairs,
        function(x) {
            i <- gsub("^([xyz])([xyz])", "\\1", x)
            j <- gsub("^([xyz])([xyz])", "\\2", x)
            .pairwiseDistance(r[[i]], r[[j]], totals[[x]] * amp)
        },
        simplify = FALSE
    )

    ## Modify distance for horizontal relationships
    d <- sapply(
        names(d),
        function(x) min(d[[x]], sum(unlist(d)) - d[[x]]),
        simplify = FALSE
    )

    ## Angle calculations
    a <- sapply(
        reps,
        function(x) {
            inc <- grepl(x, names(d))
            num <- sum(vapply(d[inc], .sq, numeric(1))) - .sq(unlist(d[!inc]))
            denom <- 2 * prod(unlist(d[inc]))
            .arccos(num / denom)
        },
        simplify = FALSE
    )
    yz <- list(x = d$xz * sin(a$z), y = d$xy * cos(a$y))

    ## PPU
    h <- list()
    h$w <- sum(
        max(r$y + yz$y, r$x, r$z - d$yz + yz$y),
        max(r$x, r$y - yz$y, r$z + d$yz - yz$y)
    )
    h$ppu <- p$w / h$w
    v <- list()
    v$w <- sum(max(r$x + yz$x, r$y, r$z), max(r$y, r$z, r$x - yz$x))
    v$ppu <- p$h / v$w
    ppu <- min(h$ppu, v$ppu)

    ## Centres
    centres <- list(
        x = c(
            max(r$x, r$y + yz$y, r$x - d$yz + yz$y),
            max(r$x, r$y - yz$x, r$z - yz$x)
        ),
        y = c(
            max(r$x - yz$y, r$y, r$z - d$yz),
            max(r$x + yz$x, r$y, r$z)
        ),
        z = c(
            max(r$x + d$yz - yz$y, r$y + d$yz, r$z),
            max(r$x + yz$x, r$y, r$z)
        )
    )
    ## Now scale for plotting
    centres <- lapply(centres, function(x) x * ppu)
    r <- lapply(r, function(x) x * ppu)
    ## The original parameterisation
    # h$x <- max(r$x, r$y + yz$y, r$x - d$yz + yz$y)
    # v$x <- max(r$x, r$y - yz$x, r$z - yz$x)
    # h$y <- max(r$x - yz$y, r$y, r$z - d$yz)
    # v$y <- max(r$x + yz$x, r$y, r$z)
    # h$z <- max(r$x + d$yz - yz$y, r$y + d$yz, r$z)
    # v$z <- max(r$x + yz$x, r$y, r$z)
    browser()

    ## Intersection Points
    # # Calculate intersection points X-Y (first inner, then outer)
    # xy_i_h_part1 = (x_h + y_h) / 2 + ((y_h - x_h) * (x_r^2 - y_r^2)) / (2 * xy_d^2)
    # xy_i_v_part1 = (x_v + y_v) / 2 + ((y_v - x_v) * (x_r^2 - y_r^2)) / (2 * sq(xy_d))
    # xy_i_h_part2 = 2 * ((x_v - y_v) / sq(xy_d)) * sqr((xy_d + x_r + y_r) * (xy_d + x_r - y_r) * (xy_d - x_r + y_r) * (-xy_d + x_r + y_r)) / 4
    # xy_i_v_part2 = 2 * ((x_h - y_h) / sq(xy_d)) * sqr((xy_d + x_r + y_r) * (xy_d + x_r - y_r) * (xy_d - x_r + y_r) * (-xy_d + x_r + y_r)) / 4
    # xy_i1_h = xy_i_h_part1 - xy_i_h_part2
    # xy_i1_v = xy_i_v_part1 + xy_i_v_part2
    # xy_i2_h = xy_i_h_part1 + xy_i_h_part2
    # xy_i2_v = xy_i_v_part1 - xy_i_v_part2

    ## Special Cases with missing overlaps

    ## Draw circles
    grid.newpage()
    lapply(
        reps,
        function(x) grid.circle(
            centres[[x]][1], centres[[x]][2], r[[x]],
            gp = gpar(fill = col[[x]], alpha = 0.5)
        )
    )


    ## Add group labels

    ## Add counts

}
ex <- list(
    x = letters[1:10],
    y = letters[6:15],
    z = c(letters[5], letters[10], letters[15:26])
)
plotOverlaps(ex)

#' @importFrom methods is
.calculateAreas <- function(data) {

    ## Limit this to 3way interactions. Anything bigger is better as an UpSet plot
    stopifnot(is(data, "list"))
    stopifnot(length(data) <= 3)
    stopifnot(all(vapply(data, is.character, logical(1))))

    n <- length(data)
    areas <- lapply(data, length)
    names(areas) <- c("x", "y", "z")[seq_along(data)]
    i <- "x"
    if (n == 2) {
        areas$xy <- length(intersect(data[[1]], data[[2]]))
        areas$x <- areas$x - areas$xy
        areas$y <- areas$y - areas$xy
        i <- c("x", "y", "xy")
    }
    if (n == 3) {
        areas$xyz <- length(intersect(intersect(data[[1]], data[[2]]), data[[3]]))
        areas$xy <- length(intersect(data[[1]], data[[2]])) - areas$xyz
        areas$xz <- length(intersect(data[[1]], data[[3]])) - areas$xyz
        areas$yz <- length(intersect(data[[2]], data[[3]])) - areas$xyz
        areas$x <- areas$x - areas$xy - areas$xz - areas$xyz
        areas$y <- areas$y - areas$xy - areas$yz - areas$xyz
        areas$z <- areas$z - areas$xz - areas$yz - areas$xyz
        i <- c("x", "y", "z", "xy", "xz", "yz", "xyz")
    }
    areas[i]

}
ex <- list(
    x = letters[1:10],
    y = letters[6:15],
    z = letters[10:26]
)
.calculateAreas(ex)


.pairwiseDistance <- function(r1, r2, ol, .err = 1e-3) {
    .getVals <- function(r1, r2, d) {
        s1 <- .sq(r1) * .arccos((.sq(d) + .sq(r1) - .sq(r2))/(2 * d * r1))
        s2 <- .sq(r2) * .arccos((.sq(d) + .sq(r2) - .sq(r1))/(2 * d * r2))
        s3 <- 0.5 * .sqr(
            round(
                (d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (-d + r1 + r2), 5
            )
        )
        s1 + s2 - s3
    }
    d <- r1 + r2
    thresh <- ol * (1 - .err)
    while(
        thresh > .getVals(r1, r2, d)
    ) {
        d <- d - min(r1, r2) * .err
    }
    d
}



#' Utility functions
.sq <- function(x) x^2
.sqr <- function(x) sqrt(round(abs(x)))
.arccos <- function(x, .digits = 5) {
    acos(
        min(
            max(-1, round(x, .digits)),
            1
        )
    )
}
